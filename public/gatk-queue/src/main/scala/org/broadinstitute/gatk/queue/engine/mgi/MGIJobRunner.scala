/*
* Copyright 2016 The McDonnell Genome Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.engine.mgi

import org.broadinstitute.gatk.queue.QException
import org.broadinstitute.gatk.queue.util.{Logging,Retry}
import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.engine.{RunnerStatus, CommandLineJobRunner}
import org.ggf.drmaa._

import scala.sys.process._
import java.io._
import java.util.{Date, Collections}

/**
 * Runs jobs using DRMAA.
 */
class MGIJobRunner(val session: Session, val function: CommandLineFunction) extends CommandLineJobRunner with Logging {
  /** Job Id of the currently executing job. */
  var jobId: String = _
  override def jobIdString = jobId

  // Don't cleanup for now
  override def cleanup() {
  }

  // Set the display name to < 512 characters of the description
  // NOTE: Not sure if this is configuration specific?
  protected val jobNameLength = 500
  protected val jobNameFilter = """[^A-Za-z0-9_]"""

  protected def functionNativeSpec() : String = {
    var nativeSpec: String = ""

    val jobName = function.jobRunnerJobName.take(jobNameLength).replaceAll(jobNameFilter, "_")
    nativeSpec += " -J " + jobName

    // If the job queue is set specify the job queue
    if (function.jobQueue != null)
      nativeSpec += " -q " + function.jobQueue

    // Set the current working directory
    if (function.commandDirectory.getPath != null)
      nativeSpec += " -cwd " + function.commandDirectory.getPath

    // Set the output file for stdout
    nativeSpec += " -oo " + function.jobOutputFile.getPath

    // If the error file is set specify the separate output for stderr
    // Otherwise the error is joined with stdout
    if (function.jobErrorFile != null) 
      nativeSpec += " -eo " + function.jobErrorFile.getPath

    // Override any command line given jobNativeArgs
    nativeSpec += " " + function.jobNativeArgs.mkString(" ")

    nativeSpec = nativeSpec.trim()
    logger.info("Native spec is: %s".format(nativeSpec))
    nativeSpec
  }

  def start() {
    session.synchronized {
      // Instead of running the function.commandLine, run "/bin/bash <jobScript>"
      val baseCmd : String = "/bin/bash"
      val script : String = jobScript.toString
      val cmd : String = Array(baseCmd, script).mkString(" ")
      logger.info("Base command is: %s".format(cmd))

      val (exitCode, stdOut, stdErr) = bsub(cmd)

      if (exitCode != 0) {
        throw new QException("Unable to submit job (BSUB error): " + stdErr)
      } else {
        logger.info("Successfully Invoked BSUB: " + stdOut)

        // store the id so it can be killed in tryStop
        val regex = "Job <([0-9]+)> .+".r
        val regex(jobNum) = stdOut.trim()
        jobId = jobNum

        logger.info("Submitted job id: " + jobId)
      }

      updateStatus(RunnerStatus.RUNNING)
    }
  }

  def updateJobStatus() = {
    session.synchronized {
      var returnStatus: RunnerStatus.Value = null

      try {
        val jobStatus = session.getJobProgramStatus(jobId);
        jobStatus match {
          case Session.QUEUED_ACTIVE => returnStatus = RunnerStatus.RUNNING
          case Session.DONE =>
            val lsfStatus = getLSFCacheJobStatus(jobId)
            lsfStatus match {
              case "DONE" => returnStatus = RunnerStatus.DONE
              case "EXIT" => returnStatus = RunnerStatus.FAILED
              case _ => returnStatus = RunnerStatus.FAILED  /* when in an unknown state, just fail */
            }
          case Session.FAILED => returnStatus = RunnerStatus.FAILED
          case Session.UNDETERMINED => logger.warn("Unable to determine status of job id " + jobId)
          case _ => returnStatus = RunnerStatus.RUNNING
        }
      } catch {
        // getJobProgramStatus will throw an exception once wait has run, as the
        // job will be reaped.  If the status is currently DONE or FAILED, return
        // the status.
        case de: DrmaaException =>
          if (lastStatus == RunnerStatus.DONE || lastStatus == RunnerStatus.FAILED)
            returnStatus = lastStatus
          else
            logger.warn("Unable to determine status of job id " + jobId, de)
      }

      if (returnStatus != null) {
        updateStatus(returnStatus)
        true
      } else {
        false
      }
    }
  }

  def tryStop() {
    session.synchronized {
      // Assumes that after being set the job may be
      // reassigned but will not be reset back to null
      if (jobId != null) {
        try {
          // Stop runners. SIGTERM(15) is preferred to SIGKILL(9).
          // Only way to send SIGTERM is for the Sys Admin set the terminate_method
          // resource of the designated queue to SIGTERM
          session.control(jobId, Session.TERMINATE)
        } catch {
          case e: Exception =>
            logger.error("Unable to kill job " + jobId, e)
        }
      }
    }
  }

  def getLSFCacheJobStatus(jobID: String): String = {
    val (exitValue, stdOut, stdErr) = bjobs(jobID)
    var status: String = ""
    if (stdErr != "") {
      logger.info("LSF Job: %s -- status not in LSF cache".format(jobID))
      status = "unknown"
    } else {
      status = stdOut.toString.split("\n")(1).split("\\s+")(2)
    }
    status
  }

  def bjobs(jobID: String): (Int, String, String) = {
    val bjobsCmd = Array("bjobs", jobID).mkString(" ")
    logger.info("(bjobs): %s".format(bjobsCmd))
    val (exitValue, stdOut, stdErr) = runCmd(bjobsCmd)
    (exitValue, stdOut.trim(), stdErr.trim())
  }

  def bsub(cmd: String): (Int, String, String) = {
    val bsubCmd = Array("bsub", functionNativeSpec(), cmd).mkString(" ")
    logger.info("(bsub): %s".format(bsubCmd))
    val (exitValue, stdOut, stdErr) = runCmd(bsubCmd)
    (exitValue, stdOut, stdErr)
  }

  def runCmd(cmd: String): (Int, String, String) = {
    val stdout = new ByteArrayOutputStream
    val stderr = new ByteArrayOutputStream
    val stdoutWriter = new PrintWriter(stdout)
    val stderrWriter = new PrintWriter(stderr)
    val exitValue = cmd.!(ProcessLogger(stdoutWriter.println, stderrWriter.println))
    stdoutWriter.close()
    stderrWriter.close()
    (exitValue, stdout.toString, stderr.toString)
  }
}
