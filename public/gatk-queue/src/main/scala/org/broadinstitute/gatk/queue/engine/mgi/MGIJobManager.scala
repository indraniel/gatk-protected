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

import org.broadinstitute.gatk.queue.function.CommandLineFunction
import org.broadinstitute.gatk.queue.engine.CommandLineJobManager
import org.broadinstitute.gatk.utils.jna.drmaa.v1_0.JnaSessionFactory
import org.ggf.drmaa.Session

/**
 * Runs job using bsub command line, but monitor with existing DRMAA
 */
class MGIJobManager extends CommandLineJobManager[MGIJobRunner] {
  protected var session: Session = _

  protected def newSession() = new JnaSessionFactory().getSession
  protected def contact = null

  override def init() {
    session = newSession()
    session.init(contact)
  }

  override def exit() {
    session.exit()
  }

  def runnerType = classOf[MGIJobRunner]
  def create(function: CommandLineFunction) = new MGIJobRunner(session, function)

  override def updateStatus(runners: Set[MGIJobRunner]) = {
    var updatedRunners = Set.empty[MGIJobRunner]
    runners.foreach(runner => if (runner.updateJobStatus()) {updatedRunners += runner})
    updatedRunners
  }
  override def tryStop(runners: Set[MGIJobRunner]) {
    runners.foreach(_.tryStop())
  }
}
