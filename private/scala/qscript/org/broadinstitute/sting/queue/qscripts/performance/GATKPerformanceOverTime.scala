package org.broadinstitute.sting.queue.qscripts.performance

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.extensions.gatk._
import java.lang.Math
import org.broadinstitute.sting.utils.baq.BAQ.CalculationMode

class GATKPerformanceOverTime extends QScript {
  val STD_RESULTS_DIR = "/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/gatkPerformanceOverTime"

  @Argument(shortName = "results", doc="results", required=false)
  val resultsDir: File = new File("runResults");

  @Argument(shortName = "test", doc="test", required=false)
  val TEST: Boolean = false;

  @Argument(shortName = "resources", doc="resources", required=true)
  val resourcesDir: String = ""

  @Argument(shortName = "myJarFile", doc="Path to the current GATK jar file", required=false)
  val myJarFile: String = "/home/unix/depristo/dev/GenomeAnalysisTK/projects/gatkPerformance/dist/GenomeAnalysisTK.jar"

  @Argument(shortName = "iterations", doc="it", required=false)
  val iterations: Int = 3;

  val nIterationsForSingleTestsPerIteration: Int = 1;

  @Argument(shortName = "maxThreads", doc="maxThreads", required=false)
  val maxThreads: Int = 6;

  @Argument(shortName = "steps", doc="steps", required=false)
  val steps: Int = 10;

  @Argument(shortName = "maxNSamples", doc="maxNSamples", required=false)
  val maxNSamples: Int = 1000000;

  val RECAL_BAM_FILENAME = "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.bam"
  val dbSNP_FILENAME = "dbsnp_132.b37.vcf"
  val RECAL_FILENAME = "NA12878.HiSeq.WGS.bwa.cleaned.recal.hg19.20.csv"
  val b37_FILENAME = "human_g1k_v37.fasta"

  def makeResource(x: String): File = new File("%s/%s".format(resourcesDir, x))
  def makeChunk(x: Int): File = makeResource("chunk_%d.vcf".format(x))
  def COMBINE_FILES: List[File] = Range(1,10).map(makeChunk).toList

  class AssessmentParameters(val name: String, val bamList: File, val intervals: File, val nSamples: Int, val dcov: Int)

  // TODO -- count the number of lines in the bam.list file
  val WGSAssessment = new AssessmentParameters("WGS", "wgs.bam.list.local.list", "wgs.bam.list.intervals", 1103, 50)
  val WExAssessment = new AssessmentParameters("WEx", "wex.bam.list.local.list", "wex.bam.list.intervals", 140, 500)

  val assessments = List(WGSAssessment, WExAssessment)

  val GATK12 = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-1.2-65-ge4a583a/GenomeAnalysisTK.jar"
  val GATK13 = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-1.3-23-g13905c0/GenomeAnalysisTK.jar"
  val GATKs = Map(
    "v1.2" -> GATK12,
    "v1.3" -> GATK13,
    "v1.cur" -> myJarFile) // TODO -- how do I get this value?

  trait UNIVERSAL_GATK_ARGS extends CommandLineGATK {
    this.logging_level = "INFO";
    this.reference_sequence = makeResource(b37_FILENAME);
    this.memoryLimit = 4
  }

  def script = {
    if ( ! resultsDir.exists ) resultsDir.mkdirs()

    for ( iteration <- 0 until iterations ) {
      for ( assess <- assessments ) {
        for (nSamples <- divideSamples(assess.nSamples) ) {
          val sublist = new SliceList(assess.name, nSamples, makeResource(assess.bamList))
          if ( iteration == 0 ) add(sublist) // todo - remove condition when Queue bug is fixed
          for ( (gatkName, gatkJar) <- GATKs ) {
            val name: String = "assess.%s_gatk.%s_iter.%d".format(assess.name, gatkName, iteration)

            trait VersionOverrides extends CommandLineGATK {
              this.jarFile = gatkJar
              this.dcov = assess.dcov

              // special handling of test intervals
              if ( TEST )
                this.intervalsString :+= "20:10,000,000-10,001,000"
              else
                this.intervals :+= assess.intervals

              this.configureJobReport(Map(
                "iteration" -> iteration,
                "gatk" -> gatkName,
                "nSamples" -> nSamples,
                "assessment" -> assess.name))
            }

            // SNP calling
            add(new Call(sublist.list, nSamples, name) with VersionOverrides)

            // CountLoci
            add(new MyCountLoci(sublist.list, nSamples, name) with VersionOverrides)

            if ( nSamples == assess.nSamples )
              addMultiThreadedTest(() => new Call(sublist.list, nSamples, name) with VersionOverrides)

            //add(new CallWithSamtools(sublist.list, nSamples, name, assess.intervals))
          }
        }
      }

      for ( subiteration <- 0 until nIterationsForSingleTestsPerIteration ) {
        for ( (gatkName, gatkJar) <- GATKs ) {
          if ( gatkName != "1.2" ) { // todo generalize

            { // BQSR
              trait VersionOverrides extends CommandLineGATK {
                this.jarFile = gatkJar
                this.dcov = 50
                this.intervalsString :+= (if ( TEST ) "20:10,000,000-10,001,000" else "20:10,000,000-20,000,000")
                this.configureJobReport(Map( "iteration" ->
                  (subiteration + (iteration * nIterationsForSingleTestsPerIteration)).toString,
                  "gatk" -> gatkName))
              }

              add(new CountCov(makeResource(RECAL_BAM_FILENAME)) with VersionOverrides)
              add(new Recal(makeResource(RECAL_BAM_FILENAME)) with VersionOverrides)

              addMultiThreadedTest(() => new CountCov(makeResource(RECAL_BAM_FILENAME)) with VersionOverrides)
            }

            { // Standard VCF tools
              trait VersionOverrides extends CommandLineGATK {
                this.jarFile = gatkJar
                this.intervalsString = (if ( TEST ) List("1:10,000,000-10,010,000") else List("1", "2", "3", "4", "5"))
                this.configureJobReport(Map( "iteration" -> iteration, "gatk" -> gatkName))
              }

              val CV = new CombineVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
              CV.variant = COMBINE_FILES
              CV.out = new File("/dev/null")
              add(CV)

              val SV = new SelectVariants with UNIVERSAL_GATK_ARGS with VersionOverrides
              SV.variant = makeResource("chunk_1.vcf")
              SV.sample_name = List("HG00096") // IMPORTANT THAT THIS SAMPLE BE IN CHUNK ONE
              SV.out = new File("/dev/null")
              add(SV)

              def makeVE(): CommandLineGATK = {
                val VE = new VariantEval with UNIVERSAL_GATK_ARGS with VersionOverrides
                VE.eval :+= makeResource("chunk_1.vcf")
                VE.out = new File("/dev/null")
                VE.comp :+= new TaggedFile(makeResource(dbSNP_FILENAME), "dbSNP")
                VE
              }

              add(makeVE())
              addMultiThreadedTest(makeVE)
            }
          }
        }
      }
    }
  }

  def addMultiThreadedTest(makeCommand: () => CommandLineGATK) {
    if ( maxThreads > 1 ) {
      for ( nt <- Range(1,maxThreads+1) ) {
        val cmd = makeCommand()
        cmd.nt = nt
        cmd.addJobReportBinding("nt", nt)
        cmd.analysisName = cmd.analysisName + ".nt"
        add(cmd)
      }
    }
  }

  def divideSamples(nTotalSamples: Int): List[Int] = {
    val maxLog10: Double = Math.log10(Math.min(maxNSamples, nTotalSamples))
    val stepSize: Double = maxLog10 / steps
    val ten: Double = 10.0
    def deLog(x: Int): Int = Math.round(Math.pow(ten, stepSize * x)).toInt
    dedupe(Range(0, steps+1).map(deLog).toList)
  }

  class Call(@Input(doc="foo") bamList: File, n: Int, name: String) extends UnifiedGenotyper with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outVCF: File = new File("/dev/null")
    this.input_file :+= bamList
    this.stand_call_conf = 10.0
    this.o = outVCF
    this.baq = org.broadinstitute.sting.utils.baq.BAQ.CalculationMode.RECALCULATE
  }

  class MyCountLoci(@Input(doc="foo") bamList: File, n: Int, name: String) extends CountLoci with UNIVERSAL_GATK_ARGS {
    @Output(doc="foo") var outFile: File = new File("/dev/null")
    this.input_file :+= bamList
    this.o = outFile
  }

  class SliceList(prefix: String, n: Int, @Input bamList: File) extends CommandLineFunction {
    this.analysisName = "SliceList"
    @Output(doc="foo") var list: File = new File("%s/%s.bams.%d.list".format(resultsDir.getPath, prefix, n))
    def commandLine = "head -n %d %s | awk '{print \"%s/\" $1}' > %s".format(n, bamList, resourcesDir, list)
  }

  class CallWithSamtools(@Input(doc="foo") bamList: File, n: Int, name: String, intervals: File) extends CommandLineFunction {
    this.analysisName = "samtools"
    @Output(doc="foo") var outVCF: File = new File("/dev/null")
    def commandLine = "./samtools mpileup -ub %s -f %s -l %s | ./bcftools view -vcg -l %s - > %s".format(bamList, b37_FILENAME, intervals, intervals, outVCF)
  }

  class CountCov(inBam: File) extends CountCovariates with UNIVERSAL_GATK_ARGS {
    this.knownSites :+= makeResource(dbSNP_FILENAME)
    this.covariate ++= List("ReadGroupCovariate", "QualityScoreCovariate", "CycleCovariate", "DinucCovariate")
    this.input_file :+= inBam
    this.recal_file = new File("/dev/null")
  }

  class Recal(inBam: File) extends TableRecalibration with UNIVERSAL_GATK_ARGS {
    this.input_file :+= inBam
    this.recal_file = makeResource(RECAL_FILENAME)
    this.baq = CalculationMode.CALCULATE_AS_NECESSARY
    this.out = new File("/dev/null")
  }

  def dedupe(elements:List[Int]):List[Int] = {
    if (elements.isEmpty)
      elements
    else
      elements.head :: dedupe(for (x <- elements.tail if x != elements.head) yield x)
  }
}
