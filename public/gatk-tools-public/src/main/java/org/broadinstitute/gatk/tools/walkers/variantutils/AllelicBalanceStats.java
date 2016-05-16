package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.engine.walkers.RodWalker;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;

import java.io.*;
import java.util.Collection;
import java.util.HashMap;

/**
 * Created by idas on 5/10/16.
 */

public class AllelicBalanceStats extends RodWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Input(fullName="snps-file", shortName = "SNP", doc="The file name to hold the SNP data", required=true)
    protected String snpsFile;

    // file streams
    PrintStream snpStream = null;

    // number of total variants counter
    private Integer nProcessedLoci;

    // number of SNP variants counter
    private Integer nSNPs;

    // number of INDEL variants counter
    private Integer nInsertions;
    private Integer nDeletions;
    private Integer nComplex;

    // number of OTHER variants counter
    private Integer nMNPs;
    private Integer nMixed;
    private Integer nSymbolic;
    private Integer nNoVariation;


    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker == null) {
            return 1;
        }
        Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());

        if ( vcs == null || vcs.isEmpty()) {
            return 0;
        }

        for (VariantContext vc : vcs) {
            nProcessedLoci++;
            processVariantSite(vc);
        }

        return null;
    }

    @Override
    public Integer reduceInit() {
        return null;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        out.println("Processed " + nProcessedLoci + " variants.");
        out.println("\t SNPs: " + nSNPs);
        out.println("\t INDELs: " + (nInsertions + nDeletions + nComplex));
        out.println("\t\t (INDEL) Simple Insertions: " + nInsertions);
        out.println("\t\t (INDEL) Simple Deletions: " + nDeletions);
        out.println("\t\t (INDEL) Complex: " + nComplex);
        out.println("\t MNPs: " + nMNPs);
        out.println("\t Mixed: " + nMixed);
        out.println("\t Symbolic: " + nSymbolic);
        out.println("\t No Variation: " + nNoVariation);

        snpStream.close();
    }

    public void initialize() {
        // initialize counters;
        nProcessedLoci = 0;
        nSNPs = 0;
        nInsertions = 0;
        nDeletions = 0;
        nComplex = 0;
        nMNPs = 0;
        nMixed = 0;
        nSymbolic = 0;
        nNoVariation = 0;

        // initialize file streams
        snpStream = openStream(snpsFile);

        // print out SNP/INDEL file headers
        snpStream.printf("#%s\n", String.join("\t", SNPRecordHeaders()));
    }

    private String[] SNPRecordHeaders() {
        return new String[]{
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "FILTER",
                "VQSLOD",
                "TRANSITION",
                "TRANSVERSION",

                "GENOTYPE_00_SAMPLE_CNT",
                "GENOTYPE_00_REF_ALLELE",
                "GENOTYPE_00_ALT_ALLELE",
                "GENOTYPE_00_TOTAL_DEPTH",
                "GENOTYPE_00_ALLELE_BALANCE",

                "GENOTYPE_01_SAMPLE_CNT",
                "GENOTYPE_01_REF_ALLELE",
                "GENOTYPE_01_ALT_ALLELE",
                "GENOTYPE_01_TOTAL_DEPTH",
                "GENOTYPE_01_ALLELE_BALANCE",

                "GENOTYPE_11_SAMPLE_CNT",
                "GENOTYPE_11_REF_ALLELE",
                "GENOTYPE_11_ALT_ALLELE",
                "GENOTYPE_11_TOTAL_DEPTH",
                "GENOTYPE_11_ALLELE_BALANCE",

                "GENOTYPE_NOCALL_SAMPLE_CNT",
                "GENOTYPE_MIXED_SAMPLE_CNT"
        };
    }

    public PrintStream openStream(String filename) {
        File file = new File(filename);
        PrintStream stream = null;
        try {
            stream = new PrintStream(file);
        } catch (FileNotFoundException e) {
            throw new GATKException(String.format("Unable to open file '%s' for writing", file), e);
        }

        return stream;
    }

    public void processVariantSite(VariantContext vc) {
        // for now, we're going to assume that the VCF file has already been decomposed
        // and normalized, making everything biallelic
        switch ( vc.getType() ) {
            case SNP:
                nSNPs++;
                recordSNP(vc);
                break;
            case MNP:
                nMNPs++;
                break;
            case INDEL:
                deriveIndelType(vc);
                break;
            case MIXED:
                nMixed++;
                break;
            case SYMBOLIC:
                nSymbolic++;
                break;
            case NO_VARIATION:
                nNoVariation++;
                break;
            default:
                throw new ReviewedGATKException("Unexpected VariantContext type " + vc.getType());
        }
    }

    public void recordSNP(VariantContext vc) {
        int transition, transversion;

        if ( GATKVariantContextUtils.isTransition(vc) ) {
            transition = 1;
            transversion = 0;
        }
        else {
            transition = 0;
            transversion = 1;
        }

        String site_filter = getVariantFilterString(vc);
        String vqslod = getVQSLODString(vc);

        HashMap<String, String> gt_counts = collectGenotypeMetrics(vc);

        String[] row = new String[] {
                vc.getContig(),
                String.format("%d", vc.getStart()),
                vc.getID(),
                vc.getReference().getDisplayString(),
                vc.getAlternateAllele(0).toString(),
                site_filter,
                vqslod,
                String.format("%d", transition),
                String.format("%d", transversion),

                gt_counts.get("gt_00_sample_cnt"),
                gt_counts.get("gt_00_ref_allele"),
                gt_counts.get("gt_00_alt_allele"),
                gt_counts.get("gt_00_total_depth"),
                gt_counts.get("gt_00_allele_bal"),

                gt_counts.get("gt_01_sample_cnt"),
                gt_counts.get("gt_01_ref_allele"),
                gt_counts.get("gt_01_alt_allele"),
                gt_counts.get("gt_01_total_depth"),
                gt_counts.get("gt_01_allele_bal"),

                gt_counts.get("gt_11_sample_cnt"),
                gt_counts.get("gt_11_ref_allele"),
                gt_counts.get("gt_11_alt_allele"),
                gt_counts.get("gt_11_total_depth"),
                gt_counts.get("gt_11_allele_bal"),

                gt_counts.get("gt_nocall_sample_cnt"),
                gt_counts.get("gt_mixed_sample_cnt"),

        };

        snpStream.printf( "%s\n", String.join("\t", row) );
    }

    public HashMap<String, String> collectGenotypeMetrics(VariantContext vc) {
        int gt_00_sample_cnt = 0;
        int gt_00_ref_allele = 0;
        int gt_00_alt_allele = 0;
        float gt_00_allele_bal = (float) 0.0;

        int gt_01_sample_cnt = 0;
        int gt_01_ref_allele = 0;
        int gt_01_alt_allele = 0;
        float gt_01_allele_bal = (float) 0.0;

        int gt_11_sample_cnt = 0;
        int gt_11_ref_allele = 0;
        int gt_11_alt_allele = 0;
        float gt_11_allele_bal = (float) 0.0;

        int gt_nocall_sample_cnt = 0;

        int gt_mixed_sample_cnt = 0;

        if ( !vc.hasGenotypes() ) {
            errorVariant(vc, "Unable to find genotypes");
        }

        GenotypesContext genotypes = vc.getGenotypes();

        for ( Genotype genotype : genotypes ) {
            String gt = genotype.getGenotypeString();

            if ( genotype.isHomRef() ) {
                gt_00_sample_cnt++;
                int[] AD = genotype.getAD();
                gt_00_ref_allele += AD[0];
                gt_00_alt_allele += AD[1];
            }
            else if ( genotype.isHomVar() ) {
                gt_11_sample_cnt++;
                int[] AD = genotype.getAD();
                gt_11_ref_allele += AD[0];
                gt_11_alt_allele += AD[1];
            }
            else if ( genotype.isHet() ) {
                gt_01_sample_cnt++;
                int[] AD = genotype.getAD();
                gt_01_ref_allele += AD[0];
                gt_01_alt_allele += AD[1];
            }
            else if ( genotype.isNoCall() ) {
                gt_nocall_sample_cnt++;
            }
            else if ( genotype.isMixed() ) {
                gt_mixed_sample_cnt++;
            }
            else {
                String msg = String.format(
                        "Unable to ascertain genotype %s on (%s)",
                        genotype.getGenotypeString(),
                        genotype.getSampleName()
                );
                errorVariant(vc, msg);
            }
        }

        if ( (gt_00_ref_allele + gt_00_alt_allele) != 0 ) {
            gt_00_allele_bal = (float) gt_00_alt_allele / (gt_00_ref_allele + gt_00_alt_allele);
        }

        if ( (gt_01_ref_allele + gt_01_alt_allele) != 0 ) {
            gt_01_allele_bal = (float) gt_01_alt_allele / (gt_01_ref_allele + gt_01_alt_allele);
        }

        if ( (gt_11_ref_allele + gt_11_alt_allele) != 0 ) {
            gt_11_allele_bal = (float) gt_11_alt_allele / (gt_11_ref_allele + gt_11_alt_allele);
        }

        HashMap<String, String> hmap = new HashMap<String, String>();

        hmap.put("gt_00_sample_cnt", String.format("%d", gt_00_sample_cnt));
        hmap.put("gt_00_ref_allele", String.format("%d", gt_00_ref_allele));
        hmap.put("gt_00_alt_allele", String.format("%d", gt_00_alt_allele));
        hmap.put("gt_00_allele_bal", String.format("%.4f", gt_00_allele_bal));
        hmap.put("gt_00_total_depth", String.format("%d", gt_00_ref_allele + gt_00_alt_allele));

        hmap.put("gt_01_sample_cnt", String.format("%d", gt_01_sample_cnt));
        hmap.put("gt_01_ref_allele", String.format("%d", gt_01_ref_allele));
        hmap.put("gt_01_alt_allele", String.format("%d", gt_01_alt_allele));
        hmap.put("gt_01_allele_bal", String.format("%.4f", gt_01_allele_bal));
        hmap.put("gt_01_total_depth", String.format("%d", gt_01_ref_allele + gt_01_alt_allele));

        hmap.put("gt_11_sample_cnt", String.format("%d", gt_11_sample_cnt));
        hmap.put("gt_11_ref_allele", String.format("%d", gt_11_ref_allele));
        hmap.put("gt_11_alt_allele", String.format("%d", gt_11_alt_allele));
        hmap.put("gt_11_allele_bal", String.format("%.4f", gt_11_allele_bal));
        hmap.put("gt_11_total_depth", String.format("%d", gt_11_ref_allele + gt_11_alt_allele));

        hmap.put("gt_nocall_sample_cnt", String.format("%d", gt_nocall_sample_cnt));
        hmap.put("gt_mixed_sample_cnt", String.format("%d", gt_mixed_sample_cnt));

        return hmap;
    }

    public void errorVariant(VariantContext vc, String msg) {
        String chrom = vc.getContig();
        String position = String.format("%d", vc.getStart());
        String id = vc.getID();
        throw new GATKException(
                String.format("On chrom: %s | position: %s | ID: %s -- %s", chrom, position, id, msg)
        );
    }

    public String deriveIndelType(VariantContext vc) {
        String indelType = null;

        if (vc.isSimpleInsertion()) {
            indelType = "I";
            nInsertions++;
        }
        else if (vc.isSimpleDeletion()) {
            indelType = "D";
            nDeletions++;
        }
        else {
            indelType = "C";
            nComplex++;
        }

        return indelType;
    }

    public String getVariantFilterString(VariantContext vc) {
        String site_filter;

        if ( vc.getFilters().isEmpty() ) {
            site_filter = "PASS";
        }
        else {
            int numElements = vc.getFilters().size();
            String[] site_filters = vc.getFilters().toArray(new String[numElements]);
            site_filter = site_filters[0];
        }

        return site_filter;
    }

    public String getVQSLODString(VariantContext vc) {
        String vqslod = null;

        if ( vc.hasAttribute("VQSLOD") ) {
            vqslod = vc.getAttribute("VQSLOD").toString();
        }
        else {
            vqslod = "NA";
        }

        return vqslod;
    }
}
