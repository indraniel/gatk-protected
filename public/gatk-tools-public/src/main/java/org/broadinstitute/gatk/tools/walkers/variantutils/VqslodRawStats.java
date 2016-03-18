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

import java.io.*;
import java.util.Collection;

/**
 * Created by idas on 3/17/16.
 */

public class VqslodRawStats extends RodWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Input(fullName="snps-file", shortName = "SNP", doc="The file name to hold the SNP data", required=true)
    protected String snpsFile;

    @Input(fullName="indels-file", shortName = "INDEL", doc="The file name to hold the INDEL data", required=true)
    protected String indelFile;

    // file streams
    PrintStream snpStream = null;
    PrintStream indelStream = null;

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
        out.println("\t (INDEL) Simple Insertions: " + nInsertions);
        out.println("\t (INDEL) Simple Deletions: " + nDeletions);
        out.println("\t (INDEL) Complex INDEL: " + nComplex);
        out.println("\t MNPs: " + nMNPs);
        out.println("\t Mixed: " + nMixed);
        out.println("\t Symbolic: " + nSymbolic);

        snpStream.close();
        indelStream.close();
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

        // initialize file streams
        snpStream = openStream(snpsFile);
        indelStream = openStream(indelFile);

        // print out SNP/INDEL file headers
        snpStream.printf("#%s\n", String.join("\t", SNPRecordHeaders()));
        indelStream.printf("#%s\n", String.join("\t", INDELRecordHeaders()));
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
                "TRANSVERSION"
        };
    }

    private String[] INDELRecordHeaders() {
        return new String[]{
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "FILTER",
                "VQSLOD",
                "TYPE",
                "IndelLength"
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
                recordINDEL(vc);
                break;
            case MIXED:
                nMixed++;
                break;
            case SYMBOLIC:
                nSymbolic++;
                break;
            case NO_VARIATION:
                String ref = vc.getReference().getDisplayString();
                String alt = vc.getAlternateAllele(0).toString();
                int pos = vc.getStart();
                String msg = String.format("(REF: %s | ALT: %s | POS: %d)", ref, alt, pos);
                out.println("Unexpected Variant type " + vc.getType() + "shouldn't get here! " + msg);
                //throw new ReviewedGATKException("Unexpected Variant type " + vc.getType() + "shouldn't get here! " + msg);
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

        String[] row = new String[] {
                vc.getContig(),
                String.format("%d", vc.getStart()),
                vc.getID(),
                vc.getReference().getDisplayString(),
                vc.getAlternateAllele(0).toString(),
                site_filter,
                vqslod,
                String.format("%d", transition),
                String.format("%d", transversion)
        };

        snpStream.printf( "%s\n", String.join("\t", row) );
    }

    public void recordINDEL(VariantContext vc) {
        String indelType = deriveIndelType(vc);
        int indelLength = calculateIndelLength(vc);
        String site_filter = getVariantFilterString(vc);
        String vqslod = getVQSLODString(vc);

        String[] row = new String[] {
                vc.getContig(),
                String.format("%d", vc.getStart()),
                vc.getID(),
                vc.getReference().getDisplayString(),
                vc.getAlternateAllele(0).toString(),
                site_filter,
                vqslod,
                indelType,
                String.format("%d", indelLength)
        };

        indelStream.printf( "%s\n", String.join("\t", row) );
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

    public int calculateIndelLength(VariantContext vc) {
        Allele alt = vc.getAlternateAllele(0);
        Allele ref = vc.getReference();

        final int alleleSize = alt.length() - ref.length();

        if ( alleleSize == 0 ) throw new ReviewedGATKException(
            "Allele size not expected to be zero for indel: alt = " + alt +
            " ref = " + ref.toString() +
            " POS = " + vc.getStart()
        );

        return alleleSize;
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
