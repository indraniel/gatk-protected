/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

@Analysis(description = "Evaluation summary for indels")
public class IndelSummary extends VariantEvaluator implements StandardEval {
    final protected static Logger logger = Logger.getLogger(IndelSummary.class);

    //
    // counts of snps and indels
    //
    @DataPoint(description = "Number of SNPs", format = "%d")
    public int n_SNPs = 0;

    @DataPoint(description = "Number of singleton SNPs", format = "%d")
    public int n_singleton_SNPs = 0;

    @DataPoint(description = "Number of indels", format = "%d")
    public int n_indels = 0;

    @DataPoint(description = "Number of singleton indels", format = "%d")
    public int n_singleton_indels = 0;

    //
    // gold standard
    //
    @DataPoint(description = "Number of Indels overlapping gold standard sites", format = "%d")
    public int n_indels_matching_gold_standard = 0;

    @DataPoint(description = "Percent of indels overlapping gold standard sites")
    public String gold_standard_matching_rate;

    //
    // multi-allelics
    //
    // Number of Indels Sites (counts one for any number of alleles at site)
    @DataPoint(description = "Number of Indels Sites [counts one for any number of alleles at site]")
    public int nIndelSites = 0;

    @DataPoint(description = "Number of sites with where the number of alleles is greater than 2")
    public int n_multiallelic_indel_sites = 0;

    @DataPoint(description = "Percent of indel sites that are multi-allelic")
    public String percent_of_sites_with_more_than_2_alleles;

    //
    // snp : indel ratios
    //
    @DataPoint(description = "SNP to indel ratio")
    public String SNP_to_indel_ratio;

    @DataPoint(description = "Singleton SNP to indel ratio")
    public String SNP_to_indel_ratio_for_singletons;

    //
    // novelty
    //
    @DataPoint(description = "Number of novel indels", format = "%d")
    public int n_novel_indels = 0;

    @DataPoint(description = "Indel novelty rate")
    public String indel_novelty_rate;

    //
    // insertions to deletions
    //
    @DataPoint(description = "Number of insertion indels")
    public int n_insertions = 0;

    @DataPoint(description = "Number of deletion indels")
    public int n_deletions = 0;

    @DataPoint(description = "Insertion to deletion ratio")
    public String insertion_to_deletion_ratio;

    @DataPoint(description = "Number of large (>10 bp) deletions")
    public int n_large_deletions = 0;

    @DataPoint(description = "Number of large (>10 bp) insertions")
    public int n_large_insertions = 0;

    @DataPoint(description = "Ratio of large (>10 bp) insertions to deletions")
    public String insertion_to_deletion_ratio_for_large_indels;

    //
    // Frameshifts
    //
    @DataPoint(description = "Number of indels in protein-coding regions labeled as frameshift")
    public int n_coding_indels_frameshifting = 0;

    @DataPoint(description = "Number of indels in protein-coding regions not labeled as frameshift")
    public int n_coding_indels_in_frame = 0;

    @DataPoint(description = "Frameshift percent")
    public String frameshift_rate_for_coding_indels;

    //
    // Het : hom ratios
    //
    @DataPoint(description = "Het to hom ratio for SNPs")
    public String SNP_het_to_hom_ratio;

    @DataPoint(description = "Het to hom ratio for indels")
    public String indel_het_to_hom_ratio;
    
    @DataPoint(description = "Number of SNP Heterozygous")
    public int nSNPHets = 0;

    @DataPoint(description = "Number of SNP Homozygous")
    public int nSNPHoms = 0;

    @DataPoint(description = "Number of Indel Heterozygous")
    public int nIndelHets = 0;

    @DataPoint(description = "Number of Indel Homozygous")
    public int nIndelHoms = 0;

    @DataPoint(description = "Insertion Count By Length 1")
    public int nInsertionLength_1_bp;

    @DataPoint(description = "Insertion Count By Length 2")
    public int nInsertionLength_2_bp;

    @DataPoint(description = "Insertion Count By Length 3")
    public int nInsertionLength_3_bp;

    @DataPoint(description = "Deletion Count By Length 1")
    public int nDeletionLength_1_bp;

    @DataPoint(description = "Deletion Count By Length 2")
    public int nDeletionLength_2_bp;

    @DataPoint(description = "Deletion Count By Length 3")
    public int nDeletionLength_3_bp;

    int[] insertionCountByLength = new int[]{0, 0, 0, 0}; // note that the first element isn't used
    int[] deletionCountByLength = new int[]{0, 0, 0, 0}; // note that the first element isn't used

    // - Since 1 & 2 bp insertions and 1 & 2 bp deletions are equally likely to cause a
    // downstream frameshift, if we make the simplifying assumptions that 3 bp ins
    // and 3bp del (adding/subtracting 1 AA in general) are roughly comparably
    // selected against, we should see a consistent 1+2 : 3 bp ratio for insertions
    // as for deletions, and certainly would expect consistency between in/dels that
    // multiple methods find and in/dels that are unique to one method  (since deletions
    // are more common and the artifacts differ, it is probably worth looking at the totals,
    // overlaps and ratios for insertions and deletions separately in the methods
    // comparison and in this case don't even need to make the simplifying in = del functional assumption

    @DataPoint(description = "ratio of 1 and 2 bp insertions to 3 bp insertions")
    public String ratio_of_1_and_2_to_3_bp_insertions;

    @DataPoint(description = "ratio of 1 and 2 bp deletions to 3 bp deletions")
    public String ratio_of_1_and_2_to_3_bp_deletions;

    public final static int LARGE_INDEL_SIZE_THRESHOLD = 10;

    @Override public int getComparisonOrder() { return 2; }

    public void update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval == null || (getWalker().ignoreAC0Sites() && eval.isMonomorphicInSamples()) )
            return;

        // update counts
        switch ( eval.getType() ) {
            case SNP:
                n_SNPs += eval.getNAlleles() - 1; // -1 for ref
                if ( variantWasSingleton(eval) ) n_singleton_SNPs++;

                // collect information about het / hom ratio
                for ( final Genotype g : eval.getGenotypes() ) {
                    if ( g.isHet() ) nSNPHets++;
                    if ( g.isHomVar() ) nSNPHoms++;
                }
                break;
            case INDEL:
                final VariantContext gold = getWalker().goldStandard == null ? null : tracker.getFirstValue(getWalker().goldStandard);

                nIndelSites++;
                if ( ! eval.isBiallelic() ) n_multiallelic_indel_sites++;

                // collect information about het / hom ratio
                for ( final Genotype g : eval.getGenotypes() ) {
                    if ( g.isHet() ) nIndelHets++;
                    if ( g.isHomVar() ) nIndelHoms++;
                }

                for ( Allele alt : eval.getAlternateAlleles() ) {
                    n_indels++; // +1 for each alt allele
                    if ( variantWasSingleton(eval) ) n_singleton_indels++;
                    if ( comp == null ) n_novel_indels++; // TODO -- make this test allele specific?
                    if ( gold != null ) n_indels_matching_gold_standard++;

                    // ins : del ratios
                    final int alleleSize = alt.length() - eval.getReference().length();
                    if ( alleleSize == 0 ) throw new ReviewedGATKException("Allele size not expected to be zero for indel: alt = " + alt + " ref = " + eval.getReference());
                    if ( alleleSize > 0 ) n_insertions++;
                    if ( alleleSize < 0 ) n_deletions++;

                    // requires snpEFF annotations
                    if ( eval.getAttributeAsString("SNPEFF_GENE_BIOTYPE", "missing").equals("protein_coding") ) {
                        final String effect = eval.getAttributeAsString("SNPEFF_EFFECT", "missing");
                        if ( effect.equals("missing") ) 
                            throw new ReviewedGATKException("Saw SNPEFF_GENE_BIOTYPE but unexpected no SNPEFF_EFFECT at " + eval);
                        if ( effect.equals("FRAME_SHIFT") )
                            n_coding_indels_frameshifting++;
                        else if ( effect.startsWith("CODON") )
                            n_coding_indels_in_frame++;
                        else
                            ; // lots of protein coding effects that shouldn't be counted, such as INTRON
                    }

                    if ( alleleSize > LARGE_INDEL_SIZE_THRESHOLD )
                        n_large_insertions++;
                    else if ( alleleSize < -LARGE_INDEL_SIZE_THRESHOLD )
                        n_large_deletions++;
                    
                    // update the baby histogram
                    final int[] countByLength = alleleSize < 0 ? deletionCountByLength : insertionCountByLength;
                    final int absSize = Math.abs(alleleSize);
                    if ( absSize < countByLength.length ) countByLength[absSize]++;

                }

                break;
            default:
                // TODO - MIXED, SYMBOLIC, and MNP records are skipped over
                //throw new UserException.BadInput("Unexpected variant context type: " + eval);
                break;
        }

        return;
    }

    public void finalizeEvaluation() {
        percent_of_sites_with_more_than_2_alleles = Utils.formattedPercent(n_multiallelic_indel_sites, nIndelSites);
        SNP_to_indel_ratio = Utils.formattedRatio(n_SNPs, n_indels);
        SNP_to_indel_ratio_for_singletons = Utils.formattedRatio(n_singleton_SNPs, n_singleton_indels);

        gold_standard_matching_rate = Utils.formattedPercent(n_indels_matching_gold_standard, n_indels);
        indel_novelty_rate = Utils.formattedNoveltyRate(n_indels - n_novel_indels, n_indels);
        frameshift_rate_for_coding_indels = Utils.formattedPercent(n_coding_indels_frameshifting, n_coding_indels_in_frame + n_coding_indels_frameshifting);

        nDeletionLength_1_bp = deletionCountByLength[1];
        nDeletionLength_2_bp = deletionCountByLength[2];
        nDeletionLength_3_bp = deletionCountByLength[3];
        ratio_of_1_and_2_to_3_bp_deletions = Utils.formattedRatio(deletionCountByLength[1] + deletionCountByLength[2], deletionCountByLength[3]);

        nInsertionLength_1_bp = insertionCountByLength[1];
        nInsertionLength_2_bp = insertionCountByLength[2];
        nInsertionLength_3_bp = insertionCountByLength[3];
        ratio_of_1_and_2_to_3_bp_insertions = Utils.formattedRatio(insertionCountByLength[1] + insertionCountByLength[2], insertionCountByLength[3]);

        SNP_het_to_hom_ratio = Utils.formattedRatio(nSNPHets, nSNPHoms);
        indel_het_to_hom_ratio = Utils.formattedRatio(nIndelHets, nIndelHoms);

        insertion_to_deletion_ratio = Utils.formattedRatio(n_insertions, n_deletions);
        insertion_to_deletion_ratio_for_large_indels = Utils.formattedRatio(n_large_insertions, n_large_deletions);

    }
}
