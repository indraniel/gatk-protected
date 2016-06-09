package org.broadinstitute.gatk.tools.walkers.annotator;

/**
 * Created by idas on 6/9/16.
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SimpleABHet extends InfoFieldAnnotation {
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        //if ( stratifiedContexts.size() == 0 )
        //    return null;

        if ( !vc.isBiallelic() )
            return null;
        final GenotypesContext genotypes = vc.getGenotypes();
        if ( !vc.hasGenotypes() )
            return null;

        Map<String, Object> gt_counts = collectGenotypeMetrics(vc);

        int sample_count     = (int) gt_counts.get("gt_01_sample_cnt");
        int ref_allele_count = (int) gt_counts.get("gt_01_ref_allele");
        int alt_allele_count = (int) gt_counts.get("gt_01_alt_allele");
        int total_depth      = (int) gt_counts.get("gt_01_total_depth");
        float allele_balance = (float) gt_counts.get("gt_01_allele_bal");

        Map<String, Object> map = new HashMap<>();
        if ( allele_balance > 0.0 ) {
            map.put(GATKVCFConstants.AB_SNP_HET_RATIO, allele_balance);
        }

        map.put(GATKVCFConstants.AB_SNP_HET_REF_ALLELE_COUNT, ref_allele_count);
        map.put(GATKVCFConstants.AB_SNP_HET_ALT_ALLELE_COUNT, alt_allele_count);
        map.put(GATKVCFConstants.AB_SNP_HET_SAMPLE_COUNT, sample_count);
        map.put(GATKVCFConstants.AB_SNP_HET_TOTAL_ALLELE_COUNT, total_depth);

        return map;
    }

    public Map<String, Object> collectGenotypeMetrics(VariantContext vc) {
        int gt_01_sample_cnt = 0;
        int gt_01_ref_allele = 0;
        int gt_01_alt_allele = 0;
        int gt_01_total_depth = 0;
        float gt_01_allele_bal = (float) -1.0;

        if ( !vc.hasGenotypes() ) {
            errorVariant(vc, "Unable to find genotypes");
        }

        GenotypesContext genotypes = vc.getGenotypes();

        for ( Genotype genotype : genotypes ) {
            String gt = genotype.getGenotypeString();

            if ( genotype.isHomRef() ) {
                continue;
            }
            else if ( genotype.isHomVar() ) {
                continue;
            }
            else if ( genotype.isHet() ) {
                gt_01_sample_cnt++;
                int[] AD = genotype.getAD();
                gt_01_ref_allele += AD[0];
                gt_01_alt_allele += AD[1];
            }
            else if ( genotype.isNoCall() ) {
                continue;
            }
            else if ( genotype.isMixed() ) {
                continue;
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

        gt_01_total_depth = gt_01_ref_allele + gt_01_alt_allele ;
        if ( (gt_01_total_depth) != 0 ) {
            gt_01_allele_bal = (float) gt_01_alt_allele / gt_01_total_depth;
        }

        Map<String, Object> hmap = new HashMap<>();

        hmap.put("gt_01_sample_cnt",   gt_01_sample_cnt);
        hmap.put("gt_01_ref_allele",   gt_01_ref_allele);
        hmap.put("gt_01_alt_allele",   gt_01_alt_allele);
        hmap.put("gt_01_total_depth", gt_01_total_depth);
        hmap.put("gt_01_allele_bal",   gt_01_allele_bal);

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

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(
                GATKVCFConstants.AB_SNP_HET_REF_ALLELE_COUNT,
                GATKVCFConstants.AB_SNP_HET_ALT_ALLELE_COUNT,
                GATKVCFConstants.AB_SNP_HET_RATIO,
                GATKVCFConstants.AB_SNP_HET_SAMPLE_COUNT,
                GATKVCFConstants.AB_SNP_HET_TOTAL_ALLELE_COUNT
        );
    }
}

