/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.haplotypecaller.graphs;

import org.broadinstitute.sting.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

public class SeqGraphUnitTest extends BaseTest {
    private final static boolean DEBUG = false;

    private class MergeNodesWithNoVariationTestProvider extends TestDataProvider {
        public byte[] sequence;
        public int KMER_LENGTH;

        public MergeNodesWithNoVariationTestProvider(String seq, int kmer) {
            super(MergeNodesWithNoVariationTestProvider.class, String.format("Merge nodes with no variation test. kmer = %d, seq = %s", kmer, seq));
            sequence = seq.getBytes();
            KMER_LENGTH = kmer;
        }

        public SeqGraph calcGraph() {
            final DeBruijnGraph deBruijnGraph = new DeBruijnGraph();
            final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
            for (int i = 0; i < kmersInSequence - 1; i++) {
                // get the kmers
                final byte[] kmer1 = new byte[KMER_LENGTH];
                System.arraycopy(sequence, i, kmer1, 0, KMER_LENGTH);
                final byte[] kmer2 = new byte[KMER_LENGTH];
                System.arraycopy(sequence, i+1, kmer2, 0, KMER_LENGTH);

                deBruijnGraph.addKmersToGraph(kmer1, kmer2, false, 1);
            }
            final SeqGraph seqGraph = deBruijnGraph.convertToSequenceGraph();
            seqGraph.simplifyGraph();
            return seqGraph;
        }
    }

    @DataProvider(name = "MergeNodesWithNoVariationTestProvider")
    public Object[][] makeMergeNodesWithNoVariationTests() {
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 3);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 4);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 5);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 6);
        new MergeNodesWithNoVariationTestProvider("GGTTAACC", 7);
        new MergeNodesWithNoVariationTestProvider("GGTTAACCATGCAGACGGGAGGCTGAGCGAGAGTTTT", 6);
        new MergeNodesWithNoVariationTestProvider("AATACCATTGGAGTTTTTTTCCAGGTTAAGATGGTGCATTGAATCCACCCATCTACTTTTGCTCCTCCCAAAACTCACTAAAACTATTATAAAGGGATTTTGTTTAAAGACACAAACTCATGAGGACAGAGAGAACAGAGTAGACAATAGTGGGGGAAAAATAAGTTGGAAGATAGAAAACAGATGGGTGAGTGGTAATCGACTCAGCAGCCCCAAGAAAGCTGAAACCCAGGGAAAGTTAAGAGTAGCCCTATTTTCATGGCAAAATCCAAGGGGGGGTGGGGAAAGAAAGAAAAACAGAAAAAAAAATGGGAATTGGCAGTCCTAGATATCTCTGGTACTGGGCAAGCCAAAGAATCAGGATAACTGGGTGAAAGGTGATTGGGAAGCAGTTAAAATCTTAGTTCCCCTCTTCCACTCTCCGAGCAGCAGGTTTCTCTCTCTCATCAGGCAGAGGGCTGGAGAT", 66);
        new MergeNodesWithNoVariationTestProvider("AATACCATTGGAGTTTTTTTCCAGGTTAAGATGGTGCATTGAATCCACCCATCTACTTTTGCTCCTCCCAAAACTCACTAAAACTATTATAAAGGGATTTTGTTTAAAGACACAAACTCATGAGGACAGAGAGAACAGAGTAGACAATAGTGGGGGAAAAATAAGTTGGAAGATAGAAAACAGATGGGTGAGTGGTAATCGACTCAGCAGCCCCAAGAAAGCTGAAACCCAGGGAAAGTTAAGAGTAGCCCTATTTTCATGGCAAAATCCAAGGGGGGGTGGGGAAAGAAAGAAAAACAGAAAAAAAAATGGGAATTGGCAGTCCTAGATATCTCTGGTACTGGGCAAGCCAAAGAATCAGGATAACTGGGTGAAAGGTGATTGGGAAGCAGTTAAAATCTTAGTTCCCCTCTTCCACTCTCCGAGCAGCAGGTTTCTCTCTCTCATCAGGCAGAGGGCTGGAGAT", 76);

        return MergeNodesWithNoVariationTestProvider.getTests(MergeNodesWithNoVariationTestProvider.class);
    }

    @Test(dataProvider = "MergeNodesWithNoVariationTestProvider", enabled = !DEBUG)
    public void testMergeNodesWithNoVariation(MergeNodesWithNoVariationTestProvider cfg) {
        logger.warn(String.format("Test: %s", cfg.toString()));

        final SeqGraph actual = cfg.calcGraph();
        Assert.assertEquals(actual.vertexSet().size(), 1);
        final SeqVertex actualV = actual.vertexSet().iterator().next();
        Assert.assertEquals(actualV.getSequence(), cfg.sequence);
    }

    @DataProvider(name = "IsDiamondData")
    public Object[][] makeIsDiamondData() throws Exception {
        List<Object[]> tests = new ArrayList<Object[]>();

        SeqGraph graph;
        SeqVertex pre1, pre2, top, middle1, middle2, middle3, bottom, tail1, tail2;

        graph = new SeqGraph();

        pre1 = new SeqVertex("ACT");
        pre2 = new SeqVertex("AGT");
        top = new SeqVertex("A");
        middle1 = new SeqVertex("CT");
        middle2 = new SeqVertex("CG");
        middle3 = new SeqVertex("CA");
        bottom = new SeqVertex("AA");
        tail1 = new SeqVertex("GC");
        tail2 = new SeqVertex("GC");

        graph.addVertices(pre1, pre2, top, middle1, middle2, middle3, bottom, tail1, tail2);
        graph.addEdges(pre1, top, middle1, bottom, tail1);
        graph.addEdges(pre2, top, middle2, bottom, tail1);
        graph.addEdges(top, middle3, bottom);
        graph.addEdges(bottom, tail2);

        for ( final SeqVertex no : Arrays.asList(pre1, pre2, middle1, middle2, middle3, bottom, tail1, tail2)) {
            tests.add(new Object[]{graph, no, false});
        }
        tests.add(new Object[]{graph, top, true});

        final SeqGraph danglingMiddleGraph = (SeqGraph)graph.clone();
        final SeqVertex danglingMiddle = new SeqVertex("A");
        danglingMiddleGraph.addVertex(danglingMiddle);
        danglingMiddleGraph.addEdge(top, danglingMiddle);
        tests.add(new Object[]{danglingMiddleGraph, top, false});

        final SeqGraph strangerToBottom = (SeqGraph)graph.clone();
        final SeqVertex notAttachedToTop = new SeqVertex("A");
        strangerToBottom.addVertex(notAttachedToTop);
        strangerToBottom.addEdge(notAttachedToTop, bottom);
        tests.add(new Object[]{strangerToBottom, top, false});

        final SeqGraph strangerToMiddle = (SeqGraph)graph.clone();
        final SeqVertex attachedToMiddle = new SeqVertex("A");
        strangerToMiddle.addVertex(attachedToMiddle);
        strangerToMiddle.addEdge(attachedToMiddle, middle1);
        tests.add(new Object[]{strangerToMiddle, top, false});

        // middle1 has outgoing edge to non-bottom
        final SeqGraph middleExtraOut = (SeqGraph)graph.clone();
        final SeqVertex fromMiddle = new SeqVertex("A");
        middleExtraOut.addVertex(fromMiddle);
        middleExtraOut.addEdge(middle1, fromMiddle);
        tests.add(new Object[]{middleExtraOut, top, false});

        // top connects to bottom directly as well
        {
            final SeqGraph topConnectsToBottomToo = new SeqGraph();
            final SeqVertex top2 = new SeqVertex("A");
            final SeqVertex middle4 = new SeqVertex("C");
            final SeqVertex bottom2 = new SeqVertex("G");
            topConnectsToBottomToo.addVertices(top2, middle4, bottom2);
            topConnectsToBottomToo.addEdges(top2, middle4, bottom2);
            topConnectsToBottomToo.addEdges(top2, bottom2);
            tests.add(new Object[]{topConnectsToBottomToo, top2, false});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "IsDiamondData", enabled = !DEBUG)
    public void testIsDiamond(final SeqGraph graph, final SeqVertex v, final boolean isRootOfDiamond) {
        final SeqGraph.MergeDiamonds merger = graph.new MergeDiamonds();
        merger.setDontModifyGraphEvenIfPossible();
        Assert.assertEquals(merger.tryToTransform(v), isRootOfDiamond);
    }

    @DataProvider(name = "MergingData")
    public Object[][] makeMergingData() throws Exception {
        List<Object[]> tests = new ArrayList<Object[]>();

        final SeqGraph graph = new SeqGraph();

        SeqVertex pre1 = new SeqVertex("ACT");
        SeqVertex pre2 = new SeqVertex("AGT");
        SeqVertex top = new SeqVertex("A");
        SeqVertex middle1 = new SeqVertex("GC");
        SeqVertex middle2 = new SeqVertex("TC");
        SeqVertex middle3 = new SeqVertex("AC");
        SeqVertex middle4 = new SeqVertex("GCAC");
        SeqVertex bottom = new SeqVertex("AA");
        SeqVertex tail1 = new SeqVertex("GC");
        SeqVertex tail2 = new SeqVertex("GC");

        // just a single vertex
        graph.addVertices(pre1);
        tests.add(new Object[]{graph.clone(), graph.clone()});

        // pre1 -> top = pre1 + top
        {
            graph.addVertices(top);
            graph.addEdges(pre1, top);
            final SeqVertex pre1_top = new SeqVertex(pre1.getSequenceString() + top.getSequenceString());
            final SeqGraph expected = new SeqGraph();
            expected.addVertex(pre1_top);
            tests.add(new Object[]{graph.clone(), expected.clone()});
        }

        // pre1 -> top -> middle1 = pre1 + top + middle1
        {
            graph.addVertices(middle1);
            graph.addEdges(top, middle1);
            final SeqGraph expected = new SeqGraph();
            final SeqVertex pre1_top_middle1 = new SeqVertex(pre1.getSequenceString() + top.getSequenceString() + middle1.getSequenceString());
            expected.addVertex(pre1_top_middle1);
            tests.add(new Object[]{graph.clone(), expected});
        }

        // pre1 -> top -> middle1 & top -> middle2 = pre1 + top -> middle1 & -> middle2
        {
            graph.addVertices(middle2);
            graph.addEdges(top, middle2);
            final SeqGraph expected = new SeqGraph();
            final SeqVertex pre1_top = new SeqVertex(pre1.getSequenceString() + top.getSequenceString());
            expected.addVertices(pre1_top, middle1, middle2);
            expected.addEdges(pre1_top, middle1);
            expected.addEdges(pre1_top, middle2);
            tests.add(new Object[]{graph.clone(), expected});
        }

        // An actual diamond event to merge!
        {
            graph.addVertices(bottom);
            graph.addEdges(middle1, bottom);
            graph.addEdges(middle2, bottom);
            final SeqGraph expected = new SeqGraph();
            final SeqVertex pre1_top = new SeqVertex(pre1.getSequenceString() + top.getSequenceString());
            final SeqVertex newMiddle1 = new SeqVertex("G");
            final SeqVertex newMiddle2 = new SeqVertex("T");
            final SeqVertex newBottom = new SeqVertex("C" + bottom.getSequenceString());
            expected.addVertices(pre1_top, newMiddle1, newMiddle2, newBottom);
            expected.addEdges(pre1_top, newMiddle1, newBottom);
            expected.addEdges(pre1_top, newMiddle2, newBottom);
            tests.add(new Object[]{graph.clone(), expected.clone()});

            graph.addVertices(middle3);
            graph.addEdges(top, middle3, bottom);
            final SeqVertex newMiddle3 = new SeqVertex("A");
            expected.addVertices(newMiddle3);
            expected.addEdges(pre1_top, newMiddle3, newBottom);
            tests.add(new Object[]{graph.clone(), expected.clone()});

            graph.addVertices(middle4);
            graph.addEdges(top, middle4, bottom);
            final SeqVertex newMiddle4 = new SeqVertex("GCA");
            expected.addVertices(newMiddle4);
            expected.addEdges(pre1_top, newMiddle4, newBottom);
            tests.add(new Object[]{graph.clone(), expected.clone()});
        }

        { // all the nodes -> lots of merging and motion of nodes
            final SeqGraph all = new SeqGraph();
            all.addVertices(pre1, pre2, top, middle1, middle2, bottom, tail1, tail2);
            all.addEdges(pre1, top, middle1, bottom, tail1);
            all.addEdges(pre2, top, middle2, bottom, tail2);

            final SeqGraph expected = new SeqGraph();
            final SeqVertex newMiddle1 = new SeqVertex("G");
            final SeqVertex newMiddle2 = new SeqVertex("T");
            final SeqVertex newBottom = new SeqVertex("C" + bottom.getSequenceString());
            final SeqVertex newTop = new SeqVertex("A");
            final SeqVertex newTopDown1 = new SeqVertex("G");
            final SeqVertex newTopDown2 = new SeqVertex("C");
            final SeqVertex newTopBottomMerged = new SeqVertex("TA");
            expected.addVertices(newTop, newTopDown1, newTopDown2, newTopBottomMerged, newMiddle1, newMiddle2, newBottom, tail1, tail2);
            expected.addEdges(newTop, newTopDown1, newTopBottomMerged, newMiddle1, newBottom, tail1);
            expected.addEdges(newTop, newTopDown2, newTopBottomMerged, newMiddle2, newBottom, tail2);
            tests.add(new Object[]{all.clone(), expected.clone()});
        }

        // test the case where we delete a middle node away because the common sequence is all of its sequence
        {
            final SeqGraph graph2 = new SeqGraph();
            final SeqVertex mytop = new SeqVertex("A");
            final SeqVertex mid1 = new SeqVertex("AC");
            final SeqVertex mid2 = new SeqVertex("C");
            final SeqVertex bot = new SeqVertex("G");
            graph2.addVertices(mytop, mid1, mid2, bot);
            graph2.addEdges(mytop, mid1, bot);
            graph2.addEdges(mytop, mid2, bot);

            final SeqGraph expected = new SeqGraph();
            final SeqVertex newMid1 = new SeqVertex("A");
            final SeqVertex newBottom = new SeqVertex("CG");
            expected.addVertices(mytop, newMid1, newBottom);
            expected.addEdges(mytop, newMid1, newBottom);
            expected.addEdges(mytop, newBottom);
            tests.add(new Object[]{graph2, expected});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "MergingData", enabled = !DEBUG)
    public void testMerging(final SeqGraph graph, final SeqGraph expected) {
        final SeqGraph merged = (SeqGraph)graph.clone();
        merged.simplifyGraph(1);
        try {
            Assert.assertTrue(SeqGraph.graphEquals(merged, expected));
        } catch (AssertionError e) {
//            if ( ! SeqGraph.graphEquals(merged, expected) ) {
//                graph.printGraph(new File("graph.dot"), 0);
//                merged.printGraph(new File("merged.dot"), 0);
//                expected.printGraph(new File("expected.dot"), 0);
//            }
            throw e;
        }
    }

    // A -> ACT -> C [non-ref]
    // A -> ACT -> C [non-ref]
    // A -> ACT -> C [ref]
    //
    // Should become A -> ACT -> C [ref and non-ref edges]
    //
    @Test(enabled = !DEBUG)
    public void testBubbleSameBasesWithRef() {
        final SeqGraph graph = new SeqGraph();
        final SeqVertex top = new SeqVertex("A");
        final SeqVertex mid1 = new SeqVertex("ACT");
        final SeqVertex mid2 = new SeqVertex("ACT");
        final SeqVertex bot = new SeqVertex("C");
        graph.addVertices(top, mid1, mid2, bot);
        graph.addEdges(top, mid2, bot);
        graph.addEdge(top, mid1, new BaseEdge(true, 1));
        graph.addEdge(mid1, bot, new BaseEdge(true, 1));

        final SeqGraph expected = new SeqGraph();
        expected.addVertex(new SeqVertex("AACTC"));
        final SeqGraph actual = ((SeqGraph)graph.clone());
        actual.simplifyGraph();
        Assert.assertTrue(BaseGraph.graphEquals(actual, expected), "Wrong merging result after complete merging");
    }

    @DataProvider(name = "LinearZipData")
    public Object[][] makeLinearZipData() throws Exception {
        List<Object[]> tests = new ArrayList<Object[]>();

        SeqGraph graph = new SeqGraph();
        SeqGraph expected = new SeqGraph();

        // empty graph => empty graph
        tests.add(new Object[]{graph.clone(), expected.clone()});

        SeqVertex a1 = new SeqVertex("A");
        SeqVertex c1 = new SeqVertex("C");
        SeqVertex ac1 = new SeqVertex("AC");

        // just a single vertex
        graph.addVertices(a1, c1);
        expected.addVertices(a1, c1);

        tests.add(new Object[]{graph.clone(), expected.clone()});

        graph.addEdges(a1, c1);
        expected = new SeqGraph();
        expected.addVertices(ac1);
        tests.add(new Object[]{graph.clone(), expected.clone()});

        // three long chain merged corrected
        SeqVertex g1 = new SeqVertex("G");
        graph.addVertices(g1);
        graph.addEdges(c1, g1);
        expected = new SeqGraph();
        expected.addVertex(new SeqVertex("ACG"));
        tests.add(new Object[]{graph.clone(), expected.clone()});

        // adding something that isn't connected isn't a problem
        SeqVertex t1 = new SeqVertex("T");
        graph.addVertices(t1);
        expected = new SeqGraph();
        expected.addVertices(new SeqVertex("ACG"), new SeqVertex("T"));
        tests.add(new Object[]{graph.clone(), expected.clone()});

        // splitting chain with branch produces the correct zipped subgraphs
        final SeqVertex a2 = new SeqVertex("A");
        final SeqVertex c2 = new SeqVertex("C");
        graph = new SeqGraph();
        graph.addVertices(a1, c1, g1, t1, a2, c2);
        graph.addEdges(a1, c1, g1, t1, a2);
        graph.addEdges(g1, c2);
        expected = new SeqGraph();
        SeqVertex acg = new SeqVertex("ACG");
        SeqVertex ta = new SeqVertex("TA");
        expected.addVertices(acg, ta, c2);
        expected.addEdges(acg, ta);
        expected.addEdges(acg, c2);
        tests.add(new Object[]{graph.clone(), expected.clone()});

        // Can merge chains with loops in them
        {
            graph = new SeqGraph();
            graph.addVertices(a1, c1, g1);
            graph.addEdges(a1, c1, g1);
            graph.addEdges(a1, a1);
            expected = new SeqGraph();

            SeqVertex ac = new SeqVertex("AC");
            SeqVertex cg = new SeqVertex("CG");

            expected.addVertices(a1, cg);
            expected.addEdges(a1, cg);
            expected.addEdges(a1, a1);
            tests.add(new Object[]{graph.clone(), expected.clone()});

            graph.removeEdge(a1, a1);
            graph.addEdges(c1, c1);
            tests.add(new Object[]{graph.clone(), graph.clone()});

            graph.removeEdge(c1, c1);
            graph.addEdges(g1, g1);
            expected = new SeqGraph();
            expected.addVertices(ac, g1);
            expected.addEdges(ac, g1, g1);
            tests.add(new Object[]{graph.clone(), expected.clone()});
        }

        // check building n element long chains
        {
            final List<String> bases = Arrays.asList("A", "C", "G", "T", "TT", "GG", "CC", "AA");
            for ( final int len : Arrays.asList(1, 2, 10, 100, 1000)) {
                graph = new SeqGraph();
                expected = new SeqGraph();
                SeqVertex last = null;
                String expectedBases = "";
                for ( int i = 0; i < len; i++ ) {
                    final String seq = bases.get(i % bases.size());
                    expectedBases += seq;
                    SeqVertex a = new SeqVertex(seq);
                    graph.addVertex(a);
                    if ( last != null ) graph.addEdge(last, a);
                    last = a;
                }
                expected.addVertex(new SeqVertex(expectedBases));
                tests.add(new Object[]{graph.clone(), expected.clone()});
            }
        }

        // check that edge connections are properly maintained
        {
            int edgeWeight = 1;
            for ( final int nIncoming : Arrays.asList(0, 2, 5, 10) ) {
                for ( final int nOutgoing : Arrays.asList(0, 2, 5, 10) ) {
                    graph = new SeqGraph();
                    expected = new SeqGraph();

                    graph.addVertices(a1, c1, g1);
                    graph.addEdges(a1, c1, g1);
                    expected.addVertex(acg);

                    for ( final SeqVertex v : makeVertices(nIncoming) ) {
                        final BaseEdge e = new BaseEdge(false, edgeWeight++);
                        graph.addVertices(v);
                        graph.addEdge(v, a1, e);
                        expected.addVertex(v);
                        expected.addEdge(v, acg, e);
                    }

                    for ( final SeqVertex v : makeVertices(nOutgoing) ) {
                        final BaseEdge e = new BaseEdge(false, edgeWeight++);
                        graph.addVertices(v);
                        graph.addEdge(g1, v, e);
                        expected.addVertex(v);
                        expected.addEdge(acg, v, e);
                    }

                    tests.add(new Object[]{graph, expected});
                }
            }
        }

        return tests.toArray(new Object[][]{});
    }

    private List<SeqVertex> makeVertices(final int n) {
        final List<SeqVertex> vs = new LinkedList<SeqVertex>();
        final List<String> bases = Arrays.asList("A", "C", "G", "T", "TT", "GG", "CC", "AA");

        for ( int i = 0; i < n; i++ )
            vs.add(new SeqVertex(bases.get(i % bases.size())));
        return vs;
    }

    @Test(dataProvider = "LinearZipData", enabled = true)
    public void testLinearZip(final SeqGraph graph, final SeqGraph expected) {
        final SeqGraph merged = (SeqGraph)graph.clone();
        merged.zipLinearChains();
        try {
            Assert.assertTrue(SeqGraph.graphEquals(merged, expected));
        } catch (AssertionError e) {
            if ( ! SeqGraph.graphEquals(merged, expected) ) {
                graph.printGraph(new File("graph.dot"), 0);
                merged.printGraph(new File("merged.dot"), 0);
                expected.printGraph(new File("expected.dot"), 0);
            }
            throw e;
        }
    }
}
