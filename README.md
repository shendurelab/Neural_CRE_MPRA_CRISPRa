# Diversified, miniaturized and ancestral parts for mammalian genome engineering and molecular recording

Troy A. McDiarmid1,2,†, Nicholas Page3,4,5,†, Florence M. Chardon1,2,†, Riza Daza1,2, George T. Chen6,17, Michael Kosicki8 , Lucas M. James, Hannah C. Nourie, Dianne Cintron-Laboy, Arthur S. Lee9,10,11,12, Paula Vij, Diego Calderon1, Jean-Benoît Lalanne1, Beth K. Martin1,2, Kyle Fink, Michael E. Talkowski10,11,12, Alysson Muotri24, Ben Philpot, Len A. Pennacchio8,15,16, Daniel H. Geschwind6,7,17,18,19, Stephan Sanders3,5,20, Nadav Ahituv4,5*, and Jay Shendure1,2,21,22,23*

Affiliations:
1. Department of Genome Sciences, University of Washington, Seattle, WA, USA
2. Seattle Hub for Synthetic Biology, Seattle, WA, USA
3. Department of Psychiatry and Behavioral Sciences, University of California, San Francisco, San Francisco, CA, USA
4. Department of Bioengineering and Therapeutic Sciences, University of California, San Francisco, San Francisco, CA, USA
5. Institute for Human Genetics, University of California, San Francisco, San Francisco, CA, USA
6. Program in Neurogenetics, Department of Neurology, David Geffen School of Medicine, UCLA, Los Angeles, CA, USA
7. Department of Psychiatry and Behavioral Sciences; David Geffen School of Medicine, UCLA, Los Angeles, CA, USA
8. Environmental Genomics & System Biology Division, Lawrence Berkeley National Laboratory, Berkeley, CA, USA
9. Department of Neurology, Boston Children’s Hospital and Harvard Medical School, Boston, MA, USA
10. Kirby Neurobiology Center, Boston Children’s Hospital, Boston, MA, USA
11. Manton Center for Orphan Disease Research, Boston Children’s Hospital, Boston, MA, USA
12. Program in Medical and Population Genetics, Broad Institute of MIT and Harvard, Cambridge, MA, USA
13. Center for Genomic Medicine, Massachusetts General Hospital, Boston, MA, USA
14. Department of Neurology, Massachusetts General Hospital and Harvard Medical School, Boston, MA, USA
15. U.S. Department of Energy Joint Genome Institute, One Cyclotron Road, Berkeley, CA, USA
16. Comparative Biochemistry Program, University of California, Berkeley, CA, USA
17. Center for Autism Research and Treatment, Semel Institute, UCLA, Los Angeles, CA, USA
18. Department of Human Genetics, David Geffen School of Medicine, UCLA, Los Angeles, CA, USA
19. Institute of Precision Health, David Geffen School of Medicine, UCLA, Los Angeles, CA, USA
20. Institute of Developmental and Regenerative Medicine, Department of Paediatrics, University of Oxford, Oxford, UK
21. Brotman Baty Institute for Precision Medicine, Seattle, WA, USA
22. Howard Hughes Medical Institute, Seattle, WA, USA
23. Allen Discovery Center for Cell Lineage Tracing, Seattle, WA, USA
24. Departments of Pediatrics & Cellular Molecular Medicine; School of Medicine, UCSD, La Jolla, CA 92093, USA

†. These authors contributed equally to this work
          
* Correspondence to N.A. (Nadav.Ahituv@ucsf.edu) or J.S. (shendure@uw.edu)

Abstract
CRISPR-based gene activation (CRISPRa) has emerged as a promising therapeutic approach for neurodevelopmental disorders (NDD) caused by haploinsufficiency. However, scaling this cis-regulatory therapy (CRT) paradigm requires pinpointing which candidate cis-regulatory elements (cCREs) are active in human neurons, and which can be targeted with CRISPRa to yield specific and therapeutic levels of target gene upregulation. Here, we combined Massively Parallel Reporter Assays (MPRAs) and multiplex single-cell CRISPRa screening to discover functional human neural enhancers whose CRISPRa targeting yields specific upregulation of NDD risk genes. First, we tested 5,425 cCREs with MPRA, identifying 2,422 enhancer sequences that are active in human neurons. Selected cCREs also displayed specific, autonomous in vivo activity in the developing mouse central nervous system. Next, we applied multiplex single-cell CRISPRa screening with 15,643 gRNAs to test all MPRA-prioritized cCREs and 761 promoters of NDD genes in their endogenous genomic contexts. We identified hundreds of promoter- and enhancer-targeting CRISPRa gRNAs that upregulated 200 of the 337 NDD genes in human neurons, including 91 novel enhancer-gene pairs. Finally, we confirmed that several of the CRISPRa gRNAs identified here demonstrated selective and therapeutically relevant levels of upregulation of SCN2A, CHD8, CTNND2, and TCF4 when delivered virally to patient cell lines, human cerebral organoids, and a Tcf4 humanized mouse model. Our results provide a comprehensive set of active and target-linked human neural enhancers for NDD genes that can be used to develop CRTs, and contribute to our understanding of the functional circuitry of neural gene regulation. Finally, we anticipate that our strategy for discovering gRNA candidates for CRT will generalize to other cell types and classes of haploinsufficient disorders.





# Repository Contents


This repository contains design, analysis, and visualization code for "Diversified, miniaturized and ancestral parts for mammalian genome engineering and molecular recording". 

Raw sequencing data, custom sequencing amplicons, and processed data files generated in this study have been deposited on the IGVF portal and are freely available with accession numbers XX (MPRA) and IGVFDS2290SSEF (single cell CRISPRa screen). 

The visualization scripts labeled "Fig_1_Viz.Rmd", "Fig_2_Viz.Rmd" etc. can then be used to recreate all visualizations in the manuscript. All processed data and associated metadata required to recreate visualizations are provided in the final figure data set folders labeled "Fig1_Final_Figure_Datasets", "Fig2_Final_Figure_Datasets" etc. More formal results tables including relevant test sequences for all experiments are also availale in Tables S1-S7 in the manuscript. 
