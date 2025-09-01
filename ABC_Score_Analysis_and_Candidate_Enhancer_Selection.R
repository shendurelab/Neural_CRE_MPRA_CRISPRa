########################################################################################################################
#
#   02_ABC_Score_Analysis_and_Candidate_Enhancer_Selection.R
#
#   Analyze ABC Score results from WTC11 NGN2-iPSC derived excitatory neurons and Human GW18 prefrontal cortex (PFC), 
#   create metadata results table, and integrate with previously published data from Cicero, Domcke et al. 2020, 
#   (Science), Trevino et al. 2021 (Cell), and Ziffra et al. 2021 (Nature).
#
#   Nicholas Page, August 2021
#   Ahituv Lab, Dept. of Bioengineering and Therapeutic Sciences, University of California, San Francisco
#   Sanders Lab, Dept. of Psychiatry, University of California, San Francisco
#   Neuroscience Graduate Program, University of California, San Francisco
#
########################################################################################################################

rm(list = ls())
options(java.parameters = "-Xmx8000m")
options(stringsAsFactors = FALSE)
options(scipen = 5)

library(ggpointdensity)
library(ComplexHeatmap)
library(tidyverse)
library(UpSetR)
library(sjmisc)
library(readxl)
library(GGally)
library(Cairo)
library(xlsx)
library(ape)

setwd('~/Dropbox/Neurohub_CRISPRa_CRT/')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Create supplemental tables for GW18 PFC and WTC11 NGN2-iPSC ABC score analyses
#   (2) Match promoter contacts and generate supplemental tables for midfetal cerebrum ExN cicero data
#   (3) Get Trevino and Ziffra datasets and create bedpe files
#   (4) Create hg19 and hg38 bigInteract files for each dataset
#   (5) Get candidate NDD enhancers for each dataset and merge them into a single bed file and make upset plot
#   (6) Get genome wide background set of candidate enhancers and make upset plot
#   (7) Plot pairwise enhancer score correlations for all 5 prioritization datasets
#   (8) Get merged bedpe tracks for all 5 prioritization datasets
#
########################################################################################################################

### (1) Create supplemental tables for GW18 PFC and WTC11 NGN2-iPSC ABC score analyses #################################

# Get gencode annotations for gene names
gencode_hg38 <- read.gff(file = "~/genomes/hg38/gencode.v38.basic.annotation.gff3", GFF3 = TRUE)
gencode_hg38_genes <- gencode_hg38[which(gencode_hg38$type == "gene"),]
annotations <- data.frame("gene_name" = substring(sapply(strsplit(gencode_hg38_genes$attributes, ";"), `[`, 4), 11),
                          "gene_id" = substring(sapply(strsplit(gencode_hg38_genes$attributes, ";"), `[`, 2), 9, 23),
                          "hgnc_id" = substring(sapply(strsplit(gencode_hg38_genes$attributes, ";"), `[`, 6), 9),
                          "chrom" = gencode_hg38_genes$seqid,
                          "start_position_hg38" = gencode_hg38_genes$start,
                          "end_position_hg38" = gencode_hg38_genes$end,
                          "strand" = gencode_hg38_genes$strand)


# Get promoter annotation file used in ABC score analysis
promoter_annotations <- read.table(file = "./bed/annotations/hg19/gencode.v38lift37.2000bp_promoter.matched_gene_bounds.bed", header = FALSE, sep = "\t",
                                   col.names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "geneIDs"))

# Get GW18 PFC ABC score output
human_18gw_pfc_abc_score <- read.table(file = "./data/ABC_Score/18gw_pfc_ABC_Score.txt", header = TRUE, sep = "\t")

# Convert promoter IDs in ABC score output to gene names 
human_18gw_pfc_abc_score$matched_transcripts <- promoter_annotations$geneIDs[match(human_18gw_pfc_abc_score$TargetGene, promoter_annotations$name)]
human_18gw_pfc_abc_score$geneName <- sapply(strsplit(human_18gw_pfc_abc_score$matched_transcripts, ";"), `[`, 1)
human_18gw_pfc_abc_score$transcripts <- sapply(strsplit(human_18gw_pfc_abc_score$matched_transcripts, ";"), `[`, -1)
human_18gw_pfc_abc_score$transcripts <- vapply(human_18gw_pfc_abc_score$transcripts, function(.x) {
  y <- unlist(.x, use.names = TRUE) 
  paste0(names(y), y, collapse = ';')
}, FUN.VALUE = character(1))

# Plot the enhancer-TSS distance for GW18
fold_plot_1 <- ggplot(human_18gw_pfc_abc_score, aes(x = distance, y = ABC.Score)) +
  geom_point(size = 0.1, color='#27AAE1') + 
  # geom_density_2d(bins = 5, color = "white") + 
  scale_x_continuous(breaks = c(0 , 1000000, 2000000, 3000000, 4000000, 5000000),labels = c("", "1 Mb", "2 Mb", "3 Mb", "4 Mb", "5 Mb")) +
  ggtitle("GW18_PFC\n") +
  theme_classic() +
  xlab("\nTSS-Enhancer Distance") + 
  ylab("ABC Score\n") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(1, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(15, 5, 15, 15))
plot(fold_plot_1)


# Get WTC11 NGN2-iPSC ExN ABC score output
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score <- read.table(file = "./data/ABC_Score/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score.txt", header = TRUE, sep = "\t")

# Convert promoter IDs in ABC score output to gene names 
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$matched_transcripts <- promoter_annotations$geneIDs[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TargetGene, promoter_annotations$name)]
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName <- sapply(strsplit(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$matched_transcripts, ";"), `[`, 1)
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$transcripts <- sapply(strsplit(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$matched_transcripts, ";"), `[`, -1)
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$transcripts <- vapply(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$transcripts, function(.x) {
  y <- unlist(.x, use.names = TRUE) 
  paste0(names(y), y, collapse = ';')
}, FUN.VALUE = character(1))

# Plot the enhancer-TSS distance for NGN2-iNeurons
fold_plot_1 <- ggplot(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score, aes(x = distance, y = ABC.Score)) +
  geom_point(size = 0.1, color='#27AAE1') + 
  # geom_density_2d(bins = 5, color = "white") + 
  scale_x_continuous(breaks = c(0 , 1000000, 2000000, 3000000, 4000000, 5000000),labels = c("", "1 Mb", "2 Mb", "3 Mb", "4 Mb", "5 Mb")) +
  ggtitle("NGN2-iNeurons\n") +
  theme_classic() +
  xlab("\nTSS-Enhancer Distance") + 
  ylab("ABC Score\n") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(1, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(15, 5, 15, 15))
plot(fold_plot_1)


# Create hg19 bedpe file for human 18GW PFC dataset
human_18gw_pfc_abc_score_bedpe <- data.frame("chrom1" = human_18gw_pfc_abc_score$chr,
                                             "start1" = human_18gw_pfc_abc_score$start,
                                             "end1" = human_18gw_pfc_abc_score$end,
                                             "chrom2" = human_18gw_pfc_abc_score$chr,
                                             "start2" = human_18gw_pfc_abc_score$TargetGeneTSS,
                                             "end2" = human_18gw_pfc_abc_score$TargetGeneTSS + 1,
                                             "name" = paste0(human_18gw_pfc_abc_score$geneName, "|", human_18gw_pfc_abc_score$chr, ":", human_18gw_pfc_abc_score$start, "-", human_18gw_pfc_abc_score$end),
                                             "score" = human_18gw_pfc_abc_score$ABC.Score,
                                             "strand1" = ".",
                                             "strand2" = ".")

write.table(human_18gw_pfc_abc_score_bedpe, file = "./bed/bedpe/hg19/gw18_pfc_ABC_Score.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Create hg19 bedpe file for WTC11 NGN2-iPSC ExN dataset
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe <- data.frame("chrom1" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chr,
                                                        "start1" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$start,
                                                        "end1" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$end,
                                                        "chrom2" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chr,
                                                        "start2" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TargetGeneTSS,
                                                        "end2" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TargetGeneTSS + 1,
                                                        "name" = paste0(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, "|", WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chr, ":", WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$start, "-", WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$end),
                                                        "score" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$ABC.Score,
                                                        "strand1" = ".",
                                                        "strand2" = ".")

write.table(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe, file = "./bed/bedpe/hg19/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# LiftOver hg19 bedpe files for each dataset to hg38
# NOTE: the MUC12 entry causes a liftOver error for the WTC11 NGN2-iPSC ExN dataset so that entry is removed
# Need to fix python2 vs python3 error in liftOverBedpe.py!!!
system("python3 ./scripts/UCSC/liftOverBedpe/liftOverBedpe.py --lift ./scripts/UCSC/liftOver --chain ./scripts/UCSC/hg19ToHg38.over.chain --i ./bed/bedpe/hg19/gw18_pfc_ABC_Score.bedpe --o ./bed/bedpe/hg38/gw18_pfc_ABC_Score_liftOver_hg38.bedpe")

# remove MUC12 because of liftOver error
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe_no_MUC12 <- WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe[which(sapply(strsplit(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$name, "|", fixed = TRUE), `[`, 1) != "MUC12"),]
write.table(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe_no_MUC12, file = "./bed/bedpe/hg19/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_no_MUC12.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
system("python3 ./scripts/UCSC/liftOverBedpe/liftOverBedpe.py --lift ./scripts/UCSC/liftOver --chain ./scripts/UCSC/hg19ToHg38.over.chain --i ./bed/bedpe/hg19/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_no_MUC12.bedpe --o ./bed/bedpe/hg38/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_liftOver_hg38.bedpe")

# Remove temporary file with no MUC12
system("rm ./bed/bedpe/hg19/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_no_MUC12.bedpe")

# Load hg38 bedpe files for each dataset
human_18gw_pfc_abc_score_bedpe <- read.table(file = "./bed/bedpe/hg38/gw18_pfc_ABC_Score_liftOver_hg38.bedpe",
                                            header = FALSE, sep = "\t",
                                            col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))
WTC11_NGN2_iPSC_abc_score_bedpe <- read.table(file = "./bed/bedpe/hg38/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_liftOver_hg38.bedpe",
                                             header = FALSE, sep = "\t",
                                             col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))

# add a new name column to the hg19 human 18GW PFC dataset
human_18gw_pfc_abc_score$new_name <- paste0(human_18gw_pfc_abc_score$geneName, "|", human_18gw_pfc_abc_score$chr, ":", human_18gw_pfc_abc_score$start, "-", human_18gw_pfc_abc_score$end)

# add hg38 coordinates to human 18GW PFC dataset
head(human_18gw_pfc_abc_score)
human_18gw_pfc_abc_score$chromStart_hg38 <- human_18gw_pfc_abc_score_bedpe$start1[match(human_18gw_pfc_abc_score$new_name, human_18gw_pfc_abc_score_bedpe$name)]
human_18gw_pfc_abc_score$chromEnd_hg38 <- human_18gw_pfc_abc_score_bedpe$end1[match(human_18gw_pfc_abc_score$new_name, human_18gw_pfc_abc_score_bedpe$name)]
human_18gw_pfc_abc_score$TSSchromStart_hg38 <- human_18gw_pfc_abc_score_bedpe$start2[match(human_18gw_pfc_abc_score$new_name, human_18gw_pfc_abc_score_bedpe$name)]
human_18gw_pfc_abc_score$TSSchromEnd_hg38 <- human_18gw_pfc_abc_score_bedpe$end2[match(human_18gw_pfc_abc_score$new_name, human_18gw_pfc_abc_score_bedpe$name)]

# add hg38 distance metric to human 18GW PFC dataset
human_18gw_pfc_abc_score$distance_hg38 <- ifelse(human_18gw_pfc_abc_score$chromEnd_hg38 < human_18gw_pfc_abc_score$TSSchromStart_hg38,
                                                 abs(human_18gw_pfc_abc_score$TSSchromStart_hg38 - human_18gw_pfc_abc_score$chromEnd_hg38),
                                                 abs(human_18gw_pfc_abc_score$chromStart_hg38 - human_18gw_pfc_abc_score$TSSchromEnd_hg38))


# add a new name column to the hg19 WTC11 NGN2-iPSC ExN dataset
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$new_name <- paste0(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, "|", WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chr, ":", WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$start, "-", WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$end)

# add hg38 coordinates to WTC11 NGN2-iPSC ExN dataset
head(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score)
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chromStart_hg38 <- WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$start1[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$new_name, WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$name)]
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chromEnd_hg38 <- WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$end1[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$new_name, WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$name)]
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TSSchromStart_hg38 <- WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$start2[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$new_name, WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$name)]
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TSSchromEnd_hg38 <- WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$end2[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$new_name, WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_bedpe$name)]

# add hg38 distance metric to WTC11 NGN2-iPSC ExN dataset
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$distance_hg38 <- ifelse(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chromEnd_hg38 < WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TSSchromStart_hg38,
                                                            abs(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TSSchromStart_hg38 - WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chromEnd_hg38),
                                                            abs(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chromStart_hg38 - WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TSSchromEnd_hg38))

# Get BrainVar gene expression data
brainvar <- data.frame(read_xlsx(path = "./data/published_datasets/1-s2.0-S2211124720303673-mmc3.xlsx", sheet = "ExpressedGenes_23782"))

# Get Gnomad LOEUF data
gnomad <- read.table(file = "./data/published_datasets/gnomad.v4.1.constraint_metrics.tsv", header = TRUE, sep = "\t")
gnomad_mane <- gnomad[which(gnomad$mane_select == "true" & gnomad$transcript_type == "protein_coding"),]

# Create final supplemental table for human 18GW PFC dataset
human_18gw_pfc_abc_score_table <- data.frame("target_gene_symbol" = human_18gw_pfc_abc_score$geneName,
                                             "ensembl_id" = annotations$gene_id[match(human_18gw_pfc_abc_score$geneName, annotations$gene_name)],
                                             "hgnc_id" = annotations$hgnc_id[match(human_18gw_pfc_abc_score$geneName, annotations$gene_name)],
                                             "dataset" = "GW18_PFC",
                                             "chrom" = human_18gw_pfc_abc_score$chr,
                                             "enhancer_start_hg19" = human_18gw_pfc_abc_score$start,
                                             "enhancer_end_hg19" = human_18gw_pfc_abc_score$end,
                                             "enhancer_start_hg38" = human_18gw_pfc_abc_score$chromStart_hg38,
                                             "enhancer_end_hg38" = human_18gw_pfc_abc_score$chromEnd_hg38,
                                             "enhancer_class" = human_18gw_pfc_abc_score$class,
                                             "TSS_start_hg19" = human_18gw_pfc_abc_score$TargetGeneTSS,
                                             "TSS_end_hg19" = human_18gw_pfc_abc_score$TargetGeneTSS + 1,
                                             "TSS_start_hg38" = human_18gw_pfc_abc_score$TSSchromStart_hg38,
                                             "TSS_end_hg38" = human_18gw_pfc_abc_score$TSSchromEnd_hg38,
                                             "median_brain_expression" = brainvar$MedianExpr[match(human_18gw_pfc_abc_score$geneName, brainvar$GeneSymbol)],
                                             "midfetal_brain_expression" = brainvar$MedianExpr_E1[match(human_18gw_pfc_abc_score$geneName, brainvar$GeneSymbol)],
                                             "postnatal_brain_expression" = brainvar$MedianExpr_E3[match(human_18gw_pfc_abc_score$geneName, brainvar$GeneSymbol)],
                                             "loeuf" = gnomad_mane$lof.oe_ci.upper[match(human_18gw_pfc_abc_score$geneName, gnomad_mane$gene)],
                                             "distance_hg19" = human_18gw_pfc_abc_score$distance,
                                             "distance_hg38" = human_18gw_pfc_abc_score$distance_hg38,
                                             "enhancer_activity" = human_18gw_pfc_abc_score$activity_base,
                                             "hic_contact" = human_18gw_pfc_abc_score$hic_contact,
                                             "hic_contact_scaled" = human_18gw_pfc_abc_score$hic_contact_pl_scaled_adj,
                                             "abc_score_unnormalized" = human_18gw_pfc_abc_score$ABC.Score.Numerator,
                                             "abc_score_normalized" = human_18gw_pfc_abc_score$ABC.Score)


# Create final supplemental table for WTC11 NGN2-iPSC ExN dataset
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_table <- data.frame("target_gene_symbol" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName,
                                                        "ensembl_id" = annotations$gene_id[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, annotations$gene_name)],
                                                        "hgnc_id" = annotations$hgnc_id[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, annotations$gene_name)],
                                                        "dataset" = "NGN2-iNeuron",
                                                        "chrom" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chr,
                                                        "enhancer_start_hg19" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$start,
                                                        "enhancer_end_hg19" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$end,
                                                        "enhancer_start_hg38" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chromStart_hg38,
                                                        "enhancer_end_hg38" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$chromEnd_hg38,
                                                        "enhancer_class" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$class,
                                                        "TSS_start_hg19" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TargetGeneTSS,
                                                        "TSS_end_hg19" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TargetGeneTSS + 1,
                                                        "TSS_start_hg38" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TSSchromStart_hg38,
                                                        "TSS_end_hg38" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$TSSchromEnd_hg38,
                                                        "median_brain_expression" = brainvar$MedianExpr[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, brainvar$GeneSymbol)],
                                                        "midfetal_brain_expression" = brainvar$MedianExpr_E1[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, brainvar$GeneSymbol)],
                                                        "postnatal_brain_expression" = brainvar$MedianExpr_E3[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, brainvar$GeneSymbol)],
                                                        "loeuf" = gnomad_mane$lof.oe_ci.upper[match(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$geneName, gnomad_mane$gene)],
                                                        "distance_hg19" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$distance,
                                                        "distance_hg38" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$distance_hg38,
                                                        "enhancer_activity" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$activity_base,
                                                        "hic_contact" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$hic_contact,
                                                        "hic_contact_scaled" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$hic_contact_pl_scaled_adj,
                                                        "abc_score_unnormalized" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$ABC.Score.Numerator,
                                                        "abc_score_normalized" = WTC11_NGN2_iPSC_7_8wk_ExN_abc_score$ABC.Score)

# combined results into final supplemental table and sort by coordinate
abc_score_table <- rbind(human_18gw_pfc_abc_score_table,
                         WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_table)
abc_score_table <- abc_score_table[order(abc_score_table$chrom, abc_score_table$enhancer_start_hg38, decreasing = FALSE),]

# Double check final table
head(abc_score_table)
dim(abc_score_table)

# Convert supplemental table for ABC Score data into an excel file
write.xlsx2(abc_score_table, file = "./tables/Table_S2_Activity_by_Contact_Analysis_in_GW18_PFC_and_NGN2_iNeurons.xlsx", sheetName = "ABC_Score", 
            col.names = TRUE, row.names = FALSE, append = FALSE)


# Get the brain expression decile for all target genes from human 18GW PFC ABC Score dataset
human_18gw_pfc_abc_score_table$decile <- ntile(human_18gw_pfc_abc_score_table$median_brain_expression, 10)

# Get the number of candidate enhancers per gene from human 18GW PFC ABC Score dataset
human_18gw_pfc_enhancers_per_gene <- human_18gw_pfc_abc_score_table %>% dplyr::count(target_gene_symbol, median_brain_expression, decile)

# Get the mean number of enhancers vs the target gene expression decile for the human 18GW PFC ABC Score dataset
human_18gw_pfc_enhancers_by_expression_decile <- human_18gw_pfc_enhancers_per_gene %>% group_by(decile) %>% summarise(mean = mean(n))

# Get the brain expression decile for all target genes from human 18GW PFC ABC Score dataset
WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_table$decile <- ntile(WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_table$median_brain_expression, 10)

# Get the number of candidate enhancers per gene from human 18GW PFC ABC Score dataset
WTC11_NGN2_iPSC_7_8wk_ExN_enhancers_per_gene <- WTC11_NGN2_iPSC_7_8wk_ExN_abc_score_table %>% dplyr::count(target_gene_symbol, median_brain_expression, decile)

# Get the mean number of enhancers vs the target gene expression decile for the human 18GW PFC ABC Score dataset
WTC11_NGN2_iPSC_7_8wk_ExN_enhancers_by_expression_decile <- WTC11_NGN2_iPSC_7_8wk_ExN_enhancers_per_gene %>% group_by(decile) %>% summarise(mean = mean(n))

# merge datasets and plot final figure
human_18gw_pfc_enhancers_by_expression_decile$dataset <- "GW18_PFC"
WTC11_NGN2_iPSC_7_8wk_ExN_enhancers_by_expression_decile$dataset <- "NGN2_iNeuron"
enhancers_by_expression_decile <- rbind(human_18gw_pfc_enhancers_by_expression_decile,
                                        WTC11_NGN2_iPSC_7_8wk_ExN_enhancers_by_expression_decile)

# change the order of factors for the legend
enhancers_by_expression_decile$dataset <- factor(enhancers_by_expression_decile$dataset, levels = c("NGN2_iNeuron", "GW18_PFC"))

ggplot(data = enhancers_by_expression_decile, aes(x = decile, y = mean, color = dataset)) +
  geom_line(linetype = "dashed", size = 1.5) +
  geom_point(size = 4) + 
  scale_x_continuous(breaks = seq(1, 10, by = 1)) +
  scale_colour_manual(values = c("#276FBF", "#183059")) +
  xlab("\nHuman brain expression decile") +
  ylab("enhancers per gene\n") +
  ylim(c(0, 11)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))


# Load the list of CRISPRa CRT candidate genes
CRISPRa_CRT_candidates <- data.frame(read_xlsx(path = "./tables/Table_S1_Cis_Regulation_Therapy_Candidate_Genes.xlsx", sheet = "Candidate_Genes"))
CRISPRa_CRT_candidates <- CRISPRa_CRT_candidates$gene_symbol

# Label the ABC Score supplemental table by CRT candidate status
abc_score_table$CRISPRa_CRT_candidates <- ifelse(abc_score_table$target_gene_symbol %in% CRISPRa_CRT_candidates, "CRT_candidates", "other_genes")

# Get the number of candidate enhancers per gene from human 18GW PFC ABC Score dataset
enhancers_per_gene <- abc_score_table %>% dplyr::count(target_gene_symbol, dataset, CRISPRa_CRT_candidates)

# Get the mean number of enhancers vs the target gene expression decile for the human 18GW PFC ABC Score dataset
enhancers_by_gene_list <- enhancers_per_gene %>% group_by(dataset, CRISPRa_CRT_candidates) %>% summarise(mean = mean(n))
enhancers_by_gene_list$dataset[which(enhancers_by_gene_list$dataset == "NGN2-iNeuron")] <- "NGN2_iNeuron"

# switch the x-axis order for plotting
enhancers_by_gene_list$CRISPRa_CRT_candidates <- factor(enhancers_by_gene_list$CRISPRa_CRT_candidates, levels = c("other_genes", "CRT_candidates"))

# change the order of factors for the legend
enhancers_by_gene_list$dataset <- factor(enhancers_by_gene_list$dataset, levels = c("NGN2_iNeuron", "GW18_PFC"))

# Plot the number of enhancers by gene list and dataset
ggplot(data = enhancers_by_gene_list, aes(x = CRISPRa_CRT_candidates, y = mean, color = dataset, group = dataset)) +
  geom_line(linetype = "dashed", size = 1.5) +
  geom_point(size = 4) + 
  scale_colour_manual(values = c("#276FBF", "#183059")) +
  ylab("enhancers per gene\n") +
  ylim(c(0, 11)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
        axis.text.y = element_text(size = 16))



### (2) Match promoter contacts and generate supplemental tables for midfetal cerebrum ExN cicero data #################

# Import midfetal ExN cicero data from Domcke et al. 2020 
cicero_coaccessability <- read.table(file = "./data/published_datasets/cicero_coaccess_scores_by_cell_type_ExN.csv", header = TRUE, sep = "\t")

# Convert cicero data frame to a bedpe format
cicero_coaccessability_bedpe <- data.frame("chrom1" = unlist(lapply(strsplit(cicero_coaccessability$Peak1, "-"), `[`, 1)),
                                           "start1" = unlist(lapply(strsplit(cicero_coaccessability$Peak1, "-"), `[`, 2)),
                                           "end1" = unlist(lapply(strsplit(cicero_coaccessability$Peak1, "-"), `[`, 3)),
                                           "chrom2" = unlist(lapply(strsplit(cicero_coaccessability$Peak2, "-"), `[`, 1)),
                                           "start2" = unlist(lapply(strsplit(cicero_coaccessability$Peak2, "-"), `[`, 2)),
                                           "end2" = unlist(lapply(strsplit(cicero_coaccessability$Peak2, "-"), `[`, 3)),
                                           "name" = "NA",
                                           "score" = cicero_coaccessability$ExN,
                                           "strand1" = ".",
                                           "strand2" = ".")

# Create a bedpe file for midfetal cerebrum ExN cicero data
write.table(cicero_coaccessability_bedpe, file = "./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Get the intersection of midfetal cerebrum ExN cicero peaks and 2000bp merged promoters
# This code will output a cicero peak if EITHER end of the pair overlaps with a 2000bp promoter
system("bedtools pairtobed -a ./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.bedpe -b ./bed/annotations/hg19/gencode.v38lift37.basic.coding.transcripts.2000bp_promoters.bed -type either > ./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.2000bp_promoter_overlap.bedpe")

# Get promoter overlapped midfetal cerebrum ExN cicero peaks
cicero_coaccessability_bedpe <- read.table(file = "./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.2000bp_promoter_overlap.bedpe", header = TRUE, sep = "\t",
                                           col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2", 
                                                         "promoter_chrom", "promoter_chromStart", "promoter_chromEnd", "promoter_name", "promoter_score", "promoter_strand"))

# Determine which midfetal cerebrum ExN cicero peak overlaps an annotated promoter
cicero_coaccessability_bedpe$promoter_overlapping_peak <- ifelse(abs(rowMeans(data.frame(cicero_coaccessability_bedpe$promoter_chromStart, cicero_coaccessability_bedpe$promoter_chromEnd)) - rowMeans(data.frame(cicero_coaccessability_bedpe$start1, cicero_coaccessability_bedpe$end1))) < abs(rowMeans(data.frame(cicero_coaccessability_bedpe$promoter_chromStart, cicero_coaccessability_bedpe$promoter_chromEnd)) - rowMeans(data.frame(cicero_coaccessability_bedpe$start2, cicero_coaccessability_bedpe$end2))), "peak_1", "peak_2")

# Create final bedpe file for midfetal cerebrum ExN cicero data
# In the bedpe file that is created, the candidate enhancer peaks comes 1st and the annoated promoter overlapping peaks comes 2nd
cicero_coaccessability_bedpe_FINAL <- data.frame("chrom1" = ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$chrom1, cicero_coaccessability_bedpe$chrom2),
                                                 "start1" = ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$start1, cicero_coaccessability_bedpe$start2),
                                                 "end1" = ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$end1, cicero_coaccessability_bedpe$end2),
                                                 "chrom2" = ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$chrom2, cicero_coaccessability_bedpe$chrom1),
                                                 "start2" = ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$start2, cicero_coaccessability_bedpe$start1),
                                                 "end2" = ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$end2, cicero_coaccessability_bedpe$end1),
                                                 "name" = paste0(sapply(strsplit(cicero_coaccessability_bedpe$promoter_name, ";"), `[`, 1), "|", ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$chrom1, cicero_coaccessability_bedpe$chrom2), ":", ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$start1, cicero_coaccessability_bedpe$start2), "-", ifelse(cicero_coaccessability_bedpe$promoter_overlapping_peak == "peak_2", cicero_coaccessability_bedpe$end1, cicero_coaccessability_bedpe$end2)),
                                                 "score" = cicero_coaccessability_bedpe$score,
                                                 "strand1" = ".",
                                                 "strand2" = ".")

# Write final bedpe file for midfetal cerebrum ExN cicero data
write.table(cicero_coaccessability_bedpe_FINAL, file = "./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Remove temporary files
system("rm ./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.2000bp_promoter_overlap.bedpe")

# LiftOver hg19 bedpe files for Cicero dataset to hg38
system("python3 ./scripts/UCSC/liftOverBedpe/liftOverBedpe.py --lift ./scripts/UCSC/liftOver --chain ./scripts/UCSC/hg19ToHg38.over.chain --i ./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.bedpe --o ./bed/bedpe/hg38/midfetal_cerebrum_ExN_Cicero_liftOver_hg38.bedpe")

# load hg38 dataset for Cicero
cicero_coaccessability_bedpe_hg38 <- read.table(file = "./bed/bedpe/hg38/midfetal_cerebrum_ExN_Cicero_liftOver_hg38.bedpe",
                                                header = FALSE, sep = "\t",
                                                col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))



### (3) Get Trevino and Ziffra datasets and create bedpe files #########################################################

# Create final bedpe file for Trevino dataset
trevino_enhancers_hg38 <- data.frame(read_xlsx(path = "./data/published_datasets/1-s2.0-S0092867421009429-mmc2.xlsx", sheet = "A", skip = 1))
trevino_enhancers_hg38_bedpe <- data.frame("chrom1" = trevino_enhancers_hg38$Peak.chromosome,
                                           "start1" = trevino_enhancers_hg38$Peak.start,
                                           "end1" = trevino_enhancers_hg38$Peak.end,
                                           "chrom2" = trevino_enhancers_hg38$Peak.chromosome,
                                           "start2" = trevino_enhancers_hg38$Gene.TSS.position,
                                           "end2" = trevino_enhancers_hg38$Gene.TSS.position + 1,
                                           "name" = paste0(trevino_enhancers_hg38$Gene.symbol, "|", trevino_enhancers_hg38$Peak.chromosome, ":", trevino_enhancers_hg38$Peak.start, "-", trevino_enhancers_hg38$Peak.end),
                                           "score" = trevino_enhancers_hg38$correlation,
                                           "strand1" = ".",
                                           "strand2" = ".")

write.table(trevino_enhancers_hg38_bedpe, file = "./bed/bedpe/hg38/trevino_fetal_cortex_enhancer_predictions.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Get ABC contacts from Ziffra dataset
ziffra_abc_score_peaks_hg38 <- data.frame(read_xlsx(path = "./data/published_datasets/41586_2021_3209_MOESM4_ESM.xlsx", sheet = "ST2 AllPrimaryPeaks"))[,c(1:4)]
ziffra_abc_score_contacts_hg38 <- data.frame(read_xlsx(path = "./data/published_datasets/41586_2021_3209_MOESM4_ESM.xlsx", sheet = "ST5 EnhancerPeaks"))
ziffra_abc_score_contacts_hg38 <- ziffra_abc_score_contacts_hg38[,c(1, 13, 14, 21)]
ziffra_abc_score_contacts_hg38$seqnames <- ziffra_abc_score_peaks_hg38$seqnames[match(ziffra_abc_score_contacts_hg38$peak_name, ziffra_abc_score_peaks_hg38$peak_name)]
ziffra_abc_score_contacts_hg38$start <- ziffra_abc_score_peaks_hg38$start[match(ziffra_abc_score_contacts_hg38$peak_name, ziffra_abc_score_peaks_hg38$peak_name)]
ziffra_abc_score_contacts_hg38$end <- ziffra_abc_score_peaks_hg38$end[match(ziffra_abc_score_contacts_hg38$peak_name, ziffra_abc_score_peaks_hg38$peak_name)]
ziffra_abc_score <- ziffra_abc_score_contacts_hg38[,c(5:7, 1:4)]
head(ziffra_abc_score)

ziffra_abc_score$length <- ziffra_abc_score$end - ziffra_abc_score$start
summary(ziffra_abc_score$length) # max peak length is 4,497 bp

# Convert bigBed contact files to bedpe files for Ziffra dataset
system("./scripts/UCSC/bigBedToBed ./data/published_datasets/bigBed/dlEN.inter.bb ./data/published_datasets/bedpe/dlEN.inter.bedpe")
system("./scripts/UCSC/bigBedToBed ./data/published_datasets/bigBed/earlyEN.inter.bb ./data/published_datasets/bedpe/earlyEN.inter.bedpe")
system("./scripts/UCSC/bigBedToBed ./data/published_datasets/bigBed/ulEN.inter.bb ./data/published_datasets/bedpe/ulEN.inter.bedpe")

system("head ./data/published_datasets/bedpe/dlEN.inter.bedpe")


ziffra_dlEN_hg38 <- read.table(file = "./data/published_datasets/bedpe/dlEN.inter.bedpe", header = FALSE, sep = "\t",
                               col.names = c("chromWindow", "chromStartWindow", "chromEndWindow", "nameWindow", "null", "scoreWindow", "strandWindow", "color",
                                             "chromPeak", "chromStartPeak", "chromEndPeak", "namePeak", "strandPeak", "chromTSS", "chromStartTSS", "chromEndTSS", "nameTSS", "strandTSS"))
head(ziffra_dlEN_hg38)
ziffra_dlEN_hg38$length_peak <- ziffra_dlEN_hg38$chromEndPeak - ziffra_dlEN_hg38$chromStartPeak
summary(ziffra_dlEN_hg38$length_peak)

ziffra_dlEN_hg38_bedpe <- data.frame("chrom1" = ziffra_dlEN_hg38$chromPeak,
                                     "start1" = ziffra_dlEN_hg38$chromStartPeak,
                                     "end1" = ziffra_dlEN_hg38$chromEndPeak,
                                     "chrom2" = ziffra_dlEN_hg38$chromTSS,
                                     "start2" = ziffra_dlEN_hg38$chromStartTSS,
                                     "end2" = ziffra_dlEN_hg38$chromEndTSS + 1,
                                     "name" = ziffra_dlEN_hg38$nameWindow,
                                     "score" = ziffra_dlEN_hg38$scoreWindow,
                                     "strand1" = ".",
                                     "strand2" = ".")

ziffra_earlyEN_hg38 <- read.table(file = "./data/published_datasets/bedpe/earlyEN.inter.bedpe", header = FALSE, sep = "\t",
                                  col.names = c("chromWindow", "chromStartWindow", "chromEndWindow", "nameWindow", "null", "scoreWindow", "strandWindow", "color",
                                                "chromPeak", "chromStartPeak", "chromEndPeak", "namePeak", "strandPeak", "chromTSS", "chromStartTSS", "chromEndTSS", "nameTSS", "strandTSS"))
ziffra_earlyEN_hg38_bedpe <- data.frame("chrom1" = ziffra_earlyEN_hg38$chromPeak,
                                        "start1" = ziffra_earlyEN_hg38$chromStartPeak,
                                        "end1" = ziffra_earlyEN_hg38$chromEndPeak,
                                        "chrom2" = ziffra_earlyEN_hg38$chromTSS,
                                        "start2" = ziffra_earlyEN_hg38$chromStartTSS,
                                        "end2" = ziffra_earlyEN_hg38$chromEndTSS + 1,
                                        "name" = ziffra_earlyEN_hg38$nameWindow,
                                        "score" = ziffra_earlyEN_hg38$scoreWindow,
                                        "strand1" = ".",
                                        "strand2" = ".")

ziffra_ulEN_hg38 <- read.table(file = "./data/published_datasets/bedpe/ulEN.inter.bedpe", header = FALSE, sep = "\t",
                               col.names = c("chromWindow", "chromStartWindow", "chromEndWindow", "nameWindow", "null", "scoreWindow", "strandWindow", "color",
                                             "chromPeak", "chromStartPeak", "chromEndPeak", "namePeak", "strandPeak", "chromTSS", "chromStartTSS", "chromEndTSS", "nameTSS", "strandTSS"))
ziffra_ulEN_hg38_bedpe <- data.frame("chrom1" = ziffra_ulEN_hg38$chromPeak,
                                     "start1" = ziffra_ulEN_hg38$chromStartPeak,
                                     "end1" = ziffra_ulEN_hg38$chromEndPeak,
                                     "chrom2" = ziffra_ulEN_hg38$chromTSS,
                                     "start2" = ziffra_ulEN_hg38$chromStartTSS,
                                     "end2" = ziffra_ulEN_hg38$chromEndTSS + 1,
                                     "name" = ziffra_ulEN_hg38$nameWindow,
                                     "score" = ziffra_ulEN_hg38$scoreWindow,
                                     "strand1" = ".",
                                     "strand2" = ".")

# Merge bedpe files from each cell-type in the Ziffra dataset
ziffra_hg38_bedpe <- rbind(ziffra_dlEN_hg38_bedpe,
                           ziffra_earlyEN_hg38_bedpe,
                           ziffra_ulEN_hg38_bedpe)
write.table(ziffra_hg38_bedpe, file = "./bed/bedpe/hg38/ziffra_fetal_human_cortex_ABC_Score.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



### (4) Create hg19 and hg38 bigInteract files for each dataset ########################################################

# Create hg19 bigInteract file for each dataset
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg19/gw18_pfc_ABC_Score.bedpe --genome hg19 --outdir ./UCSC_browser_tracks/hg19/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg19/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score.bedpe --genome hg19 --outdir ./UCSC_browser_tracks/hg19/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg19/midfetal_cerebrum_ExN_Cicero.bedpe --genome hg19 --outdir ./UCSC_browser_tracks/hg19/bigInteract/")

# Create hg38 bigInteract file for each dataset
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/gw18_pfc_ABC_Score_liftOver_hg38.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_liftOver_hg38.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/midfetal_cerebrum_ExN_Cicero_liftOver_hg38.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/trevino_fetal_cortex_enhancer_predictions.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/ziffra_fetal_human_cortex_ABC_Score.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")



### (5) Get candidate NDD enhancers for each dataset and merge them into a single bed file and make upset plot #########

# Load the list of CRISPRa CRT candidate genes
CRISPRa_CRT_candidates <- data.frame(read_xlsx(path = "./tables/Table_S1_Cis_Regulation_Therapy_Candidate_Genes.xlsx", sheet = "Candidate_Genes"))
CRISPRa_CRT_candidates <- CRISPRa_CRT_candidates$gene_symbol

# Get candidate NDD enhancers for each dataset
human_18gw_pfc_abc_score_bedpe$gene_symbol <- sapply(strsplit(human_18gw_pfc_abc_score_bedpe$name, "|", fixed = TRUE), `[`, 1)
human_18gw_pfc_abc_score_bedpe_NDD <- human_18gw_pfc_abc_score_bedpe[which(human_18gw_pfc_abc_score_bedpe$gene_symbol %in% CRISPRa_CRT_candidates),]
WTC11_NGN2_iPSC_abc_score_bedpe$gene_symbol <- sapply(strsplit(WTC11_NGN2_iPSC_abc_score_bedpe$name, "|", fixed = TRUE), `[`, 1)
WTC11_NGN2_iPSC_abc_score_bedpe_NDD <- WTC11_NGN2_iPSC_abc_score_bedpe[which(WTC11_NGN2_iPSC_abc_score_bedpe$gene_symbol %in% CRISPRa_CRT_candidates),]
cicero_coaccessability_bedpe_hg38$gene_symbol <- sapply(strsplit(cicero_coaccessability_bedpe_hg38$name, "|", fixed = TRUE), `[`, 1)
cicero_coaccessability_bedpe_NDD <- cicero_coaccessability_bedpe_hg38[which(cicero_coaccessability_bedpe_hg38$gene_symbol %in% CRISPRa_CRT_candidates),]
trevino_enhancers_hg38_bedpe$gene_symbol <- sapply(strsplit(trevino_enhancers_hg38_bedpe$name, "|", fixed = TRUE), `[`, 1)
trevino_enhancers_bedpe_NDD <- trevino_enhancers_hg38_bedpe[which(trevino_enhancers_hg38_bedpe$gene_symbol %in% CRISPRa_CRT_candidates),]
ziffra_abc_score_NDD <- ziffra_abc_score[which(ziffra_abc_score$dlEN_targetgene %in% CRISPRa_CRT_candidates |
                                               ziffra_abc_score$earlyEN_targetgene %in% CRISPRa_CRT_candidates |
                                               ziffra_abc_score$ulEN_targetgene %in% CRISPRa_CRT_candidates),]

dim(human_18gw_pfc_abc_score_bedpe_NDD) # 2,483
dim(WTC11_NGN2_iPSC_abc_score_bedpe_NDD) # 2,479
dim(cicero_coaccessability_bedpe_NDD) # 3,890
dim(trevino_enhancers_bedpe_NDD) # 2,841
dim(ziffra_abc_score_NDD) # 598


# merge ziffra abc score gene symbols into a single list
gene_list <- c()
for(i in 1:dim(ziffra_abc_score_NDD)[1]){
  gene_symbols <- paste0(unique(as.list(ziffra_abc_score_NDD[i, c(5:7)][c(ziffra_abc_score_NDD$dlEN_targetgene[i],
                                                                          ziffra_abc_score_NDD$earlyEN_targetgene[i],
                                                                          ziffra_abc_score_NDD$ulEN_targetgene[i]) %in% CRISPRa_CRT_candidates])), collapse = ",")
  gene_list <- c(gene_list, gene_symbols)
}
ziffra_abc_score_NDD$gene_symbol <- gene_list
head(ziffra_abc_score_NDD)

# create separate rows for candidate enhancers linked to multiple CRT candidate genes
ziffra_abc_score_NDD_split <- ziffra_abc_score_NDD %>%
  separate_rows(gene_symbol, sep = ",")

# get TSS coordinates from gencode
ziffra_abc_score_NDD_split$TSS <- annotations$start_position_hg38[match(ziffra_abc_score_NDD_split$gene_symbol, annotations$gene_name)]

# create bedpe file for NDD specific ziffra abc score dataset
ziffra_abc_score_bedpe_NDD <- data.frame("chrom1" = ziffra_abc_score_NDD_split$seqnames,
                                         "start1" = ziffra_abc_score_NDD_split$start,
                                         "end1" = ziffra_abc_score_NDD_split$end,
                                         "chrom2" = ziffra_abc_score_NDD_split$seqnames,
                                         "start2" = ziffra_abc_score_NDD_split$TSS,
                                         "end2" = ziffra_abc_score_NDD_split$TSS + 1,
                                         "name" = paste0(ziffra_abc_score_NDD_split$gene_symbol, "|", ziffra_abc_score_NDD_split$seqnames, ":", ziffra_abc_score_NDD_split$start, "-", ziffra_abc_score_NDD_split$end),
                                         "score" = "1",
                                         "strand1" = ".",
                                         "strand2" = ".")

# Write NDD specific bedpe files
write.table(human_18gw_pfc_abc_score_bedpe_NDD, file = "./bed/bedpe/hg38/gw18_pfc_ABC_Score_liftOver_hg38_NDD.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(WTC11_NGN2_iPSC_abc_score_bedpe_NDD, file = "./bed/bedpe/hg38/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_liftOver_hg38_NDD.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(cicero_coaccessability_bedpe_NDD, file = "./bed/bedpe/hg38/midfetal_cerebrum_ExN_Cicero_liftOver_hg38_NDD.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(trevino_enhancers_bedpe_NDD, file = "./bed/bedpe/hg38/trevino_fetal_cortex_enhancer_predictions_NDD.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(ziffra_abc_score_bedpe_NDD, file = "./bed/bedpe/hg38/ziffra_fetal_human_cortex_ABC_Score_NDD.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Create hg38 bigInteract file for NDD specific datasets
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/gw18_pfc_ABC_Score_liftOver_hg38_NDD.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_liftOver_hg38_NDD.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/midfetal_cerebrum_ExN_Cicero_liftOver_hg38_NDD.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/trevino_fetal_cortex_enhancer_predictions_NDD.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")
system("python3 ./scripts/python/get_bigInteract_from_bedpe.py --bedpe ./bed/bedpe/hg38/ziffra_fetal_human_cortex_ABC_Score_NDD.bedpe --genome hg38 --outdir ./UCSC_browser_tracks/hg38/bigInteract/")


# Merge all candidate NDD enhancers into a single bed file
all_canididate_peaks <- data.frame("chrom" = human_18gw_pfc_abc_score_bedpe_NDD$chrom1,
                                   "chromStart" = human_18gw_pfc_abc_score_bedpe_NDD$start1,
                                   "chromEnd" = human_18gw_pfc_abc_score_bedpe_NDD$end1,
                                   "name" = paste0(human_18gw_pfc_abc_score_bedpe_NDD$gene_symbol, ".GW18_PFC_ABC"))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = WTC11_NGN2_iPSC_abc_score_bedpe_NDD$chrom1,
                                         "chromStart" = WTC11_NGN2_iPSC_abc_score_bedpe_NDD$start1, 
                                         "chromEnd" = WTC11_NGN2_iPSC_abc_score_bedpe_NDD$end1,
                                         "name" = paste0(WTC11_NGN2_iPSC_abc_score_bedpe_NDD$gene_symbol, ".NGN2_iPSC_ABC")))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = cicero_coaccessability_bedpe_NDD$chrom1,
                                         "chromStart" = cicero_coaccessability_bedpe_NDD$start1, 
                                         "chromEnd" = cicero_coaccessability_bedpe_NDD$end1,
                                         "name" = paste0(cicero_coaccessability_bedpe_NDD$gene_symbol, ".Fetal_Cerebrum_Cicero")))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = trevino_enhancers_bedpe_NDD$chrom1,
                                         "chromStart" = trevino_enhancers_bedpe_NDD$start1, 
                                         "chromEnd" = trevino_enhancers_bedpe_NDD$end1,
                                         "name" = paste0(trevino_enhancers_bedpe_NDD$gene_symbol, ".Midfetal_Cortex_Trevino")))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = ziffra_abc_score_NDD$seqnames,
                                         "chromStart" = ziffra_abc_score_NDD$start, 
                                         "chromEnd" = ziffra_abc_score_NDD$end,
                                         "name" = paste0(ziffra_abc_score_NDD$gene_symbol, ".Midfetal_Cortex_Ziffra_ABC")))

# Extend all peaks to at least 500bp
for (i in 1:dim(all_canididate_peaks)[1]){
  if (all_canididate_peaks$chromEnd[i] - all_canididate_peaks$chromStart[i] < 500){
    length <- 500 - (all_canididate_peaks$chromEnd[i] - all_canididate_peaks$chromStart[i])
    all_canididate_peaks$chromStart[i] <- all_canididate_peaks$chromStart[i] - floor(length/2)
    all_canididate_peaks$chromEnd[i] <- all_canididate_peaks$chromEnd[i] + ceiling(length/2)
  }
}


write.table(all_canididate_peaks, file = "./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

system("bedtools sort -i ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed > ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.bed")
system("bedtools merge -i ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.bed -c 4 -o distinct > ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.bed")
system("wc -l ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.bed") # 6,765 total candidate enhancer regions

merged_candidate_enhancer_peaks <- read.table(file = "./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.bed", header = FALSE, sep = "\t",
                                              col.names = c("chrom", "chromStart", "chromEnd", "name"))

# Remove peaks that overlap with a promoter or exon of a protein coding gene
system("bedtools intersect -a ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.bed -b ./bed/annotations/hg38/gencode.v38.basic.coding.transcripts.2000bp_promoters.bed ./bed/annotations/hg38/gencode.v38.basic.annotation.coding.exons.bed -wa -v > ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.promoters_exons_removed.bed")
system("wc -l ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.promoters_exons_removed.bed") # 5,425 final candidate enhancer regions

# Remove temporary files and rename final enhancer list
system("rm ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed")
system("rm ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.bed")
system("rm ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.bed")
system("mv ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.sorted.merged.promoters_exons_removed.bed ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed")

candidate_enhancer_regions_for_mpra <- read.table(file = "./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed", header = FALSE, sep = "\t",
                                                     col.names = c("chrom", "chromStart", "chromEnd", "name"))

# Get bigBed file of candidate enhancers
system("./scripts/UCSC/bedToBigBed ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed ~/genomes/hg38/hg38.chrom.sizes ./UCSC_browser_tracks/hg38/bigBed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bigBed")

head(candidate_enhancer_regions_for_mpra)
dim(candidate_enhancer_regions_for_mpra) # 5,425 total candidate enhancer regions

# annotate candidate enhancer regions for all genes by dataset
datasets <- data.frame("GW18_PFC_ABC" = c(),
                       "NGN2_iNeuron_ABC" = c(),
                       "Midfetal_Cortex_Ziffra_ABC" = c(),
                       "Fetal_Cerebrum_Cicero" = c(),
                       "Midfetal_Cortex_Trevino" = c())

for(i in 1:dim(candidate_enhancer_regions_for_mpra)[1]){
  which_datasets <- data.frame("GW18_PFC_ABC" = c(0),
                               "NGN2_iNeuron_ABC" = c(0),
                               "Midfetal_Cortex_Ziffra_ABC" = c(0),
                               "Fetal_Cerebrum_Cicero" = c(0),
                               "Midfetal_Cortex_Trevino" = c(0))
  
  if(str_contains(candidate_enhancer_regions_for_mpra$name[i], "GW18_PFC_ABC")){
    which_datasets$GW18_PFC_ABC <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_mpra$name[i], "NGN2_iPSC_ABC")){
    which_datasets$NGN2_iNeuron_ABC <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_mpra$name[i], "Midfetal_Cortex_Ziffra_ABC")){
    which_datasets$Midfetal_Cortex_Ziffra_ABC <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_mpra$name[i], "Fetal_Cerebrum_Cicero")){
    which_datasets$Fetal_Cerebrum_Cicero <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_mpra$name[i], "Midfetal_Cortex_Trevino")){
    which_datasets$Midfetal_Cortex_Trevino <- 1
  }
  
  datasets <- rbind(datasets, which_datasets)
  
}

head(datasets)
dim(datasets) # 5,425 combined elements

# Combine dataset list with candidate enhancer table
candidate_enhancer_regions_for_mpra <- cbind(candidate_enhancer_regions_for_mpra, datasets)
head(candidate_enhancer_regions_for_mpra)
dim(candidate_enhancer_regions_for_mpra)
colnames(candidate_enhancer_regions_for_mpra)

# Make combination matrix for candidate enhancer table
m1 = make_comb_mat(candidate_enhancer_regions_for_mpra[,c(5:9)])
m1

# Write candidate enhancer table to a pdf
pdf('./plots/Candidate_Enhancers/candidate_enhancer_regions_for_mpra_upset_plot.pdf', height = 3, width = 8)
UpSet(m1, set_order = c("Midfetal_Cortex_Trevino", "NGN2_iNeuron_ABC", "GW18_PFC_ABC", "Fetal_Cerebrum_Cicero", "Midfetal_Cortex_Ziffra_ABC"), comb_order = order(comb_size(m1), decreasing = TRUE), pt_size = unit(2, "mm"), row_names_side = "left", left_annotation = NULL, right_annotation = NULL,
      top_annotation = upset_top_annotation(m1, add_numbers = TRUE), pt_size = unit(5, "mm"), lwd = 3)
dev.off()



### (6) Get genome wide background set of candidate enhancers and make upset plot ######################################

# Load hg38 bedpe files for each dataset
human_18gw_pfc_abc_score_hg38 <- read.table(file = "./bed/bedpe/hg38/gw18_pfc_ABC_Score_liftOver_hg38.bedpe",
                                       header = FALSE, sep = "\t",
                                       col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))
WTC11_NGN2_iPSC_abc_score_hg38 <- read.table(file = "./bed/bedpe/hg38/WTC11_NGN2_iPSC_7-8wk_ExN_ABC_Score_liftOver_hg38.bedpe",
                                        header = FALSE, sep = "\t",
                                        col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))
cicero_coaccessability_hg38 <- read.table(file = "./bed/bedpe/hg38/midfetal_cerebrum_ExN_Cicero_liftOver_hg38.bedpe",
                                          header = FALSE, sep = "\t",
                                          col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))
trevino_enhancers_hg38 <- read.table(file = "./bed/bedpe/hg38/trevino_fetal_cortex_enhancer_predictions.bedpe",
                                     header = FALSE, sep = "\t",
                                     col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))
ziffra_abc_score_hg38 <- read.table(file = "./bed/bedpe/hg38/ziffra_fetal_human_cortex_ABC_Score.bedpe",
                                    header = FALSE, sep = "\t",
                                    col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))

human_18gw_pfc_abc_score_hg38$gene_symbol <- sapply(strsplit(human_18gw_pfc_abc_score_hg38$name, "|", fixed = TRUE), `[`, 1)
human_18gw_pfc_abc_score_hg38_NDD <- human_18gw_pfc_abc_score_hg38[which(human_18gw_pfc_abc_score_hg38$gene_symbol %in% CRISPRa_CRT_candidates),]
WTC11_NGN2_iPSC_abc_score_hg38$gene_symbol <- sapply(strsplit(WTC11_NGN2_iPSC_abc_score_hg38$name, "|", fixed = TRUE), `[`, 1)
WTC11_NGN2_iPSC_abc_score_hg38_NDD <- WTC11_NGN2_iPSC_abc_score_hg38[which(WTC11_NGN2_iPSC_abc_score_hg38$gene_symbol %in% CRISPRa_CRT_candidates),]
cicero_coaccessability_hg38$gene_symbol <- sapply(strsplit(cicero_coaccessability_hg38$name, "|", fixed = TRUE), `[`, 1)
cicero_coaccessability_hg38_NDD <- cicero_coaccessability_hg38[which(cicero_coaccessability_hg38$gene_symbol %in% CRISPRa_CRT_candidates),]
trevino_enhancers_hg38$gene_symbol <- sapply(strsplit(trevino_enhancers_hg38$name, "|", fixed = TRUE), `[`, 1)
trevino_enhancers_hg38_NDD <- trevino_enhancers_hg38[which(trevino_enhancers_hg38$gene_symbol %in% CRISPRa_CRT_candidates),]
ziffra_abc_score_hg38$gene_symbol <- sapply(strsplit(ziffra_abc_score_hg38$name, "_", fixed = TRUE), `[`, 1)
ziffra_abc_score_hg38_NDD <- ziffra_abc_score_hg38[which(ziffra_abc_score_hg38$gene_symbol %in% CRISPRa_CRT_candidates),]

# Make an NDD specific version of ziffra_abc_score for plotting only
write.table(ziffra_abc_score_hg38_NDD, file = "./bed/bedpe/hg38/ziffra_fetal_human_cortex_ABC_Score_NDD_for_gviz_figures_only.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# Merge all candidate NDD enhancers into a single bed file
all_canididate_peaks <- data.frame("chrom" = human_18gw_pfc_abc_score_hg38$chrom1,
                                   "chromStart" = human_18gw_pfc_abc_score_hg38$start1,
                                   "chromEnd" = human_18gw_pfc_abc_score_hg38$end1,
                                   "name" = paste0(human_18gw_pfc_abc_score_hg38$score, "|", human_18gw_pfc_abc_score_hg38$gene_symbol, "|GW18_PFC_ABC"))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = WTC11_NGN2_iPSC_abc_score_hg38$chrom1,
                                         "chromStart" = WTC11_NGN2_iPSC_abc_score_hg38$start1, 
                                         "chromEnd" = WTC11_NGN2_iPSC_abc_score_hg38$end1,
                                         "name" = paste0(WTC11_NGN2_iPSC_abc_score_hg38$score, "|", WTC11_NGN2_iPSC_abc_score_hg38$gene_symbol, "|NGN2_iPSC_ABC")))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = cicero_coaccessability_hg38$chrom1,
                                         "chromStart" = cicero_coaccessability_hg38$start1, 
                                         "chromEnd" = cicero_coaccessability_hg38$end1,
                                         "name" = paste0(cicero_coaccessability_hg38$score, "|", cicero_coaccessability_hg38$gene_symbol, "|Fetal_Cerebrum_Cicero")))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = trevino_enhancers_hg38$chrom1,
                                         "chromStart" = trevino_enhancers_hg38$start1, 
                                         "chromEnd" = trevino_enhancers_hg38$end1,
                                         "name" = paste0(trevino_enhancers_hg38$score, "|", trevino_enhancers_hg38$gene_symbol, "|Midfetal_Cortex_Trevino")))
all_canididate_peaks <- rbind(all_canididate_peaks,
                              data.frame("chrom" = ziffra_abc_score_hg38$chrom1, 
                                         "chromStart" = ziffra_abc_score_hg38$start1, 
                                         "chromEnd" = ziffra_abc_score_hg38$end1,
                                         "name" = paste0(ziffra_abc_score_hg38$score, "|", ziffra_abc_score_hg38$gene_symbol, "|Midfetal_Cortex_Ziffra_ABC")))

# Extend all peaks to at least 500bp
for (i in 1:dim(all_canididate_peaks)[1]){
  if (all_canididate_peaks$chromEnd[i] - all_canididate_peaks$chromStart[i] < 500){
    length <- 500 - (all_canididate_peaks$chromEnd[i] - all_canididate_peaks$chromStart[i])
    all_canididate_peaks$chromStart[i] <- all_canididate_peaks$chromStart[i] - floor(length/2)
    all_canididate_peaks$chromEnd[i] <- all_canididate_peaks$chromEnd[i] + ceiling(length/2)
  }
}


head(all_canididate_peaks)
dim(all_canididate_peaks)

all_canididate_peaks$length <- all_canididate_peaks$chromEnd - all_canididate_peaks$chromStart
summary(all_canididate_peaks$length)

write.table(all_canididate_peaks, file = "./bed/candidate_enhancer_regions_for_all_genes.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

system("bedtools sort -i ./bed/candidate_enhancer_regions_for_all_genes.bed > ./bed/candidate_enhancer_regions_for_all_genes.sorted.bed")
system("bedtools merge -i ./bed/candidate_enhancer_regions_for_all_genes.sorted.bed -c 4 -o distinct > ./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.bed")
system("wc -l ./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.bed") # 98,360 total candidate enhancer regions

merged_candidate_enhancer_peaks <- read.table(file = "./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.bed", header = FALSE, sep = "\t",
                                              col.names = c("chrom", "chromStart", "chromEnd", "name"))

# Remove peaks that overlap with a promoter or exon of a protein coding gene
system("bedtools intersect -a ./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.bed -b ./bed/annotations/hg38/gencode.v38.basic.coding.transcripts.2000bp_promoters.bed ./bed/annotations/hg38/gencode.v38.basic.annotation.coding.exons.bed -wa -v > ./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.promoters_exons_removed.bed")
system("wc -l ./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.promoters_exons_removed.bed") # 80,533 final candidate enhancer regions

# Remove temporary files and rename final enhancer list
system("rm ./bed/candidate_enhancer_regions_for_all_genes.bed")
system("rm ./bed/candidate_enhancer_regions_for_all_genes.sorted.bed")
system("rm ./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.bed")
system("mv ./bed/candidate_enhancer_regions_for_all_genes.sorted.merged.promoters_exons_removed.bed ./bed/candidate_enhancer_regions_for_all_genes.bed")

candidate_enhancer_regions_for_all_genes <- read.table(file = "./bed/candidate_enhancer_regions_for_all_genes.bed", header = FALSE, sep = "\t",
                                                       col.names = c("chrom", "chromStart", "chromEnd", "name"))

head(candidate_enhancer_regions_for_all_genes)
dim(candidate_enhancer_regions_for_all_genes) # 80,533 candidate enhancers across all 5 datasets


# annotate candidate enhancer regions for all genes by dataset
datasets <- data.frame("GW18_PFC_ABC" = c(),
                       "NGN2_iNeuron_ABC" = c(),
                       "Midfetal_Cortex_Ziffra_ABC" = c(),
                       "Fetal_Cerebrum_Cicero" = c(),
                       "Midfetal_Cortex_Trevino" = c())

for(i in 1:dim(candidate_enhancer_regions_for_all_genes)[1]){
  which_datasets <- data.frame("GW18_PFC_ABC" = c(0),
                               "NGN2_iNeuron_ABC" = c(0),
                               "Midfetal_Cortex_Ziffra_ABC" = c(0),
                               "Fetal_Cerebrum_Cicero" = c(0),
                               "Midfetal_Cortex_Trevino" = c(0))
  
  if(str_contains(candidate_enhancer_regions_for_all_genes$name[i], "GW18_PFC_ABC")){
    which_datasets$GW18_PFC_ABC <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_all_genes$name[i], "NGN2_iPSC_ABC")){
    which_datasets$NGN2_iNeuron_ABC <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_all_genes$name[i], "Midfetal_Cortex_Ziffra_ABC")){
    which_datasets$Midfetal_Cortex_Ziffra_ABC <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_all_genes$name[i], "Fetal_Cerebrum_Cicero")){
    which_datasets$Fetal_Cerebrum_Cicero <- 1
  }
  
  if(str_contains(candidate_enhancer_regions_for_all_genes$name[i], "Midfetal_Cortex_Trevino")){
    which_datasets$Midfetal_Cortex_Trevino <- 1
  }
  
  datasets <- rbind(datasets, which_datasets)
  
}

head(datasets)
dim(datasets) # 80,533 combined elements

# Combine dataset list with candidate enhancer table
candidate_enhancer_regions_for_all_genes <- cbind(candidate_enhancer_regions_for_all_genes, datasets)
head(candidate_enhancer_regions_for_all_genes)
dim(candidate_enhancer_regions_for_all_genes)
colnames(candidate_enhancer_regions_for_all_genes)

# Make combination matrix for candidate enhancer table
m1 = make_comb_mat(candidate_enhancer_regions_for_all_genes[,c(5:9)])
m1

# Write candidate enhancer table to a pdf
pdf('./plots/Candidate_Enhancers/candidate_enhancer_regions_for_all_genes_upset_plot.pdf', height = 3, width = 8)
UpSet(m1, set_order = c("Midfetal_Cortex_Trevino", "NGN2_iNeuron_ABC", "Fetal_Cerebrum_Cicero", "GW18_PFC_ABC", "Midfetal_Cortex_Ziffra_ABC"), comb_order = order(comb_size(m1), decreasing = TRUE), row_names_side = "left", left_annotation = NULL, right_annotation = NULL,
      top_annotation = upset_top_annotation(m1, add_numbers = TRUE), pt_size = unit(4, "mm"), lwd = 3)
dev.off()



### (7) Plot pairwise enhancer score correlations for all 5 prioritization datasets ####################################

# create a name for every element in final candidate enhancer list
candidate_enhancer_regions_for_all_genes$coords <- paste0(candidate_enhancer_regions_for_all_genes$chrom, ":", candidate_enhancer_regions_for_all_genes$chromStart, "-", candidate_enhancer_regions_for_all_genes$chromEnd)

# now separate comma lists into distinct elements
candidate_enhancer_regions_for_all_genes <- candidate_enhancer_regions_for_all_genes %>% separate_rows(name, sep = ",")
candidate_enhancer_regions_for_all_genes <- data.frame(candidate_enhancer_regions_for_all_genes)

# get the dataset, gene name, and score for each enhancer
candidate_enhancer_regions_for_all_genes$score <- sapply(strsplit(candidate_enhancer_regions_for_all_genes$name, "|", fixed = TRUE), `[`, 1)
candidate_enhancer_regions_for_all_genes$gene <- sapply(strsplit(candidate_enhancer_regions_for_all_genes$name, "|", fixed = TRUE), `[`, 2)
candidate_enhancer_regions_for_all_genes$dataset <- sapply(strsplit(candidate_enhancer_regions_for_all_genes$name, "|", fixed = TRUE), `[`, 3)

# get a list of gene enhancer pairs for all candidate enhancers
gene_enhancer_pairs <- unique(candidate_enhancer_regions_for_all_genes[,c(10, 12)])

# get the score information for each gene enhancer pair across multiple prioritization methods
score_table <- data.frame("GW18" = c(), 
                          "WTC11" = c(), 
                          "Cicero" = c(), 
                          "Trevino" = c(), 
                          "Ziffra" = c())

### Takes a long time to run ###
for(i in 1:dim(gene_enhancer_pairs)[1]){
  GW18 <- candidate_enhancer_regions_for_all_genes$score[intersect(intersect(which(candidate_enhancer_regions_for_all_genes$coords == gene_enhancer_pairs$coords[i]), which(candidate_enhancer_regions_for_all_genes$gene == gene_enhancer_pairs$gene[i])), which(candidate_enhancer_regions_for_all_genes$dataset == "GW18_PFC_ABC"))]
  WTC11 <- candidate_enhancer_regions_for_all_genes$score[intersect(intersect(which(candidate_enhancer_regions_for_all_genes$coords == gene_enhancer_pairs$coords[i]), which(candidate_enhancer_regions_for_all_genes$gene == gene_enhancer_pairs$gene[i])), which(candidate_enhancer_regions_for_all_genes$dataset == "NGN2_iPSC_ABC"))]
  Cicero <- candidate_enhancer_regions_for_all_genes$score[intersect(intersect(which(candidate_enhancer_regions_for_all_genes$coords == gene_enhancer_pairs$coords[i]), which(candidate_enhancer_regions_for_all_genes$gene == gene_enhancer_pairs$gene[i])), which(candidate_enhancer_regions_for_all_genes$dataset == "Fetal_Cerebrum_Cicero"))]
  Trevino <- candidate_enhancer_regions_for_all_genes$score[intersect(intersect(which(candidate_enhancer_regions_for_all_genes$coords == gene_enhancer_pairs$coords[i]), which(candidate_enhancer_regions_for_all_genes$gene == gene_enhancer_pairs$gene[i])), which(candidate_enhancer_regions_for_all_genes$dataset == "Midfetal_Cortex_Trevino"))]
  Ziffra <- candidate_enhancer_regions_for_all_genes$score[intersect(intersect(which(candidate_enhancer_regions_for_all_genes$coords == gene_enhancer_pairs$coords[i]), which(candidate_enhancer_regions_for_all_genes$gene == gene_enhancer_pairs$gene[i])), which(candidate_enhancer_regions_for_all_genes$dataset == "Midfetal_Cortex_Ziffra_ABC"))]
  
  GW18 <- ifelse(length(GW18) == 0, NA, GW18)
  WTC11 <- ifelse(length(WTC11) == 0, NA, WTC11)
  Cicero <- ifelse(length(Cicero) == 0, NA, Cicero)
  Trevino <- ifelse(length(Trevino) == 0, NA, Trevino)
  Ziffra <- ifelse(length(Ziffra) == 0, NA, Ziffra)
  
  score_table <- rbind(score_table,
                         data.frame("GW18" = GW18, 
                                    "WTC11" = WTC11, 
                                    "Cicero" = Cicero, 
                                    "Trevino" = Trevino, 
                                    "Ziffra" = Ziffra))
}

# Merge score table with list of gene enhancer pairs
gene_enhancer_pairs_final <- cbind(gene_enhancer_pairs, score_table)
write.table(gene_enhancer_pairs_final, file = "./gene_enhancer_pairs_by_dataset.txt", # save the dataset to skip rerunning code
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
gene_enhancer_pairs_final <- read.table(file = "./gene_enhancer_pairs_by_dataset.txt", header = TRUE, sep = "\t")
head(gene_enhancer_pairs_final)
dim(gene_enhancer_pairs_final)

# Replace negative values for Trevino dataset with NA
gene_enhancer_pairs_final$Trevino <- ifelse(gene_enhancer_pairs_final$Trevino > 0, gene_enhancer_pairs_final$Trevino, NA)

# Make the ggpairs plot for enhancer score correlations
ggpairs(gene_enhancer_pairs_final, upper = list(continuous = wrap("cor", method = "spearman")),
        # diag = list(continuous = "blankDiag"),
        columns = c("GW18", "WTC11", "Cicero", "Trevino", "Ziffra"), 
        columnLabels = c("GW18", "WTC11", "Cicero", "Trevino", "Ziffra")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



### (8) Get merged bedpe tracks for all 5 prioritization datasets ######################################################

# create merged bedpe file from all datasets
all_datasets_combined_genome_wide_hg38 <- rbind(human_18gw_pfc_abc_score_hg38,
                                                WTC11_NGN2_iPSC_abc_score_hg38,
                                                cicero_coaccessability_hg38,
                                                trevino_enhancers_hg38,
                                                ziffra_abc_score_hg38)
head(all_datasets_combined_genome_wide_hg38)
dim(all_datasets_combined_genome_wide_hg38)

# write combined bedpe file to text 
write.table(all_datasets_combined_genome_wide_hg38, file = "./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# now intersect combined datasets with candidate regions for CRISPRa
system("bedtools intersect -a ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38.bedpe -b ./bed/CRISPRa_candidate_cCREs.bed -wa > ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_CRISPRa_candidate_cCREs_only.bedpe")
system("head ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_CRISPRa_candidate_cCREs_only.bedpe")
system("wc -l ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_CRISPRa_candidate_cCREs_only.bedpe") # 15,697 candidate enhancers


# now filter this dataset for only our validated CRISPRa gene-enhancer contacts
all_datasets_combined_genome_wide_hg38 <- rbind(human_18gw_pfc_abc_score_hg38_NDD,
                                                WTC11_NGN2_iPSC_abc_score_hg38_NDD,
                                                cicero_coaccessability_hg38_NDD,
                                                trevino_enhancers_hg38_NDD,
                                                ziffra_abc_score_hg38_NDD)
head(all_datasets_combined_genome_wide_hg38_NDD)
dim(all_datasets_combined_genome_wide_hg38_NDD)

# write the filtered bedpe dataset to a text file
write.table(all_datasets_combined_genome_wide_hg38_NDD[,c(1:10)], file = "./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_NDD.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# filter for MPRA and CRISPRa validated pairs gene-enhancer pairs
system("bedtools intersect -a ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_NDD.bedpe -b ./MPRA_CRISPRa_Validated_Enhancers_hg38.bed -wa > ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_NDD_CRISPRa_validated_cCREs_only.bedpe")
system("head ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_NDD_CRISPRa_validated_cCREs_only.bedpe")
system("wc -l ./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_NDD_CRISPRa_validated_cCREs_only.bedpe") # 256 candidate enhancers

all_datasets_combined_genome_wide_hg38_NDD_hits <- read.table(file = "./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_NDD_CRISPRa_validated_cCREs_only.bedpe", header = TRUE, sep = "\t",
                                                              col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2"))

# get list genes with validated CRISPRa activity
CRISPRa_hit_genes <- data.frame(read_excel("./MPRA_CRISPRa_Validated_Enhancers.xlsx"))$Target_Gene
all_datasets_combined_genome_wide_hg38_NDD_hits <- all_datasets_combined_genome_wide_hg38_NDD_hits[which(sapply(strsplit(all_datasets_combined_genome_wide_hg38_NDD_hits$name, "|", fixed = TRUE), `[`, 1) %in% CRISPRa_hit_genes),]
dim(all_datasets_combined_genome_wide_hg38_NDD_hits)
dim(unique(all_datasets_combined_genome_wide_hg38_NDD_hits[,c(1:3)]))

# write the filtered bedpe dataset to a text file
write.table(all_datasets_combined_genome_wide_hg38_NDD_hits, file = "./bed/bedpe/hg38/all_datasets_combined_genome_wide_hg38_NDD_CRISPRa_validated_cCREs_contacts_only.bedpe",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



###################################################### END SCRIPT ######################################################




