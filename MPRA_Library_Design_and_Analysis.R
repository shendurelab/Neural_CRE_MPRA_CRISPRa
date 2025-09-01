########################################################################################################################
#
#   03_MPRA_Library_Design_and_Analysis.R
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
library(data.table)
library(phylotools)
library(tidyverse)
library(homerkit)
library(ggrepel)
library(htmltab)
library(DESeq2)
library(UpSetR)
library(sjmisc)
library(readxl)
library(GGally)
library(Cairo)
library(xlsx)
library(grid)
library(ape)

setwd('~/Dropbox/Neurohub_CRISPRa_CRT/')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Get candidate NDD enhancers for each dataset
#   (2) Add positive and negative control sequences from Gerkis et al 2021 (Neuron) MPRA of fetal brain HARs
#   (3) Convert final candidate enhancer list to overlapping 270bp regions for an MPRA
#   (4) Convert MPRA region bed file into fasta sequences and add 15bp adaptors to each end
#   (5) Add positive and negative control regions from Inoue et al 2019 (Cell Stem Cell)
#   (6) Convert MPRA region bed file into fasta sequences and add 15bp adaptors to each end
#   (7) Load results from MPRA and combine into a single dataframe
#   (8) Identify active oligos by comparing DNA/RNA barcode counts with DESeq2
#   (9) Create combined MPRA metadata results table
#   (10) Make figures for replicate correlations and sequence class activity
#   (11) Get active and inactive MPRA oligos and perform Homer motif analysis
#   (12) Perform 90bp sub region analysis of tiling oligos
#   (13) Overlap MPRA results with human brain chip-seq experiments
#   (14) Comparison of MPRA results to Vista enhancer atlas 
#   (15) Select final candidate enhancer regions for CRISPRa screen
#
########################################################################################################################

### (1) Get candidate NDD enhancers for each dataset ###################################################################

candidate_enhancer_regions_for_mpra <- read.table(file = "./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed", header = FALSE, sep = "\t",
                                                  col.names = c("chrom", "chromStart", "chromEnd", "name"))

head(candidate_enhancer_regions_for_mpra)
dim(candidate_enhancer_regions_for_mpra) # 5,425 total candidate enhancer regions

# Add enhancer boundaries to final candidate enhancer list
candidate_enhancer_regions_for_mpra$name <- paste0(candidate_enhancer_regions_for_mpra$chrom, ":", candidate_enhancer_regions_for_mpra$chromStart, "-", candidate_enhancer_regions_for_mpra$chromEnd, "|", candidate_enhancer_regions_for_mpra$name)
head(candidate_enhancer_regions_for_mpra)



### (2) Add positive and negative control sequences from Gerkis et al 2021 (Neuron) MPRA of fetal brain HARs ###########

# Get results from Gerkis et al 2021 (Neuron) MPRA of fetal brain HARs
gerkis_HARs <- data.frame(read_xlsx(path = "./data/published_datasets/1-s2.0-S0896627321005808-mmc7.xlsx", sheet = "HARDatabase", skip = 1))

# Get high confidence enhancers HARs that show MPRA activity and fetal brain enhancer signatures
gerkis_HARs_top <- gerkis_HARs[which(gerkis_HARs$caMPRA.enhancer.activity == "YES" &
                                       gerkis_HARs$fetal.brain.DHS..fbDNA. == "YES" &
                                       gerkis_HARs$Neuronal.ChIP.seq.enhancer.signal..nChIPSignal. == "YES" &
                                       gerkis_HARs$caMPRA...fbDHS == "YES" &
                                       gerkis_HARs$caMPRA...fbDHS...nChIPSignal == "YES"),]

# Get negative control HARs that lack MPRA activity and fetal brain enhancer signatures
gerkis_HARs_bottom <- gerkis_HARs[which(gerkis_HARs$caMPRA.enhancer.activity == "NO" &
                                          gerkis_HARs$fetal.brain.DHS..fbDNA. == "NO" &
                                          gerkis_HARs$Neuronal.ChIP.seq.enhancer.signal..nChIPSignal. == "NO" &
                                          gerkis_HARs$caMPRA...fbDHS == "NO" &
                                          gerkis_HARs$caMPRA...fbDHS...nChIPSignal == "NO"),]

# Write positive and negative control HAR sequences to a bed file
write.table(gerkis_HARs_top[,c(2:4,1)], file = "./bed/gerkis_HARs_top.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(gerkis_HARs_bottom[,c(2:4,1)], file = "./bed/gerkis_HARs_bottom.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Remove positive and negative control sequences that overlap with one of the candidate enhancers to be tested
system("bedtools intersect -a ./bed/gerkis_HARs_top.bed -b ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed -wa -v > ./bed/gerkis_HARs_top.overlap_removed.bed")
system("bedtools intersect -a ./bed/gerkis_HARs_bottom.bed -b ./bed/candidate_enhancer_regions_for_CRISPRa_CRT_genes.bed -wa -v > ./bed/gerkis_HARs_bottom.overlap_removed.bed")

# Select 100 random HARs as positive control regions
gerkis_HARs_top_no_overlap <- read.table(file = "./bed/gerkis_HARs_top.overlap_removed.bed", header = FALSE, sep = "\t",
                                         col.names = c("chrom", "chromStart", "chromEnd", "name"))
set.seed(0) # Always remember to set seed!!
gerkis_HARs_positive_controls <- gerkis_HARs_top_no_overlap[which(gerkis_HARs_top_no_overlap$name %in% sample(gerkis_HARs_top_no_overlap$name, 150)),]
gerkis_HARs_positive_controls$name <- paste0(gerkis_HARs_positive_controls$name, ";positive_control")

# Select 100 random HARs as negative control regions
gerkis_HARs_bottom_no_overlap <- read.table(file = "./bed/gerkis_HARs_bottom.overlap_removed.bed", header = FALSE, sep = "\t",
                                            col.names = c("chrom", "chromStart", "chromEnd", "name"))
set.seed(0) # Always remember to set seed!!
gerkis_HARs_negative_controls <- gerkis_HARs_bottom_no_overlap[which(gerkis_HARs_bottom_no_overlap$name %in% sample(gerkis_HARs_bottom_no_overlap$name, 150)),]
gerkis_HARs_negative_controls$name <- paste0(gerkis_HARs_negative_controls$name, ";negative_control")

# Add genome coordinates to HAR positive control names
gerkis_HARs_positive_controls$name <- paste0(gerkis_HARs_positive_controls$chrom, ":", gerkis_HARs_positive_controls$chromStart, "-", gerkis_HARs_positive_controls$chromEnd, "|", gerkis_HARs_positive_controls$name)
head(gerkis_HARs_positive_controls)

# Add genome coordinates to HAR negative control names
gerkis_HARs_negative_controls$name <- paste0(gerkis_HARs_negative_controls$chrom, ":", gerkis_HARs_negative_controls$chromStart, "-", gerkis_HARs_negative_controls$chromEnd, "|", gerkis_HARs_negative_controls$name)
head(gerkis_HARs_negative_controls)

# Add positive and negative control regions to the final candidate enhancer list
final_candidate_enhancer_list <- rbind(candidate_enhancer_regions_for_mpra,
                                       gerkis_HARs_positive_controls,
                                       gerkis_HARs_negative_controls)

final_positive_negative_control_list <- rbind(gerkis_HARs_positive_controls,
                                              gerkis_HARs_negative_controls)

write.table(final_positive_negative_control_list, file = "./bed/controls/gerkis_positive_and_negative_control_HARs.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


dim(final_candidate_enhancer_list) # 5,725 total enhancer regions to test

# Remove temporary files
system("rm ./bed/gerkis_HARs_top.bed")
system("rm ./bed/gerkis_HARs_bottom.bed")
system("rm ./bed/gerkis_HARs_top.overlap_removed.bed")
system("rm ./bed/gerkis_HARs_bottom.overlap_removed.bed")



### (4) Convert final candidate enhancer list to overlapping 270bp regions for an MPRA #################################

# Get the length of each candidate enhancer region
final_candidate_enhancer_list$length <- final_candidate_enhancer_list$chromEnd - final_candidate_enhancer_list$chromStart

# Get 270bp MPRA candidate regions
mpra_candidates <- data.frame("chrom" = c(), "chromStart" = c(), "chromEnd" = c(), "name" = c())
for (i in 1:dim(final_candidate_enhancer_list)[1]){
  num_regions <- ceiling(final_candidate_enhancer_list$length[i] / 270)
  for (j in 1:num_regions) {
    new_candidates <- data.frame("chrom" = final_candidate_enhancer_list$chrom[i],
                                 "chromStart" = final_candidate_enhancer_list$chromStart[i] + (j-1)*270,
                                 "chromEnd" = final_candidate_enhancer_list$chromStart[i] + (j)*270,
                                 "name" = final_candidate_enhancer_list$name[i])
    colnames(new_candidates)[4] <- "name"
    mpra_candidates <- rbind(mpra_candidates, new_candidates)
    
    new_candidates <- data.frame("chrom" = final_candidate_enhancer_list$chrom[i],
                                 "chromStart" = final_candidate_enhancer_list$chromStart[i] + (j-1)*270 - 90,
                                 "chromEnd" = final_candidate_enhancer_list$chromStart[i] + (j)*270 - 90,
                                 "name" = final_candidate_enhancer_list$name[i])
    colnames(new_candidates)[4] <- "name"
    mpra_candidates <- rbind(mpra_candidates, new_candidates)
    
    new_candidates <- data.frame("chrom" = final_candidate_enhancer_list$chrom[i],
                                 "chromStart" = final_candidate_enhancer_list$chromStart[i] + (j-1)*270 + 90,
                                 "chromEnd" = final_candidate_enhancer_list$chromStart[i] + (j)*270 + 90,
                                 "name" = final_candidate_enhancer_list$name[i])
    colnames(new_candidates)[4] <- "name"
    mpra_candidates <- rbind(mpra_candidates, new_candidates)
  }
}

head(mpra_candidates)
dim(mpra_candidates) # 45,642 total MPRA regions



### (5) Add positive and negative control regions from Inoue et al 2019 (Cell Stem Cell) ###############################

inoue_mpra_results <- data.frame(read_xlsx(path = "./data/published_datasets/1-s2.0-S1934590919304217-mmc4.xlsx", sheet = "Designed_sequences"))

# Order MPRA results from largest to smallest and remove scrambled sequences
inoue_mpra_results <- inoue_mpra_results[order(inoue_mpra_results$ratio.72hr, decreasing = TRUE),]
inoue_mpra_results <- inoue_mpra_results[which(inoue_mpra_results$Criteria.for.choosing != "Scrambles"),]

# Get the top and bottom 100 sequences from Inoue et al's MPRA as positive and negative controls respectively
inoue_mpra_controls <- rbind(data.frame("chrom" = sapply(strsplit(inoue_mpra_results$Sequence.coordinates, ":"), `[`, 1)[1:100],
                                        "chromStart" = sapply(strsplit(sapply(strsplit(inoue_mpra_results$Sequence.coordinates, ":"), `[`, 2), "-"), `[`, 1)[1:100],
                                        "chromEnd" = sapply(strsplit(sapply(strsplit(inoue_mpra_results$Sequence.coordinates, ":"), `[`, 2), "-"), `[`, 2)[1:100],
                                        "name" = "inoue;positive_control"),
                             data.frame("chrom" = sapply(strsplit(inoue_mpra_results$Sequence.coordinates, ":"), `[`, 1)[(dim(inoue_mpra_results)[1]-99):dim(inoue_mpra_results)[1]],
                                        "chromStart" = sapply(strsplit(sapply(strsplit(inoue_mpra_results$Sequence.coordinates, ":"), `[`, 2), "-"), `[`, 1)[(dim(inoue_mpra_results)[1]-99):dim(inoue_mpra_results)[1]],
                                        "chromEnd" = sapply(strsplit(sapply(strsplit(inoue_mpra_results$Sequence.coordinates, ":"), `[`, 2), "-"), `[`, 2)[(dim(inoue_mpra_results)[1]-99):dim(inoue_mpra_results)[1]],
                                        "name" = "inoue;negative_control"))

# Write Inoue control sequences to a bed file and then lift over coordinates to hg38
write.table(inoue_mpra_controls, file = "./bed/controls/inoue_positive_and_negative_controls.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
system("./scripts/UCSC/liftover ./bed/controls/inoue_positive_and_negative_controls.bed ./scripts/UCSC/hg19ToHg38.over.chain ./bed/controls/inoue_positive_and_negative_controls_liftover_hg38.bed ./bed/controls/inoue_positive_and_negative_controls.unmapped.bed")
system("rm ./bed/controls/inoue_positive_and_negative_controls.unmapped.bed")

# Get hg38 Inoue control sequences and extend every sequence to 270bp
inoue_mpra_controls_hg38 <- read.table(file = "./bed/controls/inoue_positive_and_negative_controls_liftover_hg38.bed", header = FALSE, sep = "\t",
                                       col.names = c("chrom", "chromStart", "chromEnd", "name"))

# Extend all Inoue et al MPRA regions to at least 270bp
for (i in 1:dim(inoue_mpra_controls_hg38)[1]){
  if (inoue_mpra_controls_hg38$chromEnd[i] - inoue_mpra_controls_hg38$chromStart[i] < 270){
    length <- 270 - (inoue_mpra_controls_hg38$chromEnd[i] - inoue_mpra_controls_hg38$chromStart[i])
    inoue_mpra_controls_hg38$chromStart[i] <- inoue_mpra_controls_hg38$chromStart[i] - floor(length/2)
    inoue_mpra_controls_hg38$chromEnd[i] <- inoue_mpra_controls_hg38$chromEnd[i] + ceiling(length/2)
  }
}

# add coordinates to the name for inoue positive controls
inoue_mpra_controls_hg38$name <- paste0(inoue_mpra_controls_hg38$chrom, ":", inoue_mpra_controls_hg38$chromStart, "-", inoue_mpra_controls_hg38$chromEnd, "|", inoue_mpra_controls_hg38$name)

# Add Inoue et al controls to the final list of MPRA candidate regions
mpra_candidates <- rbind(mpra_candidates,
                         inoue_mpra_controls_hg38)

# Write bed file of final MPRA candidate regions
write.table(mpra_candidates, file = "./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Get bigBed files for enhancer and MPRA candidate regions
system("bedtools sort -i ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed > ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.sorted.bed")
system("rm ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed")
system("mv ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.sorted.bed ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed")
system("./scripts/UCSC/bedToBigBed ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed ~/genomes/hg38/hg38.chrom.sizes ./UCSC_browser_tracks/hg38/bigBed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bigBed")
system("wc -l ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed") # 45,842 final regions for an MPRA



### (6) Convert MPRA region bed file into fasta sequences and add 15bp adaptors to each end ############################

# Get fasta sequences for MPRA regions in bed file
system("bedtools getfasta -fi ~/genomes/hg38/fasta/hg38.fasta -bed ./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed -name > ./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.fasta")

# Add 500 scrambled controls to MPRA library 
# system("python3 ./scripts/python/get_scrambled_controls.py --fasta ./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.fasta --num 500")

# Add adaptor sequences to the MPRA library
system("python3 ./scripts/python/add_adaptors_to_fasta.py --fasta ./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.controls_added.fasta --fivePrimeAdaptor AGGACCGGATCAACT --threePrimeAdaptor CATTGCGTGAACCGA")

# Double check for SceI restriction enzyme sites
system("grep TAGGGATAACAGGGTAAT ./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.controls_added.adaptors_added.fasta")
system("grep ATTACCCTGTTATCCCTA ./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.controls_added.adaptors_added.fasta")

system("head ./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.controls_added.adaptors_added.fasta")
system("wc -l ./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.controls_added.adaptors_added.fasta")
# 46,342 total sequences tested via MPRA



### (7) Load results from MPRA and combine into a single dataframe #####################################################

# Load results from MPRA sequencing run 1
mpra_run1 <- fread(file = "./MPRA/results/WTC11_NGN2_iPSC_ExN_14DIV_final_labeled_counts_run1.txt",
                            select = c(1, 4:9))
colnames(mpra_run1) <- c("CRS", "DNA_S1_B1", "DNA_S2_B1", "DNA_S3_B1", "RNA_S1_B1", "RNA_S2_B1", "RNA_S3_B1")

# Load results from MPRA sequencing run 2
mpra_run2 <- fread(file = "./MPRA/results/WTC11_NGN2_iPSC_ExN_14DIV_final_labeled_counts_run2.txt",
                            select = c(1, 4:9))
colnames(mpra_run2) <- c("CRS", "DNA_S1_B2", "DNA_S2_B2", "DNA_S3_B2", "RNA_S1_B2", "RNA_S2_B2", "RNA_S3_B2")

# Replace NA values with 0
mpra_run1[is.na(mpra_run1)] <- 0
mpra_run2[is.na(mpra_run2)] <- 0

# Remove CRSs with no barcode
mpra_run1 <- mpra_run1[which(mpra_run1$CRS != "no_BC"),]
mpra_run2 <- mpra_run2[which(mpra_run2$CRS != "no_BC"),]

# get the number of DNA counts per barcode
get_barcode_count <- rbind(mpra_run1, mpra_run2, use.names = FALSE)
dim(get_barcode_count) # 15,040,363 barcodes
DNA_counts_per_barcode <- data.frame("count" = rowSums(get_barcode_count[,c(2:4)]))
summary(DNA_counts_per_barcode$count)

ggplot(DNA_counts_per_barcode, aes(x = count)) + 
  geom_histogram(binwidth = 1, color = "#56B4E9", fill = "#56B4E9") +
  scale_y_continuous(breaks = c(0 , 2000000, 4000000), labels = c("0", "2e+6", "4e+6")) +
  xlab("# counts") + ylab("# DNA barcodes") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))

# get the number of RNA counts per barcode
RNA_counts_per_barcode <- data.frame("count" = rowSums(get_barcode_count[,c(5:7)]))
summary(RNA_counts_per_barcode$count)

ggplot(RNA_counts_per_barcode, aes(x = count)) + 
  geom_histogram(binwidth = 1, color = "#56B4E9", fill = "#56B4E9") +
  scale_y_continuous(breaks = c(0 , 100000, 300000), labels = c("0", "1e+5", "3e+5")) +
  xlab("# counts") + ylab("# RNA barcodes") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))

# Get the number of DNA barcodes per CRS
barcodes_per_CRS <- get_barcode_count %>% dplyr::count(CRS) 

ggplot(barcodes_per_CRS, aes(x = n)) + 
  geom_histogram(binwidth = 1, color = "#56B4E9", fill = "#56B4E9") +
  xlab("# barcodes") + ylab("# CRSs") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y = element_text(size = 14))

# merge barcodes to get CRS levels results for sequencing run 1
summary_run1 <- mpra_run1 %>% group_by(CRS) %>% 
  summarise_if(is.numeric, sum)
summary_run1 <- data.frame(summary_run1)

# create entry for CRS that is missing in sequencing batch 1
summary_run1 <- rbind(summary_run1,
                      data.frame("CRS" = "GW18_PFC_ABC::chr3:125141315-125141585",
                                 "DNA_S1_B1" = 0,
                                 "DNA_S2_B1" = 0,
                                 "DNA_S3_B1" = 0,
                                 "RNA_S1_B1" = 0, 
                                 "RNA_S2_B1" = 0,
                                 "RNA_S3_B1" = 0))

# merge barcodes to get CRS levels results for sequencing run 2
summary_run2 <- mpra_run2 %>% group_by(CRS) %>% 
  summarise_if(is.numeric, sum)
summary_run2 <- data.frame(summary_run2)

# order results from sequencing run 1 and 2
summary_run1 <- summary_run1[order(summary_run1$CRS),]
summary_run2 <- summary_run2[order(summary_run2$CRS),]

# Double check that row names are the same for sequencing run 1 and 2
for(i in 1:dim(summary_run1)[1]){
  if(summary_run1$CRS[i] != summary_run2$CRS[i]){
    print(i)
    break
  }
} # looks okay

# merge the sequencing runs into a single dataframe
mpra_results <- cbind(summary_run1, summary_run2[,2:7])

# add columns for the DNA and RNA count totals
mpra_results$DNA_S1_Total <- mpra_results$DNA_S1_B1 + mpra_results$DNA_S1_B2
mpra_results$DNA_S2_Total <- mpra_results$DNA_S2_B1 + mpra_results$DNA_S2_B2
mpra_results$DNA_S3_Total <- mpra_results$DNA_S3_B1 + mpra_results$DNA_S3_B2
mpra_results$RNA_S1_Total <- mpra_results$RNA_S1_B1 + mpra_results$RNA_S1_B2
mpra_results$RNA_S2_Total <- mpra_results$RNA_S2_B1 + mpra_results$RNA_S2_B2
mpra_results$RNA_S3_Total <- mpra_results$RNA_S3_B1 + mpra_results$RNA_S3_B2

# add columns for the average number of DNA and RNA counts
mpra_results$DNA_Average <- rowSums(mpra_results[,c(14:16)]) / 3
mpra_results$RNA_Average <- rowSums(mpra_results[,c(17:19)]) / 3

# calculate the log2 count ratios for each replicate
mpra_results$Rep1_Log2_Count_Ratio <- log2((mpra_results$RNA_S1_Total + 1) / (mpra_results$DNA_S1_Total + 1))
mpra_results$Rep2_Log2_Count_Ratio <- log2((mpra_results$RNA_S2_Total + 1) / (mpra_results$DNA_S2_Total + 1))
mpra_results$Rep3_Log2_Count_Ratio <- log2((mpra_results$RNA_S3_Total + 1) / (mpra_results$DNA_S3_Total + 1))

# get the average log2 count ratio for each CRS and order the results table
mpra_results$Average_Log2_Count_Ratio <- rowSums(mpra_results[,c(22:24)]) / 3
mpra_results <- mpra_results[order(mpra_results$Average_Log2_Count_Ratio, decreasing = TRUE),]



### (8) Identify active oligos by comparing DNA/RNA barcode counts with DESeq2 #########################################

# filter MPRA results table to supply as input for DESeq2
mpra_result_for_deseq2 <- as.matrix(mpra_results[,c(2:13)])
rownames(mpra_result_for_deseq2) <- mpra_results$CRS

# filter MPRA results table for CRSs with a DNA barcode count >= 5
mpra_result_for_deseq2 <- mpra_result_for_deseq2[which(rowSums(mpra_result_for_deseq2[,c(1:3, 7:9)]) >= 5),]

# convert rownames to be compatible with DESeq2 naming conventions
rownames_conversion <- data.frame("old_names" = rownames(mpra_result_for_deseq2),
                                  "new_names" = sapply(strsplit(rownames(mpra_result_for_deseq2), "::"), `[`, 2))
rownames(mpra_result_for_deseq2) <- sapply(strsplit(rownames(mpra_result_for_deseq2), "::"), `[`, 2)
rownames_conversion$new_names <- paste0("crs", 1:dim(mpra_result_for_deseq2)[1], "_", rownames(mpra_result_for_deseq2))
rownames(mpra_result_for_deseq2) <- paste0("crs", 1:dim(mpra_result_for_deseq2)[1], "_", rownames(mpra_result_for_deseq2))

# annotate DNA vs RNA fraction and sequencing batch
fraction <- data.frame("DNA_or_RNA" = c("DNA", "DNA", "DNA", "RNA", "RNA", "RNA",
                                        "DNA", "DNA", "DNA", "RNA", "RNA", "RNA"),
                       "Batch" = c("B1", "B1", "B1", "B1", "B1", "B1",
                                   "B2", "B2", "B2", "B2", "B2", "B2"))
rownames(fraction) <- colnames(mpra_result_for_deseq2)
fraction$DNA_or_RNA <- factor(fraction$DNA_or_RNA, levels = c("DNA", "RNA"))
fraction$Batch <- factor(fraction$Batch)

# run DESeq2 with the model Counts ~ Batch + DNA_or_RNA
dds <- DESeqDataSetFromMatrix(countData = mpra_result_for_deseq2,
                              colData = fraction,
                              design = ~ Batch + DNA_or_RNA)

# run DESeq2 and get the results
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj, decreasing = FALSE),]
res <- data.frame(res)

# determine which oligos are active
res$is_active <- ifelse(res$log2FoldChange > 0 & res$padj < 0.01, "Active", "Inactive")

# plot the ratio of RNA to DNA counts
ggplot(data = res, aes(x = log2FoldChange, y = -log10(pvalue), color = is_active)) +
  geom_point() +
  scale_color_manual(values = c("#56B4E9", "#AAAAAA")) + 
  geom_hline(yintercept = -log10(0.01), color = "red", linetype = "dashed") + 
  ggtitle("DESeq2 RNA vs DNA counts") + xlab("log2FoldChange") + ylab("-log10(p-value)") +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# covert row names back to the original annotations
rownames(res) <- rownames_conversion$old_names[match(rownames(res), rownames_conversion$new_names)]

# match the results from DESeq2 the the combined mpra results table
DESeq2_Log2FoldChange_List <- c()
DESeq2_P_Value_List <- c()
DESeq2_FDR_List <- c()
for(i in 1:dim(mpra_results)[1]){
  DESeq2_Log2FoldChange_List[i] <- ifelse(mpra_results$CRS[i] %in% rownames(res),
                                          res$log2FoldChange[match(mpra_results$CRS[i], rownames(res))],
                                          NA)
  DESeq2_P_Value_List[i] <- ifelse(mpra_results$CRS[i] %in% rownames(res),
                                   res$pvalue[match(mpra_results$CRS[i], rownames(res))],
                                   NA)
  DESeq2_FDR_List[i] <- ifelse(mpra_results$CRS[i] %in% rownames(res),
                               res$padj[match(mpra_results$CRS[i], rownames(res))],
                               NA)
  
}

# add DESeq2 results to the MPRA results table
mpra_results$DESeq2_Log2FoldChange <- DESeq2_Log2FoldChange_List
mpra_results$DESeq2_P_Value <- DESeq2_P_Value_List
mpra_results$DESeq2_FDR <- DESeq2_FDR_List

# define active CRSs as those with an adjusted pvalue < 0.01 and a log2FC > 0
mpra_results$Oligo_Is_Active <- ifelse(mpra_results$DESeq2_Log2FoldChange > 0 & mpra_results$DESeq2_FDR < 0.01, TRUE, FALSE)



### (9) Create combined MPRA metadata results table ####################################################################

# read in MPRA library design fasta file
mpra_library = read.fasta("./MPRA/sequences/MPRA_sequences_for_CRISPRa_CRT_candidate_enhancers.controls_added.adaptors_added_from_Twist.fasta")

# create columns for 5' and 3' library adaptor sequences
mpra_library$Five_Primer_Adaptor <- substr(mpra_library$seq.text, 1, 15)
mpra_library$Three_Primer_Adaptor <- substr(mpra_library$seq.text, 286, 300)
mpra_library$Oligo_Sequence <- mpra_library$seq.text

# remove adaptors from test sequence and rename columns
mpra_library$seq.name <- gsub(",", ";", mpra_library$seq.name)
mpra_library$seq.text <- substr(mpra_library$seq.text, 16, 285)
colnames(mpra_library)[1:2] <- c("CRS", "Test_Sequence")

# add CRS categories to the mpra library design dataframe
categories <- c()
for(i in 1:dim(mpra_library)[1]){
  if(str_contains(mpra_library$CRS[i], "negative") & str_contains(mpra_library$CRS[i], "HAR")){
    categories <- c(categories, "HAR_Negative")
    next
  }
  if(str_contains(mpra_library$CRS[i], "positive") & str_contains(mpra_library$CRS[i], "HAR")){
    categories <- c(categories, "HAR_Positive")
    next
  }
  if(str_contains(mpra_library$CRS[i], "negative") & str_contains(mpra_library$CRS[i], "inoue")){
    categories <- c(categories, "Inoue_Negative")
    next
  }
  if(str_contains(mpra_library$CRS[i], "positive") & str_contains(mpra_library$CRS[i], "inoue")){
    categories <- c(categories, "Inoue_Positive")
    next
  }
  if(str_contains(mpra_library$CRS[i], "scrambled_control")){
    categories <- c(categories, "Scrambled")
    next
  }
  categories <- c(categories, "Test")
}

# add CRS categories to the MPRA library design dataframe
mpra_library$Sequence_Class <- categories

# create column for oligo coordinates in the MPRA library design dataframe
mpra_library$Oligo_Coordinates_hg38 <- sapply(strsplit(mpra_library$CRS, "::"), `[`, 2)

# add annotation metadata to the library design table
annotations <- read.table("./bed/MPRA_regions_for_CRISPRa_CRT_candidate_enhancers.bed",
                          header = FALSE, sep = "\t", col.names = c("chrom", "chromStartOligo", "chromEndOligo", "name"))

# get the CRS coordinates for each candidate enhancer tested
annotations$CRS_Coordinates_hg38 <- sapply(strsplit(annotations$name, "|", fixed = TRUE), `[`, 1)
length(unique(annotations$CRS_Coordinates_hg38))

# get the annotations by dataset and gene name
annotations$annotation <- sapply(strsplit(annotations$name, "|", fixed = TRUE), `[`, 2)

datasets <- data.frame("GW18_PFC_ABC" = c(),
                       "NGN2_iPSC_ABC" = c(),
                       "Midfetal_Cortex_Ziffra_ABC" = c(),
                       "Fetal_Cerebrum_Cicero" = c(),
                       "Midfetal_Cortex_Trevino" = c())

for(i in 1:dim(annotations)[1]){
  which_datasets <- data.frame("GW18_PFC_ABC" = c(0),
                               "NGN2_iPSC_ABC" = c(0),
                               "Midfetal_Cortex_Ziffra_ABC" = c(0),
                               "Fetal_Cerebrum_Cicero" = c(0),
                               "Midfetal_Cortex_Trevino" = c(0))
  
  if(str_contains(annotations$annotation[i], "GW18_PFC_ABC")){
    which_datasets$GW18_PFC_ABC <- 1
  }
  
  if(str_contains(annotations$annotation[i], "NGN2_iPSC_ABC")){
    which_datasets$NGN2_iPSC_ABC <- 1
  }
  
  if(str_contains(annotations$annotation[i], "Midfetal_Cortex_Ziffra_ABC")){
    which_datasets$Midfetal_Cortex_Ziffra_ABC <- 1
  }
  
  if(str_contains(annotations$annotation[i], "Fetal_Cerebrum_Cicero")){
    which_datasets$Fetal_Cerebrum_Cicero <- 1
  }
  
  if(str_contains(annotations$annotation[i], "Midfetal_Cortex_Trevino")){
    which_datasets$Midfetal_Cortex_Trevino <- 1
  }
  
  datasets <- rbind(datasets, which_datasets)
  
}

# add dataset list to the annotation dataframe
annotations <- cbind(annotations, datasets)

# add genes from GW18_PFC dataset to the metadata table
GW18_PFC_ABC_gene_list <- c()
for(i in 1:dim(annotations)[1]){
  if(annotations$GW18_PFC_ABC[i] == 1){
    all_annotations <- strsplit(annotations$annotation[i], ",")[[1]]
    select_annotations <- all_annotations[grep("GW18_PFC_ABC", all_annotations, fixed = TRUE)]
    gene_list <- unique(sapply(strsplit(select_annotations, ".", fixed = TRUE), `[`, 1))
    gene_list_final <- paste0(gene_list, collapse = ",")
    GW18_PFC_ABC_gene_list <- c(GW18_PFC_ABC_gene_list, gene_list_final)
    next
  }
  GW18_PFC_ABC_gene_list <- c(GW18_PFC_ABC_gene_list, "NA")
}

annotations$GW18_PFC_ABC_Gene_List <- GW18_PFC_ABC_gene_list

# add genes from NGN2_iPSC dataset to the metadata table
NGN2_iPSC_ABC_gene_list <- c()
for(i in 1:dim(annotations)[1]){
  if(annotations$NGN2_iPSC_ABC[i] == 1){
    all_annotations <- strsplit(annotations$annotation[i], ",")[[1]]
    select_annotations <- all_annotations[grep("NGN2_iPSC_ABC", all_annotations, fixed = TRUE)]
    gene_list <- unique(sapply(strsplit(select_annotations, ".", fixed = TRUE), `[`, 1))
    gene_list_final <- paste0(gene_list, collapse = ",")
    NGN2_iPSC_ABC_gene_list <- c(NGN2_iPSC_ABC_gene_list, gene_list_final)
    next
  }
  NGN2_iPSC_ABC_gene_list <- c(NGN2_iPSC_ABC_gene_list, "NA")
}

annotations$NGN2_iPSC_ABC_Gene_List <- NGN2_iPSC_ABC_gene_list

# add genes from Ziffra dataset to the metadata table
Midfetal_Cortex_Ziffra_ABC_gene_list <- c()
for(i in 1:dim(annotations)[1]){
  if(annotations$Midfetal_Cortex_Ziffra_ABC[i] == 1){
    all_annotations <- strsplit(annotations$annotation[i], ",")[[1]]
    select_annotations <- all_annotations[grep("Midfetal_Cortex_Ziffra_ABC", all_annotations, fixed = TRUE)]
    gene_list <- unique(sapply(strsplit(select_annotations, ".", fixed = TRUE), `[`, 1))
    gene_list_final <- paste0(gene_list, collapse = ",")
    Midfetal_Cortex_Ziffra_ABC_gene_list <- c(Midfetal_Cortex_Ziffra_ABC_gene_list, gene_list_final)
    next
  }
  Midfetal_Cortex_Ziffra_ABC_gene_list <- c(Midfetal_Cortex_Ziffra_ABC_gene_list, "NA")
}

annotations$Midfetal_Cortex_Ziffra_ABC_Gene_List <- Midfetal_Cortex_Ziffra_ABC_gene_list

# add genes from Cicero datasets to the metadata table
Fetal_Cerebrum_Cicero_gene_list <- c()
for(i in 1:dim(annotations)[1]){
  if(annotations$Fetal_Cerebrum_Cicero[i] == 1){
    all_annotations <- strsplit(annotations$annotation[i], ",")[[1]]
    select_annotations <- all_annotations[grep("Fetal_Cerebrum_Cicero", all_annotations, fixed = TRUE)]
    gene_list <- unique(sapply(strsplit(select_annotations, ".", fixed = TRUE), `[`, 1))
    gene_list_final <- paste0(gene_list, collapse = ",")
    Fetal_Cerebrum_Cicero_gene_list <- c(Fetal_Cerebrum_Cicero_gene_list, gene_list_final)
    next
  }
  Fetal_Cerebrum_Cicero_gene_list <- c(Fetal_Cerebrum_Cicero_gene_list, "NA")
}

annotations$Fetal_Cerebrum_Cicero_Gene_List <- Fetal_Cerebrum_Cicero_gene_list

### add genes from Trevino datasets to the metadata table
Midfetal_Cortex_Trevino_gene_list <- c()
for(i in 1:dim(annotations)[1]){
  if(annotations$Midfetal_Cortex_Trevino[i] == 1){
    all_annotations <- strsplit(annotations$annotation[i], ",")[[1]]
    select_annotations <- all_annotations[grep("Midfetal_Cortex_Trevino", all_annotations, fixed = TRUE)]
    gene_list <- unique(sapply(strsplit(select_annotations, ".", fixed = TRUE), `[`, 1))
    gene_list_final <- paste0(gene_list, collapse = ",")
    Midfetal_Cortex_Trevino_gene_list <- c(Midfetal_Cortex_Trevino_gene_list, gene_list_final)
    next
  }
  Midfetal_Cortex_Trevino_gene_list <- c(Midfetal_Cortex_Trevino_gene_list, "NA")
}

annotations$Midfetal_Cortex_Trevino_Gene_List <- Midfetal_Cortex_Trevino_gene_list

# create a combined gene list across all annotation datasets
combined_gene_list <- c()
for(i in 1:dim(annotations)[1]){
  gene_list <- unique(unlist(list(annotations[i, 12:16])))
  gene_list <- unique(unlist(strsplit(gene_list, ",")))
  gene_list_final <- gene_list[which(gene_list != "NA")]
  gene_list_final <- paste0(gene_list_final, collapse = ",")
  combined_gene_list <- c(combined_gene_list, gene_list_final)
}

# add combined gene list to the annotation dataframe
annotations$Merged_Gene_List <- combined_gene_list
annotations$Merged_Gene_List <- ifelse(annotations$Merged_Gene_List == "", "NA", annotations$Merged_Gene_List)

# add the oligo coordinates to the combined annotation table
annotations$Oligo_Coordinates_hg38 <- paste0(annotations$chrom, ":", annotations$chromStartOligo, "-", annotations$chromEndOligo)

# trim white space
annotations$Oligo_Coordinates_hg38 <- gsub(" ", "", annotations$Oligo_Coordinates_hg38, fixed = TRUE)
annotations$CRS_Coordinates_hg38 <- gsub(" ", "", annotations$CRS_Coordinates_hg38, fixed = TRUE)

# reorder columns in the annotation dataframe
annotations <- annotations[,c(18, 5, 7:17)]

# combine annotation dataframe to oligo dataframe
mpra_library <- dplyr::left_join(mpra_library, annotations, by = "Oligo_Coordinates_hg38")

# combine MPRA library metadata with the results table
mpra_library_and_results <- dplyr::left_join(mpra_library, mpra_results, by = "CRS")

# change the name of the first column
colnames(mpra_library_and_results)[1] <- "Oligo_Name"

# get the z-score of the average log2 count ratio relative to the scrambled negative controls for each oligo
scrambled_mean <- mean(mpra_library_and_results$Average_Log2_Count_Ratio[which(mpra_library_and_results$Sequence_Class == "Scrambled")], na.rm = TRUE)
scrambled_sd <- sd(mpra_library_and_results$Average_Log2_Count_Ratio[which(mpra_library_and_results$Sequence_Class == "Scrambled")], na.rm = TRUE)
mpra_library_and_results <- add_column(mpra_library_and_results, data.frame("Average_Log2_Count_Ratio_Z_Score" = (mpra_library_and_results$Average_Log2_Count_Ratio - scrambled_mean) / scrambled_sd), .after = 43)

# get the list of CRSs that contain at least one oligo in the top 10%
quantile(mpra_library_and_results[which(mpra_library_and_results$Sequence_Class == "Test"),]$Average_Log2_Count_Ratio, probs = seq(.1, .95, by = .05), na.rm = TRUE)
mpra_library_and_results$CRS_Contains_Oligo_in_Top_10_Percent <- ifelse(mpra_library_and_results$CRS_Coordinates_hg38 %in% mpra_library_and_results$CRS_Coordinates_hg38[which(mpra_library_and_results$Average_Log2_Count_Ratio > 4.278279 & mpra_library_and_results$Sequence_Class == "Test")], TRUE, FALSE)

# double check final table
colnames(mpra_library_and_results)
dim(unique(mpra_library_and_results[,c(8, 49)]))
table(unique(mpra_library_and_results[,c(8, 49)])$CRS_Contains_Oligo_in_Top_10_Percent) # looks good

# convert NA values to write supplemental table
mpra_library_and_results_for_xlsx <- mpra_library_and_results
mpra_library_and_results_for_xlsx[is.na(mpra_library_and_results_for_xlsx)] <- "NA"

# Write the final MPRA results to an excel file as a supplementary table
write.xlsx2(mpra_library_and_results_for_xlsx, file = "./tables/Table_S3_MPRA_Design_Annotation_and_Results.xlsx", sheetName = "MPRA_Results", 
            col.names = TRUE, row.names = FALSE, append = FALSE, showNA = TRUE)



### (10) Make figures for replicate correlations and sequence class activity ###########################################

# make MPRA activity bar plot for test sequences vs scrambled negative controls
to_plot <- mpra_library_and_results[which(mpra_library_and_results$Sequence_Class %in% c("Test", "Scrambled") & is.na(mpra_library_and_results$Average_Log2_Count_Ratio) == FALSE),]
p <- ggplot(data = to_plot, aes(x = Sequence_Class, y = Average_Log2_Count_Ratio, fill = Sequence_Class)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#AAAAAA", "#56B4E9")) +
  guides(fill = "none") +
  ylab("log2(RNA/DNA)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
p


# make MPRA activity bar plot for all sequence classes
to_plot <- mpra_library_and_results[which(is.na(mpra_library_and_results$Average_Log2_Count_Ratio) == FALSE),]
to_plot$Sequence_Class <- factor(to_plot$Sequence_Class, levels = c("Scrambled", "Test", "HAR_Negative", "HAR_Positive", "Inoue_Negative", "Inoue_Positive"))
p <- ggplot(data = to_plot, aes(x = Sequence_Class, y = Average_Log2_Count_Ratio, fill = Sequence_Class)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#AAAAAA", "#56B4E9", "#E69F00", "#E69F00", "#E69F00", "#E69F00")) +
  guides(fill = "none") +
  ylab("log2(RNA/DNA)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14))
p


# get the off-center distance for every MPRA oligo relative to its CRS
mpra_library_and_results_off_center <- mpra_library_and_results[which(is.na(mpra_library_and_results$Oligo_Is_Active) == FALSE), c(7, 8, 38, 39, 48)]
mpra_library_and_results_off_center$Count_Ratio <- mpra_library_and_results_off_center$RNA_Average / mpra_library_and_results_off_center$DNA_Average
mpra_library_and_results_off_center$oligoStart <- sapply(strsplit(sapply(strsplit(mpra_library_and_results_off_center$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 1)
mpra_library_and_results_off_center$oligoEnd <- sapply(strsplit(sapply(strsplit(mpra_library_and_results_off_center$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 2)
mpra_library_and_results_off_center$oligoMidpoint <- round((as.numeric(mpra_library_and_results_off_center$oligoEnd) - as.numeric(mpra_library_and_results_off_center$oligoStart))/2) + as.numeric(mpra_library_and_results_off_center$oligoStart)
mpra_library_and_results_off_center$crsStart <- sapply(strsplit(sapply(strsplit(mpra_library_and_results_off_center$CRS_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 1)
mpra_library_and_results_off_center$crsEnd <- sapply(strsplit(sapply(strsplit(mpra_library_and_results_off_center$CRS_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 2)
mpra_library_and_results_off_center$crsMidpoint <- round((as.numeric(mpra_library_and_results_off_center$crsEnd) - as.numeric(mpra_library_and_results_off_center$crsStart))/2) + as.numeric(mpra_library_and_results_off_center$crsStart)
mpra_library_and_results_off_center$off_center_distance <- mpra_library_and_results_off_center$crsMidpoint - mpra_library_and_results_off_center$oligoMidpoint
head(mpra_library_and_results_off_center)


# plot the off-center distance for every MPRA oligo relative to its CRS
mpra_library_and_results_off_center$Oligo_Is_Active <- factor(mpra_library_and_results_off_center$Oligo_Is_Active, levels = c(TRUE, FALSE))
p <- ggplot(mpra_library_and_results_off_center, aes(x = off_center_distance, y = Count_Ratio, color = Oligo_Is_Active)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c("#56B4E9", "#AAAAAA"), labels = c("Active", "Inactive")) +
  scale_x_continuous(breaks = c(-2000, 0, 2000)) +
  xlab("Distance to peak center") + ylab("RNA/DNA Counts") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
p


# remove NA values from the MPRA results table
mpra_results_no_NA <- mpra_library_and_results[which(is.na(mpra_library_and_results$Oligo_Is_Active) == FALSE),]


# get correlation for rep1 vs rep2
pearson <- round(cor(x = mpra_results_no_NA$Rep1_Log2_Count_Ratio, y = mpra_results_no_NA$Rep2_Log2_Count_Ratio, method = "pearson"), 2)
spearman <- round(cor(x = mpra_results_no_NA$Rep1_Log2_Count_Ratio, y = mpra_results_no_NA$Rep2_Log2_Count_Ratio, method = "spearman"), 2)
correlation_text <- grobTree(textGrob(paste0("r = ", pearson, "\n", "rho = ", spearman), x = 0.05,  y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 12, fontface = "italic")))
p <- ggplot(mpra_results_no_NA, aes(x = Rep1_Log2_Count_Ratio, y = Rep2_Log2_Count_Ratio)) +
  geom_point(size = 0.1, alpha = 0.1) +
  annotation_custom(correlation_text) +
  xlab("Rep1 log2(RNA/DNA") + ylab("Rep2 log2(RNA/DNA") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
p


# get correlation for rep1 vs rep3
pearson <- round(cor(x = mpra_results_no_NA$Rep1_Log2_Count_Ratio, y = mpra_results_no_NA$Rep3_Log2_Count_Ratio, method = "pearson"), 2)
spearman <- round(cor(x = mpra_results_no_NA$Rep1_Log2_Count_Ratio, y = mpra_results_no_NA$Rep3_Log2_Count_Ratio, method = "spearman"), 2)
correlation_text <- grobTree(textGrob(paste0("r = ", pearson, "\n", "rho = ", spearman), x = 0.05,  y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 12, fontface = "italic")))
p <- ggplot(mpra_results_no_NA, aes(x = Rep1_Log2_Count_Ratio, y = Rep3_Log2_Count_Ratio)) +
  geom_point(size = 0.1, alpha = 0.1) +
  annotation_custom(correlation_text) +
  xlab("Rep1 log2(RNA/DNA") + ylab("Rep3 log2(RNA/DNA") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
p


# get correlation for rep2 vs rep3
pearson <- round(cor(x = mpra_results_no_NA$Rep2_Log2_Count_Ratio, y = mpra_results_no_NA$Rep3_Log2_Count_Ratio, method = "pearson"), 2)
spearman <- round(cor(x = mpra_results_no_NA$Rep2_Log2_Count_Ratio, y = mpra_results_no_NA$Rep3_Log2_Count_Ratio, method = "spearman"), 2)
correlation_text <- grobTree(textGrob(paste0("r = ", pearson, "\n", "rho = ", spearman), x = 0.05,  y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 12, fontface = "italic")))
p <- ggplot(mpra_results_no_NA, aes(x = Rep2_Log2_Count_Ratio, y = Rep3_Log2_Count_Ratio)) +
  geom_point(size = 0.1, alpha = 0.1) +
  annotation_custom(correlation_text) +
  xlab("Rep2 log2(RNA/DNA") + ylab("Rep3 log2(RNA/DNA") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
p



### (11) Get active and inactive MPRA oligos and perform Homer motif analysis ##########################################

# get active and MPRA oligos and write them to a bed file
active_oligos <- mpra_library_and_results[which(mpra_library_and_results$Oligo_Is_Active == TRUE & mpra_library_and_results$Sequence_Class == "Test"),]
active_oligos_bed <- data.frame("chrom" = sapply(strsplit(active_oligos$Oligo_Coordinates_hg38, ":"), `[`, 1),
                                "chromStart" = sapply(strsplit(sapply(strsplit(active_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 1),
                                "chromEnd" = sapply(strsplit(sapply(strsplit(active_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 2))
write.table(active_oligos_bed , file = "./bed/active_oligos_for_homer_motif_analysis.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# get inactive and MPRA oligos and write them to a bed file
inactive_oligos <- mpra_library_and_results[which(mpra_library_and_results$Oligo_Is_Active == FALSE & mpra_library_and_results$Sequence_Class == "Test"),]
inactive_oligos_bed <- data.frame("chrom" = sapply(strsplit(inactive_oligos$Oligo_Coordinates_hg38, ":"), `[`, 1),
                                "chromStart" = sapply(strsplit(sapply(strsplit(inactive_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 1),
                                "chromEnd" = sapply(strsplit(sapply(strsplit(inactive_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 2))
write.table(inactive_oligos_bed , file = "./bed/inactive_oligos_for_homer_motif_analysis.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# get all and MPRA oligos and write them to a bed file
all_oligos <- mpra_library_and_results[which(mpra_library_and_results$Sequence_Class == "Test"),]
all_oligos_bed <- data.frame("chrom" = sapply(strsplit(all_oligos$Oligo_Coordinates_hg38, ":"), `[`, 1),
                                  "chromStart" = sapply(strsplit(sapply(strsplit(all_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 1),
                                  "chromEnd" = sapply(strsplit(sapply(strsplit(all_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 2))
write.table(all_oligos_bed , file = "./bed/all_oligos_for_homer_motif_analysis.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# read in results from Homer with homerkit package
homer_results <- read_homer_output("./mpra_homer_motif_analysis", max_files = 200)
homer_motif_enrichment <- data.frame(homer_results$known_motif_table)
homer_motif_peaks <- data.frame(homer_results$known_motif_peaks)

# get a bed file of homer motif peaks
motif_table <- homer_motif_peaks[,c(2:4, 26)]

# adjust motif start coordinates to match MPRA oligo coordinates
motif_table$start <- motif_table$start - 1
motif_table$Oligo_Coordinates_hg38 <- paste0(motif_table$chr, ":", motif_table$start, "-", motif_table$end)

# add MPRA activity to the homer motif table
motif_table$mpra_activity <- mpra_library_and_results$Average_Log2_Count_Ratio[match(motif_table$Oligo_Coordinates_hg38, mpra_library_and_results$Oligo_Coordinates_hg38)]

# group MPRA activity results by motif
motif_results <- motif_table %>% group_by(best_guess) %>%
  summarize(count = n(),
            mean_activity = mean(as.numeric(mpra_activity), na.rm = TRUE))

# sort motif results by average MRPA activity
motif_results <- motif_results[order(motif_results$mean_activity, decreasing = TRUE),]

# add p-values for motif enrichment to the motif table
motif_results$p_value <- homer_motif_enrichment$p_value[match(motif_results$best_guess, homer_motif_enrichment$motif_name)]

# get table of abbreviated motif names for plotting
motif_names_table <- read.table(file = "./mpra_homer_motif_analysis/motif_names_table.txt", header = TRUE, sep = "\t")

# add abbreviated motif names to the motif table
motif_results$short_name <- motif_names_table$motif_name_short[match(motif_names_table$motif_name_full, motif_results$best_guess)]

# filter motif names for plotting with ggplot2
motif_results$short_name <- ifelse(-log10(motif_results$p_value) > 20 | motif_results$mean_activity > 4, motif_results$short_name, "")

# plot homer motif enrichment results with ggplot2
p <- ggplot(motif_results, aes(x = mean_activity, y = -log10(p_value), label = short_name)) + 
  geom_point() +
  ggtitle("Active vs Inactive") +
  xlab("Mean MPRA Activity") + ylab("-log10(p-value)") +
  xlim(c(3.85, 4.15)) + ylim(c(0, 40)) +
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
p



### (12) Perform 90bp sub region analysis of tiling oligos #############################################################

# filter MPRA results for test sequences
mpra_results_filtered <- mpra_library_and_results[which(mpra_library_and_results$Sequence_Class == "Test"), c(7, 8, 38, 39)]

# split the filtered MPRA results into 90bp sub-regions
mpra_results_90bp <- data.frame()
for(i in 1:dim(mpra_results_filtered)[1]){
  oligo_start <- as.numeric(strsplit(strsplit(mpra_results_filtered$Oligo_Coordinates_hg38[i], ":")[[1]][2], "-")[[1]][1])
  mpra_results_90bp <- rbind(mpra_results_90bp,
                              data.frame(mpra_results_filtered[i,],
                                         "X90bp_start" = oligo_start,
                                         "X90bp_end" = oligo_start + 90),
                              data.frame(mpra_results_filtered[i,],
                                         "X90bp_start" = oligo_start + 90,
                                         "X90bp_end" = oligo_start + 180),
                              data.frame(mpra_results_filtered[i,],
                                         "X90bp_start" = oligo_start + 180,
                                         "X90bp_end" = oligo_start + 270))
}

# write 90 bp results to a text file
write.table(mpra_results_90bp, file = "./tables/MPRA_Results_90bp.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# skip to this line to load 90 bp results text file while running script
mpra_results_90bp <- read.table(file = "./tables/MPRA_Results_90bp.txt", header = TRUE, sep = "\t")


# make a table to divide all CRSs into two groups, this is necessary so we don't accidentally merge
# oligos from neighboring CRSs into a single bedgraph element
mpra_results_filtered_divide <- data.frame("crs" = c(), "oligo" = c(), "split" = c())
for(i in 1:dim(mpra_results_filtered)[1]){
  if(i < 2){
    mpra_results_filtered_divide <- rbind(mpra_results_filtered_divide,
                                          data.frame("crs" = mpra_results_filtered$CRS_Coordinates_hg38[i],
                                                     "oligo" = mpra_results_filtered$Oligo_Coordinates_hg38[i],
                                                     "split" = 1))
    next
  }
  if(mpra_results_filtered_divide$crs[i-1] == mpra_results_filtered$CRS_Coordinates_hg38[i]){
    mpra_results_filtered_divide <- rbind(mpra_results_filtered_divide,
                                          data.frame("crs" = mpra_results_filtered$CRS_Coordinates_hg38[i],
                                                     "oligo" = mpra_results_filtered$Oligo_Coordinates_hg38[i],
                                                     "split" = mpra_results_filtered_divide$split[i-1]))
    next
  } else {
    mpra_results_filtered_divide <- rbind(mpra_results_filtered_divide,
                                          data.frame("crs" = mpra_results_filtered$CRS_Coordinates_hg38[i],
                                                     "oligo" = mpra_results_filtered$Oligo_Coordinates_hg38[i],
                                                     "split" = mpra_results_filtered_divide$split[i-1] + 1))
  }
}

# write bedtools splits for each CRS to a text file
write.table(mpra_results_filtered_divide, file = "./tables/MPRA_Results_90bp_Splits_for_Bedtools.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# skip to this line to load CRS splits for bedtools
mpra_results_filtered_divide <- read.table(file = "./tables/MPRA_Results_90bp_Splits_for_Bedtools.txt", header = TRUE, sep = "\t")

# split the MRPA results into even and odd CRSs
is.even <- function(x) x %% 2 == 0
mpra_results_filtered_divide$is_even <- is.even(mpra_results_filtered_divide$split)
table(mpra_results_filtered_divide$is_even)


# get first split of the 90bp MPRA results
mpra_results_90bp_v1 <- mpra_results_90bp[which(mpra_results_90bp$CRS_Coordinates_hg38 %in% mpra_results_filtered_divide$crs[which(mpra_results_filtered_divide$is_even == TRUE)]),]

# remove oligos where the DNA barcodes were not detected and write to a bed file
mpra_results_90bp_v1 <- mpra_results_90bp_v1[which(mpra_results_90bp_v1$DNA_Average != "NA" & mpra_results_90bp_v1$RNA_Average != "NA"),]

# write 90 bp MPRA results for the first split to a bed file
mpra_results_90bp_v1_bed <- data.frame("chrom" = sapply(strsplit(mpra_results_90bp_v1$Oligo_Coordinates_hg38, ":"), `[`, 1),
                                        "chromStart" = mpra_results_90bp_v1$X90bp_start,
                                        "chromEnd" = mpra_results_90bp_v1$X90bp_end - 1,
                                        "name" = (mpra_results_90bp_v1$RNA_Average + 1) / (mpra_results_90bp_v1$DNA_Average + 1))
write.table(mpra_results_90bp_v1_bed, file = "~/mpra_results_90bp_v1_bed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# sort and merge the bed file to get the 90 bp sub-region analysis results
system("bedtools sort -i ~/mpra_results_90bp_v1_bed.bed > ~/mpra_results_90bp_v1_bed.sorted.bed")
system("bedtools merge -i ~/mpra_results_90bp_v1_bed.sorted.bed -c 4 -o mean > ~/mpra_results_90bp_v1_bed.sorted.merged.bed")

# load the sub-region analysis results for the first MPRA split
mpra_results_90bp_merged_v1 <- read.table(file = "~/mpra_results_90bp_v1_bed.sorted.merged.bed", header = FALSE, sep = "\t",
                                           col.names = c("chrom", "chromStart", "chromEnd", "RNA_DNA_Ratio"))

# add +1 bp back to the chromEnd for the first split of the MPRA sub-region analysis
mpra_results_90bp_merged_v1$chromEnd <- mpra_results_90bp_merged_v1$chromEnd + 1


# get second split of the 90bp MPRA results
mpra_results_90bp_v2 <- mpra_results_90bp[which(mpra_results_90bp$CRS_Coordinates_hg38 %in% mpra_results_filtered_divide$crs[which(mpra_results_filtered_divide$is_even == TRUE)]),]

# remove oligos where the DNA barcodes were not detected and write to a bed file
mpra_results_90bp_v2 <- mpra_results_90bp_v2[which(mpra_results_90bp_v2$DNA_Average != "NA" & mpra_results_90bp_v2$RNA_Average != "NA"),]

# write 90 bp MPRA results for the second split to a bed file
mpra_results_90bp_v2_bed <- data.frame("chrom" = sapply(strsplit(mpra_results_90bp_v2$Oligo_Coordinates_hg38, ":"), `[`, 1),
                                       "chromStart" = mpra_results_90bp_v2$X90bp_start,
                                       "chromEnd" = mpra_results_90bp_v2$X90bp_end - 1,
                                       "name" = (mpra_results_90bp_v2$RNA_Average + 1) / (mpra_results_90bp_v2$DNA_Average + 1))
write.table(mpra_results_90bp_v2_bed, file = "~/mpra_results_90bp_v2_bed.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# sort and merge the bed file to get the 90 bp sub-region analysis results
system("bedtools sort -i ~/mpra_results_90bp_v2_bed.bed > ~/mpra_results_90bp_v2_bed.sorted.bed")
system("bedtools merge -i ~/mpra_results_90bp_v2_bed.sorted.bed -c 4 -o mean > ~/mpra_results_90bp_v2_bed.sorted.merged.bed")

# load the sub-region analysis results for the second MPRA split
mpra_results_90bp_merged_v2 <- read.table(file = "~/mpra_results_90bp_v2_bed.sorted.merged.bed", header = FALSE, sep = "\t",
                                          col.names = c("chrom", "chromStart", "chromEnd", "RNA_DNA_Ratio"))

# add +1 bp back to the chromEnd for the second split of the MPRA sub-region analysis
mpra_results_90bp_merged_v2$chromEnd <- mpra_results_90bp_merged_v2$chromEnd + 1


# merge the splits for the 90 bp sub-region analysis
mpra_results_90bp_results <- rbind(mpra_results_90bp_merged_v1,
                                   mpra_results_90bp_merged_v2)

# write the results of the sub-region analysis to a bedgraph file
write.table(mpra_results_90bp_results, file = "./mpra_results_90bp.bedgraph",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



### (13) Overlap MPRA results with human brain chip-seq experiments ####################################################

# get a bed file of all oligo regions for bigWigAverageOverBed analysis and biochemical model
all_oligos <- mpra_library_and_results[which(mpra_library_and_results$Sequence_Class == "Test"),]
all_oligos_bed <- data.frame("chrom" = sapply(strsplit(all_oligos$Oligo_Coordinates_hg38, ":"), `[`, 1),
                             "chromStart" = sapply(strsplit(sapply(strsplit(all_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 1),
                             "chromEnd" = sapply(strsplit(sapply(strsplit(all_oligos$Oligo_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 2),
                             "name" = all_oligos$Oligo_Coordinates_hg38)
all_oligos_bed_no_duplicates <- unique(all_oligos_bed)
write.table(all_oligos_bed_no_duplicates, file = "./bed/all_oligos_no_duplicates_for_bigWigAverageOverBed_analysis.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# load neuron specific and cell-line data
get_coords <- read.table(file = "./bigwig_activity_neuron/ASCL1_DLPFC-NEUN_Output_basename_prefix.txt", header = FALSE, sep = "\t",
                         col.names = c("name", "size", "covered", "sum", "mean0", "mean"))
bigwig_activities_neuron <- data.frame("Oligo_Coordinates_hg38" = get_coords$name)

# rename datasets by target
files <- list.files(path = "/Users/npage/Dropbox/Neurohub_CRISPRa_CRT/bigwig_activity_neuron", pattern = "*.txt", full.names = TRUE, recursive = FALSE)
datasets <- str_replace_all(files, "/Users/npage/Dropbox/Neurohub_CRISPRa_CRT/bigwig_activity_neuron/", "")
datasets <- str_replace_all(datasets, "_Output_basename_prefix.txt", "")
datasets <- str_replace_all(datasets, "_basename_prefix.txt", "")
datasets <- str_replace_all(datasets, "_Output_sorted_basename_prefix.txt", "")
datasets <- str_replace_all(datasets, "_Output_sorted", "")
datasets <- str_replace_all(datasets, "-", "_")
datasets <- str_replace_all(datasets, "1224_1230_", "")
datasets <- str_replace_all(datasets, "1238_1242_", "")
datasets <- str_replace_all(datasets, "WTC11_NGN2_iPSC_7_8wk_ExN_ATAC_seq_Song_2019_hg38.txt", "ATAC_WTC11_NGN2_iPSC_7_8wk_ExN")
datasets <- gsub("WTC11_NGN2_iPSC_7_8wk_ExN_H3K27ac_CUT+RUN_Song_2019_hg38.txt", "H3K27ac_WTC11_NGN2_iPSC_7_8wk_ExN", datasets, fixed = TRUE)
datasets

# loop through all files in analysis to get their average bigwig activity scores
for(i in 1:length(files)){
  bigwig_activity_neuron <- read.table(file = files[i], header = FALSE, sep = "\t",
                                       col.names = c("name", "size", "covered", "sum", "mean0", "mean"))
  bigwig_activities_neuron[, datasets[i]] <- log10(bigwig_activity_neuron$mean + 0.1)
  
}

# add MPRA activity to bigwig results dataframe
bigwig_activities_neuron$mpra_activity <- mpra_library_and_results$Average_Log2_Count_Ratio[match(bigwig_activities_neuron$Oligo_Coordinates_hg38, mpra_library_and_results$Oligo_Coordinates_hg38)]

# filter for non-NA MPRA activities
bigwig_activities_neuron_no_NA <- bigwig_activities_neuron[which(is.na(bigwig_activities_neuron$mpra_activity) == FALSE),]

# scale the data
bigwig_activities_neuron_no_NA[,c(2:259)] <- scale(bigwig_activities_neuron_no_NA[,c(2:259)])

# get the list of all unique targets
target_list <- unique(sapply(strsplit(colnames(bigwig_activities_neuron_no_NA)[2:258], "_"), `[`, 1))
length(target_list) # 74 unique targets

# merge the bigwig activities tracks by target
bigwig_activities_neuron_no_NA_target_merged <- bigwig_activities_neuron_no_NA[, c(1, 259)]
for(target in target_list){
  bigwig_activities_neuron_no_NA_target_merged[, target] <- rowMeans(select(bigwig_activities_neuron_no_NA, contains(target)))
}


# create empty correlation table to fill with chip-seq results
correlation_table_neuron <- data.frame("dataset" = c(),
                                       "spearman" = c())

# get the spearman correlation between every target track and the MPRA activity
for(i in 3:76){
  correlation <- data.frame("dataset" = colnames(bigwig_activities_neuron_no_NA_target_merged)[i],
                            "spearman" = cor(x = bigwig_activities_neuron_no_NA_target_merged[,i], y = bigwig_activities_neuron_no_NA_target_merged$mpra_activity, method = "spearman"))
  correlation_table_neuron <- rbind(correlation_table_neuron,
                                    correlation)
}

# order the correlation table results by factor
correlation_table_neuron <- correlation_table_neuron[order(correlation_table_neuron$spearman, decreasing = FALSE),]
correlation_table_neuron$dataset <- factor(correlation_table_neuron$dataset, levels = correlation_table_neuron$dataset)

# plot the spearman correlation between every bigwig track and MPRA activity
ggplot(correlation_table_neuron, aes(x = dataset, y = spearman)) +
  geom_point() + 
  geom_segment(aes(x = dataset, xend = dataset, y = 0, yend = spearman)) +
  ggtitle("Neuron-Specific") +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 6))


# load the bulk bigwig track data
get_coords <- read.table(file = "./bigwig_activity_bulk/ARID1B_CB_Output_30M.txt", header = FALSE, sep = "\t",
                         col.names = c("name", "size", "covered", "sum", "mean0", "mean"))
bigwig_activities_bulk <- data.frame("Oligo_Coordinates_hg38" = get_coords$name)

# rename datasets by target
files <- list.files(path = "/Users/npage/Dropbox/Neurohub_CRISPRa_CRT/bigwig_activity_bulk", pattern = "*.txt", full.names = TRUE, recursive = FALSE)
datasets <- str_replace_all(files, "/Users/npage/Dropbox/Neurohub_CRISPRa_CRT/bigwig_activity_bulk/", "")
datasets <- str_replace_all(datasets, "_Output_25M.txt", "")
datasets <- str_replace_all(datasets, "_Output_30M.txt", "")
datasets <- str_replace_all(datasets, "_Output_40M.txt", "")
datasets <- str_replace_all(datasets, "_Output_45M.txt", "")
datasets <- str_replace_all(datasets, "_dedup_50M_Output_50M.txt", "")
datasets <- str_replace_all(datasets, "_30M.txt", "")
datasets <- str_replace_all(datasets, "GW23_hs-PFC_ARID1B_S47_L007_R1_001.txt", "ARID1B_DLPFC_GW23")
datasets <- str_replace_all(datasets, "GW23_hs-PFC_Bcl11a_S41_L006_R1_001.txt", "BCL11A_DLPFC_GW23")
datasets <- str_replace_all(datasets, "GW23_hs-PFC_FoxP1_S44_L006_R1_001.txt", "FOXP1_DLPFC_GW23")
datasets <- str_replace_all(datasets, "GW23_hs-PFC_Tbr1_S38_L006_R1_001.txt", "TBR1_DLPFC_GW23")
datasets <- str_replace_all(datasets, "GW23_hs-PFC_Tcf4_S50_L007_R1_001.txt", "TCF4_DLPFC_GW23")
datasets <- str_replace_all(datasets, "ATAC_1224_1230_DLPFC_Bulk_basename_prefix.txt", "ATAC_DLPFC")
datasets <- str_replace_all(datasets, "GSM4495201_atac-14gw_cge-rep1-nodup-fold_change-hg38.txt", "ATAC_CGE_GW14")
datasets <- str_replace_all(datasets, "GSM4495202_atac-14gw_lge-rep1-nodup-fold_change-hg38.txt", "ATAC_LGE_GW14")
datasets <- str_replace_all(datasets, "GSM4495203_atac-14gw_mge-rep1-nodup-fold_change-hg38.txt", "ATAC_MGE_GW14")
datasets <- str_replace_all(datasets, "GSM4495204_atac-14gw_motor-rep1-nodup-fold_change-hg38.txt", "ATAC_Motor_GW14")
datasets <- str_replace_all(datasets, "GSM4495205_atac-14gw_pfc-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_GW14")
datasets <- str_replace_all(datasets, "GSM4495206_atac-14gw_s1-rep1-nodup-fold_change-hg38.txt", "ATAC_S1_GW14")
datasets <- str_replace_all(datasets, "GSM4495207_atac-14gw_v1-rep1-nodup-fold_change-hg38.txt", "ATAC_V1_GW14")
datasets <- str_replace_all(datasets, "GSM4495208_atac-17gw_cge-rep1-nodup-fold_change-hg38.txt", "ATAC_CGE_GW17")
datasets <- str_replace_all(datasets, "GSM4495209_atac-17gw_motor-rep1-nodup-fold_change-hg38.txt", "ATAC_Motor_GW17")
datasets <- str_replace_all(datasets, "GSM4495210_atac-17gw_pfc-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_GW17")
datasets <- str_replace_all(datasets, "GSM4495211_atac-17gw_s1-rep1-nodup-fold_change-hg38.txt", "ATAC_S1_GW17")
datasets <- str_replace_all(datasets, "GSM4495212_atac-17gw_v1-rep1-nodup-fold_change-hg38.txt", "ATAC_V1_GW17")
datasets <- str_replace_all(datasets, "GSM4495213_atac-18gw_cge-rep1-nodup-fold_change-hg38.txt", "ATAC_CGE_GW18")
datasets <- str_replace_all(datasets, "GSM4495214_atac-18gw_lge-rep1-nodup-fold_change-hg38.txt", "ATAC_LGE_GW18")
datasets <- str_replace_all(datasets, "GSM4495215_atac-18gw_mge-rep1-nodup-fold_change-hg38.txt", "ATAC_MGE_GW18")
datasets <- str_replace_all(datasets, "GSM4495216_atac-18gw_pfc-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_GW18")
datasets <- str_replace_all(datasets, "GSM4495217_atac-18gw_pfc_dl-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_DL_GW18")
datasets <- str_replace_all(datasets, "GSM4495218_atac-18gw_pfc_ul-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_UL_GW18")
datasets <- str_replace_all(datasets, "GSM4495219_atac-19gw_cge-rep1-nodup-fold_change-hg38.txt", "ATAC_CGE_GW19")
datasets <- str_replace_all(datasets, "GSM4495220_atac-19gw_lge-rep1-nodup-fold_change-hg38.txt", "ATAC_LGE_GW19")
datasets <- str_replace_all(datasets, "GSM4495221_atac-19gw_mge-rep1-nodup-fold_change-hg38.txt", "ATAC_MGE_GW19")
datasets <- str_replace_all(datasets, "GSM4495222_atac-19gw_motor-rep1-nodup-fold_change-hg38.txt", "ATAC_Motor_GW19")
datasets <- str_replace_all(datasets, "GSM4495223_atac-19gw_parietal-rep1-nodup-fold_change-hg38.txt", "ATAC_Parietal_GW19")
datasets <- str_replace_all(datasets, "GSM4495224_atac-19gw_pfc-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_GW19")
datasets <- str_replace_all(datasets, "GSM4495225_atac-19gw_pfc_dl-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_DL_GW19")
datasets <- str_replace_all(datasets, "GSM4495226_atac-19gw_pfc_ul-rep1-nodup-fold_change-hg38.txt", "ATAC_PFC_UL_GW19")
datasets <- str_replace_all(datasets, "GSM4495227_atac-19gw_s1-rep1-nodup-fold_change-hg38.txt", "ATAC_S1_GW19")
datasets <- str_replace_all(datasets, "GSM4495228_atac-19gw_temporal-rep1-nodup-fold_change-hg38.txt", "ATAC_Temporal_GW19")
datasets <- str_replace_all(datasets, "GSM4495229_atac-19gw_v1-rep1-nodup-fold_change-hg38.txt", "ATAC_V1_GW19")
datasets <- str_replace_all(datasets, "GSM4495230_chip-15gw_H3K27ac-nodup-hg38.txt", "H3K27ac_GW15")
datasets <- str_replace_all(datasets, "GSM4495231_chip-15gw_H3K4me1-nodup-hg38.txt", "H3K4me1_GW15")
datasets <- str_replace_all(datasets, "GSM4495232_chip-15gw_H4K20me3-nodup-hg38.txt", "H4K20me3_GW15")
datasets <- str_replace_all(datasets, "GSM4495233_chip-17gw_H3K27ac-nodup-hg38.txt", "H3K27ac_GW17")
datasets <- str_replace_all(datasets, "GSM4495234_chip-17gw_H3K27me3-nodup-hg38.txt", "H3K27me3_GW17")
datasets <- str_replace_all(datasets, "GSM4495235_chip-17gw_H3K4me1-nodup-hg38.txt", "H3K4me1_GW17")
datasets <- str_replace_all(datasets, "GSM4495236_chip-18gw_H3K27ac-nodup-hg38.txt", "H3K27ac_GW18")
datasets <- str_replace_all(datasets, "GSM4495237_chip-18gw_H4K20me3-nodup-hg38.txt", "H4K20me3_GW18")
datasets <- str_replace_all(datasets, "GSM4495238_chip-22gw_H3K20me3_CP-nodup-hg38.txt", "H3K20me3_CP_GW22")
datasets <- str_replace_all(datasets, "GSM4495239_chip-22gw_H3K27me3-nodup-hg38.txt", "H3K27me3_GW22")
datasets <- str_replace_all(datasets, "GSM4495240_chip-22gw_H3K9me3_CP-nodup-hg38.txt", "H3K9me3_CP_GW22")
datasets <- str_replace_all(datasets, ".txt", "")
datasets

# loop through all files in analysis to get their average bigwig activity scores
for(i in 1:length(files)){
  bigwig_activity_bulk <- read.table(file = files[i], header = FALSE, sep = "\t",
                                       col.names = c("name", "size", "covered", "sum", "mean0", "mean"))
  bigwig_activities_bulk[, datasets[i]] <- log10(bigwig_activity_bulk$mean + 0.1)
  
}

# add MPRA activity to bigwig results dataframe
bigwig_activities_bulk$mpra_activity <- mpra_library_and_results$Average_Log2_Count_Ratio[match(bigwig_activities_bulk$Oligo_Coordinates_hg38, mpra_library_and_results$Oligo_Coordinates_hg38)]

# filter for non-NA MPRA activities
bigwig_activities_bulk_no_NA <- bigwig_activities_bulk[which(is.na(bigwig_activities_bulk$mpra_activity) == FALSE),]

# scale the data
bigwig_activities_bulk_no_NA[,c(2:580)] <- scale(bigwig_activities_bulk_no_NA[,c(2:580)])

# get the list of all unique targets
target_list <- unique(sapply(strsplit(colnames(bigwig_activities_bulk_no_NA)[2:580], "_"), `[`, 1))
length(target_list) # 109 unique targets

# merge the bigwig activities tracks by target
bigwig_activities_bulk_no_NA_target_merged <- bigwig_activities_bulk_no_NA[, c(1, 581)]
for(target in target_list){
  bigwig_activities_bulk_no_NA_target_merged[, target] <- rowMeans(select(bigwig_activities_bulk_no_NA, contains(target)))
}


# create empty correlation table to fill with chip-seq results
correlation_table_bulk <- data.frame("dataset" = c(),
                                     "spearman" = c())

# get the spearman correlation between every target track and the MPRA activity
for(i in 3:111){
  correlation <- data.frame("dataset" = colnames(bigwig_activities_bulk_no_NA_target_merged)[i],
                            "spearman" = cor(x = bigwig_activities_bulk_no_NA_target_merged[,i], y = bigwig_activities_bulk_no_NA_target_merged$mpra_activity, method = "spearman"))
  correlation_table_bulk <- rbind(correlation_table_bulk,
                                    correlation)
}

# order the correlation table results by factor
correlation_table_bulk <- correlation_table_bulk[order(correlation_table_bulk$spearman, decreasing = FALSE),]
correlation_table_bulk$dataset <- factor(correlation_table_bulk$dataset, levels = correlation_table_bulk$dataset)

# plot the spearman correlation between every bigwig track and MPRA activity
ggplot(correlation_table_bulk, aes(x = dataset, y = spearman)) +
  geom_point() + 
  geom_segment(aes(x = dataset, xend = dataset, y = 0, yend = spearman)) +
  ggtitle("Bulk-Tissue") +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 5))



### (14) Comparison of MPRA results to Vista enhancer atlas ############################################################

# get vista validated enhancer elements
vista_elements <- read.table(file = "./data/published_datasets/vista/experiments.tsv", header = TRUE, sep = "\t")
head(vista_elements[,c(1:18)])
dim(vista_elements)
colnames(vista_elements)

# get a bed file of the test CRSs to overlap with vista elements
test_crs <- unique(mpra_library_and_results[which(mpra_library_and_results$Sequence_Class == "Test"), c("CRS_Coordinates_hg38", "CRS_Contains_Oligo_in_Top_10_Percent")])
test_crs_bed <- data.frame("chrom" = sapply(strsplit(test_crs$CRS_Coordinates_hg38, ":"), `[`, 1),
                           "chromStart" = sapply(strsplit(sapply(strsplit(test_crs$CRS_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 1),
                           "chromEnd" = sapply(strsplit(sapply(strsplit(test_crs$CRS_Coordinates_hg38, ":"), `[`, 2), "-"), `[`, 2), 
                           "name" = test_crs$CRS_Contains_Oligo_in_Top_10_Percent)
write.table(test_crs_bed, file = "./bed/CRS_Contains_Oligo_in_Top_10_Percent.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# get a bed file of human positive vista elements
vista_elements_filtered <- vista_elements[which(vista_elements$organism == "Human" & vista_elements$curation_status != "allelic"),]
vista_elements_bed <- data.frame("chrom" = sapply(strsplit(vista_elements_filtered$coordinate_hg38, ":"), `[`, 1),
                                 "chromStart" = sapply(strsplit(sapply(strsplit(vista_elements_filtered$coordinate_hg38, ":"), `[`, 2), "-"), `[`, 1),
                                 "chromEnd" = sapply(strsplit(sapply(strsplit(vista_elements_filtered$coordinate_hg38, ":"), `[`, 2), "-"), `[`, 2), 
                                 "name" = vista_elements_filtered$tissue)
write.table(vista_elements_bed, file = "./bed/Human_Vista_Elements.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


# bedtools intersect vista positive enhancers with CRSs containing at least one oligo in the top 10%
system("bedtools intersect -a ./bed/CRS_Contains_Oligo_in_Top_10_Percent.bed -b ./bed/Human_Vista_Elements.bed -wa -wb > ./bed/CRS_Overlap_Vista_Elements.txt")
crs_overlap_vista_elements <- read.table(file = "./bed/CRS_Overlap_Vista_Elements.txt", header = FALSE, sep = "\t",
                                         col.names = c("chromCRS", "chromStartCRS", "chromEndCRS", "contains_oligo_in_top_10_Percent", 
                                                       "chromVista", "chromStartVista", "chromEndVista", "tissue_positive"))

# get a list of tissues showing vista activity
tissue_list <-  table(strsplit(paste(crs_overlap_vista_elements$tissue_positive, collapse = ";"), ";")[[1]])
tissue_list <- tissue_list[which(tissue_list != "")]
tissue_list

# create a table of tissue activity for each CRS tested via MPRA
for(i in 2:length(tissue_list)){
  crs_overlap_vista_elements[, names(tissue_list[i])] <- ifelse(grepl(names(tissue_list[i]), crs_overlap_vista_elements$tissue_positive, fixed = TRUE), 1, 0)
}

# create dataframe to track activity-tissue overlaps
activity_tissue_overlap <- data.frame("tissue" = c(),
                                      "OR" = c(),
                                      "lower_conf" = c(),
                                      "upper_conf" = c(),
                                      "p-value" = c())

# perform Fisher's exact test between MPRA active CRSs and each vista positive tissue
for (i in 9:26){
  
  # Skip tissues with no positive vista enhancers
  if (table(crs_overlap_vista_elements[i])[1] == dim(crs_overlap_vista_elements)[1]) {
    next
  }
  
  # Skip tissues with less than 20 total positive vista enhancers
  if (table(crs_overlap_vista_elements[i])[2] < 20) {
    next
  }
  
  # Get the number of vista tissue positive enhancers with MPRA activity
  enhancer_with_activity <- dim(crs_overlap_vista_elements[which(crs_overlap_vista_elements$contains_oligo_in_top_10_Percent == TRUE & crs_overlap_vista_elements[,i] == 1),])[1]
  
  # Get the number of vista tissue positive enhancers with no MPRA activity
  enhancer_no_activity <- dim(crs_overlap_vista_elements[which(crs_overlap_vista_elements$contains_oligo_in_top_10_Percent == FALSE & crs_overlap_vista_elements[,i] == 1),])[1]
  
  # Get the number of vista tissue negative enhancers with MPRA activity
  not_enhancer_with_activity <- dim(test_crs_bed[which(test_crs_bed$name == TRUE),])[1] - enhancer_with_activity
  
  # Get the number of vista tissue negative enhancers with no MPRA activity
  not_enhancer_no_activity <- dim(test_crs_bed[which(test_crs_bed$name == FALSE),])[1] - enhancer_no_activity
  
  # Generate a contingency table for computing the odds ratio
  enhancer_contingency <- data.frame(x1 = c(enhancer_with_activity, enhancer_no_activity), 
                                     x2 = c(not_enhancer_with_activity, not_enhancer_no_activity))
  colnames(enhancer_contingency) <- c("enhancer", "not_enhancer")
  rownames(enhancer_contingency) <- c("activity", "no_activity")
  enhancer_contingency
  
  # Calculate the odds of being a tissue specific vista enhancer
  OR_enhancer <- fisher.test(enhancer_contingency)
  OR_enhancer
  
  # Add tissue specific results to results list
  activity_tissue_overlap <- rbind(activity_tissue_overlap,
                               data.frame("tissue" = colnames(crs_overlap_vista_elements)[i],
                                          "OR" = OR_enhancer$estimate,
                                          "lower_conf" = OR_enhancer$conf.int[1],
                                          "upper_conf" = OR_enhancer$conf.int[2],
                                          "p-value" = OR_enhancer$p.value))
  
}

# Calculate FDR for all tests performed
activity_tissue_overlap$fdr <- p.adjust(activity_tissue_overlap$p.value, method = "fdr")

# rename tissues in results table
activity_tissue_overlap$tissue <- c("forebrain", "hindbrain", "midbrain", "neural_tube")
activity_tissue_overlap

# Plot results
p1 <- ggplot(activity_tissue_overlap, aes(x = tissue, y = OR)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower_conf, ymax = upper_conf), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label = ifelse(fdr < 0.0001, "****", ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", "")))), y = upper_conf), vjust = 0, color = "red", size = 8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  ggtitle("cCREs with top 10% active tile") +
  ylim(0, 13) +
  theme_classic() +
  theme(plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14))
print(p1)



### (15) Select final candidate enhancer regions for CRISPRa screen ####################################################

# final candidates for CRISPRa screen are cCREs with top 10% active tile
CRISPRa_candidate_cCREs <- unique(mpra_library_and_results$CRS_Coordinates_hg38[which(mpra_library_and_results$CRS_Contains_Oligo_in_Top_10_Percent)])

# write candidate enhancer regions for CRISPRa screen to a bed file
CRISPRa_candidate_cCREs_bed <- data.frame("chrom" = sapply(strsplit(CRISPRa_candidate_cCREs, ":"), `[`, 1),
                                          "chromStart" = sapply(strsplit(sapply(strsplit(CRISPRa_candidate_cCREs, ":"), `[`, 2), "-"), `[`, 1),
                                          "chromEnd" = sapply(strsplit(sapply(strsplit(CRISPRa_candidate_cCREs, ":"), `[`, 2), "-"), `[`, 2),
                                          CRISPRa_candidate_cCREs)
write.table(CRISPRa_candidate_cCREs_bed, file = "./bed/CRISPRa_candidate_cCREs.bed",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# double check the number of cCREs for CRISPRa screen
system("head ./bed/CRISPRa_candidate_cCREs.bed") 
system("wc -l ./bed/CRISPRa_candidate_cCREs.bed") # looks good



###################################################### END SCRIPT ######################################################




