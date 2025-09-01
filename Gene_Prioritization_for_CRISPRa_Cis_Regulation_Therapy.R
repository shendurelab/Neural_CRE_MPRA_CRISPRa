########################################################################################################################
#
#   01_Gene_Prioritization_for_CRISPRa_Cis_Regulation_Therapy.R
#
#   Determine the list of genes whose predicted developmental enhancers will be tested via CRISPR-QTL for the discovery 
#   of candidates to be used in CRT therapy for heterozygous loss of function mutations. This list is based on whole 
#   exome sequencing based gene discovery in autism (Satterstrom et al. 2020; Cell), developmental disorders (Kaplanis 
#   et al. 2020; Nature), schizophrenia (Singh et al. 2022; Nature), and Tourette's (Wang et al. 2018; Cell Reports),
#   obsessive compulsive disorder (Halvorsen et al. 2021; Nature Neuroscience), and clinician curated developmental 
#   disorder gene lists (DDG2P).
#
#   Nicholas Page, January 2021
#   Sanders Lab, Dept. of Psychiatry, University of California, San Francisco
#   Neuroscience Graduate Program, University of California, San Francisco
#
########################################################################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(scipen = 20)

library(denovolyzeR)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ggvenn)
library(Cairo)
library(xlsx)
library(ape)

setwd('~/Dropbox/Neurohub_CRISPRa_CRT/')

########################################################################################################################
#
#   Pipeline:
#
#   (1) Download gencode hg38 annotations
#   (2) Download de novo mutations from whole exome sequencing studies
#   (3) Download DDG2P database, filter for genes causing monoallelic brain disorders, and compile candidate genes
#   (4) Combine de novo mutation counts from each whole exome sequencing study
#   (5) Use denovolyzer to test for an excess of de novo mutations across all genes
#   (6) Visualize results and output final gene list
#
########################################################################################################################

### (1) Download gencode hg38 annotations ##############################################################################

# Download hg38 gencode annotation
if(file.exists("~/genomes/hg38/gencode.v38.basic.annotation.gff3") == FALSE){
  setwd("~/genomes/hg38/")
  system("wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.basic.annotation.gff3.gz")
  system("gunzip ~/genomes/hg38/gencode.v38.basic.annotation.gff3.gz")
  setwd('~/Dropbox/Neurohub_CRISPRa_CRT/')
}
gencode_hg38 <- read.gff(file = "~/genomes/hg38/gencode.v38.basic.annotation.gff3", GFF3 = TRUE)
gencode_hg38_genes <- gencode_hg38[which(gencode_hg38$type == "gene"),]
annotations <- data.frame("gene_name" = substring(sapply(strsplit(gencode_hg38_genes$attributes, ";"), `[`, 4), 11),
                          "gene_id" = substring(sapply(strsplit(gencode_hg38_genes$attributes, ";"), `[`, 2), 9, 23),
                          "hgnc_id" = substring(sapply(strsplit(gencode_hg38_genes$attributes, ";"), `[`, 6), 9),
                          "chrom" = gencode_hg38_genes$seqid,
                          "start_position_hg38" = gencode_hg38_genes$start,
                          "end_position_hg38" = gencode_hg38_genes$end,
                          "strand" = gencode_hg38_genes$strand)

transcript_lengths <- read.table(file = "./data/annotations/GRCh38.p13_ensembl_biomart_transcript_lengths.txt", header = TRUE, sep = "\t")
cds_lengths <- read.table(file = "./data/annotations/GRCh38.p13_ensembl_biomart_cds_lengths.txt", header = TRUE, sep = "\t")

new_annotations <- data.frame("transcript_ids" = c(),
                              "transcript_lengths" = c(),
                              "cds_lengths" = c())

for(i in 1:dim(annotations)[1]){
  new_annotations <- rbind(new_annotations,
                           data.frame("transcript_ids" = c(paste(transcript_lengths$Transcript.stable.ID[which(transcript_lengths$Gene.stable.ID %in% annotations$gene_id[i])], collapse = ";")),
                                      "transcript_lengths" = c(paste(transcript_lengths$Transcript.length..including.UTRs.and.CDS.[which(transcript_lengths$Gene.stable.ID %in% annotations$gene_id[i])], collapse = ";")),
                                      "cds_lengths" = c(paste(cds_lengths$CDS.Length[which(cds_lengths$Gene.stable.ID %in% annotations$gene_id[i])], collapse = ";"))))
}

annotations <- cbind(annotations, new_annotations)

write.table(annotations, file = "./bed/annotations/hg38/gencode.v38.basic.annotation.coding.transcripts.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
annotations <- read.table(file = "./bed/annotations/hg38/gencode.v38.basic.annotation.coding.transcripts.txt", header = TRUE, sep = "\t")



### (2) Download de novo mutations from whole exome sequencing studies #################################################

# Contains list of the de novo variants called in the 6,430 individuals with ASD and 2,179 individuals without ASD 
# (a total of 15,789 de novo variants were identified covering 9,345 coding variants and 6,444 non-coding variants)
satterstrom_ASD <- read_xlsx(path = "./genetics/1-s2.0-S0092867419313984-mmc1.xlsx", sheet = "De novo variants")

# Split the satterstrom_ASD dataset into and affected vs control subset
satterstrom_ASD_affected <- satterstrom_ASD[which(satterstrom_ASD$Affected_Status == 2),]
satterstrom_ASD_control <- satterstrom_ASD[which(satterstrom_ASD$Affected_Status == 1),]

# Results from analysis on the full cohort and on the undiagnosed subset, along with gene-level DNM counts
# (a total of 44,901 variants from 31,058 trios diagnosed with a developmental disorder) (Kaplanis et al. 2020)
kaplanis_DDD <- read_xlsx(path = "./genetics/41586_2020_2832_MOESM4_ESM.xlsx", sheet = "kaplanis_samocha_denovoWEST_res")

# Get results from Singh et al. 2022 (Nature) (Schema) which includes 3,402 trios diagnosed with schizophrenia
singh_SCZ <- data.frame(read_xlsx(path = "./genetics/41586_2022_4556_MOESM3_ESM.xlsx", sheet = "Table S5 - Gene Results"))
singh_SCZ[singh_SCZ == "NA"] <- 0
singh_SCZ[is.na(singh_SCZ)] <- 0

# Get results from Wang et al. 2018 (Cell Reports) which includes 802 trios diagnosed with Tourette's
wang_TD <- data.frame(read_xlsx(path = "./genetics/mmc4.xlsx", sheet = "Non-multiplex families"))

# Get results from Halvorsen et al. 2021 (Nat Neuro) which includes 671 de novo variants from 587 trios and 41
# quartets diagnosed with obsessive compulsive disorder (OCD)
halvorsen_OCD_trios <- data.frame(read_xlsx(path = "./genetics/41593_2021_876_MOESM3_ESM.xlsx", sheet = "Table S8", skip = 1))
halvorsen_OCD_quartets <- data.frame(read_xlsx(path = "./genetics/41593_2021_876_MOESM3_ESM.xlsx", sheet = "Table S9", skip = 1))
halvorsen_OCD <- rbind(halvorsen_OCD_trios, halvorsen_OCD_quartets)



### (3) Download DDG2P database, filter for genes causing monoallelic brain disorders, and compile candidate genes #####

# Download DDG2P database
decipher_DDG2P <- data.frame(read_csv(file = "./genetics/DDG2P.csv"))

# Filter for genes with hemizygous or monoallelic mutations that result in brain related disorders
decipher_DDG2P_brain <- decipher_DDG2P[which(grepl("Brain", decipher_DDG2P$organ.specificity.list, fixed = TRUE)),]
decipher_DDG2P_brain_het_lof <- decipher_DDG2P_brain[which((decipher_DDG2P_brain$allelic.requirement == "monoallelic") & 
                                                             decipher_DDG2P_brain$mutation.consequence == "loss of function"),]
decipher_DDG2P_brain_het_lof_gene_list <- decipher_DDG2P_brain_het_lof$gene.symbol

# Get 102 ASD candidate genes from Satterstrom et al. 2020
satterstrom_ASD_gene_list <- read_xlsx(path = "./genetics/1-s2.0-S0092867419313984-mmc2.xlsx", sheet = "102_ASD", n_max = 102)
satterstrom_ASD_gene_list <- satterstrom_ASD_gene_list$gene
satterstrom_ASD_gene_list <- replace(satterstrom_ASD_gene_list, which(satterstrom_ASD_gene_list == "SUV420H1"), "KMT5B")

# Get 285 genome wide significant developmental disorder genes from Kaplanis et al. 2020
kaplanis_DDD_gene_list <- read_xlsx(path = "./genetics/41586_2020_2832_MOESM4_ESM.xlsx", sheet = "kaplanis_samocha_denovoWEST_res")
kaplanis_DDD_gene_list <- kaplanis_DDD_gene_list$symbol[which(kaplanis_DDD_gene_list$significant == TRUE)]

# Get 32 candidate schizophrenia genes from Singh et al. 2020 (Schema)
singh_SCZ_gene_list <- data.frame(read_xlsx(path = "./genetics/41586_2022_4556_MOESM3_ESM.xlsx", sheet = "Table S5 - Gene Results"))
singh_SCZ_gene_list <- singh_SCZ_gene_list$Gene.Symbol[which(singh_SCZ_gene_list$Q.meta < 0.05)]

# Get 2 high-confidence Tourette's associated genes from Wang et al. 2018
wang_TD_gene_list <- wang_TD$gene.id[which(wang_TD$qval < 0.1)]

# Merge candidate gene lists
candidate_genes <- unique(c(decipher_DDG2P_brain_het_lof_gene_list,
                            satterstrom_ASD_gene_list,
                            kaplanis_DDD_gene_list,
                            singh_SCZ_gene_list,
                            wang_TD_gene_list))

length(candidate_genes) # 466 curated candidate neurodevelopmental and neuropsychiatric disorder genes



### (4) Combine de novo mutation counts from each whole exome sequencing study #########################################

satterstrom_ASD_affected_DN <- data.frame("gene" = c(), 
                                          "dn.non" = c(), 
                                          "dn.mis" = c(),
                                          "dn.frameshift" = c(),
                                          "dn.splice" = c(),
                                          "dn.syn" = c())

for(i in 1:dim(satterstrom_ASD_affected)[1]){
  
  if(satterstrom_ASD_affected$VEP_functional_class_canonical_simplified[i] == "stop_gained"){
    satterstrom_ASD_affected_DN <- rbind(satterstrom_ASD_affected_DN,
                                data.frame("gene" = c(satterstrom_ASD_affected$GENE_NAME[i]), 
                                           "dn.non" = c(1), 
                                           "dn.mis" = c(0),
                                           "dn.frameshift" = c(0),
                                           "dn.splice" = c(0),
                                           "dn.syn" = c(0)))
  }
  
  if(satterstrom_ASD_affected$VEP_functional_class_canonical_simplified[i] == "missense_variant"){
    satterstrom_ASD_affected_DN <- rbind(satterstrom_ASD_affected_DN,
                                data.frame("gene" = c(satterstrom_ASD_affected$GENE_NAME[i]), 
                                           "dn.non" = c(0), 
                                           "dn.mis" = c(1),
                                           "dn.frameshift" = c(0),
                                           "dn.splice" = c(0),
                                           "dn.syn" = c(0))) 
  }
  
  if(satterstrom_ASD_affected$VEP_functional_class_canonical_simplified[i] == "frameshift_variant"){
    satterstrom_ASD_affected_DN <- rbind(satterstrom_ASD_affected_DN,
                                data.frame("gene" = c(satterstrom_ASD_affected$GENE_NAME[i]), 
                                           "dn.non" = c(0), 
                                           "dn.mis" = c(0),
                                           "dn.frameshift" = c(1),
                                           "dn.splice" = c(0),
                                           "dn.syn" = c(0)))
  }
  
  if(satterstrom_ASD_affected$VEP_functional_class_canonical_simplified[i] == "splice_donor_variant" |
     satterstrom_ASD_affected$VEP_functional_class_canonical_simplified[i] == "splice_acceptor_variant"){
    satterstrom_ASD_affected_DN <- rbind(satterstrom_ASD_affected_DN,
                                data.frame("gene" = c(satterstrom_ASD_affected$GENE_NAME[i]), 
                                           "dn.non" = c(0), 
                                           "dn.mis" = c(0),
                                           "dn.frameshift" = c(0),
                                           "dn.splice" = c(1),
                                           "dn.syn" = c(0)))
  }
  
  if(satterstrom_ASD_affected$VEP_functional_class_canonical_simplified[i] == "synonymous_variant"){
    satterstrom_ASD_affected_DN <- rbind(satterstrom_ASD_affected_DN,
                                data.frame("gene" = c(satterstrom_ASD_affected$GENE_NAME[i]), 
                                           "dn.non" = c(0), 
                                           "dn.mis" = c(0),
                                           "dn.frameshift" = c(0),
                                           "dn.splice" = c(0),
                                           "dn.syn" = c(1)))
  }
}

# Summarize the de novo mutation counts from Satterstrom et al. 2020 by gene
satterstrom_ASD_affected_DN <- satterstrom_ASD_affected_DN %>% 
                        group_by(gene) %>% 
                        summarise_all(funs(sum))
satterstrom_ASD_affected_DN <- data.frame(satterstrom_ASD_affected_DN)
satterstrom_ASD_affected_DN$gene <- replace(satterstrom_ASD_affected_DN$gene, which(satterstrom_ASD_affected_DN$gene == "SUV420H1"), "KMT5B")


# Get the mutation counts by gene from Kaplanis et al. 2020
kaplanis_DDD_DN <- data.frame("gene" = kaplanis_DDD$symbol, 
                              "dn.non" = kaplanis_DDD$stop_gained, 
                              "dn.mis" = kaplanis_DDD$missense_variant,
                              "dn.frameshift" = kaplanis_DDD$frameshift_variant,
                              "dn.splice" = kaplanis_DDD$splice_acceptor_variant + kaplanis_DDD$splice_donor_variant,
                              "dn.syn" = kaplanis_DDD$synonymous_variant)

# Fix gene name formatting errors in kaplanis_DDD_DN dataframe
kaplanis_DDD_DN <- kaplanis_DDD_DN[order(kaplanis_DDD_DN$gene, decreasing = FALSE),]
kaplanis_DDD_DN$gene[1:28] <- c("MTARC1", "MARCHF1", "SEPTIN1", "DELEC1", "MTARC2", "MARCHF2", "SEPTIN2",
                                "MARCHF3", "SEPTIN3", "MARCHF4", "SEPTIN4", "MARCHF5", "SEPTIN5", "MARCHF6",
                                "SEPTIN6", "MARCHF7", "SEPTIN7", "MARCHF8", "SEPTIN8", "MARCHF9", "SEPTIN9",
                                "MARCHF10", "SEPTIN10", "MARCHF11", "SEPTIN11", "SEPTIN12", "SEPTIN14", "SEPTIN15")


# Get the mutation counts by gene from Singh et al. 2022
singh_SCZ_DN <- data.frame("gene" = singh_SCZ$Gene.Symbol, 
                           "dn.non" = as.integer(singh_SCZ$De.novo.PTV), 
                           "dn.mis" = as.integer(singh_SCZ$De.novo.mis3) + as.integer(singh_SCZ$De.novo.mis2),
                           "dn.frameshift" = 0,
                           "dn.splice" = 0,
                           "dn.syn" = 0)


# Get Tourette's data
wang_TD_DN <- data.frame("gene" = wang_TD$gene.id, 
                         "dn.non" = as.integer(wang_TD$dn.lof), 
                         "dn.mis" = as.integer(wang_TD$dn.mis3),
                         "dn.frameshift" = 0,
                         "dn.splice" = 0,
                         "dn.syn" = 0)


# Get OCD data
halvorsen_OCD_DN <- data.frame("gene" = c(), 
                               "dn.non" = c(), 
                               "dn.mis" = c(),
                               "dn.frameshift" = c(),
                               "dn.splice" = c(),
                               "dn.syn" = c())


for(i in 1:dim(halvorsen_OCD)[1]){
  
  if(halvorsen_OCD$EffectShort[i] == "non"){
    halvorsen_OCD_DN <- rbind(halvorsen_OCD_DN,
                                         data.frame("gene" = c(halvorsen_OCD$Gene[i]), 
                                                    "dn.non" = c(1), 
                                                    "dn.mis" = c(0),
                                                    "dn.frameshift" = c(0),
                                                    "dn.splice" = c(0),
                                                    "dn.syn" = c(0)))
  }
  
  if(halvorsen_OCD$EffectShort[i] == "misD"){
    halvorsen_OCD_DN <- rbind(halvorsen_OCD_DN,
                                         data.frame("gene" = c(halvorsen_OCD$Gene[i]), 
                                                    "dn.non" = c(0), 
                                                    "dn.mis" = c(1),
                                                    "dn.frameshift" = c(0),
                                                    "dn.splice" = c(0),
                                                    "dn.syn" = c(0))) 
  }
  
  if(halvorsen_OCD$EffectShort[i] == "frameshift"){
    halvorsen_OCD_DN <- rbind(halvorsen_OCD_DN,
                                         data.frame("gene" = c(halvorsen_OCD$Gene[i]), 
                                                    "dn.non" = c(0), 
                                                    "dn.mis" = c(0),
                                                    "dn.frameshift" = c(1),
                                                    "dn.splice" = c(0),
                                                    "dn.syn" = c(0)))
  }
  
  if(halvorsen_OCD$EffectShort[i] == "splice"){
    halvorsen_OCD_DN <- rbind(halvorsen_OCD_DN,
                                         data.frame("gene" = c(halvorsen_OCD$Gene[i]), 
                                                    "dn.non" = c(0), 
                                                    "dn.mis" = c(0),
                                                    "dn.frameshift" = c(0),
                                                    "dn.splice" = c(1),
                                                    "dn.syn" = c(0)))
  }
  
  if(halvorsen_OCD$EffectShort[i] == "syn"){
    halvorsen_OCD_DN <- rbind(halvorsen_OCD_DN,
                                         data.frame("gene" = c(halvorsen_OCD$Gene[i]), 
                                                    "dn.non" = c(0), 
                                                    "dn.mis" = c(0),
                                                    "dn.frameshift" = c(0),
                                                    "dn.splice" = c(0),
                                                    "dn.syn" = c(1)))
  }
}


# Now get the combined number of de novo mutations per gene in each dataset
combined_DN <- rbind(satterstrom_ASD_affected_DN, kaplanis_DDD_DN, wang_TD_DN, singh_SCZ_DN, halvorsen_OCD_DN) %>% 
  group_by(gene) %>% 
  summarise_all(funs(sum))
combined_DN <- data.frame(combined_DN)



### (5) Use denovolyzer to test for an excess of de novo mutations across all genes ####################################

# Convert dataframe of de novo mutations to the format required by denovolyzeR
combined_DeNovos <- data.frame("gene" = c(), "class" = c())
for (i in 1:dim(combined_DN)[1]){
  combined_DeNovos <- rbind(combined_DeNovos,
                            data.frame("gene" = rep(combined_DN$gene[i], combined_DN$dn.non[i]), 
                                       "class" = rep("non", combined_DN$dn.non[i])),
                            data.frame("gene" = rep(combined_DN$gene[i], combined_DN$dn.mis[i]), 
                                       "class" = rep("mis", combined_DN$dn.mis[i])),
                            data.frame("gene" = rep(combined_DN$gene[i], combined_DN$dn.frameshift[i]), 
                                       "class" = rep("frameshift", combined_DN$dn.frameshift[i])),
                            data.frame("gene" = rep(combined_DN$gene[i], combined_DN$dn.splice[i]), 
                                       "class" = rep("splice", combined_DN$dn.splice[i])),
                            data.frame("gene" = rep(combined_DN$gene[i], combined_DN$dn.syn[i]), 
                                       "class" = rep("syn", combined_DN$dn.syn[i])))
}

table(combined_DeNovos$class)



########################################################################################################################
#
#   NOTE: For the following section of code nsamples = 42,320 which is the sum of the 6,430 ASD cases included in the
#         Satterstrom et al. 2020 dataset, the 31,058 samples from the Kaplanis et al. 2020 dataset, the 628 samples 
#         from the Halvorsen et al. 2021 dataset, the 802 samples from the Wang et al 2018. dataset and the 3,402
#         samples from the Singh et al. 2021 dataset.
#
########################################################################################################################

# Use denovolyzeR to determine which classes possess an excess of de novo mutations
denovo_by_class <- denovolyzeByClass(genes = combined_DeNovos$gene,
                                     classes = combined_DeNovos$class,
                                     nsamples = 42320,
                                     roundExpected = 5,
                                     includeClasses = c("lof", "mis"))

# Use denovolyzeR to determine which genes possess an excess of de novo mutations
denovo_by_gene <- denovolyzeByGene(genes = combined_DeNovos$gene,
                                   classes = combined_DeNovos$class,
                                   nsamples = 42320,
                                   roundExpected = 20,
                                   includeClasses = c("lof", "mis"))

# Filter results to include the 466 candidate genes for CRISPRa cis-regulation therapy (CRT)
denovo_by_gene_CRT_candidates <- denovo_by_gene[which(denovo_by_gene$gene %in% candidate_genes),]

# replace incorrect gene names with HUGO gene symbols
denovo_by_gene_CRT_candidates$gene <- replace(denovo_by_gene_CRT_candidates$gene, which(denovo_by_gene_CRT_candidates$gene == "HIST1H1E"), "H1-4")
denovo_by_gene_CRT_candidates$gene <- replace(denovo_by_gene_CRT_candidates$gene, which(denovo_by_gene_CRT_candidates$gene == "H3F3A"), "H3-3A")
denovo_by_gene_CRT_candidates$gene <- replace(denovo_by_gene_CRT_candidates$gene, which(denovo_by_gene_CRT_candidates$gene == "COL4A3BP"), "CERT1")
denovo_by_gene_CRT_candidates$gene <- replace(denovo_by_gene_CRT_candidates$gene, which(denovo_by_gene_CRT_candidates$gene == "SRPR"), "SRPRA")

# Calculate the fold obs/exp de novo mutations for each CRISPRa CRT candidate gene
denovo_by_gene_CRT_candidates$fold_mis <- denovo_by_gene_CRT_candidates$mis_observed / denovo_by_gene_CRT_candidates$mis_expected
denovo_by_gene_CRT_candidates$fold_lof <- denovo_by_gene_CRT_candidates$lof_observed / denovo_by_gene_CRT_candidates$lof_expected

# Remove all genes with fold_lof <= 1, these genes will be excluded from downstream experiments
denovo_by_gene_CRT_candidates_FINAL <- denovo_by_gene_CRT_candidates[which(denovo_by_gene_CRT_candidates$fold_lof > 1),]



### (6) Visualize results and output final gene list ###################################################################

# Create combined gene list from each prioritization dataset
combined_gene_list <- unique(c(decipher_DDG2P_brain_het_lof_gene_list,
                               satterstrom_ASD_gene_list,
                               kaplanis_DDD_gene_list,
                               singh_SCZ_gene_list,
                               wang_TD_gene_list))
length(combined_gene_list)

# Convert combined gene lists into ggplot2 format
combined_table <- data.frame("dataset" = "Combined",
                             "genes" = denovo_by_gene_CRT_candidates_FINAL$gene)
decipher_DDG2P_brain_table <- data.frame("dataset" = "DDG2P_Brain_Het_LoF",
                                         "genes" = decipher_DDG2P_brain_het_lof_gene_list)
satterstrom_ASD_table <- data.frame("dataset" = "Satterstrom_ASD_102",
                                    "genes" = satterstrom_ASD_gene_list)
kaplanis_DDD_table <- data.frame("dataset" = "Kaplanis_DDD_285",
                                 "genes" = kaplanis_DDD_gene_list)
singh_SCZ_table <- data.frame("dataset" = "Singh_SCZ_32",
                              "genes" = singh_SCZ_gene_list)
wang_TD_table <- data.frame("dataset" = "Wang_TD_2",
                            "genes" = wang_TD_gene_list)

# Convert combined gene lists into ggplot2 format
to_plot <- rbind(combined_table,
                 decipher_DDG2P_brain_table,
                 satterstrom_ASD_table,
                 kaplanis_DDD_table,
                 singh_SCZ_table,
                 wang_TD_table)
to_plot$CRT_Candidates <- ifelse(to_plot$genes %in% denovo_by_gene_CRT_candidates_FINAL$gene, "CRT_Candidates", "other")

# Merge counts for combined gene lists into ggplot2 format
to_plot <- to_plot %>% 
  group_by(CRT_Candidates, dataset) %>% 
  summarise(n = n())
to_plot <- data.frame(to_plot)
to_plot <- rbind(to_plot,
                 data.frame("CRT_Candidates" = "other",
                            "dataset" = "Combined",
                            "n" = length(combined_gene_list) - to_plot$n[which(to_plot$CRT_Candidates == "CRT_Candidates" & to_plot$dataset == "Combined")]))

# Order columns for plotting
positions <- c("Combined", "Kaplanis_DDD_285", "DDG2P_Brain_Het_LoF", "Satterstrom_ASD_102", "Singh_SCZ_32", "Wang_TD_2")
to_plot$CRT_Candidates <- factor(to_plot$CRT_Candidates, levels = c("other", "CRT_Candidates"))

# Create plot
ggplot(data = to_plot, aes(x = dataset, y = n, fill = CRT_Candidates)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = positions) +
  scale_fill_manual(values=c("#D1D3D4", "#27AAE1")) +
  geom_text(aes(label = n)) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  ylab("# genes\n") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12))



# Make another version of this plot for the project schematic
to_plot_CRT_only <- to_plot[which(to_plot$CRT_Candidates == "CRT_Candidates"),]
to_plot_CRT_only$dataset <- c("Combined", "DDG2P", "DDD", "ASD", "SCZ", "TD")

# Order columns for plotting
positions <- c("Combined", "DDD", "DDG2P", "ASD", "SCZ", "TD")
to_plot$CRT_Candidates <- factor(to_plot$CRT_Candidates, levels = c("other", "CRT_Candidates"))

ggplot(data = to_plot_CRT_only, aes(x = dataset, y = n, fill = CRT_Candidates)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = positions) +
  scale_fill_manual(values=c("#27AAE1")) +
  geom_text(aes(label = n)) + 
  ylab("# genes\n") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        legend.text = element_blank())


# Get the upper and lower confidence intervals for denovo_by_class
denovo_by_class$lower_conf <- denovo_by_class$enrichment - (1.96 * (denovo_by_class$enrichment / sqrt(denovo_by_class$observed)))
denovo_by_class$upper_conf <- denovo_by_class$enrichment + (1.96 * (denovo_by_class$enrichment / sqrt(denovo_by_class$observed)))

# Set significance value
sig <- 1

# Plot the de novo observed vs expected enrichment by class
p <- ggplot(denovo_by_class, aes(x = class, y = enrichment)) + 
  geom_bar(stat = "identity", color = "#27AAE1", fill = "#27AAE1", position = position_dodge()) +
  geom_hline(yintercept = sig, color = "red", linetype = "dashed") + 
  geom_errorbar(aes(ymin = lower_conf, ymax = upper_conf), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label = c("****", "****"), y = upper_conf), vjust = -0.8, color = "black", size = 8) +
  #ggtitle("De novo variant meta-analysis from 42,320 cases\n") +
  xlab("") + ylab("Observed vs Expected\n") + ylim(0, 3) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
        axis.text.y = element_text(size = 24),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))
print(p)


# Calculate the fold obs/exp de novo mutations for each CRISPRa CRT candidate gene for manhattan plots
denovo_by_gene_to_plot <- denovo_by_gene_CRT_candidates

head(denovo_by_gene_to_plot)
dim(denovo_by_gene_to_plot)

# Set genes with a -log10(p-value) > Inf to 1e-300
denovo_by_gene_to_plot$lof_pValue <- ifelse(denovo_by_gene_to_plot$lof_pValue == 0, 1e-300, denovo_by_gene_to_plot$lof_pValue)

# Get summary statistics for manhattan plot
denovo_by_gene_to_plot$chr <- annotations$chrom[match(denovo_by_gene_to_plot$gene, annotations$gene_name)]
denovo_by_gene_to_plot$bp <- annotations$start_position_hg38[match(denovo_by_gene_to_plot$gene, annotations$gene_name)]
denovo_by_gene_to_plot <- denovo_by_gene_to_plot[which(denovo_by_gene_to_plot$chr != "<NA>"),]
denovo_by_gene_to_plot$chr <- factor(denovo_by_gene_to_plot$chr, levels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                                                            "chr6", "chr7", "chr8", "chr9", "chr10",
                                                                            "chr11", "chr12", "chr13", "chr14", "chr15",
                                                                            "chr16", "chr17", "chr18", "chr19", "chr20",
                                                                            "chr21", "chr22", "chrX", "chrY"),
                                     labels = c("chr1", "chr2", "chr3", "chr4", "chr5",
                                                "chr6", "chr7", "chr8", "chr9", "chr10",
                                                "chr11", "chr12", "chr13", "chr14", "chr15",
                                                "chr16", "chr17", "chr18", "chr19", "chr20",
                                                "chr21", "chr22", "chrX", "chrY"))

# Transform data to create manhattan plot
data_cum <- denovo_by_gene_to_plot %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  dplyr::select(chr, bp_add)

# Transform data to create manhattan plot
data_cum <- data.frame(data_cum)

# Transform data to create manhattan plot
denovo_by_gene_to_plot <- denovo_by_gene_to_plot %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

# Transform data to create manhattan plot
axis_set <- denovo_by_gene_to_plot %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

# Transform data to create manhattan plot
axis_set <- data.frame(axis_set)

# Set significance threshold
sig <- 1

# Plot manhattan plot
manhplot <- ggplot(denovo_by_gene_to_plot, aes(x = bp_cum, y = -log(lof_pValue), 
                                  color = as_factor(chr), size = -log(lof_pValue))) +
  geom_hline(yintercept = 645, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 710)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  geom_text_repel(data = denovo_by_gene_to_plot[which(-log(denovo_by_gene_to_plot$lof_pValue) > 100),], aes(label = denovo_by_gene_to_plot$gene[which(-log(denovo_by_gene_to_plot$lof_pValue) > 100)]), size = 3.5) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p-value)\n") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

manhplot


# Calculate the fold obs/exp de novo mutations for each CRISPRa CRT candidate gene for manhattan plots
denovo_by_gene_to_plot <- denovo_by_gene_CRT_candidates

head(denovo_by_gene_to_plot)
dim(denovo_by_gene_to_plot)

# Set genes with a fold > Inf to 100
denovo_by_gene_to_plot$fold_lof <- ifelse(denovo_by_gene_to_plot$fold_lof > 50, 50, denovo_by_gene_to_plot$fold_lof)
denovo_by_gene_to_plot$fold_mis <- ifelse(denovo_by_gene_to_plot$fold_mis > 50, 50, denovo_by_gene_to_plot$fold_mis)

# Make plots comparing missense and loss of function observed vs expected ratios
denovo_by_gene_to_plot_FINAL <- denovo_by_gene_to_plot[which(denovo_by_gene_to_plot$gene %in% combined_gene_list),]
denovo_by_gene_to_plot_FINAL$color <- ifelse(denovo_by_gene_to_plot_FINAL$fold_lof < 1, 
                                             ifelse(denovo_by_gene_to_plot_FINAL$fold_mis > 1, "red", "#D1D3D4"), 
                                             "#27AAE1")
target_list <- c("PTPN11", "SCN2A")

# Create fold plot for main figures
fold_plot_1 <- ggplot(denovo_by_gene_to_plot_FINAL, aes(x = fold_lof, y = fold_mis, color = color)) +
  geom_point() + 
  #geom_text_repel(denovo_by_gene_CRT_candidates_TO_PLOT[1:10,], mapping = aes(label = denovo_by_gene_CRT_candidates_TO_PLOT$gene[1:10])) +
  scale_color_manual(values = c("#27AAE1", "#D1D3D4", "red")) +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed") + 
  geom_vline(xintercept = 49, color = "#D1D3D4", linetype = "dashed") + 
  geom_text_repel(data = denovo_by_gene_to_plot_FINAL[which(denovo_by_gene_to_plot_FINAL$gene %in% target_list),], aes(label = denovo_by_gene_to_plot_FINAL$gene[which(denovo_by_gene_to_plot_FINAL$gene %in% target_list)]), color = "black", size = 4) +
  #ylim(0, 50) +
  #xlim(0, 50) +
  theme_classic() +
  xlab("\nLoF (Observed/Expected)") + 
  ylab("Mis (Observed/Expected)\n") +
  theme(plot.title = element_blank(), 
        axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 2, b = 0, l = 0)),
        axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.text=element_text(size = 8),
        legend.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(size = 10),
        legend.spacing.y = unit(1, "mm"),
        legend.box.margin = margin(-10, -5, -10, -10),
        plot.margin = margin(15, 5, 15, 15))
plot(fold_plot_1)


# Get BrainVar gene expression data and make a figure
brainvar <- data.frame(read_xlsx(path = "./data/published_datasets/1-s2.0-S2211124720303673-mmc3.xlsx", sheet = "ExpressedGenes_23782"))
brainvar$CRT_candidates <- ifelse(brainvar$GeneSymbol %in% denovo_by_gene_CRT_candidates_FINAL$gene, "CRT_candidates", "other_genes")
brainvar$CRT_candidates <- factor(brainvar$CRT_candidates, levels = c("other_genes", "CRT_candidates"))

# Create plot for brainvar data
p <- ggplot(brainvar, aes(x = CRT_candidates, y = MedianExpr, fill = CRT_candidates)) + 
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values=c("#D1D3D4", "#27AAE1")) +
  ylab("Median Expr Human Brain\n") +
  ylim(c(-7, 15)) + 
  theme_classic() +
  theme(axis.title.x = element_blank(),legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        legend.position = "none")
p


# Get Gnomad LOEUF data for CRT candidates and make a figure
gnomad <- read.table(file = "./data/published_datasets/gnomad.v4.1.constraint_metrics.tsv", header = TRUE, sep = "\t")
gnomad_mane <- gnomad[which(gnomad$mane_select == "true" & gnomad$transcript_type == "protein_coding"),]
gnomad_mane$CRT_candidates <- ifelse(gnomad_mane$gene %in% denovo_by_gene_CRT_candidates_FINAL$gene, "CRT_candidates", "other_genes")
gnomad_mane$CRT_candidates <- factor(gnomad_mane$CRT_candidates, levels = c("other_genes", "CRT_candidates"))

# Create plot for Gnomad LOEUF data
p <- ggplot(gnomad_mane, aes(x = CRT_candidates, y = lof.oe_ci.upper, fill = CRT_candidates)) + 
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values=c("#D1D3D4", "#27AAE1")) +
  ylab("LOEUF Score\n") +
  ylim(c(2, 0)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        legend.position = "none")
p


### Make final supplementary table_s1 ###

# Create metadata table for the 337 genes moving on to downstream analysis
denovo_by_gene_CRT_candidates_FINAL <- denovo_by_gene_CRT_candidates_FINAL[order(denovo_by_gene_CRT_candidates_FINAL$gene, decreasing = FALSE),]
CRT_candidates_FINAL_metadata <- data.frame("gene_symbol" = denovo_by_gene_CRT_candidates_FINAL$gene,
                                            "ensembl_id" = annotations$gene_id[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "hgnc_id" = annotations$hgnc_id[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "chrom" = annotations$chrom[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "start_hg38" = annotations$start_position_hg38[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "end_hg38" = annotations$end_position_hg38[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "strand" = annotations$strand[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "transcript_id" = annotations$transcript_ids[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "transcript_length" = annotations$transcript_lengths[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "cds_length" = annotations$cds_lengths[match(denovo_by_gene_CRT_candidates_FINAL$gene, annotations$gene_name)],
                                            "median_brain_expression" = brainvar$MedianExpr[match(denovo_by_gene_CRT_candidates_FINAL$gene, brainvar$GeneSymbol)],
                                            "midfetal_brain_expression" = brainvar$MedianExpr_E1[match(denovo_by_gene_CRT_candidates_FINAL$gene, brainvar$GeneSymbol)],
                                            "postnatal_brain_expression" = brainvar$MedianExpr_E3[match(denovo_by_gene_CRT_candidates_FINAL$gene, brainvar$GeneSymbol)],
                                            "loeuf" = gnomad_mane$lof.oe_ci.upper[match(denovo_by_gene_CRT_candidates_FINAL$gene, gnomad_mane$gene)],
                                            "dn.lof.observed" = denovo_by_gene_CRT_candidates_FINAL$lof_observed,
                                            "dn.lof.expected" = denovo_by_gene_CRT_candidates_FINAL$lof_expected,
                                            "dn.lof.pval" = denovo_by_gene_CRT_candidates_FINAL$lof_pValue,
                                            "dn.mis.observed" = denovo_by_gene_CRT_candidates_FINAL$mis_observed,
                                            "dn.mis.expected" = denovo_by_gene_CRT_candidates_FINAL$mis_expected,
                                            "dn.mis.pval" = denovo_by_gene_CRT_candidates_FINAL$mis_pValue,
                                            "Kaplanis_DDD_285" = ifelse(denovo_by_gene_CRT_candidates_FINAL$gene %in% kaplanis_DDD_gene_list, "True", "False"),
                                            "Satterstrom_ASD_102" = ifelse(denovo_by_gene_CRT_candidates_FINAL$gene %in% satterstrom_ASD_gene_list, "True", "False"),
                                            "Singh_SCZ_32" = ifelse(denovo_by_gene_CRT_candidates_FINAL$gene %in% singh_SCZ_gene_list, "True", "False"),
                                            "Wang_TD_2" = ifelse(denovo_by_gene_CRT_candidates_FINAL$gene %in% wang_TD_gene_list, "True", "False"),
                                            "DDG2P_Brain_Het_Lof" = ifelse(denovo_by_gene_CRT_candidates_FINAL$gene %in% decipher_DDG2P_brain_het_lof_gene_list, "True", "False"))

write.xlsx(CRT_candidates_FINAL_metadata, file = "./tables/Table_S1_Cis_Regulation_Therapy_Candidate_Genes.xlsx", sheetName = "Candidate_Genes", 
           col.names = TRUE, row.names = FALSE, append = FALSE)



###################################################### END SCRIPT ######################################################




