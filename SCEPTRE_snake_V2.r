##Install and load SCEPTRE package
#devtools::install_github("katsevich-lab/sceptre", lib = "/net/shendure/vol10/projects/troym/ResQTL/nobackup/sceptre/")

library(sceptre, lib.loc = "/net/shendure/vol10/projects/troym/ResQTL/nobackup/sceptre/")
library(tidyverse); library(data.table)

args <- commandArgs(trailingOnly = TRUE)
part <- args[1]; # part <- 1

# Filter data for corresponding partition
gRNA_target_data_frame <- readRDS(paste0("/net/shendure/vol10/projects/troym/ResQTL/nobackup/scaled_screen/gRNA_mapping/gRNA_target_data_frame_partitioned.Rds")) %>%
    filter(partition==part | grna_target == 'non-targeting') %>% 
    select(-partition)

# Filter based on partition in target df
discovery_pairs <- fread(paste0("/net/shendure/vol10/projects/troym/ResQTL/nobackup/scaled_screen/gRNA_mapping/Predicted_Target_SCEPTRE_V09_Discovery_Pairs.csv")) %>% 
    filter(grna_target %in% gRNA_target_data_frame$grna_target)

##Load the aggregated data

gRNA_matrix <- read_rds("/net/shendure/vol10/projects/troym/ResQTL/nobackup/scaled_screen/SCEPTRE_Processed_Matrices/Scaled_Screen_gRNA_Matrix.mtx")
gene_matrix <- read_rds("/net/shendure/vol10/projects/troym/ResQTL/nobackup/scaled_screen/SCEPTRE_Processed_Matrices/Scaled_Screen_gene_Matrix.mtx")
covariate_matrix <- read_rds("/net/shendure/vol10/projects/troym/ResQTL/nobackup/scaled_screen/SCEPTRE_Processed_Matrices/Scaled_Screen_covariate_Matrix.mtx")

##Select the extra covariates 

Extra_Covariates <-  as.data.frame(covariate_matrix)

Extra_Covariates <- Extra_Covariates %>% 
  select(p_mito, lane)

##create sceptre object

sceptre_object <- import_data(
  response_matrix = gene_matrix,
  grna_matrix = gRNA_matrix,
  grna_target_data_frame = gRNA_target_data_frame,
  extra_covariates = Extra_Covariates,
  moi = "high",
)

##Filter the list to only test genes that are detected
Detected_Genes <- rownames(sceptre_object@response_matrix) 
discovery_pairs <- discovery_pairs %>%
    filter(response_id %in% Detected_Genes)

# make version for testing a subset of gRNAs    
nrespmat <- sceptre_object@response_matrix[unique(discovery_pairs$response_id),]
sot <- import_data(response_matrix=nrespmat,
    grna_matrix=sceptre_object@grna_matrix[gRNA_target_data_frame$grna_id,],
    grna_target_data_frame=gRNA_target_data_frame,
    moi='high', response_names=rownames(nrespmat))
sot@covariate_data_frame <- sceptre_object@covariate_data_frame
rm(sceptre_object); gc()

# apply the pipeline functions to the sceptre_object in order
sot <- sot |> # |> is R's base pipe, similar to %>%
  set_analysis_parameters(discovery_pairs, side = "right")

sot <- sot |> 
  assign_grnas(sceptre_object = sot, method = "thresholding", threshold = 5) |> 
  run_qc(p_mito_threshold = 0, response_n_nonzero_range = c(0, 1), response_n_umis_range = c(0, 1))

sot <- sot |> 
  run_calibration_check() |>
  run_discovery_analysis()

# write the results to disk
write_outputs_to_directory(sot, paste0("/net/shendure/vol10/projects/troym/ResQTL/nobackup/",
    "scaled_screen/SCEPTRE_PT/outputs/sceptre_parallel_part", part))
 

