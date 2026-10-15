
############################################################
# Split SingleCellExperiment by cell type and cohort
# Purpose:
#   Prepare per–cell-type, per–cohort SCE objects for downstream
#   DE analyses (MAST-RE/edgeR/pseudobulk, etc), followed by
#   meta-analysis across cohorts.
#
# Data source:
#   snRNA-seq dataset from Synapse (Sage Bionetworks):
#   https://www.synapse.org/Synapse:syn51535901
#   ~40 GB, 16,798 genes x 491,427 nuclei
#
# Notes:
#   - Output objects use the SingleCellExperiment format, including assays,
#     colData and reducedDims.
#   - Split by Celltype_inferred and cohort (McLean, MtSinai) for
#     cohort-specific differential expression and subsequent meta-analysis.

# Inputs:
#   ./data/sce_snBD_McLean_MtSinai.rds (not bundled).
# Outputs:
#   ./data/cell_type_data_McLean/snBD_sce_format_<cell_type>.rds
#   ./data/cell_type_data_Sinai/snBD_sce_format_<cell_type>.rds
# Dependencies:
#   data.table, SingleCellExperiment, readr.
# Usage:
#   From repository root: Rscript differential_expression/create_cell_type_file.R
#
# Date: 2023-3
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.5 LTS
############################################################


library(data.table)
library(SingleCellExperiment)

dir_mclean <- "./data/cell_type_data_McLean"
dir_sinai  <- "./data/cell_type_data_Sinai"
if (!dir.exists(dir_mclean)) dir.create(dir_mclean, recursive = TRUE)
if (!dir.exists(dir_sinai))  dir.create(dir_sinai,  recursive = TRUE)


# import data: rows (gene): 16798; column (cells): 491427
sce = readr::read_rds("./data/sce_snBD_McLean_MtSinai.rds")

table(sce$Celltype, useNA = "ifany")
unique(sce$Celltype_inferred)


# Split by cell type × cohort
for(cell_type_i in unique(sce$Celltype_inferred)) {
  print(cell_type_i)
  
  sce_one <- sce[, sce$Celltype_inferred == cell_type_i]
  for (cohort in c("McLean", "MtSinai")) {
    if(cohort == "McLean") {
      file = paste0("./data/cell_type_data_McLean/snBD_sce_format_" ,cell_type_i,".rds")
      file
      sce_s <- sce_one[, sce_one$Cohort == "McLean"]
      table(sce_s$Phenotype)
      readr::write_rds(sce_s, file)
    }
    
    if(cohort == "MtSinai") {
      file = paste0("./data/cell_type_data_Sinai/snBD_sce_format_" ,cell_type_i,".rds")
      file
      sce_s <- sce_one[, sce_one$Cohort == "MtSinai"]
      table(sce_s$Phenotype)
      readr::write_rds(sce_s, file)
    }
  }
}


message("\nSession info:")
print(sessionInfo())
