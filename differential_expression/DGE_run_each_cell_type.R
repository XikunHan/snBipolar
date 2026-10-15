#!/usr/bin/env Rscript

###############################################################################
# snBD — Differential expression by cell type and cohort
#
# Purpose:
#   Run differential expression for one cell type and cohort using MAST-RE,
#   NEBULA or edgeR, with three RUVr factors and cell-level QC covariates.
#   MAST-RE jobs can process a selected gene block.
#
# Inputs (produced by split script):
#   McLean: ./data/cell_type_data_McLean/snBD_sce_format_<celltype>.rds
#   Sinai : ./data/cell_type_data_Sinai/snBD_sce_format_<celltype>.rds
#
# Outputs:
#   <sub_outpath>/<DGE>_DGE__<trait>__<celltype>[__<block>].csv
#   (Plus job logs written by SLURM launcher)
#
# Command-line arguments (set by the launcher):
#   --trait        BD
#   --cell_type    Harmonized cell type name (Celltype_inferred)
#   --method       Tag for outputs (e.g., "RUV")
#   --DGE          One of: edger | nebula | MASTlmmBlock
#   --cohort       McLean | Sinai  (maps to the corresponding file tree)
#   --sub_outpath  Output directory (created if missing)
#   --block        Integer block index (MASTlmmBlock), or NA otherwise
#
# Model & covariates:
#   - RUVr factors K=3 (the helper calls RUVSeq::RUVr).
#   - Per-cell QC covariates: ngenes_z, Mito_z, Ribo_z.
#   - Cohort-specific covariates:
#       McLean: age, sex, PMI.z, [RUV W’s]
#       Sinai : Age, Ethnicity, Sex, PMI, [RUV W’s]
#   - Pseudobulk (edgeR): donor-level aggregation by SubID.
#
# Dependencies:
#   sc_helper.R in this folder; data.table, optparse, SingleCellExperiment,
#   scuttle, scater, stringr, RUVSeq and the selected model packages.
# Usage:
#   From repository root: bash differential_expression/DGE_run_submit.sh
#
# Date: 2023-3
# R version: 4.2.0 (2022-04-22) on x86_64-pc-linux-gnu
###############################################################################


# Record the repository root before changing the working directory.
project_root <- normalizePath(Sys.getenv("SNBD_ROOT", unset = "."))

# 1. Packages and command-line arguments -------------------------------
library(data.table)
library(optparse)
library(SingleCellExperiment)


option_list <- list(
  make_option(c("--trait"),       type="character"),
  make_option(c("--cell_type"),   type="character"),
  make_option(c("--method"),      type="character"),
  make_option(c("--DGE"),         type="character"),
  make_option(c("--cohort"),      type="character"),
  make_option(c("--sub_outpath"), type="character"),
  make_option(c("--block"),       type="character", default=NA)
)


opt <- parse_args(OptionParser(option_list=option_list))
print(opt)


trait       <- opt$trait
method      <- opt$method
DGE         <- opt$DGE
cohort  <- opt$cohort      # "McLean" or "Sinai" (CLI)
cell_type   <- opt$cell_type
sub_outpath <- opt$sub_outpath
block       <- suppressWarnings(as.integer(opt$block))
stopifnot(length(cohort) == 1L, cohort %in% c("McLean", "Sinai"),
          length(trait) == 1L, trait == "BD",
          length(DGE) == 1L, DGE %in% c("edger", "nebula", "MASTlmmBlock"))

if (DGE == "MASTlmmBlock") {
  if (length(block) != 1L || is.na(block) || block < 1L)
    stop("Provide a positive --block index.")
  gene_list_file <- file.path(project_root, "out/MASTlmm/data",
                             sprintf("gene_list_McLean_BD_Ex-L23_b%s.csv", block))
  if (!file.exists(gene_list_file)) stop("Missing shared gene list: ", gene_list_file)
}


source(file.path(project_root, "differential_expression", "sc_helper.R"))

dir.create(sub_outpath, showWarnings = FALSE, recursive = TRUE)
setwd(sub_outpath)



if(0) {
  block <- 1
  trait <- "BD"
  cell_type <- "Endo"
  method <- "RUV"
  cohort <- "McLean"
  DGE <- "edger"
  sub_outpath <- paste0(project_root, "/out/DGE/test/edger_McLean_RUV")
  
}


# 2. Read the selected cohort/cell type --------------------------------
if(cohort == "McLean") {
    file = paste0(project_root, "/data/cell_type_data_McLean/snBD_sce_format_", cell_type, ".rds")
    } else if(cohort == "Sinai") {
    file = paste0(project_root, "/data/cell_type_data_Sinai/snBD_sce_format_", cell_type, ".rds")
}


sce <- readRDS(file = file)
sce

v_case <- trait
v_control <- 'CON'


# use 3 RUV
factorK <- 3
print(sprintf("factorK: %s", factorK))


# Cohort-specific covariates
if(cohort == "McLean") {
  v_cov <- c("age", "sex", "PMI.z", "ngenes_z", "Mito_z", "Ribo_z", paste0("RUVW_", 1:factorK))
  } else if (cohort == "Sinai") {
    v_cov <- c("Age", "Ethnicity", "Sex", "PMI", "ngenes_z", "Mito_z", "Ribo_z", paste0("RUVW_", 1:factorK))
  }
  

message("Covariates requested: ", paste(v_cov, collapse=", "))



# for MAST-RE create gene block, n=250 genes
# 3. Prepare 250-gene blocks -------------------------------------------
requested_block <- block
# run once
if(0) {
  
  # total gene = 16798, BD for McLean and Sinai
  sce
  nrow(sce)
  
  size <- 250
  nrow(sce)/size
  
  n_block <- ceiling(nrow(sce)/size)
  n_block
  
  # block <- 1:n_block
  # block
  
  block=1
  start <- ((block-1)*size + 1)
  end <- ((block)*size)
  
  gene_list <- rownames(sce)[na.omit((1:nrow(sce))[start:end])]
  gene_list
  
  for(block in 1:n_block){
    print(block)
    
    start <- ((block-1)*size + 1)
    end <- ((block)*size)
    start
    end 
    gene_list <- rownames(sce)[na.omit((1:nrow(sce))[start:end])]
    gene_list
    
    fwrite(list(gene_list), sprintf("gene_list_%s_%s_%s_b%s.csv", cohort, trait, cell_type, block), col.names = FALSE)
  }
}



block <- requested_block

# select cell type
# 4. Select cases/controls and calculate cell QC -----------------------
df <- sce[, sce$Celltype_inferred == cell_type]
df$outcome <- ifelse(df@colData$Phenotype == "BD", trait, ifelse(df@colData$Phenotype == "CON", "CON", NA))
df$outcome <- as.character(df$outcome)
df$outcome <- factor(df$outcome)
df$outcome <- relevel(df$outcome, ref = "CON")
table(df$outcome, useNA = "ifany")


# select outcome 
df <- df[, !is.na(df$outcome)]

df


## add MT and ribosomal genes
is.mito <- grep("MT-", rownames(df))
table(is.mito)
is.ribo <- grep("RPL|RPS", rownames(df))
table(is.ribo)
df <- scuttle::addPerCellQC(df, subsets=list(Mito=is.mito, Ribo=is.ribo))
df$SampleQuality <- c(scale(log10(df$sum)))
df$ngenes_z <- c(scale(df$detected))
df$Mito_z <- c(scale(df$subsets_Mito_percent))
df$Ribo_z <- c(scale(df$subsets_Ribo_percent))
df$SubID <- stringr::str_remove_all(df$SubID, "_")


if(cohort == "McLean") {
  df$Batch <- factor(as.character(df$Batch))
  df$HTO <- factor(as.character(df$HTO))
}


# MAST-RE specific filtering (cell-level models)
if(DGE == "MASTlmm" || DGE == "MASTlmmBlock" ) {
  
  # 1) gene in more than 10 cells 
  d <- scater::nexprs(df, byrow=TRUE)
  df <- df[d > 10, ]
  df
  
  # 2) more than 1% of cells, in both cases and controls
  is.case <- df$outcome == trait
  table(is.case)
  is.control <- df$outcome !=  trait
  table(is.control)
  
  d <- scater::perFeatureQCMetrics(df, subsets=list(Case = is.case, Control = is.control))
  df <- df[d$subsets_Case_detected > 1 & d$subsets_Control_detected > 1, ]
  df
  
  # 3) keep cells with > 1000 detected genes
  d <- scuttle::perCellQCMetrics(df)
  df <- df[, d$detected > 1000 ]
}




# add RUV 
# 5. Estimate residual-based RUVr covariates ---------------------------
library(RUVSeq)

set.seed(123)
df = deg.ruvseq(df,
  sample="SubID",
  pathology="outcome",
  NRUV= factorK)



# Pseudobulk differential expression ##############

## edger ####
# 6. Run the selected differential-expression method -------------------
if(DGE %in% c("edger"))  {
  df_PB = se_make_pseudobulk(se = df, on = "SubID") 
  
  dd <- as.data.table(df@colData)
  
  dd <- dd[, lapply(.SD, mean), by= SubID, .SDcols = c("SampleQuality",'ngenes_z', "Mito_z", "Ribo_z")]
  
  dd <- dd[match(as.character(df_PB$SubID), as.character(dd$SubID))]
  stopifnot(identical(as.character(df_PB$SubID), as.character(dd$SubID)))
  for (cn in setdiff(names(dd), "SubID")) {
    SummarizedExperiment::colData(df_PB)[[cn]] <- dd[[cn]]
  }
  
   if(DGE == "edger") {
    
    print(sprintf("Runing: %s", DGE))
    
    res <- deg.edger(df_PB, 
      pathology = "outcome",
      case = v_case,
      control = v_control,
      covariates = v_cov)
  }
  
  
  
  print(table(res$FDR < 0.05))
  print(table(res$FDR < 0.05, abs(res$log2FC) > 0.2))
  
  file_name <- sprintf("%s_DGE__%s__%s.csv", DGE, trait, cell_type)
  file_name
  fwrite(res, file_name)
  

} else if(DGE == "nebula") {
  ## nebula ####
  
  print(v_cov)
  print(df)
  
  print(sprintf("Runing: %s", DGE))
  
  res <- deg.nebula(df, 
    sample = "SubID", 
    pathology = "outcome",
    case = v_case,
    control = v_control,
    covariates = v_cov,
    factorize_pathology = FALSE,
    NRUV= 0) # RUVs were provided as covariates
  
  
  head(res)
  
  print(table(res$FDR < 0.05))
  print(table(res$FDR < 0.05, abs(res$log2FC) > 0.2))
  
  file_name <- sprintf("%s_DGE__%s__%s.csv", DGE, trait, cell_type)
  file_name
  fwrite(res, file_name)


  
} else if(DGE == "MASTlmmBlock") {
  ## MASTlmm block ####
  
  print(v_cov)
  
  print(sprintf("Runing: %s %s", DGE, block))
  
  # the gene list is the same for all cell types
  gene_list <- fread(gene_list_file, header = FALSE)
  if (ncol(gene_list) != 1L || !nrow(gene_list) || anyNA(gene_list$V1) || anyDuplicated(gene_list$V1))
    stop("Shared gene list must contain one headerless column of unique gene IDs.")
  if (any(!gene_list$V1 %in% rownames(sce)))
    stop("Shared gene list contains IDs absent from this SCE; verify its header and reference gene order.")
  
  print(gene_list)
  
  res <- deg.MAST(df,
                  pathology = "outcome",
                  sample = "SubID",
                  case = v_case,
                  control = v_control,
                  covariates = v_cov,
                  detection_rate = FALSE,
                  method = "glmer", 
                  ebayes = F, 
                  gene_list = gene_list$V1,
                  filter_only_case_control=FALSE)
  
  
  print(table(res$FDR < 0.05))
  print(table(res$FDR < 0.05, abs(res$logFC) > 0.2))
  
  file_name <- sprintf("%s_DGE__%s__%s__%s.csv", DGE, trait, cell_type, block)
  file_name
  fwrite(res, file_name)
  
}

print(sessionInfo())
print(date())
