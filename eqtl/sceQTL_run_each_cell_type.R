#!/usr/bin/env Rscript

###############################################################################
# snBD — Cell-type-specific cis-eQTL mapping
#
# Purpose:
#   Run the linear MatrixEQTL analysis for one cell type using
#   prepared genotype, expression and covariate matrices.
#
# Inputs:
#   --cell_type; --sub_outpath, normally out/sceQTL430/result.
#   out/sceQTL430/data/QTL/: df_genotype_<cell_type>.bed/.bim/.fam,
#   df__<cell_type>.txt, df__<cell_type>_cov.txt, SNP_location.txt, gene_location.txt.
#   Matrix files use IDs in column 1 and matched donor columns thereafter.
#
# Outputs:
#   genotype_<cell_type>.txt in the QTL input folder;
#   out__<cell_type>__cis_eqtls.txt.gz in the result folder.
#
# Method overview:
#   Decode PLINK genotypes with genio; load SlicedData objects; fit modelLINEAR
#   with cisDist=1e6, cis P threshold=1, trans P threshold=0; attach BIM fields.
#
# Dependencies:
#   data.table, optparse, MatrixEQTL, genio.
#
# Usage:
#   From repository root: bash eqtl/sceQTL_run_submit.sh
#
#
###############################################################################

if(0) {
  devtools::install_github("andreyshabalin/MatrixEQTL")
  devtools::install_github("OchoaLab/genio", build_vignettes = TRUE)
}


# Record the repository root before changing the working directory.
project_root <- normalizePath(Sys.getenv("SNBD_ROOT", unset = "."))
qtl_dir <- file.path(project_root, "out/sceQTL430/data/QTL")

library(data.table)
library(optparse)


# 1. Arguments ---------------------------------------------------------
option_list <- list(
    make_option(c("--cell_type"), type="character", default=NULL, help="cell_type"),
    make_option(c("--sub_outpath"), type="character", default=NULL, help="sub_outpath")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)


cell_type <- opt$cell_type
cell_type

sub_outpath <- opt$sub_outpath

# Example arguments for an interactive test 
if(0) {
    cell_type <- "Oli"
    cell_type <- "Endo"
    sub_outpath <- paste0(project_root, "/out/sceQTL430/result")
    }

working_dir <- sub_outpath
dir.create(working_dir, recursive = TRUE, showWarnings = FALSE)
setwd(working_dir)



# run QTL ####

if(1) {

# pipeline
library("MatrixEQTL")

useModel = modelLINEAR
library(genio)


# 2. Read and export genotype matrix -----------------------------------
df_genotype <- read_plink(sprintf(paste0(project_root, "/out/sceQTL430/data/QTL/df_genotype_%s"), cell_type))
df_genotype <- df_genotype$X


df_genotype <- as.data.table(df_genotype, keep.rownames = TRUE)
setnames(df_genotype, "rn", "id")
head(df_genotype[, seq_len(min(10L, ncol(df_genotype))), with = FALSE])

fwrite(df_genotype, file.path(qtl_dir, sprintf("genotype_%s.txt", cell_type)), sep = "\t")

rm(df_genotype); gc()


# 3. Load genotype, expression and covariate slices --------------------
SNP_file_name = file.path(qtl_dir, sprintf("genotype_%s.txt", cell_type))
expression_file_name = file.path(qtl_dir, sprintf("df__%s.txt", cell_type))
covariates_file_name = file.path(qtl_dir, sprintf("df__%s_cov.txt", cell_type))


snps = SlicedData$new();
snps$fileDelimiter = "\t";      
snps$fileOmitCharacters = "NA"; 
snps$fileSkipRows = 1;          
snps$fileSkipColumns = 1;      
snps$fileSliceSize = 2000;     
snps$LoadFile(SNP_file_name);


gene = SlicedData$new();
gene$fileDelimiter = "\t";     
gene$fileOmitCharacters = "NA";
gene$fileSkipRows = 1;         
gene$fileSkipColumns = 1;      
gene$fileSliceSize = 2000;     
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";     
cvrt$fileOmitCharacters = "NA";
cvrt$fileSkipRows = 1;         
cvrt$fileSkipColumns = 1;  
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


for (ids in list(snps$columnNames, gene$columnNames, cvrt$columnNames)) {
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Matrix donor IDs must be nonempty, nonmissing and unique.")
}
if (!identical(snps$columnNames, gene$columnNames) ||
    !identical(snps$columnNames, cvrt$columnNames))
  stop("Genotype, expression and covariate donor IDs/order differ; align input files first.")

## Run the analysis
# 4. Fit cis-eQTL models -----------------------------------------------
df_snps_pos <- fread(paste0(project_root, "/out/sceQTL430/data/QTL/SNP_location.txt"))
df_snps_pos

df_gene_pos <- fread(paste0(project_root, "/out/sceQTL430/data/QTL/gene_location.txt"))

output_file_name = tempfile()
pvOutputThreshold = 0 
pvOutputThreshold.cis = 1
errorCovariance = numeric()


df_res = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name,
  pvOutputThreshold     = pvOutputThreshold,
  pvOutputThreshold.cis = pvOutputThreshold.cis, 
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  snpspos = as.data.frame(df_snps_pos),
  genepos = as.data.frame(df_gene_pos),
  cisDist = 1e6,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)


print(sprintf("%s: eQTL finished.", date()))

str(df_res)
head(df_res$cis$eqtls)


# 5. Attach BIM fields and export --------------------------------------
# REF and ALT store BIM columns V5 and V6; the labels do not establish
# genomic reference/alternate status or dosage-effect allele coding.
df_SNP <- fread(sprintf(paste0(project_root, "/out/sceQTL430/data/QTL/df_genotype_%s.bim"), cell_type))


df_SNP
df_SNP <- df_SNP[, list(snps = V2, CHR = V1, BP = V4, REF = V5, ALT = V6 )] 


df_res_eqtl <- merge(df_res$cis$eqtls, df_SNP, by = "snps")
fwrite(df_res_eqtl, sprintf("out__%s__cis_eqtls.txt.gz", cell_type), sep="\t")
print(sprintf("%s: All finished.", date()))

}


