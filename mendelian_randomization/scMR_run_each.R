#!/usr/bin/env Rscript

###############################################################################
# snBD — Two-sample MR for cell-type-specific eQTL instruments
#
# Purpose:
#   Estimate associations of genetically predicted expression with BD or SZ
#   using TwoSampleMR, with heterogeneity and pleiotropy analyses.
#
# Inputs:
#   --trait BD|SZ; --cell_type; --sub_outpath, normally
#   out/sceQTL430/scMR/result. Prepared clumped instruments in ../IV/clump/merge/.
#   GWAS files in out/scMR/gwas/; columns SNP, A1, A2, BETA, SE, P.
#
# Outputs:
#   df_scMR__<trait>__<cell_type>__mr_results.csv, __mr_hetero.csv,
#   __mr_pleio.csv and __mr_single.csv.
#
# Method overview:
#   Construct exposure/outcome tables; harmonize; fit Wald ratio, IVW, weighted
#   median/mode and MR-Egger; assess heterogeneity/intercept and single-SNP effects.
#
# Dependencies:
#   data.table, optparse, TwoSampleMR.
#
# Usage:
#   Called by mendelian_randomization/scMR_run_each.sh; direct Rscript calls
#   must start from the repository root or inherit SNBD_ROOT.
#
###############################################################################

# Record the repository root before changing the working directory.
project_root <- normalizePath(Sys.getenv("SNBD_ROOT", unset = "."))

library(data.table)
library(optparse)


# 1. Arguments and output directory ------------------------------------
option_list <- list(
  make_option(c("--trait"), type="character", default=NULL, help="trait"),
  make_option(c("--cell_type"), type="character", default=NULL, help="cell_type"),
  make_option(c("--sub_outpath"), type="character", default=NULL, help="sub_outpath")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)


trait <- opt$trait
cell_type <- opt$cell_type
sub_outpath <- opt$sub_outpath
stopifnot(length(trait) == 1L, trait %in% c("BD", "SZ"))

# Example arguments for an interactive test 
if(0) {
    
    trait <- "BD"
    cell_type <- "Endo"
    sub_outpath  <- paste0(project_root, "/out/sceQTL430/scMR/result")
    
} 

dir.create(sub_outpath, recursive = TRUE, showWarnings = FALSE)
working_dir <- sub_outpath 
setwd(working_dir)

# 2. Load GWAS and prepared instruments --------------------------------
if(trait == "BD") {
  
  df_gwas <- fread(paste0(project_root, "/out/scMR/gwas/pgc_bip2021_all.txt"))
  df_gwas
  df_y <- df_gwas[, list(SNP, A1, A2, BETA, SE, P)]
  
} else if(trait == "SZ") {
  df_gwas <- fread(paste0(project_root, "/out/scMR/gwas/PGC3_SCZ_wave3.primary.autosome.public.v3.txt"))
  df_gwas
  df_y <- df_gwas[, list(SNP, A1, A2, BETA, SE, P)]
}



df_SNP <- fread(file.path(project_root, "out/sceQTL430/scMR/IV/clump/merge",
                          sprintf("df_merge_clump_%s.txt", cell_type)))
df_SNP



# 3. Format and harmonize exposure/outcome data ------------------------
library(TwoSampleMR)

Ins <- df_SNP[, list(SNP = snps, CHR, BP, A1 = REF, A2 = ALT,
                     BETA = beta, SE = beta / statistic, P = pvalue, Phenotype = paste0(cell_type, "__", gene))]
Ins


# Check for empty instrument tables, duplicate SNPs and invalid effect estimates.
if (!nrow(Ins)) stop("No instruments in the prepared clump table.")
if (anyDuplicated(Ins[, list(Phenotype, SNP)])) stop("Duplicate SNPs within an exposure.")
if (any(!is.finite(Ins$BETA) | !is.finite(Ins$SE) | Ins$SE <= 0))
  stop("Non-finite effect or non-positive exposure SE; check beta/statistic.")

Ins <- format_data(as.data.frame(Ins), type = "exposure", header = TRUE,
                   phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "BETA",
                   se_col = "SE",  effect_allele_col = "A1",
                   other_allele_col = "A2", pval_col = "P")

head(Ins)

df_y <- df_y[SNP %chin% Ins$SNP]
if (!nrow(df_y)) stop("No instruments are present in the outcome GWAS.")
if (anyDuplicated(df_y$SNP)) stop("Duplicate outcome SNP IDs; resolve before MR.")
outcome_dat <- format_data(as.data.frame(df_y), type = "outcome", header = TRUE,
                           snp_col = "SNP", beta_col = "BETA",
                           se_col = "SE",  effect_allele_col = "A1",other_allele_col = "A2", pval_col = "P")

head(outcome_dat)


dat <- harmonise_data(
  exposure_dat = Ins, 
  outcome_dat = outcome_dat
)


head(dat)
if (!nrow(dat) || !any(dat$mr_keep, na.rm = TRUE))
  stop("No usable instruments remain after harmonization.")

##run the MR and sensitivity analyses 
# 4. MR estimators and sensitivity analyses ----------------------------
mr_results <- NULL
mr_hetero <- NULL
mr_pleio <- NULL
mr_single <- NULL


mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))  # main MR analysis
mr_hetero <- mr_heterogeneity(dat) # heterogeneity test across instruments
mr_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test  
mr_single <- mr_singlesnp(dat) #single SNP MR using Wald ratio

if (is.null(mr_results) || !nrow(mr_results)) stop("MR returned no estimates; inspect instrument counts.")

# Export MR estimates and sensitivity-test results. 
# 5. Export results ----------------------------------------------------
file_base <- sprintf("df_scMR__%s__%s", trait, cell_type)
file_base

fwrite(mr_results, paste0(file_base, "__mr_results.csv"))
fwrite(mr_hetero, paste0(file_base, "__mr_hetero.csv"))
fwrite(mr_pleio, paste0(file_base, "__mr_pleio.csv"))
fwrite(mr_single, paste0(file_base, "__mr_single.csv"))
