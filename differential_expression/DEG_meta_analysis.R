#!/usr/bin/env Rscript

###############################################################################
# snBD — Meta-analysis of DEG results across cohorts (cell-type–matched)
#
# Purpose:
#   Perform fixed-effect meta-analysis of differential expression (DEG) results
#   from two independent cohorts (e.g., McLean and Sinai) for the same trait
#   and harmonized cell types. Combine effect estimates by inverse-variance
#   weighting (IVW) and P values using Fisher's method for each gene and cell type.
#
# Inputs:
#   CSV files with gene-level statistics for each cohort (one row per gene ×
#   cell type). The example uses 50 genes:
#     - ./data/DEG_MASTRE_McL_gene_test.csv
#     - ./data/DEG_MASTRE_Sinai_gene_test.csv
#
#   Required columns in each file:
#     gene       : gene symbol (character)
#     cell_type  : harmonized cell-type label (character)
#     logFC      : effect size in the input table's declared units
#     ci.lo      : lower bound of 95% CI for logFC
#     ci.hi      : upper bound of 95% CI for logFC
#     P          : cohort P value (MAST exports the hurdle-test P value)
#
# Outputs:
#   - dd: meta-analysis table in memory; optional fwrite() is commented below.
#   - ./data/meta_analysis_IVW_MASTRE.csv contains example results.
#     Columns include (not exhaustive):
#       ID           : "gene_BD_celltype" identifier
#       b_meta       : IVW meta effect in the input logFC units
#       se_meta      : standard error for b_meta
#       p_hurdle     : Fisher combined p from cohort p-values (P.x, P.y)
#       FDR          : BH-adjusted p_hurdle
#       log2FC       : log2FC
#
# Dependencies:
#   R >= 4.2
#   data.table, ggplot2, stringr, metap, dplyr
#   geneticToolBox (github.com/XikunHan/geneticToolBox)
#
# Optional installation helper:
#   if (0) {
#     install.packages("data.table")
#     install.packages("metap")
#     install.packages("ggplot2")
#     install.packages("stringr")
#     if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
#     devtools::install_github("XikunHan/geneticToolBox")
#   }
#
# Usage:
#   From repository root: Rscript differential_expression/DEG_meta_analysis.R
#   Uses the example input files; set the input paths for the full analysis.
#
# Date  : 2023-3
# R version: 4.2.0 (2022-04-22) on x86_64-pc-linux-gnu
###############################################################################



if(0) {
  install.packages("data.table")
  devtools::install_github("XikunHan/geneticToolBox")
}

# 1. Packages ----------------------------------------------------------
library(data.table)
library(geneticToolBox) 
library(ggplot2)


# meta analysis ####

# Example cohort tables; use the full merged cohort tables for the study analysis.
# 2. Read cohort statistics --------------------------------------------
df1 <- fread("./data/DEG_MASTRE_McL_gene_test.csv")
df2 <- fread("./data/DEG_MASTRE_Sinai_gene_test.csv")
# df1 <- fread("./data/merge_MASTlmm_McLean_RUV.csv.gz")
# df2 <- fread("./data/merge_MASTlmm_Sinai_RUV.csv.gz")



# 3. Prepare statistics and gene/cell-type keys ------------------------
df1[, FDR:= p.adjust(P, method = "BH")] 
df2[, FDR:= p.adjust(P, method = "BH")]
  
df1[, cell_type := stringr::str_remove_all(cell_type, "Block|V1")]
df2[, cell_type := stringr::str_remove_all(cell_type, "Block|V1")]

df1 <- df1[complete.cases(df1)]
df2 <- df2[complete.cases(df2)]

# ## for logFC
df1
df1[, se := (ci.hi - ci.lo)/2/qnorm(0.975)]
df2[, se := (ci.hi - ci.lo)/2/qnorm(0.975)]

  
df1[, ID := paste0(gene, "_BD_", cell_type)]
df2[, ID := paste0(gene, "_BD_", cell_type)]
if (anyDuplicated(df1$ID) || anyDuplicated(df2$ID))
  stop("Duplicate gene/cell-type IDs after label cleaning; resolve before meta-analysis.")
if (any(!is.finite(df1$se) | df1$se <= 0) || any(!is.finite(df2$se) | df2$se <= 0))
  stop("Invalid standard errors reconstructed from confidence intervals.")




# 4. Merge cohorts and inspect effect-based statistics -----------------
df_m <- merge(df1, df2, by = "ID")
df_m <- df_m[complete.cases(df_m)]
  

df_m$Pz.x <- 2 * pnorm(-abs(df_m$logFC.x / df_m$se.x))
df_m$Pz.y <- 2 * pnorm(-abs(df_m$logFC.y / df_m$se.y))


cor(-log10(df_m$Pz.x), -log10(df_m$P.x), method = "spearman")
cor(-log10(df_m$Pz.y), -log10(df_m$P.y), method = "spearman", use = "pair")


plot(-log10(df_m$Pz.x), -log10(df_m$P.x))

df_m[-log10(Pz.x) < 20, plot(-log10(Pz.x), -log10(P.x))]


table(-log10(df_m$Pz.x) > 4,-log10(df_m$P.x) > 4)


# using inverse_variance method, equivalent to metafor::rma with method="FE"
# 5. Inverse-variance meta-analysis ------------------------------------
dd <- geneticToolBox::meta_inverse_variance_df(df_m, 
                                                 beta1 = "logFC.x",
                                                 beta2 = "logFC.y",
                                                 se1 = "se.x",
                                                 se2 = "se.y"
  )
  

dd

# BiocManager::install("multtest")

# setnames(dd, "p_meta", "p_hurdle")

# 6. Fisher combination and BH adjustment ------------------------------
dd$p_hurdle <- apply(dd[, list(P.x, P.y)], 1, function(x) metap::sumlog(x)$p)
# BH is applied to the rows in this run; example-only FDR is not study-wide FDR.
dd[, FDR := p.adjust(p_hurdle, method = "BH")]
dd[, log2FC := b_meta/log(2)]

dd
cor(-log10(dd$p_hurdle), -log10(dd$p_meta), use = "complete")
cor(-log10(dd$p_hurdle), -log10(dd$p_meta), method = "spearman")
plot(-log10(dd$p_hurdle), -log10(dd$p_meta))

table(-log10(dd$p_hurdle) > 4, -log10(dd$p_meta) > 4)


library(dplyr)
res <- dd %>%
  group_by(cell_type.x) %>%
  summarise(
    cor_xy = cor(p_hurdle, p_meta, method = "pearson",  use = "complete"),
  )



# fwrite(dd, "./data/meta_analysis_IVW_MASTRE.csv")

# 7. Session information -----------------------------------------------
print(sessionInfo())




