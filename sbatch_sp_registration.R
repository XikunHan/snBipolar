#!/usr/bin/env Rscript

print(date())

f_log <- function(mes) {
  cat("\n", sprintf('%s: %s', date(), mes), "\n")
}


library(data.table)
library(optparse)

library("spatialLIBD")
library(SingleCellExperiment, quietly = TRUE)


option_list <- list(
  make_option(c("--cohort"), type="character", default=NULL, help="cohort"),
  make_option(c("--sub_outpath"), type="character", default=NULL, help="sub_outpath")
)


opt <- parse_args(OptionParser(option_list=option_list))
print(opt)


cohort <- opt$cohort
sub_outpath <- opt$sub_outpath

setwd(sub_outpath)


library("spatialLIBD")
library("SingleCellExperiment")

layer_modeling_results <- fetch_data(type = "modeling_results")
layer_modeling_results 
layer_modeling_results$enrichment[1:5, 1:5]


dir_publication_data <- "./data"

sce <- readr::read_rds(file = sprintf(paste0(dir_publication_data, "/%s.rds"), "sce_snBD_McLean_MtSinai"))
sce$cellType <- sce$Celltype_inferred
sce$cellType <- stringr::str_replace(sce$cellType, "-", "_")
sce 

set.seed(123)

table(sce$Cohort)

sce_test <- sce[, sce$Cohort == cohort]

gc()

df_gene <- as.data.table(rowData(sce))
fwrite(df_gene, "df_gene_annotation.csv")

sce_modeling_results <- registration_wrapper(
  sce = sce_test,
  var_registration = "cellType",
  var_sample_id = "SubID",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)



# save results 
readr::write_rds(sce_modeling_results, paste0("sce_modeling_results_", cohort, ".rds"))


if(0) {
  sce_modeling_results <-  readr::read_rds(paste0("sce_modeling_results_", cohort, ".rds"))
}



if(0) {
  
  registration_t_stats <- sce_modeling_results$enrichment[, grep("^t_stat", colnames(sce_modeling_results$enrichment))]
  colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))
  registration_t_stats
  
  df_gene <- fread("df_gene_annotation.csv")
  
  table(df_gene$gene_name == rownames(registration_t_stats))
  
  registration_t_stats$rn <- rownames(registration_t_stats)
  
  registration_t_stats
  
  df_gene
  
  registration_t_stats <- merge(registration_t_stats, df_gene[, list(gene_name, gene_id)], by.x = "rn", by.y = "gene_name") 
  registration_t_stats
  
  rownames(registration_t_stats) <- registration_t_stats$gene_id
  
  registration_t_stats$rn <- NULL
  registration_t_stats$gene_id <- NULL
  
  registration_t_stats
  
  dim(registration_t_stats)
  registration_t_stats[1:5, 1:5]
  
  cor_layer <- layer_stat_cor(
    stats = registration_t_stats,
    modeling_results = layer_modeling_results,
    model_type = "enrichment",
    top_n = 100
  )
  
  v_order_cell_type <- c( "Ex-L2", "Ex-L23", "Ex-L3", "Ex-L4_MYLK","Ex-L45_MET", "Ex-L45_LRRK1", "Ex-L5b_HTR2C", "Ex-L56", "Ex-L56_CC_NTNG2",  "Ex-L6_CC_SEMA3A", "Ex-L6b_SEMA3D", "Ex-L6b_SEMA3E", "In-Rosehip_CHST9", "In-Rosehip_TRPC3","In-Reelin","In-VIP", "In-PV_Chandelier", "In-PV_Basket", "In-SST",  "Oli", "OPC","Ast",  "Mic","Endo", "Pericytes")
  v_order_cell_type
  v_order_cell_type <- stringr::str_replace_all(v_order_cell_type, "-", "_")
  
  rownames(cor_layer) %in% v_order_cell_type
  cor_layer <- cor_layer[v_order_cell_type, ]
  
  
  cor_layer
  
  readr::write_rds(sce_modeling_results, paste0("cor_layer_", cohort, ".rds"))
  cor_layer
  
  cor_layer
  
  dir_publication <- "./publication/"
  
  
  layer_matrix_plot <- function (matrix_values, matrix_labels = NULL, xlabs = NULL, 
                                 layerHeights = NULL, mypal = c("white", (grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, 
                                                                                                                               "YlOrRd")))(50)), breaks = NULL, axis.args = NULL, srt = 45, 
                                 mar = c(12, 4 + (max(nchar(rownames(matrix_values)))%/%3) * 
                                           0.5, 4, 2) + 0.1, cex = 1.2) 
  {
    if (is.null(xlabs)) {
      if (is.null(colnames(matrix_values))) {
        xlabs <- paste0("V", seq_len(ncol(matrix_values)))
      }
      else {
        xlabs <- colnames(matrix_values)
      }
    }
    if (is.null(layerHeights)) {
      layerHeights <- c(0, seq_len(nrow(matrix_values))) * 
        15
    }
    if (is.null(matrix_labels)) {
      matrix_labels <- matrix("", ncol = ncol(matrix_values), 
                              nrow = nrow(matrix_values), dimnames = dimnames(matrix_values))
    }
    stopifnot(length(layerHeights) == nrow(matrix_values) + 
                1)
    stopifnot(length(xlabs) == ncol(matrix_values))
    stopifnot(layerHeights[1] == 0)
    midpoint <- function(x) {
      x[-length(x)] + diff(x)/2
    }
    par(mar = mar)
    fields::image.plot(x = seq(0, ncol(matrix_values), by = 1), 
                       y = layerHeights, z = as.matrix(t(matrix_values)), col = mypal, 
                       xaxt = "n", yaxt = "n", xlab = "", ylab = "", breaks = breaks, 
                       nlevel = length(mypal), axis.args = axis.args)
    axis(2, rownames(matrix_labels), at = midpoint(layerHeights), 
         las = 1)
    axis(1, rep("", ncol(matrix_values)), at = seq(0.5, ncol(matrix_values) - 
                                                     0.5))
    text(x = seq(0.5, ncol(matrix_values) - 0.5), y = -1 * max(nchar(xlabs))/2, 
         xlabs, xpd = TRUE, srt = srt, cex = cex, adj = 1)
    abline(h = layerHeights, v = c(0, seq_len(ncol(matrix_values))))
    text(x = rep(seq(0.5, ncol(matrix_values) - 0.5), each = nrow(matrix_values)), 
         y = rep(midpoint(layerHeights), ncol(matrix_values)), 
         as.character(matrix_labels), cex = cex * 3/4, font = 2)
  }
  
  layer_stat_cor_plot <- function (cor_stats_layer, max = 0.81, min = -max, layerHeights = NULL, 
                                   cex = 1.2) 
  {
    theSeq <- seq(min, max, by = 0.01)
    my.col <- (grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, 
                                                                    "PRGn")))(length(theSeq))
    cor_stats_layer[cor_stats_layer <= min] <- min
    cor_stats_layer[cor_stats_layer >= max] <- max
    mat_vals <- t(cor_stats_layer)
    if (all(rownames(mat_vals) %in% c("WM", paste0("Layer", 
                                                   seq_len(6))))) {
      rownames(mat_vals) <- gsub("ayer", "", rownames(mat_vals))
      mat_vals <- mat_vals[c("WM", paste0("L", rev(seq_len(6)))), 
                           , drop = FALSE]
      if (is.null(layerHeights)) {
        layerHeights <- c(0, 40, 55, 75, 85, 110, 120, 135)
      }
    }
    midpoints <- seq(min, max, length.out = length(my.col))
    delta <- (midpoints[2] - midpoints[1])/2
    breaks <- c(midpoints[1] - delta, midpoints + delta)
    legend_cuts <- seq(-1, 1, by = 0.1)
    legend_cuts <- legend_cuts[legend_cuts >= min & legend_cuts <= 
                                 max]
    axis.args <- list(at = legend_cuts, labels = legend_cuts)
    layer_matrix_plot(matrix_values = mat_vals, matrix_labels = NULL, 
                      xlabs = NULL, layerHeights = layerHeights, mypal = my.col, 
                      breaks = breaks, axis.args = axis.args, srt = 90, cex = cex)
  }
  
  
  
  pdf(sprintf(paste0(dir_publication, "sp_registration_%s.pdf"), cohort), width = 10, height = 8)
  layer_stat_cor_plot(cor_layer, max = max(cor_layer))
  dev.off()
  
  
  
  anno <- annotate_registered_clusters(
    cor_stats_layer = cor_layer,
    confidence_threshold = 0.25,
    cutoff_merge_ratio = 0.25
  )
  
  anno
  
  
  
  
df <- fread("./publication/df_sig_gene_sp_enrichment_reorder.csv")
  
if (!exists("modeling_results")) {
 modeling_results <- fetch_data(type = "modeling_results")
   }
modeling_results

df

df_r <- gene_set_enrichment(df,
                            modeling_results = modeling_results,
                            fdr_cut = 0.05,
                            model_type = "enrichment"
                            )
  


gene_set_enrichment_plot <- function(enrichment,
         xlabs = unique(enrichment$ID),
         PThresh = 12,
         ORcut = 3,
         enrichOnly = FALSE,
         layerHeights = c(0, seq_len(length(unique(enrichment$test)))) * 15,
         mypal = c(
           "white",
           grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(50)
         ),
         cex = 1.2) {
  ## Re-order and shorten names if they match our data
  if (all(unique(enrichment$test) %in% c("WM", paste0("Layer", seq_len(6))))) {
    enrichment$test <-
      factor(gsub("ayer", "", enrichment$test), levels = rev(c(paste0(
        "L", seq_len(6)
      ), "WM")))
  }
  
  stopifnot(is(enrichment, "data.frame"))
  stopifnot(all(c("ID", "test", "OR", "Pval") %in% colnames(enrichment)))
  stopifnot(length(layerHeights) == length(unique(enrichment$test)) + 1)
  stopifnot(ORcut <= PThresh)
  stopifnot(length(xlabs) == length(unique(enrichment$ID)))
  
  ## Convert to -log10 scale and threshold the pvalues
  enrichment$log10_P_thresh <-
    round(-log10(enrichment$Pval), 2)
  enrichment$log10_P_thresh[which(enrichment$log10_P_thresh > PThresh)] <-
    PThresh
  
  if(1) {
    if (enrichOnly) {
      enrichment$log10_P_thresh[enrichment$OR < 1] <- 0
    }
    enrichment$OR_char <- as.character(round(enrichment$OR, 1))
    # enrichment$OR_char[enrichment$log10_P_thresh < ORcut] <- ""
    enrichment$fdr <- p.adjust(enrichment$Pval, method = "BH")
    enrichment$OR_char[enrichment$fdr > 0.01] <- ""
     
  } else {
    # raw
    if (enrichOnly) {
      enrichment$log10_P_thresh[enrichment$OR < 1] <- 0
    }
    enrichment$OR_char <- as.character(round(enrichment$OR, 1))
    enrichment$OR_char[enrichment$log10_P_thresh < ORcut] <- ""
    
  }
  
  
  make_wide <- function(var = "OR_char") {
    res <-
      reshape(
        enrichment,
        idvar = "ID",
        timevar = "test",
        direction = "wide",
        drop = colnames(enrichment)[!colnames(enrichment) %in% c("ID", "test", var)],
        sep = "_mypattern_"
      )[, -1, drop = FALSE]
    colnames(res) <-
      gsub(".*_mypattern_", "", colnames(res))
    rownames(res) <- unique(enrichment$ID)
    res <- res[, levels(as.factor(enrichment$test))]
    t(res)
  }
  wide_or <- make_wide("OR_char")
  wide_p <- make_wide("log10_P_thresh")
  
  layer_matrix_plot(
    matrix_values = wide_p,
    matrix_labels = wide_or,
    # matrix_labels = NULL,
    xlabs = xlabs,
    layerHeights = layerHeights,
    mypal = mypal,
    cex = cex,
    srt = 90,
    mar = c(12, 4 + (max(nchar(rownames(wide_p))) %/% 3) * 0.5, 4, 2) + 0.1
  )
}
  
pdf(sprintf(paste0(dir_publication, "sp_enrichment_%s.pdf"), cohort), width = 10, height = 8)
gene_set_enrichment_plot(df_r)
dev.off()

pdf(sprintf(paste0(dir_publication, "sp_enrichment_%s.pdf"), cohort), width = 10, height = 8)
spatialLIBD::gene_set_enrichment_plot(df_r)
dev.off()

}



