#!/usr/bin/env Rscript

###############################################################################
# snBD — Cell-cell communication within each cohort
#
# Purpose:
#   Fit CellChat networks separately for BD and CON, then compare and visualize
#   their interaction counts, weights and signaling pathways.
#
# Inputs:
#   --cohort McLean|Sinai; --cell_type major (job label; no cell-type subset).
#   --sub_outpath: job output directory.
#   data/synapse_BD/sce_snBD_McLean_MtSinai.rds, with logcounts, Cohort,
#   Phenotype and Celltype_inferred.
#
# Outputs:
#   cell_chat.RData and cell_chat_2.RData contain saveRDS() lists. 

# Method overview:
#   Create condition-specific objects; attach CellChatDB.human; identify
#   overexpressed genes/interactions; compute/filter communication (min.cells=10).
#
# Dependencies:
#   data.table, optparse, SingleCellExperiment, CellChat, patchwork, ggplot2,
#   readr, dplyr, geneticToolBox; CellChat database/package dependencies.
#
# Usage:
#   From repository root, use cell_cell_communication/cellchat_run_submit.sh
#   or run this script with Rscript and the three arguments above.
#
#
###############################################################################

# install packages
if(0) {
  devtools::install_github("jinworks/CellChat")
  devtools::install_github("immunogenomics/presto")
  BiocManager::install("BiocNeighbors")
  BiocManager::install("ComplexHeatmap")
}



print(date())

# Record the repository root before changing the working directory.
project_root <- normalizePath(Sys.getenv("SNBD_ROOT", unset = "."))

library(data.table)
library(optparse)


# 1. Command-line arguments --------------------------------------------
option_list <- list(
  make_option(c("--cell_type"), type="character", default=NULL, help="cell_type"),
  make_option(c("--cohort"), type="character", default=NULL, help="cohort"),
  make_option(c("--sub_outpath"), type="character", default=NULL, help="sub_outpath")
)


opt <- parse_args(OptionParser(option_list=option_list))
print(opt)


cell_type <- opt$cell_type
cohort <- opt$cohort
sub_outpath <- opt$sub_outpath

dir.create(sub_outpath, recursive = TRUE, showWarnings = FALSE)
setwd(sub_outpath)
dir.create("publication", showWarnings = FALSE)

################


# library('NMF')
# devtools::install_github("jokergoo/circlize")
# library(ComplexHeatmap)


library(SingleCellExperiment)
library(CellChat)
library(ggplot2) 
library(patchwork)
options(stringsAsFactors = FALSE)
print(sessionInfo())



# 2. Read data and subset cohort ---------------------------------------
dir_publication_data <- file.path(project_root, "data")
sce <- readr::read_rds(file = sprintf(paste0(dir_publication_data, "/synapse_BD/%s.rds"), "sce_snBD_McLean_MtSinai"))
sce
table(sce$Cohort)

# select the cohort 
if(cohort == "Sinai") {
    sce <- sce[,sce$Cohort == "MtSinai"]
    sce
  } else if(cohort == "McLean") {
    sce <- sce[,sce$Cohort == "McLean"]
    sce
  } 


  

# Optional random subset for a test run
if(0) {
  sce <- sce[, sample(1:ncol(sce), 1e5)]
  sce
  gc()
}


# run analysis 
if(1) {
  
# 3. Build BD and control CellChat objects -----------------------------
  table(sce$Phenotype, sce$Celltype_inferred)
  
  cellchat_BD <- createCellChat(object = assay(sce[, sce$Phenotype == "BD"], "logcounts"), 
                                meta = as.data.frame(sce[, sce$Phenotype == "BD"]@colData), 
                                group.by = "Celltype_inferred")
  cellchat_BD
  
  cellchat_con <- createCellChat(object = assay(sce[, sce$Phenotype == "CON"], "logcounts"), 
                                 meta = as.data.frame(sce[, sce$Phenotype == "CON"]@colData), 
                                 group.by = "Celltype_inferred")
  
  cellchat_con
  
  
  CellChatDB <- CellChatDB.human 
  
  showDatabaseCategory(CellChatDB)
  
  CellChatDB.use <- CellChatDB 
  
  dplyr::glimpse(CellChatDB$interaction) 
  
  head(CellChatDB.human)
  
  
  cellchat_BD@DB <- CellChatDB.use
  cellchat_BD <- subsetData(cellchat_BD) 
  cellchat_BD <- identifyOverExpressedGenes(cellchat_BD)
  cellchat_BD <- identifyOverExpressedInteractions(cellchat_BD)
  cellchat_BD
  
  cellchat_con@DB <- CellChatDB.use
  cellchat_con <- subsetData(cellchat_con) 
  cellchat_con <- identifyOverExpressedGenes(cellchat_con)
  cellchat_con <- identifyOverExpressedInteractions(cellchat_con)
  
  saveRDS(list(cellchat_BD, cellchat_con), "cell_chat.RData")
}



# Optional loading of saved CellChat objects

if(0) {
  
  df <- readRDS(file = "cell_chat.RData")
  
  str(df)
  
  cellchat_BD <- df[[1]]
  cellchat_con <- df[[2]]
  
  rm(sce)
  gc()
}


# run analysis
if(1) {
  
# 4. Infer and aggregate communication ---------------------------------
  cellchat_BD <- computeCommunProb(cellchat_BD)
  cellchat_BD <- filterCommunication(cellchat_BD, min.cells = 10)
  cellchat_BD <- computeCommunProbPathway(cellchat_BD)
  cellchat_BD <- aggregateNet(cellchat_BD)
  groupSize <- as.numeric(table(cellchat_BD@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat_BD@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat_BD@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  cellchat_con <- computeCommunProb(cellchat_con)
  cellchat_con <- filterCommunication(cellchat_con, min.cells = 10)
  cellchat_con <- computeCommunProbPathway(cellchat_con)
  cellchat_con <- aggregateNet(cellchat_con)
  groupSize <- as.numeric(table(cellchat_con@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  
  netVisual_circle(cellchat_con@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat_con@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  saveRDS(list(cellchat_BD, cellchat_con), "cell_chat_2.RData")
  
}



####### read in processed data + results  ####

if(1) {
  
  
  df <- readRDS("cell_chat_2.RData")
  str(df)
  
  cellchat_BD <- df[[1]]
  cellchat_con <- df[[2]]
  str(cellchat_BD)

}





# Visualization of cell-cell communication network ####

if(1) {
  
# 5. Explore communication networks ------------------------------------
object.list <- list(BD = cellchat_BD, CON = cellchat_con)
object.list

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 

cellchat

cellchat@options$datatype <- "RNA"

# Objects are ordered BD, CON; the default difference is CON minus BD.
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")



cellchat <- cellchat_BD
cellchat

groupSize <- as.numeric(table(cellchat@idents))
groupSize


mat <- cellchat@net$weight

par(mfrow = c(3,4), xpd=TRUE)

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


df_pathway <- cellchat@DB$interaction

library(data.table)
dd <- as.data.table(table(df_pathway$pathway_name))


df_pathway_num <- as.data.frame(table(df_pathway$pathway_name))
df_pathway_num


pathways.show <- c("ANGPT") 

pathways.show <- c("CCK") 

pathways.show

# Optional receiver-specific plot; define vertex.receiver to select receiver groups.
if (exists("vertex.receiver")) {
  netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
}


par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
cellchat

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1

gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
gg2


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
head(pairLR)


netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = pathways.show, remove.isolate = FALSE)

}




# Compare BD and control communication networks #####


if(1) {
  
  df <- readRDS("cell_chat_2.RData")
  str(df)
  
  cellchat_BD <- df[[1]]
  cellchat_con <- df[[2]]
  
  
  # for heatmap ####
   if(0) {
    
     object.list <- list(CON = cellchat_con, BD = cellchat_BD)
     object.list
     
     cellchat <- mergeCellChat(object.list, add.names = names(object.list))
     cellchat
     
     cellchat@options$datatype <- "RNA"
     gg1 <- netVisual_heatmap(cellchat)
     gg2 <- netVisual_heatmap(cellchat, measure = "weight")
     gg1 
     gg2
     
     p <- gg1 + gg2
     p
     
 
     
     # Export the comparison heatmaps 
     pdf(sprintf("./publication/ccc2_heatmap_%s.pdf", cohort), 
         width = 13, height = 10)
     gg1 
     gg2
     graphics.off()
     
     }  

  
  
# 6. BD-control comparison figures -----------------------------------------
  ### for most figures 
  
  object.list <- list(BD = cellchat_BD, CON = cellchat_con)
  object.list
  
  
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  cellchat
  
  
  # bar plot
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  
  p <- gg1 + gg2
  p
  
  
  geneticToolBox::ggsave_all(sprintf("./publication/ccc_bar1_%s", cohort), p, 
                             width = 6, height = 4)
  
  cellchat
  
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)+ theme(legend.position="none")
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)+ theme(legend.position="none")
  gg1 + gg2
  
  
  geneticToolBox::ggsave_all(sprintf("./publication/ccc_bar2_%s_weight", cohort), gg1 + gg2, 
                             width = 8, height = 8)
  
  gg1 <- rankNet(cellchat, mode = "comparison", measure = "count", stacked = T, do.stat = TRUE, return.data = FALSE) + theme(legend.position="none")
  gg2 <- rankNet(cellchat, mode = "comparison", measure = "count", stacked = F, do.stat = TRUE) + theme(legend.position="none")
  gg1 + gg2
  
  # Optional comparison using rankNet2(), when the custom function is available.
  if (exists("rankNet2", mode = "function")) {
    gg2 <- rankNet2(cellchat, mode = "comparison", measure = "count", stacked = F, do.stat = TRUE) + theme(legend.position="none")
    gg1 + gg2
  }

    geneticToolBox::ggsave_all(sprintf("./publication/ccc_bar2_%s", cohort), gg1 + gg2, 
                           width = 8, height = 8)
  
}



print("Finished!")
print(date())

