#!/usr/bin/env Rscript

print(date())

library(data.table)
library(optparse)


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


setwd(sub_outpath)


library(SingleCellExperiment)
library(ACTIONet)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# read data
if(1) {
  
  
if(1) {
  sce <- readr::read_rds("./data/sce_snBD_McLean_MtSinai.rds")
  table(sce$Cohort)
  sce <- sce[, sce$Cohort == cohort]
  sce
  }
  
  
  if(cell_type == "major") {
    table(sce$Celltype_inferred)
    sce <- sce[, !sce$Celltype_inferred %in% c("Endo", "Pericytes")]
    sce
  }
  
  v_order_cell_type <- c( "Ex-L2", "Ex-L23", "Ex-L3", "Ex-L4_MYLK","Ex-L45_MET", "Ex-L45_LRRK1", "Ex-L5b_HTR2C", "Ex-L56", "Ex-L56_CC_NTNG2",  "Ex-L6_CC_SEMA3A", "Ex-L6b_SEMA3D", "Ex-L6b_SEMA3E", "In-Rosehip_CHST9", "In-Rosehip_TRPC3","In-Reelin","In-VIP", "In-PV_Chandelier", "In-PV_Basket", "In-SST",  "Oli", "OPC","Ast",  "Mic")
  v_order_cell_type
  
  sce$Celltype_inferred <- factor(sce$Celltype_inferred, levels = v_order_cell_type)
  }



# run analysis 
if(1) {
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
  
  # save(list = c("cellchat_BD", "cellchat_con"), file = "cell_chat.RData")
}


# run analysis
if(1) {
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


if(0) {

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


if(0) {
  cellchat <- cellchat_con
}

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
}




# compare CCC


if(0) {
  df <- readRDS("cell_chat_2.RData")
  str(df)
  
  cellchat_BD <- df[[1]]
  cellchat_con <- df[[2]]
  
  # heatmap ####
   if(0) {
    
     object.list <- list(CON = cellchat_con, BD = cellchat_BD)
     object.list
     
     cellchat <- mergeCellChat(object.list, add.names = names(object.list))
     cellchat
     
     
     cellchat@options$datatype <- "RNA"
     gg1 <- netVisual_heatmap(cellchat)
     gg2 <- netVisual_heatmap(cellchat, measure = "weight")
     gg1 + gg2
     
     p <- gg1 + gg2
     p
     
 
     pdf(sprintf("ccc_heatmap_%s.pdf", cohort), 
         width = 13, height = 10)
     p
     graphics.off()
     }  

  
  
  object.list <- list(BD = cellchat_BD, CON = cellchat_con)
  object.list
  
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  cellchat
  
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  
  p <- gg1 + gg2
  geneticToolBox::ggsave_all(sprintf("ccc_bar1_%s", cohort), 
                             width = 6, height = 4)
}


# pathway between case and controls 
if(0) {
  
  cellchat@netP$pathways
  
  prob_1 <- NULL
  for(i in 1:dim(cellchat_BD@netP$prob)[3]) {
    prob_1  <- c(prob_1, sum(cellchat_BD@netP$prob[, , i]))
  }
  prob_1
  names(prob_1) <- cellchat_BD@netP$pathways
  
  prob_2 <- NULL
  for(i in 1:dim(cellchat_con@netP$prob)[3]) {
    prob_2  <- c(prob_2, sum(cellchat_con@netP$prob[, , i]))
  }
  
  prob_2
  names(prob_2) <- cellchat_con@netP$pathways
  
  setdiff(names(prob_1), names(prob_2))
  
  gplots::venn(list(names(prob_1), names(prob_2)))
  v_share <- names(prob_1)[names(prob_1) %in% names(prob_2)]
  
  prob_1 <- prob_1[names(prob_1) %in% v_share][v_share]
  prob_2[names(prob_2) %in% v_share]
  
  prob_2 <- prob_2[names(prob_2) %in% v_share][v_share]
  
  summary(prob_1/prob_2)
  
  wilcox.test(prob_1, prob_2,  paired = TRUE)
  
  v_ratio <- prob_1/prob_2
  
  v_ratio
  
  v_ratio <- v_ratio[order(v_ratio)]
  
  v_ratio
  
  barplot(v_ratio)
  df <- data.table(pathway = names(v_ratio), ratio = v_ratio)
  
  library(ggpubr)
  
  ggbarplot(df, "pathway", "ratio",
            palette = "Paired",
            position = position_dodge(0.9)) 
  
  cellchat_BD
  
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat_BD, mode = "comparison", stacked = F, do.stat = TRUE)
  gg1 + gg2
  
}
