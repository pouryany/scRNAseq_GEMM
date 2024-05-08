rm(list = ls())

#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)





cellchat.late  <- readRDS("data/05_cellchat_scGemms_late.rds")
cellchat.early <- readRDS("data/05_cellchat_scGemms_early.rds")


object.list    <- list(Early = cellchat.early, Late = cellchat.late)


my_cols <- 
    c("#a6cee3", "#1f78b4",
      "#b2df8a", "#33a02c",
      "#fb9a99", "#e31a1c", "#fdbf6f",
      "#ff7f00", "#cab2d6",
      "#6a3d9a", "#ffff99",
      "#b15928", "#8dd3c7",
      "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
      "#fdb462", "#b3de69", "#fccde5", "#d9d9d9",
      "#bc80bd","#ccebc5","#ffed6f","#ffed6f" )



cellchat <- mergeCellChat(object.list, add.names = names(object.list))


# In the colorbar, red (or blue) represents increased (or decreased)
# signaling in the second dataset compared to the first one.
colorss <- readRDS("data/03_colors.RDS")
colorss <- colorss[names(colorss) %in% cellchat.early@meta$cell_type_secondary]


col_use <- colorss


source("codes/utils.R")


cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional")









out_address <-"figures/cell_cell_communication/to_tumors_selected"
if(!dir.exists(out_address))
  dir.create(out_address)





source("codes/utils.R")



pathways.show <-  c("MIF","CXCL", "CX3C","IFN-II",
                    "CD137","LIGHT","CD80","BTLA","ITGAL-ITGB2",
                    "TNF","TWEAK","PARs","MHC-II","CALCR")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)


for(j in seq_along(levels(cellchat@meta$labels))){

  P1 <- netVisual_bubble_modified(cellchat, targets.use =  c(1:3), 
                                  sources.use = j,point_scaler = 1,
                                  comparison = c(1, 2), pairLR.use = pairLR.CXCL,
                                  angle.x = 45,thresh = 0.01,font.size = 14,
                                  return.data = TRUE, y.size = 16,x.size = 0)
  
  box_height  <- 
    length(unique(as.character((P1$communication$interaction_name_2))))
  
  cell_Tag <- levels(cellchat@meta$labels)[j]
  ggsave(paste0(out_address,"/",cell_Tag,"_to_tumors.pdf"),P1$gg.obj,
         width = 6,height = min( box_height* 0.5, 48))
}



P1 <- netVisual_bubble_modified(cellchat, targets.use =  c(1:3), 
                                sources.use = 18,
                                comparison = c(1, 2), pairLR.use = pairLR.CXCL,
                                angle.x = 45,thresh = 0.01,font.size = 14,
                                return.data = TRUE, y.size = 20,x.size = 20)



P2  <- netVisual_bubble_modified(cellchat, targets.use =  c(1:3), 
                                sources.use = 14,
                                comparison = c(1, 2), pairLR.use = pairLR.CXCL,
                                angle.x = 45,thresh = 0.01,font.size = 14,
                                return.data = TRUE, y.size = 20,x.size = 20)

P1$gg.obj / P2$gg.obj



out_address <-"figures/cell_cell_communication/from_tumors_selected"
if(!dir.exists(out_address))
  dir.create(out_address)





for(j in seq_along(levels(cellchat@meta$labels))){
  P1 <- netVisual_bubble_modified(cellchat, targets.use =  j, 
                                  sources.use = c(1:3),
                                  comparison = c(1, 2), 
                                  angle.x = 45,
                                  pairLR.use = pairLR.CXCL,
                                  thresh = 0.01,font.size = 14,
                                  return.data = TRUE,
                                  color.text.use = FALSE,
                                  y.size = 20,x.size = 20)
  
  
  box_height  <- 
    length(unique(as.character((P1$communication$interaction_name_2))))
  
  cell_Tag <- levels(cellchat@meta$labels)[j]
  ggsave(paste0(out_address,"/",cell_Tag,"_from_tumors.pdf"),P1$gg.obj,
         width = 9,height = min(4 + box_height* 0.5, 48))
}





