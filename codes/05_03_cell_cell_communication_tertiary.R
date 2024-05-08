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









out_address <-"figures/cell_cell_communication/to_tumors_selected2"
if(!dir.exists(out_address))
  dir.create(out_address)





pathways.show <-  c("SPP1")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)
spp1 <- grep( "A9|A8|B6",pairLR.CXCL$interaction_name, value = T,invert = T)
spp1 <- data.frame("interaction_name" = spp1)





pathways.show <-  c("MIF","LIGHT","ITGAL-ITGB2", "TWEAK",
                    "TNF")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)


pairLR.CXCL <- rbind(pairLR.CXCL, spp1)

levels(cellchat@meta$labels)


P1 <- netVisual_bubble_modified(cellchat, targets.use =  c(1:3), 
                                sources.use = 10:17,point_scaler = 1,
                                comparison = c(1, 2), pairLR.use = pairLR.CXCL,
                                angle.x = 45,thresh = 0.01,font.size = 14,
                                return.data = TRUE, y.size = 16,x.size = 0)




P1 <- P1$gg.obj+
      theme(axis.text.x = element_text(angle = 90,
                                       vjust = 0.5, hjust=1))


ggsave(paste0(out_address,"/","CD8_Selected","_to_tumors.pdf"),P1,
       width = 12,height = 6)





#Immune checkpoint





pathways.show <-  c("TIGIT")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)



levels(cellchat@meta$labels)


P1 <- netVisual_bubble_modified(cellchat, targets.use =  c(1:3), 
                                sources.use = c(11:14,17),point_scaler = 1,
                                comparison = c(1, 2), pairLR.use = pairLR.CXCL,
                                angle.x = 45,thresh = 0.01,font.size = 14,
                                return.data = TRUE, y.size = 16,x.size = 0)




P1 <- P1$gg.obj+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1))


ggsave(paste0(out_address,"/","Checkpoint_CD8_Selected","_to_tumors.pdf"),P1,
       width = 8,height = 4)













out_address <-"figures/cell_cell_communication/from_tumors_selected2"
if(!dir.exists(out_address))
  dir.create(out_address)




pathways.show <-  c("MHC-I","CXCL")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)
mhc1 <- grep( "T23_CD8A|Q6|M3|CL9",pairLR.CXCL$interaction_name, value = T)
mhc1 <- data.frame("interaction_name" = mhc1)


pathways.show <-  c("SPP1")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)
spp1 <- grep( "CD44|A4",pairLR.CXCL$interaction_name, value = T)
spp1 <- data.frame("interaction_name" = spp1)





pathways.show <-  c("TNF","CX3C", "FN1","IL2", "MIF")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)


pairLR.CXCL <- rbind(pairLR.CXCL,(mhc1),spp1)

levels(cellchat@meta$labels)


P1 <- netVisual_bubble_modified(cellchat,  sources.use  =  c(1:3), 
                                targets.use = 10:17,point_scaler = 1,
                                comparison = c(1, 2), pairLR.use = pairLR.CXCL,
                                angle.x = 45,thresh = 0.01,font.size = 14,
                                return.data = TRUE, y.size = 16,x.size = 0)




P1 <- P1$gg.obj+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1))


ggsave(paste0(out_address,"/","CD8_Selected","_from_tumors.pdf"),P1,
       width = 12,height = 9)






pathways.show <-  c("PVR","NECTIN")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)


mhc1 <- grep( "TIGIT",pairLR.CXCL$interaction_name, value = T)
mhc1 <- data.frame("interaction_name" = mhc1)

levels(cellchat@meta$labels)



P1 <- netVisual_bubble_modified(cellchat,  sources.use  =  c(1:3), 
                                targets.use = 13,point_scaler = 1,
                                comparison = c(1, 2), pairLR.use = mhc1,
                                angle.x = 45,thresh = 0.01,font.size = 14,
                                return.data = TRUE, y.size = 16,x.size = 0)




P1 <- P1$gg.obj+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1))


ggsave(paste0(out_address,"/","Checkpoint_CD8_Selected","_from_tumors.pdf"),P1,
       width = 4.7,height = 4)




out_address <-"figures/cell_cell_communication/from_myeloid_selected2"
if(!dir.exists(out_address))
  dir.create(out_address)




pathways.show <-  c("PD-L1","ICOS", "CD86")
pairLR.CXCL   <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)

pairLR.CXCL <- grep( "PDCD1|CTLA",pairLR.CXCL$interaction_name, value = T)
pairLR.CXCL <- data.frame("interaction_name" = pairLR.CXCL)


levels(cellchat@meta$labels)


P1 <- netVisual_bubble_modified(cellchat,  sources.use  =  c(19:20), 
                                targets.use = 10:17,point_scaler = 1,
                                comparison = c(1, 2), pairLR.use = pairLR.CXCL,
                                angle.x = 45,thresh = 0.01,font.size = 14,
                                return.data = TRUE, y.size = 16,x.size = 0)




P1 <- P1$gg.obj+
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5, hjust=1))


ggsave(paste0(out_address,"/","checkpoint_myeloid_selected","_to_cd8.pdf"),P1,
       width = 9,height = 4)




