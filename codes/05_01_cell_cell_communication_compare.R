rm(list = ls())

#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
source("codes/utils.R")


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



colorss <- readRDS("data/03_colors.RDS")
colorss <- colorss[names(colorss) %in% cellchat.early@meta$cell_type_secondary]

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 
#> 'netP','meta', 'idents', 'var.features' ,
#>  'DB', and 'LR'.





future::plan("multisession", workers = 8) # do parallel



cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space

rankSimilarity(cellchat, type = "functional")


gg1 <- rankNet(cellchat, mode = "comparison",
               stacked = T, do.stat = TRUE,measure = "weight",font.size = 11)


cairo_pdf("figures/Ex05B_cellchat_significant_pathways.pdf",
          width = 5,height = 13, onefile = T)
gg1
dev.off()

write.csv(gg1$data,"output/SourceData Ex5_b.csv")


# Compute the network centrality scores
cellchat.late  <- netAnalysis_computeCentrality(cellchat.late,
                                                slot.name = "netP") 
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat.early <- netAnalysis_computeCentrality(cellchat.early, 
                                                slot.name = "netP") 
# the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap,
# allowing ready identification of major signaling roles of cell groups

pathways.show <- as.character(unique(gg1$data$name))


if(!dir.exists("figures/cell_cell_communication"))
   dir.create("figures/cell_cell_communication",recursive = T)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

# related to Figure 4
cairo_pdf("figures/cell_cell_communication/04B_cellchat_early_vs_late_networks.pdf",
          width = 18,height = 14, onefile = T)
#par(mfrow = c(1,2), xpd=TRUE)
for(j in pathways.show){

  # Plot a signaling pathway in early and late samples
  # If the plot is nonexistent in early or late samples, just plot an empty
  # network.
  #> Plot codes are modified to visualize the label names properly.
  tryCatch(expr = {netVisual_aggregate_modified(object.list[[1]], signaling = j,
                                        vertex.receiver = vertex.receiver,
                                        signaling.name = paste("Early ", j),
                                        vertex.label.cex = 2.1,
                                        color.use = colorss,
                                        arrow.size = 1.2,
                                        remove.isolate = F, 
                                        seed = 38)
  }, 
  error = function(e){ netVisual_aggregate_modified(object.list[[2]],
                                        signaling = pathways.show,
                                        alpha.edge = 0,
                                        layout = "circle", 
                                        edge.weight.max = weight.max[1],
                                        edge.width.max = 10,
                                        vertex.label.cex = 2.1, 
                                        signaling.name = paste("Early ", j),
                                        remove.isolate = F,
                                        color.use = colorss,
                                        seed = 38)
  }
  )
  
    tryCatch(expr = {netVisual_aggregate_modified(object.list[[2]],
                                 signaling = j,
                                 vertex.receiver = vertex.receiver,
                                 signaling.name = paste( "Late ", j),
                                 vertex.label.cex = 2.1,
                                 color.use = colorss,
                                 arrow.size = 1.2,
                                 remove.isolate = F,
                                 seed = 38)
      }, 
             error = function(e){ netVisual_aggregate_modified(
                                        object.list[[1]],
                                        signaling = pathways.show,
                                        alpha.edge = 0,
                                        layout = "circle", 
                                        edge.weight.max = weight.max[1],
                                        edge.width.max = 10,
                                        vertex.label.cex = 2.1, 
                                        signaling.name = paste("Late ", j),
                                        remove.isolate = F,
                                        color.use = colorss,
                                        seed = 38)
             }
             )
}
dev.off()


i = 1
# combining all the identified signaling pathways from different datasets 



gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
top_paths <- gg1$data

top_paths <- c(tail(top_paths,15)$name,head(top_paths,15)$name)


ht1 = netAnalysis_signalingRole_heatmap_modified(object.list[[i]],
                                                 pattern = "all",
                                                 name = "Early",
                                                 signaling = top_paths,
                                                 title = names(object.list)[i],
                                                 width = 8.5, height = 12,
                                                 color.heatmap = "GnBu",
                                                 cluster.rows = T,
                                                 show_column_dend = FALSE, 
                                                 show_row_dend = FALSE,
                                                 cluster.cols = F, 
                                                 color.use = rep("grey85",21),
                                                 show_annot = FALSE,font.size = 12)

column_order(ht1)

ht2 = netAnalysis_signalingRole_heatmap_modified(object.list[[i+1]], 
                                                 pattern = "all", 
                                                 name = "Late",
                                                 signaling = top_paths,
                                                 title = names(object.list)[i+1],
                                                 width = 8.5, height = 12,
                                                 color.heatmap = "GnBu",
                                                 cluster.rows = T,
                                                 show_column_dend = FALSE,
                                                 show_row_dend = FALSE,
                                                 col_order = column_order(ht1),
                                                 color.use = rep("grey85",21),font.size = 12)




cairo_pdf("figures/04_cellchat_heatmap_selected.pdf",
          width = 20,height = 26, onefile = T)
ht_list =  ht2 + ht1

draw(ht_list, ht_gap = unit(0.5, "cm"), main_heatmap = "Late")
dev.off()









pathway.union <- union(object.list[[i]]@netP$pathways, 
                       object.list[[i+1]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap_modified(object.list[[i]],
                                           pattern = "all",
                                           name = "Early",
                                           signaling = pathway.union,
                                           title = names(object.list)[i],
                                           width = 8, height = 36,
                                           color.heatmap = "GnBu",
                                           cluster.rows = T,
                                           show_column_dend = FALSE, 
                                           show_row_dend = FALSE,
                                           cluster.cols = T, 
                                           color.use = rep("grey85",21),
                                           show_annot = FALSE,
                                           font.size = 12)

column_order(ht1)

ht2 = netAnalysis_signalingRole_heatmap_modified(object.list[[i+1]], 
                                           pattern = "all", 
                                           name = "Late",
                                           signaling = pathway.union,
                                           title = names(object.list)[i+1],
                                           width = 8, height = 36,
                                           color.heatmap = "GnBu",
                                           cluster.rows = T,
                                           show_column_dend = FALSE,
                                           show_row_dend = FALSE,
                                           col_order = column_order(ht1),
                                           color.use = rep("grey85",21),
                                           font.size = 12)




cairo_pdf("figures/Ex05A_cellchat_heatmap.pdf",
          width = 20,height = 26, onefile = T)
ht_list =  ht1 + ht2

draw(ht_list, ht_gap = unit(0.5, "cm"), main_heatmap = "Late")
dev.off()





pathways.show <- c("CD80","BST2","CD39","IFN-II","IL2","CD137","CD22","CD40",
                   "TNF","MHC-II","MHC-I","PDL2", "CD45","IL16","RANKL","NT",
                   "VCAM","NCAM", "LIGHT","MIF")



out_address <-"figures/cell_cell_communication/to_tumors"
if(!dir.exists(out_address))
  dir.create(out_address)

for(j in seq_along(levels(cellchat@meta$labels))){
P1 <- netVisual_bubble_modified(cellchat, 
                               targets.use =  c(1:3),
                               sources.use = j,
                               comparison = c(1, 2),
                               angle.x = 45,
                               thresh = 0.01,
                               font.size = 14,
                               return.data = TRUE,
                               y.size = 16,
                               x.size = 0)

box_height  <- 
  length(unique(as.character((P1$communication$interaction_name_2))))

cell_Tag <- levels(cellchat@meta$labels)[j]
ggsave(paste0(out_address,"/",cell_Tag,"_to_tumors.pdf"),P1$gg.obj,
          width = 7,height = min( box_height* 0.4, 48))
}



out_address <-"figures/cell_cell_communication/from_tumors"
if(!dir.exists(out_address))
  dir.create(out_address)


for(j in seq_along(levels(cellchat@meta$labels))){
  P1 <- netVisual_bubble_modified(cellchat, targets.use =  j, 
                                  sources.use = c(1:3),
                                  comparison = c(1, 2),
                                  angle.x = 45,thresh = 0.01,font.size = 14,
                                  return.data = TRUE,y.size = 12,x.size = 0)
  
  
  box_height  <-
    length(unique(as.character((P1$communication$interaction_name_2))))
  
  cell_Tag <- levels(cellchat@meta$labels)[j]
  ggsave(paste0(out_address,"/",cell_Tag,"_from_tumors.pdf"),P1$gg.obj,
         width = 6,height = min( box_height* 0.5, 48))
}









out_address <-"figures/cell_cell_communication/to_myeloids"
if(!dir.exists(out_address))
  dir.create(out_address)

for(j in seq_along(levels(cellchat@meta$labels))){
  P1 <- netVisual_bubble_modified(cellchat, 
                                  targets.use =  c(19:20),
                                  sources.use = j,
                                  comparison = c(1, 2),
                                  angle.x = 45,
                                  thresh = 0.01,
                                  font.size = 14,
                                  return.data = TRUE,
                                  y.size = 16,
                                  x.size = 0)
  
  box_height  <- 
    length(unique(as.character((P1$communication$interaction_name_2))))
  
  cell_Tag <- levels(cellchat@meta$labels)[j]
  ggsave(paste0(out_address,"/",cell_Tag,"_to_myeloids.pdf"),P1$gg.obj,
         width = 7,height = min( box_height* 0.4, 48))
}



out_address <-"figures/cell_cell_communication/from_myeloids"
if(!dir.exists(out_address))
  dir.create(out_address)


for(j in seq_along(levels(cellchat@meta$labels))){
  P1 <- netVisual_bubble_modified(cellchat, targets.use =  j, 
                                  sources.use = c(19:20),
                                  comparison = c(1, 2),
                                  angle.x = 45,thresh = 0.01,font.size = 14,
                                  return.data = TRUE,y.size = 12,x.size = 0)
  
  
  box_height  <-
    length(unique(as.character((P1$communication$interaction_name_2))))
  
  cell_Tag <- levels(cellchat@meta$labels)[j]
  ggsave(paste0(out_address,"/",cell_Tag,"_from_myeloids.pdf"),P1$gg.obj,
         width = 6,height = min( box_height* 0.5, 48))
}

