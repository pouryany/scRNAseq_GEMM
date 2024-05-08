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


cellchat.late  <- readRDS("data/06_cellchat_scGemms_Tcell_tumor_Late.rds")
cellchat.early <- readRDS("data/06_cellchat_scGemms_Tcell_tumor_early.rds")


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


library(wesanderson)
cel_cols <- wes_palette("Moonrise3",
                        type = "discrete")


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
               stacked = T, do.stat = TRUE,measure = "weight", font.size = 11)


# cairo_pdf("figures/S04B_cellchat_significant_pathways.pdf",
#           width = 5,height = 13, onefile = T)
# gg1
# dev.off()




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
cairo_pdf("figures/cell_cell_communication/11_cellchat_Tcell_Tumor_early_vs_late_networks.pdf",
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
                                        color.use = cel_cols,
                                        arrow.size = 1.2,
                                        remove.isolate = F, 
                                        seed = 38)
  }, 
  error = function(e){ netVisual_aggregate_modified(object.list[[2]],
                                        signaling = j,
                                        alpha.edge = 0,
                                        layout = "circle", 
                                        edge.weight.max = weight.max[1],
                                        edge.width.max = 10,
                                        vertex.label.cex = 2.1, 
                                        signaling.name = paste("Early ", j),
                                        remove.isolate = F,
                                        color.use = cel_cols,
                                        seed = 38)
  }
  )
  
    tryCatch(expr = {netVisual_aggregate_modified(object.list[[2]],
                                 signaling = j,
                                 vertex.receiver = vertex.receiver,
                                 signaling.name = paste( "Late ", j),
                                 vertex.label.cex = 2.1,
                                 color.use = cel_cols,
                                 arrow.size = 1.2,
                                 remove.isolate = F,
                                 seed = 38)
      }, 
             error = function(e){ netVisual_aggregate_modified(
                                        object.list[[1]],
                                        signaling = j,
                                        alpha.edge = 0,
                                        layout = "circle", 
                                        edge.weight.max = weight.max[1],
                                        edge.width.max = 10,
                                        vertex.label.cex = 2.1, 
                                        signaling.name = paste("Late ", j),
                                        remove.isolate = F,
                                        color.use = cel_cols,
                                        seed = 38)
             }
             )
}
dev.off()


i = 1
# combining all the identified signaling pathways from different datasets 





out_address <-"figures/cell_cell_communication/bulk/"
if(!dir.exists(out_address))
  dir.create(out_address)
######>
######>
######>
######>
######> NEEDS SOME WORK HERE
######>
######>
######>
######>
######>
######>
for(j in seq_along(levels(cellchat@meta$labels))) {
  

  
P1 <- netVisual_bubble_modified(cellchat, 
                               targets.use = 1:3,
                               sources.use = j,
                               comparison = c(1, 2),
                               angle.x = 45,
                               thresh = 0.01,
                               font.size = 14,
                               return.data = TRUE,
                               y.size = 16,
                               x.size = 0)

celltype_sender <- levels(cellchat@meta$labels)[j]

# Make an empty matrix of all ligand receptors.
temp <- P1$communication
row_tags <- unique(temp$ligand)
col_tags <- unique(temp$receptor)

interaction_mat <- matrix(nrow = length(row_tags), 
                          ncol = length(col_tags),data = 0)

rownames(interaction_mat) <- row_tags
colnames(interaction_mat) <- col_tags




# Annotate matrices one by one for each sender-receiver pair in early and late 

for(cell_type in levels(cellchat@meta$labels)) {
  
  temp <- P1$communication
  #pvalue < 0.01
  temp <- temp[temp$pval==3,]
  temp <- temp[temp$target==cell_type,]
  row_tags <- unique(temp$ligand)
  col_tags <- unique(temp$receptor)
  
  
  interaction_mat <- matrix(nrow = length(row_tags), 
                            ncol = length(col_tags),data = 0)
  
  rownames(interaction_mat) <- row_tags
  colnames(interaction_mat) <- col_tags
  
 
  
  interaction_dat_early <- temp[temp$dataset == "Early",]
  interaction_mat_early <- interaction_mat
  
  interaction_dat_late <- temp[temp$dataset == "Late",]
  interaction_mat_late <- interaction_mat
  
  for(k in seq_len(nrow(interaction_dat_early))){
    ligand   <- interaction_dat_early$ligand[k]
    receptor <- interaction_dat_early$receptor[k]
    prob     <- interaction_dat_early$prob[k]
    interaction_mat_early[ligand,receptor] <- prob
  }
  
  
  for(k in seq_len(nrow(interaction_dat_late))){
    ligand   <- interaction_dat_late$ligand[k]
    receptor <- interaction_dat_late$receptor[k]
    prob     <- interaction_dat_late$prob[k]
    interaction_mat_late[ligand,receptor] <- prob
  }
  
  
  
  info_path <- temp[,c("ligand","pathway_name")]
  info_path <- unique(info_path)

  
  info_path %<>% dplyr::group_by(.,ligand) %>%
    dplyr::summarize(pathway_name = paste0(unique(pathway_name), collapse = '/'))
  
  labs <- as.data.frame(info_path)
  colnames(labs) <- c("ligand", "Pathway")
  rownames(labs) <- labs$ligand
  
  labs <- labs[rownames(interaction_mat_early),]
  
  labs$Pathway <- factor(labs$Pathway)
 
  
  interaction_mat_early[interaction_mat_early == 0] <- NA
  
early_heat <-  Heatmap(interaction_mat_early, 
          row_split = labs$Pathway, 
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          width  = unit(0.36 * ncol(interaction_mat_late),"cm"),
          height = unit(0.4 * nrow(interaction_mat_late),"cm"),
          row_title = paste0(celltype_sender, " ligands"),
          row_title_gp = gpar(fontsize = 20),
          column_title = paste0("Early\n",cell_type, " receptors"),
          column_title_gp = gpar(fontsize = 20),
          column_names_side = "top",
          na_col = "grey95",
          rect_gp = gpar(col = "white", lwd = 2),
          col = circlize::colorRamp2(c(0,1),
                                     c("#fc9272","#67000d") ),
          border_gp = gpar(col = "black", lty = 1),
          left_annotation = 
            rowAnnotation(foo = 
                            anno_block(gp = gpar(fill = "grey20"),labels_rot =0,
                                       labels = levels(labs$Pathway), 
                                       labels_gp = gpar(col = "white", 
                                                        fontsize = 12),
                                       width = unit(3,"cm"))),
          row_names_side = "left",
          heatmap_legend_param = list(
            title = "Communication \nStrength"
          )
)


interaction_mat_late[interaction_mat_late == 0] <- NA

late_heat <-  Heatmap(interaction_mat_late,
                       row_split = labs$Pathway,
                       cluster_rows = FALSE,
                       width  = unit(0.36 * ncol(interaction_mat_late),"cm"),
                       height = unit(0.4 * nrow(interaction_mat_late),"cm"),
                       cluster_columns = FALSE,
                       row_title = paste0(celltype_sender, " ligands"),
                       row_title_gp = gpar(fontsize = 20),
                       column_title = paste0("Late\n",cell_type, " receptors"),
                       column_title_gp = gpar(fontsize = 20),
                       column_names_side = "top",
                       na_col = "grey95",
                       rect_gp = gpar(col = "white", lwd = 2),
                       col = circlize::colorRamp2(c(0,1),
                                                 c("#fc9272","#67000d") ),
                       # col = circlize::colorRamp2(c(0,1),
                       #                            c("#f7f7f7","#ca0020") ),
                       border_gp = gpar(col = "black", lty = 1),
                       left_annotation = 
                         rowAnnotation(foo = 
                                         anno_block(gp = gpar(fill = "grey20"),labels_rot =0,
                                                    labels = levels(labs$Pathway), 
                                                    labels_gp = gpar(col = "white", 
                                                                     fontsize = 12),
                                                    width = unit(3,"cm"))),
                       row_names_side = "left", 
                      heatmap_legend_param = list(
                        title = "Communication \nStrength"
                      )
)




cairo_pdf(filename = paste0(out_address,celltype_sender,"_to_",cell_type,".pdf"),
          width = 15,height = 15,onefile = T)
draw(early_heat)
draw(late_heat)
dev.off()



cairo_pdf(filename = paste0(out_address,"SIDE_BY_SIDE_",celltype_sender,"_to_",cell_type,".pdf"),
          width = 32,height = 15,onefile = T)
draw(early_heat + late_heat, ht_gap = unit(c(3, 10), "mm"), 
     auto_adjust = FALSE)
dev.off()


}

# 
# temp <- P1$communication
# 
# 
# library(circlize)
# temp <- netVisual_chord_gene_modified(object.list[[1]], sources.use = 3,
#                      targets.use = c(1), 
#                      slot.name = 'net', 
#                      lab.cex = 0.8, small.gap = 3.5, 
#                      title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
# 
# 
# temp$df
# 
# library(reshape2)
# temp <- acast(data = temp$df,formula =  ligand ~ receptor, value.var = "prob")
# 
# 
# cell_Tag <- levels(cellchat@meta$labels)[j]
# ggsave(paste0(out_address,"/",cell_Tag,"_to_tumors.pdf"),P1$gg.obj,
#           width = 7,height = min( box_height* 0.4, 48))
}

# 
# 
# out_address <-"figures/cell_cell_communication/from_tumors"
# if(!dir.exists(out_address))
#   dir.create(out_address)
# 
# 
# for(j in seq_along(levels(cellchat@meta$labels))){
#   P1 <- netVisual_bubble_modified(cellchat, targets.use =  j, 
#                                   sources.use = c(1:3),
#                                   comparison = c(1, 2),
#                                   angle.x = 45,thresh = 0.01,font.size = 14,
#                                   return.data = TRUE,y.size = 12,x.size = 0)
#   
#   
#   box_height  <-
#     length(unique(as.character((P1$communication$interaction_name_2))))
#   
#   cell_Tag <- levels(cellchat@meta$labels)[j]
#   ggsave(paste0(out_address,"/",cell_Tag,"_from_tumors.pdf"),P1$gg.obj,
#          width = 6,height = min( box_height* 0.5, 48))
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# out_address <-"figures/cell_cell_communication/to_myeloids"
# if(!dir.exists(out_address))
#   dir.create(out_address)
# 
# for(j in seq_along(levels(cellchat@meta$labels))){
#   P1 <- netVisual_bubble_modified(cellchat, 
#                                   targets.use =  c(19:20),
#                                   sources.use = j,
#                                   comparison = c(1, 2),
#                                   angle.x = 45,
#                                   thresh = 0.01,
#                                   font.size = 14,
#                                   return.data = TRUE,
#                                   y.size = 16,
#                                   x.size = 0)
#   
#   box_height  <- 
#     length(unique(as.character((P1$communication$interaction_name_2))))
#   
#   cell_Tag <- levels(cellchat@meta$labels)[j]
#   ggsave(paste0(out_address,"/",cell_Tag,"_to_myeloids.pdf"),P1$gg.obj,
#          width = 7,height = min( box_height* 0.4, 48))
# }
# 
# 
# 
# out_address <-"figures/cell_cell_communication/from_myeloids"
# if(!dir.exists(out_address))
#   dir.create(out_address)
# 
# 
# for(j in seq_along(levels(cellchat@meta$labels))){
#   P1 <- netVisual_bubble_modified(cellchat, targets.use =  j, 
#                                   sources.use = c(19:20),
#                                   comparison = c(1, 2),
#                                   angle.x = 45,thresh = 0.01,font.size = 14,
#                                   return.data = TRUE,y.size = 12,x.size = 0)
#   
#   
#   box_height  <-
#     length(unique(as.character((P1$communication$interaction_name_2))))
#   
#   cell_Tag <- levels(cellchat@meta$labels)[j]
#   ggsave(paste0(out_address,"/",cell_Tag,"_from_myeloids.pdf"),P1$gg.obj,
#          width = 6,height = min( box_height* 0.5, 48))
# }
# 
