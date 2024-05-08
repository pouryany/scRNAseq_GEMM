rm(list = ls())


library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
source("codes/utils.R")

cell_tags <-  c(
    "Alv-Basal", "Alv-Proliferating",
    "Basal",
    "Alv 1", "Alv 2","Alv 3")    


main_list <- list()
for (ij in cell_tags){
    i <- gsub("/","_",ij)
    saveTag <- 
      paste0("output/tumor_late_vs_early/DE_genes_tumorCluster_",
             i,".csv")
    cluster0.markers <- read.csv(saveTag) 
    cluster0.markers$cell <- ij
    main_list <- rbind(main_list,cluster0.markers)
  
}

cell_tags <-  c(
    "CD4+ Naïve", "CD4+ Treg",             
    "CD8+ Naïve", "CD8+ CM",
    "CD8+ ISG","CD8+ Effector",
    "CD8+ TRM",
    "CD8+ Progenitor Ex",
    "CD8+ Ex", 
    "CD8+ Proliferating",
    "NK")    


for (ij in cell_tags){
    i <- gsub("/","_",ij)
    saveTag <-
      paste0("output/immune_late_vs_early/DE_genes_immuneCluster_",
             i,".csv")
    cluster0.markers <- read.csv(saveTag) 
    cluster0.markers$cell <- ij
    main_list <- rbind(main_list,cluster0.markers)
    
}

main_list <- main_list[,c("Gene","cell")]
names(main_list) <- c("feature","ident")

main_list$my_text <- "***"

scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
scGemms$initial_clusts <- scGemms$seurat_clusters
colorss <- readRDS("data/03_colors.RDS")



t_cells <- grep("T cell: CD8", unique(scGemms$cell_type_detail),value = T)
scGemms0 <- subset(scGemms, subset = cell_type_detail %in% t_cells)


cell_tags <- unique(as.character(scGemms0$seurat_clusters))


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]


Idents(scGemms0) <- "cell_type_secondary"

order_vec <-  c(           
    "CD8+ Naïve", "CD8+ CM",
    "CD8+ ISG","CD8+ Effector",
    "CD8+ TRM",
    "CD8+ Progenitor Ex",
    "CD8+ Ex", 
    "CD8+ Proliferating")    


scGemms0@active.ident <- factor(scGemms0@active.ident, levels = order_vec)

my_features <- c("Klf2", 
                 "Tcf7",
                 "S1pr1",
                 "Klrd1","Fcer1g",
                 "Cd160", "Tigit",
                 "Ccl4","Ccl5","Cxcr6")

vln_eGFP1 <- VlnPlot_costumized(scGemms0,my_features = my_features,
                                idents_order = order_vec,
                                signif_list = main_list,
                                plot_alphabetic = F,
                                x_size = 16)


cairo_pdf("figures/S03D_t_cell_x_markers.pdf",width = 5,height = 8,onefile = T)
vln_eGFP1
dev.off()

######
#*
#*
#* Adding NK Cells and DNAJB1 to the mix 
#* 
#* Remove if not necessary
#*
#*
#####

t_cells <- grep("T cell: CD8|NK", unique(scGemms$cell_type_detail),value = T)
scGemms0 <- subset(scGemms, subset = cell_type_detail %in% t_cells)


cell_tags <- unique(as.character(scGemms0$seurat_clusters))


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]


Idents(scGemms0) <- "cell_type_secondary"

order_vec <-  c(           
  "CD8+ Naïve", "CD8+ CM",
  "CD8+ ISG","CD8+ Effector",
  "CD8+ TRM",
  "CD8+ Progenitor Ex",
  "CD8+ Ex", 
  "CD8+ Proliferating", "NK")    


scGemms0@active.ident <- factor(scGemms0@active.ident, levels = order_vec)

my_features <- c("Klf2", 
                 "Tcf7",
                 "S1pr1",
                 "Klrd1","Fcer1g",
                 "Cd160", "Tigit",
                 "Ccl4","Ccl5","Cxcr6")

vln_eGFP1 <- VlnPlot_costumized(scGemms0,my_features = my_features,
                                idents_order = order_vec,
                                signif_list = main_list,
                                plot_alphabetic = F,
                                x_size = 16)


cairo_pdf("figures/S03D_alternative_t_cellNK_x_markers.pdf",width = 5,height = 8.5,onefile = T)
vln_eGFP1
dev.off()





order_vec <-  c(
    "Alv-Proliferating", "Alv-Basal",
    "Basal",
    "NK")    

scGemms0 <- subset(scGemms, subset = cell_type_secondary %in% order_vec)
Idents(scGemms0) <- "cell_type_secondary"

scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]

scGemms0@active.ident <- factor(x = scGemms0@active.ident, levels = order_vec)

my_features <- c("Spp1","Itgav","Itgb6", 
  "Tnfsf11","Tnfrsf11a",
  "Sell","Podxl",
  "Lta","Tnfrsf1a","Ltb","Ltbr",
  "Itgb2","Icam1",
  "Gzma","F2r","F2rl1","Pard3")


vln_eGFP <- VlnPlot_costumized(scGemms0,my_features = my_features,
                               idents_order = order_vec,
                                signif_list = main_list,
                               plot_alphabetic = F)
 ggsave("figures/S05B_nk_cell_chat.pdf",
        vln_eGFP,width = 5,height = 12)

 
 my_features <- c("Ifng","Ifngr1","Ifngr2")
 
 
 vln_eGFP1 <- VlnPlot_costumized2(scGemms0,my_features = my_features,
                                idents_order = order_vec,
                                signif_list = main_list,
                                plot_alphabetic = F) &
     theme(axis.title.x = element_blank(),
           axis.text.x =element_text(size = 20, face = "plain",
                                     angle = 90, vjust = 0.5),
           axis.text.y = element_text(size = 10, face = "plain"),
           strip.text =element_text(size = 20, face = "plain",
                                    vjust = 1,hjust = 0),
           plot.title = element_text(size = 30, face = "plain"),
           axis.title.y = element_blank(),
           title = element_text(size = 20, face = "plain"),
           legend.position = "none")

 
 my_features <- c("Mif","Cd74","Cd44", "Ackr3")
 
 
 vln_eGFP2 <- VlnPlot_costumized2(scGemms0,my_features = my_features,
                                  idents_order = order_vec,
                                  signif_list = main_list,
                                  plot_alphabetic = F) &
     theme(axis.title.x = element_blank(),
           axis.text.x =element_text(size = 20, face = "plain",
                                     angle = 90, vjust = 0.5),
           axis.text.y = element_text(size = 10, face = "plain"),
           strip.text =element_text(size = 20, face = "plain",
                                    vjust = 1,hjust = 0),
           plot.title = element_text(size = 30, face = "plain"),
           axis.title.y = element_blank(),
           title = element_text(size = 20, face = "plain"),
           legend.position = "none")


  vln_eGFP1 + vln_eGFP2 & plot_layout(heights = c(3,4))
 
 
 
 ggsave("figures/04D_select_nk_cell_chat.pdf",
        width = 3,height = 10)
 

 
 scGemms0 <- subset(scGemms,
                    cells = names(grep(pattern = "Alv-|Basal",
                              scGemms$cell_type_secondary,value = T)))
 
 
 
 order_vec <-  c( "Alv-Proliferating","Alv-Basal",  "Basal")    
 
 
 
 
 cell_tags <- unique(as.character(scGemms0$cell_type_secondary))
 Idents(scGemms0) <- "cell_type_secondary"
 
 
 scGemms0@active.ident <- factor(x = scGemms0@active.ident, levels = order_vec)

 t_cells <- grep("CD|NK", unique(scGemms$cell_type_secondary),value = T)
 t_cells <- subset(scGemms, subset = cell_type_secondary %in% t_cells)
 t_cells$cell_type_secondary <- droplevels(t_cells$cell_type_secondary)
 
 
 
 t_cells[["immune_type"]] <-plyr::mapvalues(t_cells$cell_type_secondary,
                                            from = 
                                                c("CD4+ Naïve", 
                                                  "CD4+ Treg",             
                                                  "CD8+ Naïve",
                                                  "CD8+ CM",
                                                  "CD8+ Effector",
                                                  "CD8+ ISG",
                                                  "CD8+ Progenitor Ex",
                                                  "CD8+ Ex", 
                                                  "CD8+ Proliferating",
                                                  "CD8+ TRM",
                                                  "NK"),
                                            to =
                                                c("CD4+",
                                                  "CD4+",
                                                  "CD8+",
                                                  "CD8+",
                                                  "CD8+",
                                                  "CD8+",
                                                  "CD8+",
                                                  "CD8+",
                                                  "CD8+",
                                                  "CD8+",
                                                  "NK"))
 
 
 
 
 
 
 cell_tags <- c("CD4+","CD8+","NK")
 
 main_list <- list()
 for (ij in cell_tags){
     i <- gsub("/","_",ij)
     saveTag <- paste0("output/broad_immune/DE_genes_immuneCluster_",
                       i,".csv")
     cluster0.markers <- read.csv(saveTag) 
     cluster0.markers$cell <- ij
     main_list <- rbind(main_list,cluster0.markers)
     
 }
 
 main_list <- main_list[,c("Gene","cell")]
 names(main_list) <- c("feature","ident")
 
 main_list$my_text <- "***"
 
 
 
 
 cell_tags <- c("CD4+","CD8+","NK")
 Idents(t_cells) <- "immune_type"
 order_vec <- cell_tags
 
 
 my_features <- c("Cd28",
                  "Cd27","Ccr7","Klf2","Cd69",
                  "Itga1","Itga2","S1pr1",
                  "Eomes","Tcf7")
 
 vln_eGFP1 <- VlnPlot_costumized(t_cells,my_features = my_features,
                                 idents_order = order_vec,
                                 signif_list = main_list,
                                 plot_alphabetic = F) +
   coord_fixed(0.3)
 
 
 my_features <- c("Gzma","Ifng","Prf1",
   "Pdcd1","Lag3","Cd244a","Tigit",
   "Ctla4","Cd160","Ccl3", "Ccl4","Ccl5")
 
 vln_eGFP2 <- VlnPlot_costumized(t_cells,my_features = my_features,
                                 idents_order = order_vec,
                                 signif_list = main_list,
                                 plot_alphabetic = F) + 
   coord_fixed(0.3)
 
 
 cairo_pdf("figures/03C_t_cell_x_markers_broad.pdf",
           width = 3,height = 10,onefile = T)
 vln_eGFP1 
 dev.off()
 cairo_pdf("figures/03C_t_cell_x_markers_broad_1.pdf",
           width = 3,height = 12,onefile = T)
 vln_eGFP2
 dev.off()


