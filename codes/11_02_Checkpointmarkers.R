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
    "NK", "Macrophage","Neutrophil")    


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



 
 scGemms0 <- subset(scGemms,
                    cells = names(grep(pattern = "Alv-|Basal|Macrophage|Neutrophil",
                              scGemms$cell_type_secondary,value = T)))
 
 
 
 order_vec <-  c( "Alv-Proliferating","Alv-Basal",  "Basal",
                  "Macrophage", "Neutrophil")    
 
 
 
 
 cell_tags <- unique(as.character(scGemms0$cell_type_secondary))
 Idents(scGemms0) <- "cell_type_secondary"
 
 
 scGemms0@active.ident <- factor(x = scGemms0@active.ident, levels = order_vec)

 
 my_features <- c("Cd274", "Pdcd1lg2","Pvr", "Nectin2", 
                  "Lgals9","Fgl1")
 
 vln_eGFP1 <- VlnPlot_costumized(scGemms0,my_features = my_features,
                                 idents_order = order_vec,
                                 signif_list = main_list,
                                 plot_alphabetic = F) +
   coord_fixed(0.3)
 

 cairo_pdf("figures/Sup_checkpoint_ligands.pdf",
           width = 4,height = 8,onefile = T)
 vln_eGFP1 
 dev.off()
 
 