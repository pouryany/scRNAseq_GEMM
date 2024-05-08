# Plot the expression levels of pyroptosis and necroptosis associated genes

rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggthemes)


scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
scGemms$initial_clusts <- scGemms$seurat_clusters


table(scGemms$cell_type_secondary,scGemms$stage)


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

col_use <- c(rep(my_cols[1],3),
             rep(my_cols[2],3),
             rep(my_cols[11],1),
             rep(my_cols[4],2),
             rep(my_cols[6],8),
             rep(my_cols[7],1),
             rep(my_cols[8],1),
             rep(my_cols[9],1),
             rep(my_cols[10],1))


order_vec <-  c(
    "Alv-Proliferating","Alv-Basal",
    "Basal",
    "Alv 1","Alv 2","Alv 3",
    "Alv 4", 
    "Undefined Alv","Hormone Sensing",
    "CD4+ Naïve", "CD4+ Treg",             
    "CD8+ Naïve", "CD8+ CM","CD8+ ISG",
    "CD8+ Effector",
    "CD8+ TRM",
    "CD8+ Progenitor Ex",
    "CD8+ Ex", 
    "CD8+ Proliferating",
    "NK",
    "Macrophage",
    "Neutrophil","B cell",        
    "Fibroblast")    





colorss <- readRDS("data/03_colors.RDS")



colorss <- colorss[order_vec]









scGemms0 <- subset(scGemms, subset = cell_type_primary == "Immune")

# Filter mitochondrial and ribosomal genes
rm.ind <- grep("^mt-|^Rp[ls]",rownames(scGemms0))
keep   <- rownames(scGemms0)[-rm.ind]
scGemms0 <- scGemms0[keep,]

# Filter mitochondrial and ribosomal genes
keep.ind <- rowSums(scGemms0@assays$RNA@counts)
keep.ind <- names(keep.ind[keep.ind>=20])
#keep   <- rownames(scGemms0)[-rm.ind]

scGemms0 <- scGemms0[keep.ind,]


cell_tags <- unique(as.character(scGemms0$seurat_clusters))


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]




features = c("Gsdmd","Gsdme","Nlrp3","Ripk3")
my_plots <- list()
for(thisGene in features){
  
  
  temp_cells2 <- FetchData(scGemms0, vars  = thisGene)
  temp_cells2 <- temp_cells2 >0
  
  
  scGemms0[["my_marker"]] <- plyr::mapvalues((temp_cells2),
                                             c(FALSE,TRUE),c("Not expressed",
                                                             "Expressed"))
  
  data_meta  <- scGemms0@meta.data
  
  for (stage in c("Early", "Late")){
    data_tcell <- data_meta[data_meta$stage == stage,]
    
    data_tcell <- as.data.frame.matrix(table(data_tcell$cell_type_secondary,
                                             data_tcell$my_marker))
    
    order_vec <-  c(
      "CD4+ Naïve", "CD4+ Treg",             
      "CD8+ Naïve", "CD8+ CM","CD8+ ISG",
      "CD8+ Effector",
      "CD8+ TRM",
      "CD8+ Progenitor Ex",
      "CD8+ Ex", 
      "CD8+ Proliferating",
      "NK",
      "Macrophage",
      "Neutrophil","B cell")    
    data_tcell$tot <- rowSums(data_tcell)
    data_tcell <- data_tcell/data_tcell$tot
    
    data_tcell$cell_type_secondary <- factor(rownames(data_tcell))
    
    data_tcell <- data_tcell[order_vec,]
    data_tcell <- tidyr::pivot_longer(data_tcell,cols= c("Not expressed",
                                                         "Expressed"),
                                      names_to = thisGene)
    
    
    
    data_tcell$cell_type_secondary <- factor(x = data_tcell$cell_type_secondary,
                                             levels = rev(order_vec))
    
    
    stack_plot <-  ggplot(data_tcell, aes(x = (cell_type_secondary),
                                          y =100* value,  fill = !!sym((thisGene)),
    )) +
      geom_col() + scale_fill_excel_new() +
      theme_classic() +
      labs(x = "", y = "", 
           title = paste0(thisGene," in ",stage)) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
            legend.text = element_text(size = 16),
            axis.text.x = element_text(angle = 0,
                                       hjust = 0.5,
                                       size = 16),
            axis.text.y = element_text(angle = 0,
                                       hjust = 1,
                                       size = 20, face = "plain",
                                       colour = "black"),
            axis.line = element_blank(),
            axis.ticks = element_blank(), 
            title = element_text(size = 18),
            legend.title = element_blank()) +
      coord_flip()
    
    
    my_plots[[paste0(thisGene,"_",stage)]] <- stack_plot
    
    
  }

}

p1 <- patchwork::wrap_plots(my_plots, ncol = 2)
ggsave("figures/Sup_immune_gsdm.pdf",p1, width = 14, height = 20)





my_plots <- list()
for(thisGene in features){
  
  
  temp_cells2 <- FetchData(scGemms, vars  = thisGene)
  temp_cells2 <- temp_cells2 >0
  
  
  scGemms[["my_marker"]] <- plyr::mapvalues((temp_cells2),
                                             c(FALSE,TRUE),c("Not expressed",
                                                             "Expressed"))
  
  data_meta  <- scGemms@meta.data
  order_vec <-  c("Alv-Proliferating","Alv-Basal", "Basal")
  
  data_meta <- data_meta[data_meta$cell_type_secondary %in% order_vec,]
  data_meta$cell_type_secondary <- droplevels(data_meta$cell_type_secondary)
  for (stage in c("Early", "Late")){
    data_tcell <- data_meta[data_meta$stage == stage,]
    
    data_tcell <- as.data.frame.matrix(table(data_tcell$cell_type_secondary,
                                             data_tcell$my_marker))
    
    order_vec <-    c("Alv-Proliferating","Alv-Basal", "Basal")  
    data_tcell$tot <- rowSums(data_tcell)
    data_tcell <- data_tcell/data_tcell$tot
    
    data_tcell$cell_type_secondary <- factor(rownames(data_tcell))
    
    data_tcell <- data_tcell[order_vec,]
    data_tcell <- tidyr::pivot_longer(data_tcell,cols= c("Not expressed",
                                                         "Expressed"),
                                      names_to = thisGene)
    
    
    
    data_tcell$cell_type_secondary <- factor(x = data_tcell$cell_type_secondary,
                                             levels = rev(order_vec))
    
    
    stack_plot <-  ggplot(data_tcell, aes(x = (cell_type_secondary),
                                          y =100* value,  fill = !!sym((thisGene)),
    )) +
      geom_col() + scale_fill_excel_new() +
      theme_classic() +
      labs(x = "", y = "", 
           title = paste0(thisGene," in ",stage)) +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
            legend.text = element_text(size = 16),
            axis.text.x = element_text(angle = 0,
                                       hjust = 0.5,
                                       size = 16),
            axis.text.y = element_text(angle = 0,
                                       hjust = 1,
                                       size = 20, face = "plain",
                                       colour = "black"),
            axis.line = element_blank(),
            axis.ticks = element_blank(), 
            title = element_text(size = 18),
            legend.title = element_blank()) +
      coord_flip()
    
    
    my_plots[[paste0(thisGene,"_",stage)]] <- stack_plot
    
    
  }
  
}

p1 <- patchwork::wrap_plots(my_plots, ncol = 2)
ggsave("figures/Sup_tumor_gsdm.pdf",p1, width = 14, height = 10)



source("codes/utils.R")


cell_tags <-  c(
  "Alv-Basal", "Alv-Proliferating",
  "Basal")    


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
  "NK",
  "Macrophage",
  "Neutrophil","B cell")    


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



scGemms1 <- subset(scGemms,
                   cells = names(grep(pattern = "Alv-|Basal",
                                      scGemms$cell_type_secondary,value = T)))



order_vec <-  c( "Alv-Proliferating","Alv-Basal", "Basal")    




cell_tags <- unique(as.character(scGemms1$cell_type_secondary))
Idents(scGemms1) <- "cell_type_secondary"


scGemms1@active.ident <- factor(x = scGemms1@active.ident, levels = order_vec)


my_features = c("Epcam","Gsdmd","Gsdme","Nlrp3","Ripk3")


vln_eGFP1 <- VlnPlot_costumized(scGemms1,my_features = my_features,
                                idents_order = order_vec,
                                signif_list = main_list,
                                plot_alphabetic = F) +
  coord_fixed(0.3)

ggsave("figures/Sup_epithelial_gsdm_violin.pdf",vln_eGFP1, width = 14, height = 14)

order_vec <-  c(
  "CD4+ Naïve", "CD4+ Treg",             
  "CD8+ Naïve", "CD8+ CM","CD8+ ISG",
  "CD8+ Effector",
  "CD8+ TRM",
  "CD8+ Progenitor Ex",
  "CD8+ Ex", 
  "CD8+ Proliferating",
  "NK",
  "Macrophage",
  "Neutrophil","B cell")    


my_features = c("Ptprc","Gsdmd","Gsdme","Nlrp3","Ripk3")


Idents(scGemms0) <- "cell_type_secondary"
scGemms0@active.ident <- factor(x = scGemms0@active.ident, levels = order_vec)

vln_eGFP1 <- VlnPlot_costumized(scGemms0,my_features = my_features,
                                idents_order = order_vec,
                                signif_list = main_list,
                                plot_alphabetic = F) +
  coord_fixed(0.3)
  
ggsave("figures/Sup_immune_gsdm_violin.pdf",vln_eGFP1, width = 14, height = 14)




