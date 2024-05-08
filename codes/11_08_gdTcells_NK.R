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









scGemms0 <- subset(scGemms, subset = cell_type_secondary == "NK")

# Filter mitochondrial and ribosomal genes
rm.ind <- grep("^mt-|^Rp[ls]",rownames(scGemms0))
keep   <- rownames(scGemms0)[-rm.ind]
scGemms0 <- scGemms0[keep,]

# Filter mitochondrial and ribosomal genes



cell_tags <- unique(as.character(scGemms0$seurat_clusters))


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]




features = c("Itgae","Trac","Trbc1","Trbc2",
             "Trdc","Tcrg-C1","Tcrg-C2","Tcrg-C3","Tcrg-C4", "Tcrg-V1","Tcrg-V3",
             "Trgv2","Fcer1g","Ncr1")
my_plots <- list()

  
  
  temp_cells1 <- FetchData(scGemms0, vars  = features)
  
  # temp_cells2 <- temp_cells2 >0
  # colSums(temp_cells2)/nrow(temp_cells2)
  # data_tcell <- data_tcell/data_tcell$tot
  temp_cells2 <- rownames_to_column(temp_cells1,var="cell_id")
  temp_cells2 <- tidyr::pivot_longer(temp_cells2,!cell_id)
  temp_cells2$value <- temp_cells2$value > 0
  
  temp_cells2$value <- plyr::mapvalues((temp_cells2$value),
                                             c(FALSE,TRUE),c("Not expressed",
                                                             "Expressed"))
  
  
 
  data_meta  <- scGemms0@meta.data
  
  for (stage in c("Early", "Late")){
    data_tcell0 <- data_meta[data_meta$stage == stage,]
    
    data_tcell <- temp_cells2[temp_cells2$cell_id %in% rownames(data_tcell0),]
    
    

    data_tcell <- as.data.frame.matrix(table(data_tcell$name,
                                             data_tcell$value))
    
    
    
    data_tcell <- data_tcell/nrow(data_tcell0)
    data_tcell <- rownames_to_column(data_tcell,var = "Gene")
    data_tcell <- tidyr::pivot_longer(data_tcell,cols= c("Not expressed",
                                                         "Expressed"),
                                      names_to = "Status")
    
    
    
    data_tcell$Gene  <- factor(x = data_tcell$Gene, levels = rev(features))
    data_tcell$Status <- factor(data_tcell$Status, levels = c("Not expressed",
                                                              "Expressed"))
    
    stack_plot <-  ggplot(data_tcell, aes(x = (Gene),
                                          y =100* value,  fill = (Status),
    )) +
      geom_col() + scale_fill_excel_new() +
      theme_classic() +
      labs(x = "", y = "", 
           title = paste0("NK "," in ",stage)) +
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
    
    
    my_plots[[paste0("NK","_",stage)]] <- stack_plot
    
    
  }



p1 <- patchwork::wrap_plots(my_plots, ncol = 2)
ggsave("figures/Sup_gdt_nk.pdf",p1, width = 14, height = 10)





source("codes/utils.R")


cell_tags <-  c(
  "NK")    


main_list <- list()
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
                   cells = names(grep(pattern = "NK",
                                      scGemms$cell_type_secondary,value = T)))






cell_tags <- unique(as.character(scGemms1$cell_type_secondary))
Idents(scGemms1) <- "cell_type_secondary"


my_features = features


vln_eGFP1 <- VlnPlot_costumized(scGemms1,my_features = my_features,
                                idents_order = "NK",
                                signif_list = main_list,
                                plot_alphabetic = F) +
  coord_fixed(0.1)

ggsave("figures/Sup_nk_gdT_violin.pdf",vln_eGFP1, width = 10, height = 10)
