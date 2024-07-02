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




features = c("Itgae","Trdc","Trgv2","Fcer1g","Trac","Cd3d", "Cd8a","Cd8b1")
my_plots <- list()
my_tab   <- list()

for(thisGene in features){
  temp_cells2 <- FetchData(scGemms0, vars  = thisGene)
  temp_cells2 <- temp_cells2 >0
  
  scGemms0[["my_marker"]] <- plyr::mapvalues((temp_cells2),
                                             c(FALSE,TRUE),c("Not expressed",
                                                             "Expressed"))
  
  data_meta  <- scGemms0@meta.data
  order_vec <-  c("CD8+ TRM")
  
  data_meta <- data_meta[data_meta$cell_type_secondary %in% order_vec,]
  data_meta$cell_type_secondary <- droplevels(data_meta$cell_type_secondary)
  for (stage in c("Early", "Late", "Overall")){
    data_tcell <- data_meta
    if(stage != "Overall")
    data_tcell <- data_meta[data_meta$stage == stage,]
    
    data_tcell <- as.data.frame.matrix(table(data_tcell$cell_type_secondary,
                                             data_tcell$my_marker))
    
    order_vec <-  c("CD8+ TRM")    
    data_tcell$tot <- rowSums(data_tcell)
    data_tcell <- data_tcell/data_tcell$tot
    
    data_tcell$cell_type_secondary <- factor(rownames(data_tcell))
    
    data_tcell <- data_tcell[order_vec,]
    data_tcell <- tidyr::pivot_longer(data_tcell,cols= c("Not expressed",
                                                         "Expressed"),
                                      names_to = thisGene)
    
    
    
    data_tcell$cell_type_secondary <- factor(x = data_tcell$cell_type_secondary,
                                             levels = rev(order_vec))
    
    data_tcell[,thisGene] <- factor(x =  unlist(data_tcell[,thisGene]), levels = c("Not expressed",
                                                                           "Expressed"))
    
    
    
    temp_tab <- data_tcell
    temp_tab$tot <- (thisGene)
    names(temp_tab) <- c("tot","cell_type_secondary","Status","value")
    temp_tab$stage  <- stage
    my_tab <- rbind(my_tab,temp_tab)
    
  }

}

my_tab$value <- signif(my_tab$value,3)

write.csv(my_tab, "output/TRM_proportions.csv")

