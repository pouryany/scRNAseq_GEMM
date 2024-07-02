rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)


scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
scGemms$initial_clusts <- scGemms$seurat_clusters
colorss <- readRDS("data/03_colors.RDS")




t_cells <- grep("T cell|NK", unique(scGemms$cell_type_detail),value = T)
scGemms0 <- subset(scGemms, subset = cell_type_detail %in% t_cells)



cell_tags <- unique(as.character(scGemms0$seurat_clusters))

# Filter mitochondrial and ribosomal genes
rm.ind <- grep("^mt-|^Rp[ls]",rownames(scGemms0))
keep   <- rownames(scGemms0)[-rm.ind]
scGemms0 <- scGemms0[keep,]

# Filter mitochondrial and ribosomal genes
keep.ind <- rowSums(scGemms0@assays$RNA@counts)
keep.ind <- names(keep.ind[keep.ind>=20])
#keep   <- rownames(scGemms0)[-rm.ind]
scGemms0 <- scGemms0[keep.ind,]



scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]

order_vec <-  c(
        "CD4+ Naïve", "CD4+ Treg",             
        "CD8+ Naïve", "CD8+ CM",
        "CD8+ ISG","CD8+ Effector",
        "CD8+ TRM",
        "CD8+ Progenitor Ex",
        "CD8+ Ex", 
        "CD8+ Proliferating",
        "NK")    

scGemms0$cell_type_secondary <- factor(scGemms0$cell_type_secondary,
                                       levels = order_vec)
Idents(scGemms0) <- "cell_type_secondary"




cell_tags <- unique(as.character(scGemms0$cell_type_secondary))

Idents(scGemms0) <- "cell_type_secondary"




main_list <- list()
for (ij in order_vec){
        i <- gsub("/","_",ij)
        saveTag <- 
                paste0("output/immune_markers/ImmuneCluster_",
                       i,".csv")
        
        cluster0.markers <- read.csv(saveTag, row.names = 1) 
        cluster0.markers <- cluster0.markers[cluster0.markers$avg_log2FC > 0,]
        
        if (any(cluster0.markers$p_val == 0)){
                
                
        sel_inds <- max(sum(cluster0.markers$p_val == 0),10)   
        cluster0.markers <- cluster0.markers[1:sel_inds,]
        cluster0.markers <- cluster0.markers[order(cluster0.markers$avg_log2FC,
                                                   decreasing = T),]
        }
        cluster0.markers <- cluster0.markers[1:10,]
        cluster0.markers$cell <- ij
        
        main_list <- rbind(main_list,cluster0.markers)
        
}





scGemms0$cell_type_secondary <- factor(scGemms0$cell_type_secondary,
                                       levels = order_vec)
Idents(scGemms0) <- "cell_type_secondary"

p1 <- DotPlot(scGemms0, 
        features = c(unique(main_list$Gene)
                     ),
        cols = c(low  = "#f7f7f7",
                 high = "#6e016b"),dot.min = .1)+
        scale_fill_gradient2(low  = "#2166ac",
                             mid  = "#f7f7f7",
                             high = "#b2182b") &
        coord_flip() & theme(axis.text.x = element_text(angle = 45,
                                                        hjust = 1, size = 12),
                             axis.text.y = element_text(size = 12),
                             axis.title = element_blank())



ggsave("figures/Ex04_immune_dotplot.pdf", p1,width = 7, height = 16)
write.csv(main_list,"output/top_markers_Tcells.csv")



