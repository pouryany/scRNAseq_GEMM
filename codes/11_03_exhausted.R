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



scGemms0 <- subset(scGemms, subset = cell_type_primary == "Immune")

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
        "Macrophage",
        "Neutrophil","NK","B cell")    
scGemms0$cell_type_secondary <- factor(scGemms0$cell_type_secondary,
                                       levels = order_vec)
Idents(scGemms0) <- "cell_type_secondary"




cell_tags <- unique(as.character(scGemms0$cell_type_secondary))

Idents(scGemms0) <- "cell_type_secondary"

cluster0.markers <- FindMarkers(scGemms0, ident.1 = "CD8+ Ex",
                                ident.2 = "CD8+ Progenitor Ex",
                                only.pos = F,
                                test.use = "MAST")

cluster0.markers <- rownames_to_column(cluster0.markers) 
write.csv(cluster0.markers,"output/Tex_vs_TPex.csv")
