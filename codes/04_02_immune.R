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

for (ij in cell_tags){
        i <- gsub("/","_",ij)
        cluster0.markers <- FindMarkers(scGemms0, ident.1 = ij,
                                        only.pos = F,test.use = "MAST" )
        
        cluster0.markers <- tibble::rownames_to_column(cluster0.markers,var = "Gene")
        
        cluster0.markers <- cluster0.markers[cluster0.markers$p_val_adj<0.05 & 
                                                     abs(cluster0.markers$avg_log2FC) > 0.25,]
        dirOut <- paste0("output/immune_markers")
        if(!dir.exists(dirOut))
                dir.create(dirOut)
        saveTag <- paste0("output/immune_markers/ImmuneCluster_",i,".csv")
        write.csv(cluster0.markers,saveTag)
        
}


for (ij in cell_tags){
        i <- gsub("/","_",ij)
        cluster0.markers <- FindMarkers(scGemms0, ident.1 = "Late",
                                        group.by = "stage",
                                        subset.ident = ij,
                                        only.pos = F,
                                        test.use = "MAST")
        
        cluster0.markers <- cluster0.markers[cluster0.markers$p_val_adj<0.05 & 
                                                     abs(cluster0.markers$avg_log2FC) > 0.25,]
        cluster0.markers <- tibble::rownames_to_column(cluster0.markers,var = "Gene")
        dirOut <- paste0("output/immune_late_vs_early")
        if(!dir.exists(dirOut))
                dir.create(dirOut)
        saveTag <- paste0("output/immune_late_vs_early/DE_genes_immuneCluster_",i,".csv")
        write.csv(cluster0.markers,saveTag)
        
}



order_vec <-  c(
        "CD4+ Naïve", "CD4+ Treg",             
        "CD8+ Naïve", "CD8+ CM","CD8+ Effector",
        "CD8+ ISG",
        "CD8+ TRM",
        "CD8+ Progenitor Ex",
        "CD8+ Ex", 
        "CD8+ Proliferating",
        "Macrophage",
        "Neutrophil","NK","B cell")    






scGemms0$cell_type_secondary <- factor(scGemms0$cell_type_secondary,
                                       levels = order_vec)
Idents(scGemms0) <- "cell_type_secondary"

DotPlot(scGemms0, 
        features = c("Cd19","Cd79a", "Ms4a1", # B cell
                     "Ncr1","Cd160", # NK Cell
                     "S100a9","S100a8","Il1b","Cd14", # Neutrophil
                     "Cd68","Lyz2", # Macrophage
                     "Xcl1", "Cd7","Gzmb","Klra5","Trgv2","Fcer1g",#CD8 TRM
                     "Stmn1", "Top2a", "Mki67", # Cd8 Profilerating,
                     "Cd244a" ,"Pdcd1", # Cd8 Eff/Mem
                     "Isg15", "Ifit3", "Ifit1", "Irf7", # CD8+ IFN
                     "Ly6c2", "Nkg7", "Ctsw","Itgb1", # Cd8+ Eff early
                     "Cxcr6", "Cd8b1", "Cd8a","Cd44", # CD8+ CM
                     "Ctla4","Foxp3", # Treg
                     "Tcf7","Il7r","Sell","Ccr7","Cd3d","Cd4" # T cell Cd 4 
                     ),
        cols = c(low  = "#f7f7f7",
                 high = "#6e016b"),dot.min = .1)+
        scale_fill_gradient2(low  = "#2166ac",
                             mid  = "#f7f7f7",
                             high = "#b2182b") &
        coord_flip() & theme(axis.text.x = element_text(angle = 45,
                                                        hjust = 1, size = 14),
                             axis.text.y = element_text(size = 14),
                             axis.title = element_blank())


ggsave("figures/S3A_Immune_dotplot.pdf",width = 7, height = 12)

