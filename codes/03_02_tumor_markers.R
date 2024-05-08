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


scGemms0 <- subset(scGemms, subset = cell_type_primary == "Epithelial")
scGemms0 <- subset(scGemms, cells = names(grep(pattern = "HER2",
                                   scGemms$cell_type_detail,value = T)))

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
        "Alv-Basal", "Alv-Proliferating",
        "Basal",
        "Alv 1", "Alv 2","Alv 3", 
        "Alv 4")    


cell_tags <- unique(as.character(scGemms0$cell_type_secondary))

Idents(scGemms0) <- "cell_type_secondary"

for (ij in cell_tags){
        
        i <- gsub("/","_",ij)
        cluster0.markers <- FindMarkers(scGemms0, ident.1 = ij,
                                        only.pos = F,test.use = "MAST" )
        
        cluster0.markers <- tibble::rownames_to_column(cluster0.markers,var = "Gene")
        
        cluster0.markers <-
                cluster0.markers[cluster0.markers$p_val_adj<0.05 & 
                                         abs(cluster0.markers$avg_log2FC) > 0.25,]
        dirOut <- paste0("output/tumor_markers")
        if(!dir.exists(dirOut))
                dir.create(dirOut)
        saveTag <- paste0("output/tumor_markers/tumorCluster_",i,".csv")
        write.csv(cluster0.markers,saveTag)
}



for (ij in cell_tags){
        i <- gsub("/","_",ij)
        
        tryCatch({
                cluster0.markers <- FindMarkers(scGemms0, ident.1 = "Late",
                                                group.by = "stage",
                                                subset.ident = ij,
                                                only.pos = F,
                                                test.use = "MAST",
                                                latent.vars = c("nCount_RNA") )
                
                cluster0.markers <- 
                        cluster0.markers[cluster0.markers$p_val_adj<0.05 & 
                                                 abs(cluster0.markers$avg_log2FC) > 0.25,]
                
                
                cluster0.markers <- tibble::rownames_to_column(cluster0.markers,
                                                               var = "Gene")
                
                dirOut <- paste0("output/tumor_late_vs_early")
                if(!dir.exists(dirOut))
                        dir.create(dirOut)
                saveTag <- paste0("output/",
                                  "tumor_late_vs_early/",
                                  "DE_genes_tumorCluster_",i,".csv")
                write.csv(cluster0.markers,saveTag) 
        }
        ,
        error = function(e){print("passsing")})
}


write.csv(rownames(scGemms0),file = "output/backgroundHER2.csv")
