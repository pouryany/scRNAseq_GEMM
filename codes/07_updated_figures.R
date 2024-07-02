rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(EnhancedVolcano)


scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
scGemms$initial_clusts <- scGemms$seurat_clusters


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

t_cells <- grep("CD|NK", unique(scGemms0$cell_type_secondary),value = T)
t_cells <- subset(scGemms0, subset = cell_type_secondary %in% t_cells)


t_cells$cell_type_secondary <- droplevels(t_cells$cell_type_secondary)



t_cells[["immune_type"]] <-plyr::mapvalues(t_cells$cell_type_secondary, from = 
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

Idents(t_cells) <- "immune_type"
for (ij in cell_tags){
  i <- gsub("/","_",ij)
  cluster0.markers <- FindMarkers(t_cells, ident.1 = "Late",
                                  group.by = "stage",
                                  subset.ident = ij,
                                  only.pos = F,
                                  test.use = "MAST")
  
  cluster0.markers <- cluster0.markers[cluster0.markers$p_val_adj<0.05 & 
                                      abs(cluster0.markers$avg_log2FC) > 0.25,]
  cluster0.markers <- tibble::rownames_to_column(cluster0.markers,var = "Gene")
  dirOut <- paste0("output/broad_immune")
  if(!dir.exists(dirOut))
    dir.create(dirOut)
  saveTag <- paste0("output/broad_immune/DE_genes_immuneCluster_",i,".csv")
  write.csv(cluster0.markers,saveTag)
  
}


order_vec0 <-  c(
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



markers_list <- rev(c(
                  "Gzmb","Trgv2","Fcer1g",
                  "Mki67", # Cd8 Profilerating,
                  "Pdcd1", # Cd8 Eff/Mem
                  "Isg15", # CD8+
                  "Nkg7", # Cd8+ Eff early
                  "Cxcr6", # Treg
                  "Il7r","Sell","Ccr7",
                  "Cd44","Ncr1",
                  "Cd8b1", "Cd8a", # CD8+ CM
                  "Foxp3","Cd4" # T cell Cd 4 
))



Idents(scGemms) <- "cell_type_secondary"

scGemms@active.ident <- factor(scGemms@active.ident,levels = order_vec)
temp <- as.data.frame.matrix(table(scGemms$cell_type_secondary,scGemms$stage))


write.csv(x =temp[order_vec,],
          file = "output/cell_counts.csv")

cluster0.markers <-
  read.csv("output/immune_late_vs_early/DE_genes_immuneCluster_NK.csv")



keyvals <- ifelse(
  cluster0.markers$avg_log2FC < 0, 'royalblue',
  ifelse(cluster0.markers$avg_log2FC > 0, 'red3',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red3'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'


sel_labs <- (cluster0.markers$Gene)
sel_labs <- grep(paste0("Spp1|Ifng|Igflr1|Eomes|Ifi|Isg|Irf|",
                        "Cxc|Zbp|Gbp|Stat|Irgm|Irg|Ly6|Ddx|",
                        "Pd|Ccl|Gzm|Ifn|Klr|Klf"),
                 sel_labs,value = T)

p1 <-   EnhancedVolcano(cluster0.markers,
                        lab = cluster0.markers$Gene, #res1$hgnc_symbol,
                        selectLab = sel_labs,
                        x = "avg_log2FC",
                        y = "p_val_adj",
                        xlim = c(-2,3),
                        ylim = c(0, 50),
                        caption = NULL,
                        axisLabSize = 30,
                        xlab = bquote(~Log[2]~ "fold change"),
                        ylab = bquote(~-Log[10]~"adjusted"~italic(p)),
                        title = NULL,
                        subtitle = "",
                        pCutoff = 0.05,
                        pointSize = 2.0,
                        labSize = 9,
                        labCol = 'black',
                        labFace = 'plain',
                        drawConnectors = T,
                        arrowheads = F,
                        boxedLabels = F,
                        colCustom = keyvals,
                        FCcutoff = (0.25),
                        legendPosition = "none") 


ggsave("figures/S03A_nk_DE_volcano.pdf",p1,width = 12,height = 10)



