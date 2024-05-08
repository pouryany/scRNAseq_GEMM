rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(forcats)
library(viridis)
library(CytoTRACE)



scGemms  <- readRDS("data/03_umap_labled_1_sObj.RDS")
colorss  <- readRDS("data/03_colors.RDS")
scGemms0 <- subset(scGemms, subset = cell_type_primary == "Epithelial")
scGemms0 <- subset(scGemms, cells = names(grep(pattern = "HER2",
                                   scGemms$cell_type_detail,value = T)))


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]
exp_mat <- scGemms0@assays$RNA@counts




results  <- CytoTRACE(as.matrix(exp_mat),ncores = 12,enableFast = F)
scGemms0 <- AddMetaData(scGemms0,results$CytoTRACE,col.name = "cytoTRACE_Score")

Idents(scGemms0) <- "cell_type_secondary"

p00 <- FeaturePlot(scGemms0,features = "cytoTRACE_Score", label = T, pt.size = 3, label.size = 8)+ 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))&  NoAxes() &
    labs(title = "CytoTRACE maturation score") & theme(title = element_text(hjust = 1)) &
    ylim(-12,-4) & xlim(-8,10)

ggsave("figures/01G_maturation_tumor.pdf", height = 6 , width = 13)

saveRDS(scGemms0,"data/04_cytoTRACE_tumors.RDS")

