
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Data source PMID: 35444279, GSE195937
temp <- readRDS("data/07_GSE195937_Mouse_breast_Fig1_scRNA_SeuratObj.RDS")


temp <- UpdateSeuratObject(temp)



library(AnnotationDbi)
library(org.Mm.eg.db)
my_genes <- select(org.Mm.eg.db, keys = rownames(temp), keytype = 'ENSEMBL', columns = 'SYMBOL')

my_genes <- my_genes[!is.na(my_genes$SYMBOL),]

my_genes <- my_genes[(!duplicated(my_genes$SYMBOL)),]
my_genes <- my_genes[(!duplicated(my_genes$ENSEMBL)),]


temp <- temp[my_genes$ENSEMBL,]


rownames(temp@assays$RNA@counts)  <- my_genes$SYMBOL
rownames(temp@assays$RNA@data)  <- my_genes$SYMBOL
rownames(temp@assays$RNA@scale.data)  <- my_genes$SYMBOL

Idents(temp)

p00 <- DimPlot(temp,
               combine = T, label =T, pt.size = 1) & NoAxes() & theme(legend.position="bottom")


# Stage refers to time since tumor induction
temp[["Cluster_Label"]]    <- plyr::mapvalues(temp$Cluster,
                                         0:4,
                                         c("C2","C1","C3","C5",
                                           "C4"))

temp$Cluster_Label <- factor(temp$Cluster_Label,levels= (paste0("C",1:5)))

temp[["Cluster_Label"]]    <- plyr::mapvalues(temp$Cluster_Label,
                                              paste0("C",1:5),
                                              c("C1 (Naive)","C2 (Ex)",
                                                "C3 (ILTCK)","C4 (ISG)", 
                                                "C5 (Prolif)")
                                              )
Idents(temp) <- "Cluster_Label"
rownames(temp)

my_features2 <- c("Fcer1g","Xcl1","Itgae","Itga1","Klrb1a","Klrb1c",
                  "Gzma","Gzmb","Clnk","Emid1","Chn2")

p1 <- FeaturePlot(temp,features = my_features2,ncol = 2, pt.size = 1,
                  combine = T, cols = c("#d6d6d6","#6e016b"), min.cutoff = "q1")  &
    NoAxes() & theme(legend.position="bottom") & 
    guides(colourbar = guide_legend(override.aes = list(size=12, nrow = 1))) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))




p1 <- VlnPlot(temp,features = my_features2, fill.by = "ident",
              stack = T,flip = T,ncol = 3,pt.size = 0, same.y.lims = T, combine = T)&
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16),
          axis.title.y = element_blank(),  axis.text.y = element_text(size = 16),
          strip.text =element_text(size = 16, face = "plain",vjust = 1,hjust = 0),
          plot.title = element_text(size = 30, face = "plain", hjust = 0)) & coord_fixed(ratio = 0.4)

ggsave("figures/Sup_LiLab_ILTCK_violin.pdf",p1, width = 14, height = 14)


# my_features2 <- c("Itgae","Trac","Trbc1","Trbc2",
#                   "Trdc","Tcrg-C1","Tcrg-C2","Tcrg-C3","Tcrg-C4", "Tcrg-V1","Tcrg-V3",
#                   "Trgv2","Fcer1g","Ncr1")
# p1 <- VlnPlot(temp,features = my_features2, fill.by = "ident",
#               stack = T,flip = T,ncol = 3,pt.size = 0, same.y.lims = T, combine = T)&
#     theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16),
#           axis.title.y = element_blank(),  axis.text.y = element_text(size = 16),
#           strip.text =element_text(size = 16, face = "plain",vjust = 1,hjust = 0),
#           plot.title = element_text(size = 30, face = "plain", hjust = 0)) & coord_fixed(ratio = 0.4)
# 
