
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Data Source PMID: 34079055
temp <- readRDS("~/Downloads/Harmony.rds")


temp <- UpdateSeuratObject(temp)

temp$Annotation
Idents(temp) <- "Annotation"

p00 <- DimPlot(temp,
               combine = T, label =T, pt.size = 1) & NoAxes() & 
    theme(legend.position="bottom")


p1 <- VlnPlot(temp,features = c( "Epcam","Gsdmd","Gsdme", "Nlrp3", "Ripk3"),
              fill.by = "ident",
              stack = T,flip = T,ncol = 1,pt.size = 0, same.y.lims = T,
              combine = T)&
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 20),
          axis.title.y = element_blank(),  axis.text.y = element_text(size = 16),
          strip.text   = element_text(size = 20, face = "plain",vjust = 1,
                                   hjust = 0),
          plot.title = element_text(size = 30, face = "plain", hjust = 0)) &
    coord_fixed(ratio = 0.8)

ggsave("figures/Sup_Saeki_mammary_gsdm_violin.pdf",p1, width = 14, height = 14)

