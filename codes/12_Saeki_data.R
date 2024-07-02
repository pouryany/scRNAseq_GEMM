rm(list = ls())
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Data Source PMID: 34079055
# download the Harmony data from the following url 
# (https://mouse-mammary-epithelium-integrated.cells.ucsc.edu)
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


features = c("Gsdmd","Nlrp3","Ripk3")
bigTab <- list()
for(thisGene in features){
    
    
    temp_cells2 <- FetchData(temp, vars  = thisGene)
    temp_cells2 <- temp_cells2 >0
    
    
    temp[["my_marker"]] <- plyr::mapvalues((temp_cells2),
                                              c(FALSE,TRUE),c("Not expressed",
                                                              "Expressed"))
    
    data_meta  <- temp@meta.data

    
        data_tcell <- data_meta
        
        data_tcell <- as.data.frame.matrix(table(data_tcell$Annotation,
                                                 data_tcell$my_marker))
        
        data_tcell$tot <- rowSums(data_tcell)
        data_tcell <- data_tcell/data_tcell$tot
        
        data_tcell$cell_type <- factor(rownames(data_tcell))
        

 
        data_tcell$gene <- thisGene
        rownames(data_tcell) <- NULL
        
        bigTab <- rbind(bigTab,data_tcell)
        
    
    
}




