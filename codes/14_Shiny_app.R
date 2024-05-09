### Instalaation stuf
# 
# reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2",
#            "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
# newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
# if(length(newPkg)){install.packages(newPkg)}
# 
# reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r",
#            "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
# newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
# if(length(newPkg)){install.packages(newPkg)}
rm(list = ls())
#devtools::install_github("SGDDNB/ShinyCell")
library(Seurat)
library(ShinyCell)

seu = readRDS("data/03_umap_labled_1_sObj.RDS")
my_cols = readRDS("data/03_colors.RDS")
seu@meta.data$sample <- seu@meta.data$sampleID
seu@meta.data$sampleID <- colnames(seu)






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


seu@meta.data$cell_type_secondary <- factor(x = seu@meta.data$cell_type_secondary,
                                            levels = order_vec)





scConf = createConfig(seu)
colnames(seu@meta.data)






scConf = modColours(scConf, meta.to.mod = "cell_type_secondary", 
                    new.colours= my_cols[order_vec])

scales::show_col(my_cols) 

scConf = delMeta(scConf, c("orig.ident", "RNA_snn_res.2","keep",
                           "seurat_clusters","doublets",
                           "cell_type_detail","color"))


scConf = modDefault(scConf, "cell_type_secondary", "cell_type_primary")

scConf = modMetaName(scConf, meta.to.mod = c("stage","cell_type_secondary","cell_type_primary"), 
                      new.name = c("Time since induction", "cell subtype","primary cell type"))


showLegend(scConf)



scConf = modColours(scConf, meta.to.mod = "cell subtype", 
                    new.colours= my_cols[order_vec])




makeShinyApp(seu, scConf, gene.mapping = FALSE,
             shiny.title = "GEMM Single cell") 
