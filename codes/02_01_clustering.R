rm(list = ls())


library(ggplot2)
library(harmony)
library(RColorBrewer)
library(Seurat)
library(stringr)
# Color libraries
library("wesanderson")
library("palettetown")
library("rcartocolor")
library("cartography")
library("pals")


scGemms  <- readRDS("data/02_filtered_sObj.RDS")

scGemms$seurat_clusters         <- droplevels(scGemms$RNA_snn_res.2)
levels(scGemms$seurat_clusters) <- 0:39

Idents(scGemms) <- "seurat_clusters"


scGemms[["cell_type_primary"]] <- plyr::mapvalues(scGemms$seurat_clusters,
                                                  0:39,
                                                  c("Immune",
                                                    "Epithelial",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Epithelial",
                                                    "Epithelial",
                                                    "Epithelial",
                                                    "Epithelial",
                                                    "Fibroblast",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Epithelial",
                                                    "Immune",
                                                    "Epithelial",
                                                    "Epithelial",
                                                    "Fibroblast",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Epithelial",
                                                    "Immune",
                                                    "Fibroblast",
                                                    "Immune",
                                                    "Immune",
                                                    "Immune",
                                                    "Epithelial",
                                                    "Fibroblast",
                                                    "Epithelial",
                                                    "Immune",
                                                    "Epithelial"))




my_cols <- 
c("#a6cee3", "#1f78b4",
"#b2df8a", "#33a02c",
"#fb9a99", "#e31a1c", "#fdbf6f",
"#ff7f00", "#cab2d6",
"#6a3d9a", "#ffff99",
"#b15928", "#8dd3c7",
"#ffffb3", "#bebada",
"#fb8072", "#80b1d3",
"#fdb462", "#b3de69", 
"#fccde5", "#d9d9d9",
"#bc80bd","#ccebc5",
"#ffed6f","#ffed6f" )




scGemms[["cell_type_secondary"]] <-plyr::mapvalues(scGemms$seurat_clusters,
                                       0:39,
                                       c("T cell: CD4+",
                                         "Epi: Tumor HER2+ Alv 1",
                                         "T cell: CD8+ CM",
                                         "T cell: CD4+",
                                         "T cell: CD8+ Naïve",
                                         "Neutrophil",
                                         "Epi: Tumor HER2+ Alv 2",
                                         "Epi: Tumor HER2+ Alv-Basal",
                                         "Epi: Tumor HER2+ Basal",
                                         "Epi: Alv undefined",
                                         "Fibroblast",
                                         "T cell: CD4+",
                                         "T cell: CD4+",
                                         "NK Cell",
                                         "T cell: CD8+ Effector",
                                         "T cell: CD8+ Progenitor Ex",
                                         "Macrophage",
                                         "Neutrophil",
                                         "T cell: CD8+ Ex",
                                         "Epi: Hormone Sensing",
                                         "T cell: CD8+ TRM",
                                         "Epi: Tumor HER2+ Alv 3",
                                         "Epi: Tumor HER2+ Alv-Proliferating",
                                         "Fibroblast",
                                         "T cell: CD8+ Proliferating",
                                         "T cell: CD4+",
                                         "NK Cell", 
                                         "Macrophage",
                                         "Macrophage",
                                         "Epi: Tumor HER2+ Alv 4",
                                         "T cell: CD4+ Treg",
                                         "Fibroblast",
                                         "T cell: CD8+ ISG",
                                         "B cell",
                                         "T cell: CD8+ CM",
                                         "Epi: Hormone Sensing",
                                         "Fibroblast",
                                         "Epi: Alv undefined",
                                         "Macrophage",
                                         "Epi: Hormone Sensing"))


scGemms[["cell_type_detail"]] <- scGemms$cell_type_secondary

scGemms[["cell_type_secondary"]] <-plyr::mapvalues(scGemms$cell_type_secondary,
                                       from = 
                                        c( "T cell: CD4+",
                                           "T cell: CD8+ Effector",
                                           "T cell: CD8+ CM",
                                           "T cell: CD8+ Naïve",
                                           "Neutrophil",
                                           "T cell: CD8+ ISG",
                                           "B cell",
                                           "NK Cell",
                                           "Macrophage",
                                           "T cell: CD4+ Treg",
                                           "T cell: CD8+ Proliferating",
                                           "Epi: Alv undefined",
                                           "T cell: CD8+ Progenitor Ex",
                                           "T cell: CD8+ Ex",
                                           "Fibroblast",
                                           "Epi: Hormone Sensing",
                                           "T cell: CD8+ TRM",                 
                                           "Epi: Tumor HER2+ Basal",
                                           "Epi: Tumor HER2+ Alv 2",
                                           "Epi: Tumor HER2+ Alv 1",            
                                           "Epi: Tumor HER2+ Alv-Basal",
                                           "Epi: Tumor HER2+ Alv-Proliferating",
                                           "Epi: Tumor HER2+ Alv 3",
                                           "Epi: Tumor HER2+ Alv 4"), 
                                to =    c( "CD4+ Naïve",
                                           "CD8+ Effector",
                                           "CD8+ CM",
                                           "CD8+ Naïve",
                                           "Neutrophil",
                                           "CD8+ ISG",
                                           "B cell",
                                           "NK",
                                           "Macrophage",
                                           "CD4+ Treg",
                                           "CD8+ Proliferating",
                                           "Undefined Alv",
                                           "CD8+ Progenitor Ex",
                                           "CD8+ Ex",    
                                           "Fibroblast",
                                           "Hormone Sensing",
                                           "CD8+ TRM",                 
                                           "Basal",
                                           "Alv 2",
                                           "Alv 1",            
                                           "Alv-Basal",
                                           "Alv-Proliferating",
                                           "Alv 3",
                                           "Alv 4"))


scGemms[["stage"]]    <- plyr::mapvalues(scGemms$stage,
                                         c("early","late"),
                                         c("Early","Late"))





cel_cols <- wes_palette("IsleofDogs1",
                        type = "discrete")


rcartocolor::carto_pal(3, "Magenta")

new_cols <- c(rcartocolor::carto_pal(5, "Magenta")[2:4],
              rcartocolor::carto_pal(4, "Mint"),
              c("wheat4","wheat3"),
              c('#c7e9c0','palegreen4'),
              c('mistyrose1',"rosybrown1","salmon",
                'firebrick1','firebrick3',"#980043",
                '#67001f',"black"),
              "#045a8d",
              rcartocolor::carto_pal(5, "Safe")[1],
              "darkgoldenrod3","sienna3","yellow4")


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



scGemms[["color"]] <-plyr::mapvalues(as.character(scGemms$cell_type_secondary),
                                     order_vec,
                                     new_cols)

colorss <- new_cols
names(colorss) <- order_vec


#> Saving Colors, processed seurat object, and metadata
saveRDS(colorss,"data/03_colors.RDS")
saveRDS(scGemms,"data/03_umap_labled_1_sObj.RDS")
saveRDS(scGemms@meta.data,"data/03_umap_labled_1_metadata.RDS")
