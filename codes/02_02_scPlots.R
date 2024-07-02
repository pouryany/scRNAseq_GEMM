rm(list = ls())

library(Seurat)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(wesanderson)
scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
colorss <- readRDS("data/03_colors.RDS")


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

cel_cols <- wes_palette("IsleofDogs1",
                        type = "discrete")


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




p00 <-  DimPlot(scGemms,group.by = "cell_type_secondary",
                cols = colorss[order_vec],
                combine = T, label = F, pt.size = 0.5) & 
        NoAxes() &
        labs(title = "UMAP by Cell Type") &
        theme(plot.title = element_text(hjust = 0, size = 30, face = "plain"), 
              legend.text = element_text(size = 28),
              axis.title.y.right = element_text(size = 20)) &
        guides(color = guide_legend(override.aes = list(size=8), ncol=1) )


p01 <- DimPlot(scGemms,group.by = "sampleID", cols = my_cols,
               combine = T, label = F, pt.size = 0.4) &
        NoAxes() &
        ylim(-13,13) &
        labs(title = "UMAP by Sample ID") & 
        theme(plot.title = element_text(hjust = 0, size = 30, face = "plain"), 
              legend.text = element_text(size = 28),
              axis.title.y.right = element_text(size = 20)) &
        guides(color = guide_legend(override.aes = list(size=8), ncol=1) )

p02 <- DimPlot(scGemms,group.by = "cell_type_primary", cols = cel_cols,
               combine = T, label = F, pt.size = 0.4) &
        NoAxes() & 
        ylim(-13,13) &
        labs(title = "UMAP by Primary Type") & 
        theme(plot.title = element_text(hjust = 0, size = 30, face = "plain"), 
              legend.text = element_text(size = 28),
              axis.title.y.right = element_text(size = 20)) &
        guides(color = guide_legend(override.aes = list(size=8), ncol=1) )



p03 <- DimPlot(scGemms,group.by = "stage", cols = my_cols[c(2,7)],
               combine = T, label = F, pt.size = 0.4) &
        NoAxes() &
        ylim(-13,13) &
        labs(title = "Time since induction") & 
        theme(plot.title = element_text(hjust = 0, size = 30, face = "plain"), 
              legend.text = element_text(size = 28),
              axis.title.y.right = element_text(size = 20)) &
        guides(color = guide_legend(override.aes = list(size=8), ncol=1) )




((p02/p01 /p03) | p00)  +  plot_layout(widths = c(1, 3))


#> Associated figures with Figures 1B, 1C
ggsave("figures/01B_1C_umap.pdf", height = 12 , width = 20 )








### Exporting metadata


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




scMetadata <- scGemms@meta.data



colnames(scMetadata)



scales::show_col(unique(scMetadata$color))
class(scMetadata)
scMetadata <- scMetadata[,setdiff(colnames(scMetadata),c("orig.ident","keep",
                                          "seurat_clusters","doublets",
                                          "cell_type_detail"))]

new_column_names <- c("stage" = "Time since induction", 
                      "cell_type_secondary" = "cell subtype",
                      "cell_type_primary" ="primary cell type")

# Update the column names
colnames(scMetadata)[match(names(new_column_names), colnames(scMetadata))] <- new_column_names
umapCoord <- as.data.frame(Embeddings(object = scGemms[["umap"]]))

all(rownames(umapCoord) == rownames(scMetadata))
scMetadata <- cbind(scMetadata,umapCoord)

write.csv(scMetadata,"output/SourceData1.csv")
