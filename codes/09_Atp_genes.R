rm(list = ls())

library(stringr)
library(ComplexHeatmap)
library(dplyr)


## Generate Source data figure 2


temp <- read.csv("output/new_tumors/enrichment_neg_tumor_late_vs_early.csv")

temp$Phenotype <- "DOWN"
temp <- temp[,c("ID","Description", "pvalue", "p.adjust","qvalue",
                "Phenotype","geneID")]
colnames(temp) <- c("Go.ID", "Description", "p.Val", "FDR","qvalue",
                    "Phenotype" , "Genes")

temp$Genes <- gsub("\\/",",",temp$Genes)
temp
#gsub("\\",",", temp)


temp2 <- read.table("data/Tumor_late_early--clustered default  node.csv",sep = ",", header = T)

temp2 <- temp2 |> dplyr::select(name,X__mclCluster)


temp2 <- dplyr::left_join(temp,temp2 , by = c("Go.ID" = "name"))

write.csv(temp2,"output/Supplementary Table 7.csv")





enrich_de     <- read.csv("output/new_tumors/tumor_late_vs_early.csv")
enrichmentMap <-
    read.csv("data/Tumor_late_early--clustered default  node.csv")

enrichmentMap <- enrichmentMap[enrichmentMap$X__mclCluster %in% c(2),]

late_tumors <-
    read.csv("data/enrichment_neg_tumor_late_vs_early.csv")
late_tumors <- late_tumors[late_tumors$ID %in% enrichmentMap$name,]

ATP_genes <- unique(unlist(str_split(paste0(late_tumors$geneID,collapse = "/"),
                                     pattern = "/")))



shared_genes <- as.data.frame(late_tumors)
shared_genes <- shared_genes[,c(2,3,9)]

genelists  <- shared_genes$geneID
genelists0 <- lapply(genelists,
                     function(X){unlist(str_split(X, pattern = "/"))})
allgenes   <- Reduce(union,genelists0)
genelists1 <- lapply(genelists0, function(X){ as.numeric(allgenes %in% X)})

genelists1 <- as.matrix(Reduce(rbind,genelists1))

rownames(genelists1) <- shared_genes$Description
colnames(genelists1)  <- allgenes

scaled_mat <- genelists1

enrich_de0  <- enrich_de[enrich_de$Gene %in% colnames(scaled_mat),]

write.csv(enrich_de0,"output/SourceData Ex1_e_1_Ox.csv")
rownames(enrich_de0) <- paste0(enrich_de0$Gene)
enrich_de0$Gene <- NULL

# 
# mat <- as.matrix(enrich_de0)
col_fun = circlize::colorRamp2(c(0,4,10, 100), c( "#f7f7f7",
                                                  "#d6604d",
                                                  "#ca0020",
                                                  "#67001f"))


col_fun_2 <- circlize::colorRamp2(c(-1,-0.5,-0.25,0),
                                  c("dodgerblue4","dodgerblue",
                                    "lightskyblue1","grey60"))

enrich_de0$p_value <- -log10(enrich_de0$p_val)


DE0 <- enrich_de0[colnames(scaled_mat),]
DE0 <- as.matrix(DE0[,c("p_value","avg_log2FC")])
rownames(DE0) <- rownames(enrich_de0)
DE0 <- as.data.frame(DE0)

ht1 <- Heatmap(scaled_mat,
           width =  unit(6.5,"inch"), 
           height = unit(5.7,"inch"),
           row_title = "",
           column_title = "ATP and oxidative phosphorylation",
           clustering_distance_columns =  "euclidean",
           clustering_method_columns = "ward.D2",
           col = circlize::colorRamp2(c(0,1),c("#e0e0e0","grey20") ),
           show_row_names = T, 
           row_names_side = "left",
           show_column_names = T, cluster_rows = F,
           show_row_dend = F,
           top_annotation = 
               HeatmapAnnotation("-log(p-value)" = DE0$p_value,
                                 "LogFC" = DE0$avg_log2FC,
                                 col = list("-log(p-value)" = col_fun,
                                         "LogFC" = col_fun_2),
                                 show_legend = T,
                                 show_annotation_name = T),
           show_column_dend = F,
           show_heatmap_legend = F,
           border_gp = gpar(col = "black", lty = 2),
           rect_gp = gpar(col = "white", lwd = 2),
           column_title_gp = gpar(fill = "white",col = "gray35",
                                  border = "white", cex = 1.5),
           row_names_gp = gpar(cex =1.5),
           column_names_gp = gpar(cex =1.2),
           column_title_side = "top",
           column_names_rot =  45,
           heatmap_legend_param = list(direction = "horizontal"))


enrichmentMap <-
    read.csv("data/Tumor_late_early--clustered default  node.csv")

enrichmentMap <-
    enrichmentMap[enrichmentMap$X__mclCluster %in% c(11),]

late_tumors <- 
    read.csv("data/enrichment_neg_tumor_late_vs_early.csv")

late_tumors <-
    late_tumors[late_tumors$ID %in% enrichmentMap$name,]

Mito_genes <- unique(unlist(str_split(paste0(late_tumors$geneID,collapse = "/"),
                                      pattern = "/")))





shared_genes <- as.data.frame(late_tumors)
shared_genes <- shared_genes[,c(2,3,9)]

genelists  <- shared_genes$geneID
genelists0 <- lapply(genelists,
                     function(X){unlist(str_split(X, pattern = "/"))})
allgenes   <- Reduce(union,genelists0)
genelists1 <- lapply(genelists0, function(X){ as.numeric(allgenes %in% X)})

genelists1 <- as.matrix(Reduce(rbind,genelists1))

rownames(genelists1) <- shared_genes$Description
colnames(genelists1)  <- allgenes

scaled_mat <- genelists1

enrich_de0  <- enrich_de[enrich_de$Gene %in% colnames(scaled_mat),]

write.csv(enrich_de0,"output/SourceData Ex1_e_2_MT.csv")

rownames(enrich_de0) <- paste0(enrich_de0$Gene)
enrich_de0$Gene <- NULL

# 
# mat <- as.matrix(enrich_de0)
col_fun = circlize::colorRamp2(c(0,4,10, 100), c( "#f7f7f7",
                                                  "#d6604d",
                                                  "#ca0020",
                                                  "#67001f"))


col_fun_2 <- circlize::colorRamp2(c(-1,-0.5,-0.25,0),
                                  c("dodgerblue4","dodgerblue",
                                    "lightskyblue1","grey60"))

enrich_de0$p_value <- -log10(enrich_de0$p_val)


DE0 <- enrich_de0[colnames(scaled_mat),]
DE0 <- as.matrix(DE0[,c("p_value","avg_log2FC")])
rownames(DE0) <- rownames(enrich_de0)
DE0 <- as.data.frame(DE0)

ht2 <- Heatmap(scaled_mat,
               width =  unit(6.5,"inch"), 
               height = unit(2.1,"inch"),
               row_title = "", column_title = "Mitochondrion organization",
               clustering_distance_columns =  "euclidean",
               clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(0,1),c("#e0e0e0","grey20") ),
               show_row_names = T, 
               row_names_side = "left",
               show_column_names = T, cluster_rows = F,
               show_row_dend = F,
               top_annotation = HeatmapAnnotation("-log(p-value)" = DE0$p_value,
                                  "LogFC" = DE0$avg_log2FC,
                                  col = list("-log(p-value)" = col_fun,
                                             "LogFC" = col_fun_2),
                                  show_legend = T, show_annotation_name = T),
               show_column_dend = F,
               show_heatmap_legend = F,
               border_gp = gpar(col = "black", lty = 2),
               rect_gp = gpar(col = "white", lwd = 2),
               column_title_gp = gpar(fill = "white",col = "gray35",
                                      border = "white", cex = 1.5),
               row_names_gp = gpar(cex =1.5),
               column_names_gp = gpar(cex =1.2),
               column_title_side = "top",
               column_names_rot =  45,
               heatmap_legend_param = list(direction = "horizontal"))


pdf("figures/Ex01E_ATP_mitochondrial.pdf", width = 22, height = 10)
ht1
ht2
dev.off()

