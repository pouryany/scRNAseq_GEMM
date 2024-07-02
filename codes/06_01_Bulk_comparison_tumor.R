rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(tidyr)
library(VennDiagram)
library(enrichplot)
library(forcats)
library(org.Mm.eg.db)
library(clusterProfiler)
library(forcats)
library(EnhancedVolcano)

scGemms_meta = readRDS("data/03_umap_labled_1_metadata.RDS")

DE_4T1  <- read.csv("data/Table5_DEGenes_4T1.csv")




mouse_genes <- c("Tnf", "Ifnb1", "Irf7", "Ddx58" ,
                 "Il6" ,
                 "Ripk3")
#mouse_genes <- c(rownames(res1)[1:5], mouse_genes)

cluster0.markers <- DE_4T1

keyvals <- ifelse(
        cluster0.markers$Direction =="DOWN" , 'royalblue',
        ifelse(cluster0.markers$Direction =="UP", 'red3',
               'grey35'))
keyvals[is.na(keyvals)] <- 'grey35'
names(keyvals)[keyvals == 'red3'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'


sel_labs <- cluster0.markers[cluster0.markers$Direction != "NONE",]$X
sel_labs <- grep("Ifi|Isg|Irf|Cxcl10|Zbp|Gbp|Stat|Irgm|Irgm|Ddx58|Ccl|Gzm|Ifn",
                 sel_labs,value = T)

p1 <-   EnhancedVolcano(cluster0.markers,
                        lab = cluster0.markers$X, #res1$hgnc_symbol,
                        selectLab = sel_labs,
                        x = "log2FoldChange",
                        y = "padj",
                        xlim = c(-12,12),
                        ylim = c(0, 310),
                        caption = NULL,
                        axisLabSize = 30,
                        xlab = bquote(~Log[2]~ "fold change"),
                        ylab = bquote(~-Log[10]~"adjusted"~italic(p)),
                        title = NULL,
                        subtitle = "",
                        pCutoff = 1e-4,
                        pointSize = 2.0,
                        labSize = 9,
                        labCol = 'black',
                        labFace = 'plain',
                        drawConnectors = T,
                        arrowheads = F,
                        boxedLabels = F,
                        colCustom = keyvals,
                        FCcutoff = (1),
                        legendPosition = "none") 




ggsave("figures/06C_4T1_volcano.pdf",width = 8,height = 8)






sc_back <- read.csv("output/backgroundHER2.csv")
sc_back <- sc_back$x
DE_4T1  <- DE_4T1[DE_4T1$X %in% sc_back, ]
all4T1 <- DE_4T1$X
up4T1  <- DE_4T1[DE_4T1$Direction=="UP",]$X

enrich_de <- list()
enriches <- list()



saveTag <- paste0("output/new_tumors/tumor_late_vs_early.csv")
tryCatch({
        cluster0.markers <- read.csv(saveTag)
        negs <- cluster0.markers[cluster0.markers$avg_log2FC<0,]$Gene
        negs <- negs[negs %in% all4T1]
        de_both <- intersect(negs,up4T1)
        p_val   <- phyper(q = length(de_both)-1,
                          m= length(up4T1),
                          n = length(all4T1) - length(up4T1),
                          k = length(negs), lower.tail = FALSE)
        a <- venn.diagram(
                x = list(up4T1,negs),
                category.names = c(paste0("Downregulated in late tumors"),
                                   "Upregulated in DAC treated 4T1"),
                filename = NULL,
                output=F,
                inverted = T,
                # Numbers,
                col  = rev(c("red3","royalblue")),#myCol[1:3],
                fill = rev(c("red3","royalblue")),#myCol[1:3],
                cex = 6,
                fontface = "plain",
                fontfamily = "sans",
                cat.cex = 2,
                # Set names
                cat.fontface = "bold",
                cat.default.pos = "outer",
                cat.fontfamily = "sans",
                main.cex = 4.0,
                main.fontfamily = "sans",
                cat.dist = c(0.035, 0.035) ,
                cat.pos = c(-9.5, 9.5),
                scaled = F, rotation.degree = 180
        )
        
        b <- venn.diagram(
                x = list(up4T1,negs),
                category.names = c(paste0(""),
                                   ""),
                filename = NULL,
                output=F,
                inverted = T,
                # Numbers,
                col  = (c("red3","royalblue")),#myCol[1:3],
                fill = (c("red3","royalblue")),#myCol[1:3],
                cex = 6,
                fontface = "plain",
                fontfamily = "sans",
                cat.cex = 2,
                # Set names
                cat.fontface = "bold",
                cat.default.pos = "outer",
                cat.fontfamily = "sans",
                main.cex = 4.0,
                main.fontfamily = "sans",
                cat.dist = c(0.035, 0.035) ,
                cat.pos = c(-9.5, 9.5),
                scaled = T )
        ggsave(paste0("figures/06D_scRNA_tumor_4T1_venn.PDF"),
               plot =a, width = 12, height = 12)
        ggsave(paste0("figures/06D_scRNA_tumor_4T1_venn_proportiona.PDF"),
               plot =b, width = 12, height = 12)
        write.csv(matrix(t(sort(de_both)), ncol = 7),
                  "output/new_tumors/All_tumor_down_DAC_up_4T1.csv")
        write.csv(matrix(t(sort(de_both))),
                  "output/new_tumors/All_tumor_down_DAC_up_4T1_long.csv")
        
}
,
error = function(e){print("passsing")})


background_genes <- read.csv("output/backgroundHER2.csv")
background_genes <- background_genes$x

ggo.neg  <- enrichGO(gene     = de_both,
                     universe = background_genes,
                     OrgDb    = org.Mm.eg.db,
                     keyType  = "SYMBOL", 
                     ont      = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.1,
                     qvalueCutoff  = 0.05)
temp <- as.data.frame(ggo.neg)

write.csv(temp, "output/new_tumors/12_shared_DAC_scGEMM.csv")



### Redo this part


my_cols <- 
        c("#a6cee3", "#1f78b4",
          "#b2df8a", "#33a02c",
          "#fb9a99", "#e31a1c", "#fdbf6f",
          "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99",
          "#b15928", "#8dd3c7",
          "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
          "#fdb462", "#b3de69", "#fccde5", "#d9d9d9",
          "#bc80bd","#ccebc5","#ffed6f","#ffed6f" )


edox2 <- pairwise_termsim(ggo.neg)
p2 <- treeplot(edox2,showCategory = 40, nCluster = 7, group_color = my_cols[1:7])

enrich4T1 <-  p2$data

unique_4T1 <- enrich4T1


unique_4T1 <- unique_4T1[,c("label", "color", "group")]
names(unique_4T1) <- c("name", "score","group")
unique_4T1$score <- -log10(unique_4T1$score)



unique_4T1 <- unique_4T1[1:40,]


unique_4T1 <- unique_4T1[order(unique_4T1$score, decreasing = F),]
unique_4T1 <- unique_4T1[order(unique_4T1$group, decreasing = T),]
unique_4T1$name <- str_sub(unique_4T1$name,1,80)
while(any(duplicated(unique_4T1$name)))
        unique_4T1[duplicated(unique_4T1$name),]$name <- 
        paste0(unique_4T1[duplicated(unique_4T1$name),]$name,"*")


p21 <- ggplot(unique_4T1, aes( fct_inorder(name), score, color = group)) +
        geom_point(stat = "identity", size = 7) +
        scale_color_manual("",labels = c( "Interferon beta",
                                          "viral response",
                                          "Interferon gamma",
                                          "Innate immune response",
                                          "Viral genome replication",
                                          "Antigen processing, cell killing",
                                          "Type I interferon"),
                           values= my_cols[15:21]) +
        labs(title= "",y = "-log(enrichment p-value)", x= "")+
        coord_flip() + theme_bw() +
        theme(axis.title.x =  element_text(size = 16, hjust = 0.5),
              axis.text = element_text(size = 24, color = "black"),
              title = element_text(size = 16, hjust = 1),
              legend.direction="vertical",
              legend.position="none",
              legend.text = element_text(size = 16))


#ggsave("figures/06E_shared_dac_scGemms.pdf", width = 18, height = 15)



cluster0.markers <- read.csv(saveTag)
negs <- cluster0.markers[cluster0.markers$avg_log2FC<0,]$Gene
DE_4T1_only <- setdiff(up4T1,negs)


ggo.diff  <- enrichGO(gene     = DE_4T1_only,
                      universe = all4T1,
                      OrgDb    = org.Mm.eg.db,
                      keyType  = "SYMBOL", 
                      ont      = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.05)
temp.diff <- as.data.frame(ggo.diff)

write.csv(temp.diff, "output/new_tumors/13_DAC_no_scGEMM.csv")



my_cols <- 
        c("#a6cee3", "#1f78b4",
          "#b2df8a", "#33a02c",
          "#fb9a99", "#e31a1c", "#fdbf6f",
          "#ff7f00", "#cab2d6",
          "#6a3d9a", "#ffff99",
          "#b15928", "#8dd3c7",
          "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
          "#fdb462", "#b3de69", "#fccde5", "#d9d9d9",
          "#bc80bd","#ccebc5","#ffed6f","#ffed6f" )


edox2 <- pairwise_termsim(ggo.diff)
p2 <- treeplot(edox2,showCategory = 40, nCluster = 7,
               group_color = my_cols[1:7])

enrich4T1 <-  p2$data

unique_4T1 <- enrich4T1


unique_4T1 <- unique_4T1[,c("label", "color", "group")]
names(unique_4T1) <- c("name", "score","group")
unique_4T1$score <- -log10(unique_4T1$score)



unique_4T1 <- unique_4T1[1:40,]


unique_4T1 <- unique_4T1[order(unique_4T1$score, decreasing = F),]
unique_4T1 <- unique_4T1[order(unique_4T1$group, decreasing = T),]
unique_4T1$name <- str_sub(unique_4T1$name,1,80)
while(any(duplicated(unique_4T1$name)))
        unique_4T1[duplicated(unique_4T1$name),]$name <-
        paste0(unique_4T1[duplicated(unique_4T1$name),]$name,"*")


p3 <- ggplot(unique_4T1, aes( fct_inorder(name), score, color = group)) +
        geom_point(stat = "identity", size = 7) +
        scale_color_manual("",labels = c( "Immune effector process",
                                          "Inflammation",
                                          "Immune cell activation",
                                          "Response to bacterium",
                                          "Cell killing",
                                          "Adaptive immune response",
                                          "Cell activation in immune response"),
                           values= my_cols[1:7]) +
        labs(title= "",y = "-log(enrichment p-value)", x= "")+
        coord_flip() + theme_bw() + 
        theme(axis.title.x =  element_text(size = 16, hjust = 0.5),
              axis.text = element_text(size = 24, colour = "black"),
              title = element_text(size = 16, hjust = 1),
              legend.direction="vertical",
              legend.position="none",
              legend.text = element_text(size = 16))

ggsave("figures/Ex08F_dac_no_scGemms.pdf", width = 19, height = 15)


library(patchwork)
ggsave("figures/Ex08F_Full_dac_both.pdf", p21 + p3, width = 32, height = 15)


