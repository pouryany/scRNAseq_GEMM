rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(EnhancedVolcano)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(forcats)



scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
colorss <- readRDS("data/03_colors.RDS")


scGemms0 <- subset(scGemms, subset = cell_type_primary == "Epithelial")
scGemms0 <- subset(scGemms, cells = names(grep(pattern = "HER2",
                                   scGemms$cell_type_detail,value = T)))
scGemms0 <- subset(scGemms0, 
                   cells = names(grep(pattern = "Basal|Proliferating",
                                      scGemms0$cell_type_detail,value = T)))


# Filter mitochondrial and ribosomal genes
rm.ind <- grep("^mt-|^Rp[ls]",rownames(scGemms0))
keep   <- rownames(scGemms0)[-rm.ind]
scGemms0 <- scGemms0[keep,]

# Filter mitochondrial and ribosomal genes
keep.ind <- rowSums(scGemms0@assays$RNA@counts)
keep.ind <- names(keep.ind[keep.ind>=20])
#keep   <- rownames(scGemms0)[-rm.ind]

scGemms0 <- scGemms0[keep.ind,]


if(!dir.exists("output/new_tumors"))
        dir.create("output/new_tumors",recursive = T)

write.csv(keep.ind,"output/new_tumors/tumor_background.csv")


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]





cell_tags <- unique(as.character(scGemms0$cell_type_secondary))

Idents(scGemms0) <- "cell_type_secondary"



cluster0.markers0 <- FindMarkers(scGemms0, ident.1 = "Late",
                                group.by = "stage",
                                only.pos = F,
                                test.use = "MAST",
                                latent.vars = c("nCount_RNA"))

cluster0.markers <- 
        cluster0.markers0[cluster0.markers0$p_val_adj<0.05 & 
                                abs(cluster0.markers0$avg_log2FC) > 0.25,]

cluster0.markers <- tibble::rownames_to_column(cluster0.markers,var = "Gene")
saveTag <- paste0("output/new_tumors/tumor_late_vs_early.csv")
write.csv(cluster0.markers,saveTag)




keyvals <- ifelse(
        cluster0.markers$avg_log2FC < 0, 'royalblue',
        ifelse(cluster0.markers$avg_log2FC > 0, 'red3',
               'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red3'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'


sel_labs <- (cluster0.markers$Gene)[which(names(keyvals) %in% c('low'))]
sel_labs <- grep("Ifi|Isg|Irf|Cxc|Zbp|Gbp|Stat|Irgm|Irg|Ly6|Ddx|Pd",
                 sel_labs,value = T)

p1 <-   EnhancedVolcano(cluster0.markers,
                        lab = cluster0.markers$Gene, #res1$hgnc_symbol,
                        selectLab = sel_labs,
                        x = "avg_log2FC",
                        y = "p_val_adj",
                        xlim = c(-2,2),
                        ylim = c(0, 260),
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


 ggsave("figures/02A_all_tumor_volcano.pdf",
        p1,width = 9,height = 12)
# ggsave("results_v3/Figure12_4T1_volcano.png",width = 8,height = 8)

 
 
 

 
 saveTag <- paste0("output/new_tumors/tumor_late_vs_early.csv")
 
 tryCatch(cluster0.markers <- read.csv(saveTag),
          error = function(e){print("passsing")})
 
 background_genes <- read.csv("output/new_tumors/tumor_background.csv")
 background_genes <- background_genes$x
 
 negs <- cluster0.markers[cluster0.markers$avg_log2FC<0,]$Gene
 
 
 
 
 ggo.neg  <- enrichGO(gene     = negs,
                      universe = background_genes,
                      OrgDb    = org.Mm.eg.db,
                      keyType  = "SYMBOL", 
                      ont      = "BP",
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.1,
                      qvalueCutoff  = 0.05)
 
 saveTag <-
        paste0("output/new_tumors/enrichment_neg_tumor_late_vs_early.csv")
 
 write.csv(as.data.frame(ggo.neg),saveTag)
 
 poses  <- cluster0.markers[cluster0.markers$avg_log2FC>0,]$Gene
 ggo.pos    <- enrichGO(gene     = poses,
                        universe = background_genes,
                        OrgDb    = org.Mm.eg.db,
                        keyType  = "SYMBOL", 
                        ont      = "BP",
                        pAdjustMethod = "fdr",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.05)
 
 saveTag <-
         paste0("output/new_tumors/enrichment_positive_tumor_late_vs_early.csv")
 
 write.csv(as.data.frame(ggo.pos),saveTag)
 
 
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
 nclust <- 7
 p2 <- treeplot(edox2,showCategory = 50,
                nCluster = nclust,
                group_color = my_cols[1:(nclust)])
 
 enrich4T1 <-  p2$data
 unique_4T1 <- enrich4T1
 
 
 unique_4T1 <- unique_4T1[,c("label", "color", "group")]
 names(unique_4T1) <- c("name", "score","group")
 unique_4T1$score <- -log10(unique_4T1$score)
 
 
 unique_4T1 <- unique_4T1[1:50,]
 # Manually grouped top enrichment terms
 unique_4T1 <- read.csv("data/Figure2B_table.csv")
 unique_4T1 <- unique_4T1[,c(1,2,4)]
 names(unique_4T1) <- c("name", "score","group")
 
 
 unique_4T1 <- unique_4T1[order(unique_4T1$score, decreasing = F),]
 unique_4T1 <- unique_4T1[order(unique_4T1$group, decreasing = T),]
 
 unique_4T1$name <- str_sub(unique_4T1$name,1,70)
 while(any(duplicated(unique_4T1$name)))
         unique_4T1[duplicated(unique_4T1$name),]$name <-
                paste0(unique_4T1[duplicated(unique_4T1$name),]$name,"*")
 
 p2 <- ggplot(unique_4T1, aes( fct_inorder(name), score, color = group)) +
         geom_point(stat = "identity", size = 6) + ylim(0,25) +
         scale_color_manual("",
                            labels = c( "Type I Interferon Response",
                                   "Defense Response to Virus",
                                   "Type II Interferon Response",
                                   "Cytokine-mediated Signaling",
                                   "Antigen presentation and processing",
                                   "T cell immunity and cytotoxicity",
                                   "Oxidative Phosphorylation"),
                            values= my_cols[15:22]) +
         labs(title= "",y = "-log(enrichment p-value)", x= "")+
         coord_flip() + 
         theme_bw() +
         theme(axis.title.x =  element_text(size = 20, hjust = 0.5),
                 axis.text = element_text(size = 20, colour = "black"),
                 title = element_text(size = 16, hjust = 1),
                 legend.text = element_text(size = 16),
                 legend.direction="vertical",
                 legend.position="right")

 ggsave("figures/02B_enrichment_negs.pdf",
        width = 16, height = 15)
 
 
 
 
 
 
 
library(msigdbr)



h_gene_sets = msigdbr(species = "mouse", category = "H")
head(h_gene_sets)
msigdbr_t2g = h_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()



negs <- cluster0.markers[cluster0.markers$avg_log2FC<0,]$Gene



temp0 <- enricher(gene = negs,
                 universe = background_genes,
                 TERM2GENE = msigdbr_t2g,qvalueCutoff = 0.1)

temp0 <- as.data.frame(temp0)
temp0$Direction <- "Down in late"



temp <- enricher(gene = cluster0.markers[cluster0.markers$avg_log2FC>0,]$Gene,
                 universe = background_genes,
                 TERM2GENE = msigdbr_t2g,pvalueCutoff = 0.1)

temp <- as.data.frame(temp)
temp$Direction <- "Up in late"




temp0 <- rbind(temp0,temp)

saveTag <-
        paste0("output/new_tumors/Hallmarks_enrichment_late_vs_early.csv")

write.csv(temp0,saveTag)

