rm(list = ls())


T1Expression <- read.csv("Data/RNA_Seq_T41_raw.csv")
T1Expression <- T1Expression[!duplicated(T1Expression$Gene.name),]


row.names(T1Expression) <- T1Expression$Gene.name


counts <- T1Expression[,c("Control_1", "Control_2", "Control_3",
                            "DAC_1", "DAC_2", "DAC_3")]
counts <- as.matrix((round(counts)))

counts <- apply(counts,2, as.integer)

row.names(counts) <- row.names(T1Expression)
colnames(counts) <- c("Control_1", "Control_2", "Control_3",
                      "DAC_1", "DAC_2", "DAC_3")


COVARIATES <- data.frame("Condition" = c("Control", "Control", "Control",
                                        "DAC", "DAC", "DAC"))

row.names(COVARIATES) <- colnames(counts)

COVARIATES$Condition <- factor(COVARIATES$Condition)



library(DESeq2)


dds <- DESeqDataSetFromMatrix(countData = (counts),
                              colData = DataFrame(COVARIATES),
                              design = ~ Condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
res$padj <- ifelse(is.na(res$padj), 1, res$padj)

#write.csv(res,"DESeq_4T1.csv")

sigRes <- res[abs(res$log2FoldChange) > 1  & res$padj < 1E-4,]
sigResPos <- res[(res$log2FoldChange) > 1  & res$padj < 1E-4,]
sigResNeg <- res[(res$log2FoldChange) < -1  & res$padj < 1E-4,]

write_res <- as.data.frame(res)
write_res$Direction <- "NONE"

library(dplyr)

write_res[write_res$log2FoldChange  > 1  & res$padj < 1E-4,]$Direction <-  "UP" 
write_res[write_res$log2FoldChange  < -1  & res$padj < 1E-4,]$Direction <-  "DOWN" 

if(!dir.exists("results_bulk"))
        dir.create("results_bulk")
write.csv(write_res,"results_bulk/DEGenes_4T1.csv")


erv_genes <- c("Syna", "Synb","Fv1","Peg10" )
erv_vals <- write_res[erv_genes,]

write.csv(erv_vals,"results_bulk/ervs_4T1.csv")


expression_values <- counts(dds, normalized=TRUE)
erv_exp <- expression_values[erv_genes,]
erv_exp <-data.frame((erv_exp))
erv_exp$gene <- rownames(erv_exp)

erv_exp <- reshape2::melt(erv_exp,"gene")
erv_exp$group <- gsub("_.","",erv_exp$variable)

ggplot(erv_exp, aes(group,value)) +
        geom_boxplot(width = 0.3) +
        theme_bw() +
        facet_grid(cols = vars(gene)) + 
        labs(y = "Expression", x = "")+
        theme(strip.text.x = element_text(size = 30),
              axis.title.y = element_text(size = 30),
              axis.text.y = element_text(size = 20),
              axis.text.x = element_text(size = 20,angle = 45, vjust = 0.5))

ggsave("results_bulk/S08_ERV.pdf",width = 8, height = 8)



library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


ggoPos <- enrichGO(gene     = rownames(sigResPos),
                universe = rownames(res),
                OrgDb    = org.Mm.eg.db,
                keyType  = "SYMBOL", 
                ont      = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.05)

write.csv(as.data.frame(ggoPos),"results_bulk/4T1_enrichment_upregulated.csv")




temp <- as.data.frame(ggoPos)
rownames(temp) <- NULL
temp$Phenotype <- "UP"


temp <- temp[,c("ID","Description", "pvalue", "p.adjust", "Phenotype","geneID")]
colnames(temp) <- c("Go.ID", "Description", "p.Val", "FDR","Phenotype" , "Genes")

temp$Genes <- gsub("\\/",",",temp$Genes)

write.table(temp[1:150,],"results_bulk/4T1_enrichment_upregulated.txt", 
            sep = "\t" , quote = F,  row.names = F)

write.table(temp[1:150,1:5],
            "results_bulk/Trimmed_4T1_enrichment_upregulated.txt", 
            sep = "\t" , quote = F,  row.names = F)

