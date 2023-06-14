rm(list = ls())


B16Expression <- read.csv("data/RNA_Seq_B16_raw.csv")

B16Expression <- B16Expression[!duplicated(B16Expression$Gene.name),]


row.names(B16Expression) <- B16Expression$Gene.name


counts <- B16Expression[,c("Control_1", "Control_2", "Control_3",
                            "DAC_1", "DAC_2", "DAC_3")]
counts <- as.matrix((round(counts)))

counts <- apply(counts,2, as.integer)

row.names(counts) <- row.names(B16Expression)
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

sigRes <- res[abs(res$log2FoldChange) > 1  & res$padj < 1E-02,]
sigResPos <- res[(res$log2FoldChange) > 1  & res$padj < 1E-02,]
sigResNeg <- res[(res$log2FoldChange) < -1 & res$padj < 1E-02,]

write_res <- as.data.frame(res)
write_res$Direction <- "NONE"

library(dplyr)
write_res[write_res$log2FoldChange > 1  & res$padj < 1E-02,]$Direction <-  "UP" 
write_res[write_res$log2FoldChange  < -1 & res$padj < 1E-02,]$Direction <-  "DOWN" 

write.csv(write_res,"results_bulk/DEGenes_B16.csv")

