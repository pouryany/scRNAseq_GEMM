rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(DirichletReg)
library(tibble)

scGemms_meta <- readRDS("data/03_umap_labled_1_metadata.RDS")

scGemms_meta_epi <- scGemms_meta[grep(scGemms_meta$cell_type_secondary,
                                      pattern = "Alv|Basal"),]


scGemms_meta_epi$cell_type_secondary <-
    droplevels(scGemms_meta_epi$cell_type_secondary)

temp <- table(scGemms_meta_epi$sampleID,
              scGemms_meta_epi$cell_type_secondary)


temp <- as.data.frame.matrix(temp)
temp <- temp[2:8,]
temp$stage <- rowSums(temp)
temp <- temp + (temp$stage/(sum(temp$stage)))
temp <- temp/temp$stage
temp$stage <- str_sub(rownames(temp),1,1)


#temp[,1:10] <- temp[,1:10] + 0.0002
AL <- DR_data(temp[,1:8])


fit <- DirichReg(AL ~ stage  , temp)
fit2 <- DirichReg(AL ~ stage   , temp)
anova(fit,fit2)
u = summary(fit)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit$pvals = pvals

(u$coef.mat)
fitted(fit)
confint(fit)
plot_me       <- fit$fitted.values$mu
plot_me       <- as.data.frame(plot_me)
plot_me$stage <- temp$stage

plot_me <- plot_me[!duplicated(plot_me$stage),]
rownames(plot_me) <- c("Early Tumor","Late Tumor")
plot_me$stage <- NULL

plot_me <- t(plot_me)
plot_me = rownames_to_column(as.data.frame(plot_me),var = "Cell type")
plot_me$pval <- (as.data.frame.table(pvals)$Freq)
colnames(plot_me) <- c("Cell type","Early Tumor","Late Tumor", "pvals")

plot_me$p.adj<- p.adjust(plot_me$pvals,method = "fdr")

write.csv(plot_me,"output/20_dirichletTumor.csv",row.names = F)

