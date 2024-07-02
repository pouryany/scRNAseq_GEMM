rm(list = ls())

library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(forcats)
library(viridis)
library(rlang)
library(forcats)
library(DirichletReg)
library(tibble)
library(dplyr)
library(tidyr)


scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
scGemms$initial_clusts <- scGemms$seurat_clusters
colorss <- readRDS("data/03_colors.RDS")




t_cells <- grep("T cell|NK", unique(scGemms$cell_type_detail),value = T)
scGemms0 <- subset(scGemms, subset = cell_type_detail %in% t_cells)


cell_tags <- unique(as.character(scGemms0$seurat_clusters))


# Filter mitochondrial and ribosomal genes
rm.ind <- grep("^mt-|^Rp[ls]",rownames(scGemms0))
keep   <- rownames(scGemms0)[-rm.ind]
scGemms0 <- scGemms0[keep,]

# Filter mitochondrial and ribosomal genes
keep.ind <- rowSums(scGemms0@assays$RNA@counts)
keep.ind <- names(keep.ind[keep.ind>=20])
#keep   <- rownames(scGemms0)[-rm.ind]
scGemms0 <- scGemms0[keep.ind,]


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]


Idents(scGemms0) <- "cell_type_secondary"
#scGemms0$color <- droplevels(scGemms0$color)




order_vec <-  c(
        "CD4+ Naïve", "CD4+ Treg",             
        "CD8+ Naïve", "CD8+ CM",
        "CD8+ ISG","CD8+ Effector",
        "CD8+ TRM",
        "CD8+ Progenitor Ex",
        "CD8+ Ex", 
        "CD8+ Proliferating",
        "NK")    



 
p00 <- DimPlot(scGemms0,group.by = "cell_type_secondary",
               cols = colors0[order_vec],
               combine = T, label = F, pt.size = 0.6) &  NoAxes() &
        labs(title = "T and NK cells") &  
        theme(plot.title = element_text(hjust = 0, size = 30, face = "plain"), 
              legend.text = element_text(size = 28),
              axis.title.y.right = element_text(size = 20)) &
        guides(color = guide_legend(override.aes = list(size=12), ncol=1) ) &
        ylim(-4,9) & xlim(-6,10)


ggsave("figures/03A_Tcell_umap.pdf",p00, height = 9  , width = 15 )


 
my_features2 <- c( "Cd4", "Foxp3", "Cd8a",
                   "Cd8b1", "Ncr1","Cd44",
                   "Ccr7", "Sell",  "Il7r",
                   "Isg15","Pdcd1","Mki67",
                   "Fcer1g","Gzmb")



p1 <- FeaturePlot(scGemms0,features = my_features2
                  ,ncol = 3,
                  combine = T,cols = c("#d6d6d6","#6e016b"))  &
        NoAxes()  & NoLegend() &
        theme(plot.title = element_text(hjust = 0.5, size = 20,
                                        face = "plain"), 
              legend.text = element_text(size = 24)) &
        ylim(-4,9) & xlim(-6,10)

ggsave("figures/Ex03B_markers.pdf",p1,width = 6, height = 8)





stack_plotter_sc <- function(data_meta, var1 = "cell_type_secondary",
                             var2 = "stage", var2_freq = "Early",
                             colors0 = colors0){
        data_tcell <- data_meta
        df2 <- data_tcell %>% 
                group_by(!!sym(var1), !!sym(var2)) %>% 
                summarise(n = n()) %>% 
                mutate(freq = n / sum(n)) %>%
                ungroup() %>%
                dplyr::filter(!!sym(var2) == var2_freq) %>%
                arrange(.,freq)
        
        
        
        data_tcell[[var1]] <- 
                factor(data_tcell[[var1]], levels = as.character(df2[[var1]] ))
        

        stack_plot <-  ggplot(data_tcell, aes(x = !!sym(var1) ,
                                              fill = !!sym(var2))) +
                geom_bar(position = "fill") + 
                theme_classic() + labs(x = "", y = "") + 
                scale_y_continuous(labels=scales::percent) +
                scale_fill_manual(values = colors0)+
                theme(panel.border = element_rect(colour = "black",
                                                  fill=NA, size=1),
                      axis.text.x = element_text(angle = 0,
                                                 hjust = 0.5,
                                                 size = 14),
                      axis.text.y = element_text(angle = 0,
                                                 hjust = 1,
                                                 size = 20, face = "bold"),
                      axis.line = element_blank(),
                      axis.ticks = element_blank()) +
                coord_flip()
        
        return(stack_plot)
}

data_meta  <- scGemms0@meta.data
data_tcell <- data_meta

df2 <- data_tcell %>% 
        group_by(cell_type_secondary, stage) %>% 
        summarise(n = n()) %>% 
        mutate(freq = n / sum(n)) %>%
        ungroup() 




data_tcell$cell_type_secondary <-
        factor(data_tcell$cell_type_secondary,levels = order_vec)


stack_plot <-  ggplot(data_tcell, aes(x = stage ,
                                      fill = cell_type_secondary)) +
        geom_bar(position = "fill") + 
        theme_classic() + labs(x = "", y = "") + 
        scale_fill_manual(values = colors0)+
        theme(axis.text.x = element_text(angle = 0,
                                         hjust = 0.5,
                                         size = 16),
              axis.text.y = element_text(angle = 0,
                                         hjust = 1,
                                         size = 30, face = "plain", 
                                         colour = "black"),
              axis.line = element_blank(),axis.ticks = element_blank()) +
        coord_flip()

ggsave("figures/03B_tcell_proportion_stage.pdf",stack_plot,
       height = 2.5, width = 10)







temp <- table(scGemms0$sampleID, scGemms0$cell_type_secondary)

temp <- as.data.frame.matrix(temp)
temp <- temp[1:7,]
temp$stage <- rowSums(temp)
temp <- temp + (temp$stage/(sum(temp$stage)))
temp <- temp/temp$stage
temp$stage <- str_sub(rownames(temp),1,1)


#temp[,1:10] <- temp[,1:10] + 0.0002
AL <- DR_data(temp[,1:11])


fit <- DirichReg(AL ~ stage  , temp)
u = summary(fit)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
colnames(pvals) = u$varnames
fit$pvals = pvals

fit$p.adj <- p.adjust(pvals,method = "fdr")


plot_me       <- fit$fitted.values$mu
plot_me       <- as.data.frame(plot_me)
plot_me$stage <- temp$stage

plot_me <- plot_me[!duplicated(plot_me$stage),]
rownames(plot_me) <- c("Early Tumor","Late Tumor")
plot_me$stage <- NULL

plot_me <- t(plot_me)
plot_me = rownames_to_column(as.data.frame(plot_me),var = "Cell type")
plot_me$pval <- (as.data.frame.table(pvals)$Freq)
#pvals is adjusted

plot_me0 <- plot_me
plot_me0$p.adj <- p.adjust(plot_me0$pval,method = "fdr")

write.csv(plot_me0,"output/20_Dirichlet_T_NK.csv")





scGemms0$cell_type_secondary <- factor(scGemms0$cell_type_secondary,
                                       levels = order_vec)
Idents(scGemms0) <- "cell_type_secondary"



vln_eGFP <- VlnPlot(scGemms0,features = c("Sell", "Ccr7","Lef1","Tcf7", 
                                          "Bach2","Mki67","S1pr1",
                                          "Fcer1g", "Itgae",
                                          "Isg15", "Ifit1", "Irf7", 
                                          "Xcl1", "Klrk1", "Klrc1", "Klra5",
                                          "Pdcd1","Tigit","Lag3","Cd244a",
                                          "Ctla4"),
                cols =  colorss[(names(colorss))],
                same.y.lims = F,ncol = 3,pt.size = 0)&
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(size = 30, face = "plain", hjust = 0))

ggsave("figures/Ex03C_Tcell_markers_vln_0.pdf",vln_eGFP,width = 14,height = 12)




vln_eGFP <- VlnPlot(scGemms0,features = c("Sell", "Ccr7","Lef1","Tcf7", 
                                          "Bach2","Mki67","S1pr1",
                                          "Fcer1g", "Itgae",
                                          "Isg15", "Ifit1", "Irf7", 
                                          "Xcl1", "Klrk1", "Klrc1", "Klra5",
                                          "Pdcd1","Tigit","Lag3","Cd244a",
                                          "Ctla4"),
                    cols =  colorss[(names(colorss))],
                    same.y.lims = F,ncol = 3,pt.size = 0)&
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5),
              axis.title.y = element_blank(),
              plot.title = element_text(size = 30, face = "plain", hjust = 0))

ggsave("figures/Ex03C_Tcell_markers_vln_legend.pdf",vln_eGFP,width = 14,height = 16)



