rm(list = ls())

library(Seurat)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(forcats)
library(forcats)
library(ggthemes)
library(dplyr)
library(forcats)


scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
colorss <- readRDS("data/03_colors.RDS")



order_vec <-  c(
        "Alv-Proliferating","Alv-Basal",
        "Basal",
        "Alv 1","Alv 2","Alv 3",
        "Alv 4", 
        "Undefined Alv",
        "Hormone Sensing",
        "CD4+ Naïve",
        "CD4+ Treg",             
        "CD8+ Naïve", 
        "CD8+ CM",
        "CD8+ ISG",
        "CD8+ Effector",
        "CD8+ TRM",
        "CD8+ Progenitor Ex",
        "CD8+ Ex", 
        "CD8+ Proliferating",
        "NK",
        "Macrophage",
        "Neutrophil",
        "B cell",        
        "Fibroblast")    


my_features <- c( "Epcam","Ptprc","eGFP","Dcn",
                  "Csn3", "Krt14","Erbb2", "Elf5",
                  "Aldh1a3","Mki67","Prlr",
                  "Cd8b1","Cd3d","Cd4",
                  "Ccr7","Nkg7", "Cd79a",
                  "Ncr1", "Cd14", "Pdcd1",
                  "Foxp3","S100a9", "Cd68","Fcgr3",
                  "Lyz2","Cd74", "Sell", "Cd160",
                  "Il7r", "Cd44", "Emid1", "Auts2","Il1b")



#> Individual panels corresponding to broad cellular and specific 
#> epithelial markers for Figures S1 and S2C
for(item in my_features){
        
        p1 <- FeaturePlot(scGemms,features = item
                          ,ncol = 1,
                          combine = T,cols = c("#d6d6d6","#6e016b"))  &
                NoAxes()  & NoLegend() &
                theme(plot.title = element_text(hjust = 0, size = 20,
                                                face = "plain"), 
                      legend.text = element_text(size = 24))
        
        if(!dir.exists("figures/S1_S2C/"))
                dir.create("figures/S1_S2C/",recursive = T)
        
        ggsave(paste0("figures/S1_S2C/",item,".png"),
               p1,width = 4, height = 4)
        
        
}







scGemms0  <- subset(scGemms, subset = cell_type_primary == "Epithelial")
cell_tags <- unique(as.character(scGemms0$cell_type_secondary))

Idents(scGemms0) <- "cell_type_secondary"

order_vec <-  c(
        "Alv-Proliferating","Alv-Basal", 
        "Basal",
        "Alv 1", "Alv 2","Alv 3", 
        "Alv 4",  "Undefined Alv","Hormone Sensing")    

cell_tags <- unique(as.character(scGemms0$cell_type_secondary))

Idents(scGemms0) <- "cell_type_secondary"


scGemms0@active.ident <- factor(x = scGemms0@active.ident, levels = order_vec)




p1 <- FeaturePlot(scGemms,features = c( "Epcam","eGFP"),ncol = 2, pt.size = 0.1,
                  combine = T, cols = c("#d6d6d6","#6e016b"), min.cutoff = "q1")  &
        NoAxes() & theme(legend.position="bottom") & 
        guides(colourbar = guide_legend(override.aes = list(size=12, nrow = 1))) &
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

ggsave("figures/01D_1E_markers.png",p1,width = 7, height = 5)







p1 <- DotPlot(scGemms0, 
        features = c( "Prlr","Krt19","Esr1","Pgr",
                      "Egfr","Auts2","Chrm3","Pde4d","Mdga2",
                      "Il31ra","Tnc","Il1b",
                      "Emid1","Cst3","Iigp1",
                      "Csn2","Btn1a1","Elf5",
                      "Krt14","Ecrg4","Serpine2","Col18a1","Vim",
                      "Aldh1a3","Cldn4","Hspb1","Ly6a","Anxa1","Klk10",
                      "Adamts1",
                      "Top2a","Mki67","Birc5",
                      "Erbb2","Epcam","eGFP"),
        cols = c(low  = "#f7f7f7",
                 high = "#6e016b")) & coord_flip()&
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
              axis.text.y = element_text(hjust = 1, size = 20),
              axis.title = element_blank())

ggsave("figures/S2A_epithelial_dotplot.pdf",plot = p1,
       width = 7, height = 12)




scGemms0 <- subset(scGemms, subset = cell_type_primary == "Epithelial")
scGemms0 <- subset(scGemms, cells = names(grep(pattern = "HER2",
                                   scGemms$cell_type_detail,value = T)))


scGemms0$cell_type_secondary <- droplevels(scGemms0$cell_type_secondary)
colors0 <- colorss[names(colorss) %in% levels(scGemms0$cell_type_secondary)]



order_vec <-  c(
        "Alv-Proliferating",
        "Alv-Basal",
        "Basal",
        "Alv 1",
        "Alv 2",
        "Alv 3", 
        "Alv 4")    

cell_tags <- unique(as.character(scGemms0$cell_type_secondary))
Idents(scGemms0) <- "cell_type_secondary"
scGemms0@active.ident <- factor(x = scGemms0@active.ident, levels = order_vec)




p00 <- DimPlot(scGemms0,group.by = "cell_type_secondary",
               cols = colors0[order_vec],
               combine = T, label = F, pt.size = 0.6) &  NoAxes() & NoLegend()& 
        labs(title = NULL) & theme(title = element_text(hjust = 1)) &
        guides(color = guide_legend(override.aes = list(size=12), ncol=1) )

ggsave("figures/01F_tumor_umap.pdf", height = 5, width = 10)




data_meta  <- scGemms0@meta.data
data_tcell <- data_meta[data_meta$cell_type_primary == "Epithelial",]
data_tcell$cell_type_secondary <- droplevels(data_tcell$cell_type_secondary)

df2 <- data_tcell %>% 
        group_by(cell_type_secondary, stage) %>% 
        summarise(n = n()) %>% 
        mutate(freq = n / sum(n)) %>%
        ungroup() %>% dplyr::filter(stage == "Early") %>% arrange(.,freq)



data_tcell$cell_type_secondary <- 
        factor(data_tcell$cell_type_secondary,
               levels = as.character(df2$cell_type_secondary))


base_line <- as.data.frame(table(scGemms0$stage))
base_line <- base_line$Freq[1]/sum(base_line$Freq)






df2 <- data_tcell %>% 
        group_by(cell_type_secondary, stage) %>% 
        summarise(n = n()) %>% 
        mutate(freq = n / sum(n)) %>%
        ungroup() 

data_tcell$cell_type_secondary <- factor(data_tcell$cell_type_secondary,
                                         levels = order_vec)

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
                                         size = 24, face = "plain",
                                         colour = "black"),
              axis.line = element_blank(),axis.ticks = element_blank(),
              legend.text=element_text(size=24)) + coord_flip()

ggsave("figures/S2D_tumor_proportion_Stage.pdf",stack_plot,
       height = 2.5, width = 10)




temp_cells2 <- WhichCells(scGemms, expression = eGFP > 0 & Epcam > 0)
temp_cells3 <- WhichCells(scGemms, expression = Erbb2 > 0 &
                                  Epcam > 0 &
                                  eGFP > 0)

scGemms[["true_tumor"]] <- plyr::mapvalues((colnames(scGemms) %in% temp_cells2),
                                           c(FALSE,TRUE),c("Non-tumor",
                                                           "Tumor"))

scGemms[["erbb_score"]] <- plyr::mapvalues((colnames(scGemms) %in% temp_cells3),
                                           c(FALSE,TRUE),c("HER2 NEG",
                                                           "HER2 POS"))




data_tcell <- as.data.frame.matrix(table(scGemms$cell_type_secondary,
                                         scGemms$true_tumor))

data_tcell$tot <- rowSums(data_tcell)
data_tcell <- data_tcell/data_tcell$tot
data_tcell$labs <- data_tcell$Tumor

data_tcell <- data_tcell[order(data_tcell$Tumor, decreasing = T),]
data_tcell$cell_type_secondary <- factor(rownames(data_tcell))

data_tcell <- tidyr::pivot_longer(data_tcell,cols= c("Non-tumor", "Tumor"),
                                  names_to = "type")
data_tcell[data_tcell$type == "Tumor",]$labs <- NA


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

data_tcell$cell_type_secondary <- factor(x = data_tcell$cell_type_secondary,
                                         levels = rev(order_vec))


stack_plot <-  ggplot(data_tcell, aes(x = (cell_type_secondary),
                                      y =100* value,  fill = type,
                                      )) +
        geom_col() + scale_fill_excel_new() +
        theme_classic() + labs(x = "", y = "") + 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
              legend.text = element_text(size = 16),
              axis.text.x = element_text(angle = 0,
                                         hjust = 0.5,
                                         size = 16),
              axis.text.y = element_text(angle = 0,
                                         hjust = 1,
                                         size = 20, face = "plain",
                                         colour = "black"),
              axis.line = element_blank(),
              axis.ticks = element_blank()) +
        coord_flip()

ggsave("figures/S2C_tumor_cell_proportion.pdf",stack_plot,
       height = 14 * 0.5, width = 6.5)







data_tcell <- as.data.frame.matrix(table(scGemms$cell_type_secondary,
                                         scGemms$erbb_score))
data_tcell$tot <- rowSums(data_tcell)
data_tcell <- data_tcell/data_tcell$tot

data_tcell <- data_tcell[order(data_tcell$`HER2 POS`, decreasing = T),]
data_tcell$cell_type_secondary <- factor(rownames(data_tcell))

data_tcell$labs <- data_tcell$`HER2 POS`

data_tcell <- tidyr::pivot_longer(data_tcell,cols= c("HER2 NEG","HER2 POS"),
                                  names_to = "type")
data_tcell[data_tcell$type == "HER2 POS",]$labs <- NA




data_tcell$cell_type_secondary <- factor(x = data_tcell$cell_type_secondary, 
                                         levels = rev(order_vec))

stack_plot <-  ggplot(data_tcell, aes(x = (cell_type_secondary),
                                      y = 100* value,  fill = type)) +
        geom_col() + scale_fill_brewer() +
        theme_classic() + labs(x = "", y = "") + 
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
              legend.text = element_text(size = 16),
              axis.text.x = element_text(angle = 0,
                                         hjust = 0.5,
                                         size = 16),
              axis.text.y = element_text(angle = 0,
                                         hjust = 1,
                                         size = 20, face = "plain",
                                         colour = "black"),
              axis.line = element_blank(),axis.ticks = element_blank()) +
        coord_flip()


ggsave("figures/S2C_HER2_cell_proportion.pdf",stack_plot,
       height = 14 * 0.5, width = 6.5)



