rm(list = ls())


library(ggplot2)
library(harmony)
library(Seurat)
library(stringr)
library(RColorBrewer)

# aggregated results
# Directory to aggregated count matrix
data_dir <- 'filtered_feature_bc_matrix/'

expression_matrix <- Read10X(data.dir = data_dir)
keep_feature      <- Matrix::rowSums(expression_matrix)
# Filtering out genes with 30 or less associated reads
expression_matrix <- expression_matrix[keep_feature > 30,]



scGemms <- CreateSeuratObject(counts = expression_matrix,
                              project = "GEMM",
                              min.cells = 10, 
                              min.features = 200)

# Defining mitochondrial and ribosomal content
scGemms[["percent.mt"]] <- PercentageFeatureSet(scGemms, pattern = "^mt-")
scGemms[["percent.rp"]] <- PercentageFeatureSet(scGemms, pattern = "^Rpl|^Rps")



samp_tag <- str_sub(rownames(scGemms@meta.data),-1)


# Stage refers to time since tumor induction
scGemms[["stage"]]    <- plyr::mapvalues(samp_tag,
                                         1:8,
                                         c("early","early","late","late",
                                           "early","early","late","late"))

scGemms[["sampleID"]] <- plyr::mapvalues(samp_tag,
                                         1:8,
                                         c("E1","E2","L1","L2","ET1",
                                           "ET2","LT1","LT2"))
scGemms[["Batch"]]    <- plyr::mapvalues(samp_tag,
                                         1:8,
                                         c("B1","B1","B1","B1","B2",
                                           "B2","B2","B2"))

# Quality filter
scGemms <- subset(scGemms, subset = nFeature_RNA > 500 & nCount_RNA > 1000 &
                      percent.mt < 25 & percent.rp > 1)




scGemms   <- NormalizeData(scGemms)
scGemms   <- FindVariableFeatures(scGemms,
                                  selection.method = "vst",
                                  nfeatures = 3000)
all.genes <- rownames(scGemms)
scGemms   <- ScaleData(scGemms, features = all.genes)
scGemms   <- RunPCA(scGemms,
                    features = VariableFeatures(object = scGemms))
scGemms   <- RunHarmony(scGemms,"Batch")



scGemms <- FindNeighbors(scGemms, dims = 1:30, reduction = "harmony")
scGemms <- FindClusters(scGemms, resolution = 2.0)

scGemms <- RunUMAP(scGemms, dims = 1:30,
                   reduction = "harmony",
                   n.neighbors = 10,
                   min.dist = 0.3)

Idents(scGemms) <- scGemms$seurat_clusters
saveRDS(scGemms,"data/01_unfiltered_processed_sObj.RDS")

# At res 2.0 these are outlier clusters
exclude <- c("45","42","41","40","35", "29")
scGemms[["keep"]] <- plyr::mapvalues(!(scGemms$seurat_clusters %in% exclude),
                                     c(FALSE,TRUE),c("discard","keep"))

scGemms <- subset(scGemms, subset = keep == "keep")


saveRDS(scGemms,"data/02_filtered_sObj.RDS")


