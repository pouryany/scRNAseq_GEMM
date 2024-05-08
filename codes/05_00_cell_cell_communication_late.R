rm(list = ls())

#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(harmony)
library(stringr)
library(ggplot2)
library(RColorBrewer)

scGemms <- readRDS("data/03_umap_labled_1_sObj.RDS")
scGemms$initial_clusts <- scGemms$seurat_clusters

scGemms0 <- subset(scGemms,
                   cells = names(grep(pattern = "Fibroblast|Undefined|Alv 4",
                                               scGemms$cell_type_secondary,
                                               value = T,
                                               invert = T)))

colorss <- readRDS("data/03_colors.RDS")



data.input = scGemms0@assays$RNA@data # normalized data matrix


meta = scGemms0@meta.data 
# a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$stage == "Late"] 
# extract the cell names from disease data

# Prepare input data for CelChat analysis
data.input = data.input[, cell.use]
meta = meta[cell.use, ]

meta$labels = droplevels(meta$cell_type_secondary)


order_vec <-  c(
    "Alv-Proliferating","Alv-Basal",
    "Basal",
    "Alv 1","Alv 2","Alv 3",
    "Hormone Sensing",
    "CD4+ Naïve", "CD4+ Treg",             
    "CD8+ Naïve", "CD8+ CM","CD8+ ISG",
    "CD8+ Effector",
    "CD8+ TRM",
    "CD8+ Progenitor Ex",
    "CD8+ Ex", 
    "CD8+ Proliferating",
    "NK",
    "Macrophage",
    "Neutrophil","B cell")    

meta$labels <- factor(meta$labels, levels= order_vec)


cellchat <- createCellChat(object = data.input, meta = meta,
                           group.by = "labels")
#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") 
# set "labels" as default cell identity


CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use


cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)



# project gene expression data onto PPI (Optional: when running it,
# USER should set `raw.use = FALSE` in the function `computeCommunProb()` 
# in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# the slot 'netP' means the inferred intercellular communication network of
#signaling pathways



saveRDS(cellchat, file = "data/05_cellchat_scGemms_late.rds")
