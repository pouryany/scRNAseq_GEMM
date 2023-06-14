rm(list = ls())

isg <- data.table::fread("data/20220913071924_GeneSearchResults.txt")
deg <- read.csv("output/new_tumors/tumor_late_vs_early.csv",row.names = 1)

deg <- deg[deg$avg_log2FC<0,]

deg <- merge(deg, isg, by.x = "Gene", by.y = "Gene Name" )

write.csv(deg,"output/new_tumors/ISG_downregulated.csv")
