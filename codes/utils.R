VlnPlot_costumized <- function(seurat_object,my_features,
                               split_var= "stage",
                               split_levels = c("Early","Late"),
                               text_y = 6,idents_order,signif_list,
                               limits =  c(0,7),plot_alphabetic = T, x_size = 20){
    
    data_text <- expand.grid(feature = my_features,
                             ident = idents_order,
                             split = split_levels) 
    
    data_text <- dplyr::left_join(data_text,signif_list)
    data_text$y <- text_y
    if(plot_alphabetic == FALSE)
        data_text$feature <- factor(data_text$feature,levels  = my_features)
    
    vln_eGFP1 <- VlnPlot(seurat_object,split.by = split_var,features =my_features,
                         #cols =  colors0[sort(names(colors0))],
                         same.y.lims = T,stack = T,flip = T, fill.by = "ident",
                         combine = F) &
        scale_y_continuous(limits = limits) &
        theme(axis.title.x = element_blank(),
              axis.text.x =element_text(size = x_size, face = "plain", angle = 90, vjust = 0.5),
              axis.text.y = element_text(size = 10, face = "plain"),
              strip.text =element_text(size = 16, face = "plain",vjust = 1,hjust = 0),
              plot.title = element_text(size = 30, face = "plain"),
              axis.title.y = element_blank(), title = element_text(size = 20, face = "plain"))&
        geom_text(data = data_text, aes(x = ident, y = y,feature = feature, label = my_text))
    
    return(vln_eGFP1)
}

VlnPlot_costumized2 <- function(seurat_object,my_features,split_var= "stage",split_levels = c("Early","Late"),
                                text_y = 6,idents_order,signif_list,
                                limits =  c(0,7),plot_alphabetic = T){
    
    data_text <- expand.grid(feature = my_features, ident = idents_order, split = split_levels) 
    
    data_text <- left_join(data_text,signif_list)
    data_text$y <- text_y
    if(plot_alphabetic == FALSE)
        data_text$feature <- factor(data_text$feature,levels  = my_features)
    
    vln_eGFP1 <- VlnPlot(seurat_object,split.by = split_var,features =my_features,
                         #cols =  colors0[sort(names(colors0))],
                         same.y.lims = T,stack = T,flip = T, fill.by = "ident", combine = F) &
        scale_y_continuous(limits = limits)&
        geom_text(data = data_text, aes(x = ident, y = y,feature = feature, label = my_text))
    
    return(vln_eGFP1)
}



GeneID2entrez_modified <- function (gene.IDs,
                                return.Matrix = FALSE,
                                mode = "h2h") 
        {
            stopifnot(mode %in% c("h2h", "m2m", "m2h"))
            gene.IDs <- gene.IDs[!is.na(gene.IDs)]
            gene.IDs <- trimws(gene.IDs)
            if (mode == "h2h") {
                gene.IDs <- toupper(gene.IDs)
                if (all(grepl("^ENSG", gene.IDs, perl = TRUE))) {
                    ID.type <- "ENSEMBL"
                    suppressMessages(GeneNames <- biomaRt::select(org.Hs.eg.db, 
                                                                  gene.IDs, c("ENTREZID", "ENSEMBL"), keytype = "ENSEMBL"))
                }
                else {
                    ID.type <- "SYMBOL"
                    suppressMessages(GeneNames <- biomaRt::select(org.Hs.eg.db, 
                                                                  gene.IDs, c("SYMBOL", "ENTREZID"), keytype = "ALIAS"))
                }
                matched <- match(as.character(gene.IDs), GeneNames[, 
                                                                   ID.type])
                matched_2 <- match(GeneNames[, ID.type], as.character(gene.IDs))
                if (sum(duplicated(matched_2[!is.na(matched_2)])) > 0) {
                    warning("Some genes returned 1:many mapping to ENTREZ ID. ", 
                            "Genes were assigned the first ENTREZ ID match found.\n", 
                            call. = FALSE)
                }
                message("Done! ", sum(!is.na(matched)), " genes of ", 
                        length(matched), " successfully converted.\n")
                if (return.Matrix == TRUE) {
                    if (sum(is.na(matched)) > 0) {
                        message("Couldn't find Entrez IDs for ", sum(is.na(matched)), 
                                " genes (NAs returned instead).\n")
                    }
                    return(data.frame(GENE.ID = gene.IDs, ENTREZ.ID = GeneNames[matched, 
                                                                                "ENTREZID"], stringsAsFactors = FALSE))
                }
                else {
                    if (sum(is.na(matched)) > 0) {
                        message("Couldn't find Entrez IDs for ", sum(is.na(matched)), 
                                " genes.\n")
                    }
                    return(GeneNames[matched[!is.na(matched)], "ENTREZID"])
                }
            }
            else if (mode == "m2m") {
                biomart_test <- tryCatch({
                    R.utils::withTimeout({
                        tmp <- biomaRt::listMarts()
                    }, timeout = 3, onTimeout = "warning")
                }, error = function(w) {
                    return(0)
                }, warning = function(w) {
                    return(0)
                })
                if (all(biomart_test == 0)) {
                    stop(paste0("We are having trouble reaching biomaRt.\n", 
                                "Please, try again later."))
                }
                mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
                GeneNames <- biomaRt::getBM(attributes = c("ensembl_gene_id", 
                                                           "mgi_symbol", "entrezgene_id"), values = "*", mart = mouse)
                if (all(grepl("^ENSM", gene.IDs, perl = TRUE))) {
                    ids <- GeneNames[match(gene.IDs, GeneNames$ensembl_gene_id), 
                                     c(1, 3)]
                    colnames(ids) <- c("GENE.ID", "ENTREZ.ID")
                }
                else {
                    ids <- GeneNames[match(gene.IDs, GeneNames$mgi_symbol), 
                                     c(2, 3)]
                    colnames(ids) <- c("GENE.ID", "ENTREZ.ID")
                }
                message("Done! ", sum(!is.na(ids$ENTREZ.ID)), " genes of ", 
                        dim(ids)[1], " successfully converted.\n")
                if (return.Matrix == TRUE) {
                    if (sum(is.na(ids$ENTREZ.ID)) > 0) {
                        message("Couldn't find Entrez IDs for ", sum(is.na(ids$ENTREZ.ID)), 
                                " genes (NAs returned instead).\n")
                    }
                    return(ids)
                }
                else {
                    if (sum(is.na(ids$ENTREZ.ID)) > 0) {
                        message("Couldn't find human Entrez IDs for ", 
                                sum(is.na(ids$ENTREZ.ID)), " genes.\n")
                    }
                    return(ids$ENTREZ.ID[!is.na(ids$ENTREZ.ID)])
                }
            }
            else if (mode == "m2h") {
                biomart_test <- tryCatch({
                    R.utils::withTimeout({
                        tmp <- biomaRt::listMarts()
                    }, timeout = 3, onTimeout = "warning")
                }, error = function(w) {
                    return(0)
                }, warning = function(w) {
                    return(0)
                })
                if (all(biomart_test == 0)) {
                    stop(paste0("We are having trouble reaching biomaRt.\n", 
                                "Please, try again later."))
                }
                human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
                mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
                if (all(grepl("^ENSM", gene.IDs, perl = TRUE)) == TRUE) {
                    hs_ids = biomaRt::getLDS(attributes = c("ensembl_gene_id"), 
                                             filters = "ensembl_gene_id", values = gene.IDs, 
                                             mart = mouse, attributesL = c("entrezgene_id"), 
                                             martL = human, uniqueRows = TRUE)
                    hs_ids <- data.frame(mouse.gene.ID = gene.IDs, human.gene.ID = hs_ids$NCBI.gene.ID[match(gene.IDs, 
                                                                                                             hs_ids$Gene.stable.ID)], stringsAsFactors = F)
                }
                else if (all(grepl("^\\d*$", gene.IDs, perl = TRUE)) == 
                         TRUE) {
                    hs_ids = biomaRt::getLDS(attributes = c("entrezgene_id"), 
                                             filters = "entrezgene_id", values = gene.IDs, 
                                             mart = mouse, attributesL = c("entrezgene_id"), 
                                             martL = human, uniqueRows = TRUE)
                    hs_ids <- data.frame(mouse.gene.ID = gene.IDs,
                                         human.gene.ID = hs_ids$NCBI.gene.ID.1[match(gene.IDs, hs_ids$NCBI.gene.ID)],
                                         stringsAsFactors = F)
                }
                else {
                    hs_ids = biomaRt::getLDS(attributes = c("mgi_symbol"), 
                                             filters = "mgi_symbol", values = gene.IDs, mart = mouse, 
                                             attributesL = c("entrezgene_id"), martL = human, 
                                             uniqueRows = TRUE)
                    hs_ids <- data.frame(mouse.gene.ID = gene.IDs,
                                         human.gene.ID = hs_ids[match(gene.IDs,hs_ids$MGI.symbol),2],
                                         stringsAsFactors = F)
                }
                message("Done! ", sum(!is.na(hs_ids$human.gene.ID)), 
                        " genes of ", length(gene.IDs), " successfully converted.\n")
                if (return.Matrix == TRUE) {
                    if (any(is.na(hs_ids$human.gene.ID))) {
                        message("Couldn't find human Entrez IDs for ", 
                                sum(is.na(hs_ids$human.gene.ID)), " genes (NAs returned instead).\n")
                    }
                    return(hs_ids)
                }
                else {
                    if (any(is.na(hs_ids$human.gene.ID))) {
                        message("Couldn't find human Entrez IDs for ", 
                                sum(is.na(hs_ids$human.gene.ID)), " genes.\n")
                    }
                    return(hs_ids$human.gene.ID[!is.na(hs_ids$human.gene.ID)])
                }
            }
}




netClustering_modified <- 
    function (object, slot.name = "netP", type = c("functional", 
    "structural"), comparison = NULL, k = NULL, methods = "kmeans", 
          do.plot = TRUE, fig.id = NULL, do.parallel = TRUE, nCores = 4, 
          k.eigen = NULL) 
{
    type <- match.arg(type)
    if (object@options$mode == "single") {
        comparison <- "single"
        cat("Classification learning of the signaling networks for a single dataset", 
            "\n")
    }
    else if (object@options$mode == "merged") {
        if (is.null(comparison)) {
            comparison <- 1:length(unique(object@meta$datasets))
        }
        cat("Classification learning of the signaling networks for datasets", 
            as.character(comparison), "\n")
    }
    comparison.name <- paste(comparison, collapse = "-")
    Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
    pathways.ignore <- rownames( Y[rowSums(!is.finite(Y))>0, ] )
    cellchat@options$pathways.ignore = pathways.ignore
    Y <- Y[!rowSums(!is.finite(Y)),] # filter out rows with NaN, not working downstream
    methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]] <- Y
    data.use <- Y
    if (methods == "kmeans") {
        if (!is.null(k)) {
            clusters = kmeans(data.use, k, nstart = 10)$cluster
        }
        else {
            N <- nrow(data.use)
            kRange <- seq(2, min(N - 1, 10), by = 1)
            if (do.parallel) {
                future::plan("multiprocess", workers = nCores)
                options(future.globals.maxSize = 1000 * 1024^2)
            }
            my.sapply <- ifelse(test = future::nbrOfWorkers() == 
                                    1, yes = pbapply::pbsapply, no = future.apply::future_sapply)
            results = my.sapply(X = 1:length(kRange), FUN = function(x) {
                idents <- kmeans(data.use, kRange[x], nstart = 10)$cluster
                clusIndex <- idents
                adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, 
                                                           clusIndex, FUN = "==")), nrow = N, ncol = N)
                return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
            }, simplify = FALSE)
            adjMat <- lapply(results, "[[", 1)
            CM <- Reduce("+", adjMat)/length(kRange)
            res <- computeEigengap(as.matrix(CM))
            numCluster <- res$upper_bound
            clusters = kmeans(data.use, numCluster, nstart = 10)$cluster
            if (do.plot) {
                gg <- res$gg.obj
                ggsave(filename = paste0("estimationNumCluster_", 
                                         fig.id, "_", type, "_dataset_", comparison.name, 
                                         ".pdf"), plot = gg, width = 3.5, height = 3, 
                       units = "in", dpi = 300)
            }
        }
    }
    else if (methods == "spectral") {
        A <- as.matrix(data.use)
        D <- apply(A, 1, sum)
        L <- diag(D) - A
        L <- diag(D^-0.5) %*% L %*% diag(D^-0.5)
        evL <- eigen(L, symmetric = TRUE)
        plot(rev(evL$values)[1:30])
        Z <- evL$vectors[, (ncol(evL$vectors) - k.eigen + 1):ncol(evL$vectors)]
        clusters = kmeans(Z, k, nstart = 20)$cluster
    }
    if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$group)) {
        methods::slot(object, slot.name)$similarity[[type]]$group <- NULL
    }
    methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]] <- clusters
    return(object)
}



netAnalysis_signalingRole_heatmap_modified <- function(object, signaling = NULL,
                                                       pattern = c("outgoing", "incoming", 
                                                "all"), slot.name = "netP",
                                                color.use = NULL, color.heatmap = "BuGn", 
          title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
          cluster.rows = FALSE, cluster.cols = FALSE,
          show_column_dend = FALSE, show_row_dend = FALSE,
          col_order = NULL,show_annot = FALSE, name = "Heatmap1") 
{
    pattern <- match.arg(pattern)
    if (length(slot(object, slot.name)$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
        outgoing[, i] <- centr[[i]]$outdeg
        incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
        mat <- t(outgoing)
        legend.name <- "Outgoing"
    }
    else if (pattern == "incoming") {
        mat <- t(incoming)
        legend.name <- "Incoming"
    }
    else if (pattern == "all") {
        mat <- t(outgoing + incoming)
        legend.name <- "Overall"
    }
    if (is.null(title)) {
        title <- paste0(legend.name, " signaling patterns")
    }
    else {
        title <- paste0(paste0(legend.name, " signaling patterns"), 
                        " - ", title)
    }
    if (!is.null(signaling)) {
        mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
        mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
        idx <- match(rownames(mat1), signaling)
        mat[idx[!is.na(idx)], ] <- mat1
        dimnames(mat) <- list(signaling, colnames(mat1))
    }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[mat == 0] <- NA
    if (is.null(color.use)) {
        color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                               name = color.heatmap))))(100)
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    names(color.use) <- colnames(mat)
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                        which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                        simple_anno_size = grid::unit(0.2, "cm"))
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                    border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                            show_annotation_name = FALSE)
    pSum <- rowSums(mat.ori)
    pSum.original <- pSum
    pSum <- -1/log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    if (length(idx1) > 0) {
        values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                             length.out = length(idx1))
        position <- sort(pSum.original[idx1], index.return = TRUE)$ix
        pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                        show_annotation_name = FALSE)
    if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
        legend.break <- max(mat, na.rm = T)
    }
    else {
        legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                          round(max(mat, na.rm = T), digits = 1))
    }
    mat[is.na(mat)] <- 0
    col_order0 <- seq_len(ncol(mat))
    if(is.null(col_order))
        col_order <- col_order0
    if(show_annot){
        ht1 = Heatmap(mat[,col_order], col = color.heatmap.use, 
                      na_col = "white", name = name(), bottom_annotation = col_annotation, 
                      top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                      cluster_columns = cluster.cols, row_names_side = "left", 
                      show_column_dend = show_column_dend, show_row_dend = show_row_dend,
                      row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                      column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                                                                                 "cm"), height = unit(height, "cm"), column_title = title, 
                      column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                      heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                  fontface = "plain"), title_position = "leftcenter-rot", 
                                                  border = NA, at = legend.break,
                                                  legend_height = unit(20, "mm"),
                                                  labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm")))
                                                                                                                                                                
    }else{
        
    ht1 = Heatmap(mat[,col_order], col = color.heatmap.use, na_col = "white" ,
                  cluster_rows = cluster.rows, 
                  name = name,
                  cluster_columns = cluster.cols, row_names_side = "left", 
                  show_column_dend = show_column_dend, show_row_dend = show_row_dend,
                  row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                  column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"), column_title = title, 
                  column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),
                                              title_position = "leftcenter-rot", 
                                              border = NA, at = legend.break, legend_height = unit(20, "mm"),
                                              labels_gp = gpar(fontsize = 8), grid_width = unit(2,"mm")))
    }
    return(ht1)
        
}







netVisual_circle_modified <- function (net, color.use = NULL, title.name = NULL, sources.use = NULL, 
                                       targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, 
                                       top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                                       vertex.size.max = NULL, vertex.label.cex = 1, vertex.label.color = "black", 
                                       edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                                       label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
                                       edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                                       margin = 0.2, vertex.size = NULL, arrow.width = 1, arrow.size = 0.2,
                                       seed = 1,low_dist = 0,high_dist = 3) 
{
    if (!is.null(vertex.size)) {
        warning("'vertex.size' is deprecated. Use `vertex.weight`")
    }
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }
    options(warn = -1)
    thresh <- stats::quantile(net, probs = 1 - top)
    net[net < thresh] <- 0
    if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use))) {
        if (is.null(rownames(net))) {
            stop("The input weighted matrix should have rownames!")
        }
        cells.level <- rownames(net)
        df.net <- reshape2::melt(net, value.name = "value")
        colnames(df.net)[1:2] <- c("source", "target")
        if (!is.null(sources.use)) {
            if (is.numeric(sources.use)) {
                sources.use <- cells.level[sources.use]
            }
            df.net <- subset(df.net, source %in% sources.use)
        }
        if (!is.null(targets.use)) {
            if (is.numeric(targets.use)) {
                targets.use <- cells.level[targets.use]
            }
            df.net <- subset(df.net, target %in% targets.use)
        }
        if (!is.null(idents.use)) {
            if (is.numeric(idents.use)) {
                idents.use <- cells.level[idents.use]
            }
            df.net <- filter(df.net, (source %in% idents.use) | 
                                 (target %in% idents.use))
        }
        df.net$source <- factor(df.net$source, levels = cells.level)
        df.net$target <- factor(df.net$target, levels = cells.level)
        df.net$value[is.na(df.net$value)] <- 0
        net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                              df.net[["target"]]), sum)
    }
    net[is.na(net)] <- 0
    if (remove.isolate) {
        idx1 <- which(Matrix::rowSums(net) == 0)
        idx2 <- which(Matrix::colSums(net) == 0)
        idx <- intersect(idx1, idx2)
        net <- net[-idx, ]
        net <- net[, -idx]
    }
    g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                     weighted = T)
    edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
    coords <- layout_(g, layout)
    if (nrow(coords) != 1) {
        coords_scale = scale(coords)
    }
    else {
        coords_scale <- coords
    }
    if (is.null(color.use)) {
        color.use = scPalette(length(igraph::V(g)))
    }
    if (is.null(vertex.weight.max)) {
        vertex.weight.max <- max(vertex.weight)
    }
    vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
        5
    loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0,
                         -atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g), 1]),
                         pi - atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g), 1]))
    igraph::V(g)$size <- vertex.weight
    igraph::V(g)$color <- color.use[igraph::V(g)]
    igraph::V(g)$frame.color <- color.use[igraph::V(g)]
    igraph::V(g)$label.color <- vertex.label.color
    igraph::V(g)$label.cex <- vertex.label.cex
    if (label.edge) {
        igraph::E(g)$label <- igraph::E(g)$weight
        igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
    }
    if (is.null(edge.weight.max)) {
        edge.weight.max <- max(igraph::E(g)$weight)
    }
    if (weight.scale == TRUE) {
        igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
            edge.width.max
    }
    else {
        igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
    }
    igraph::E(g)$arrow.width <- arrow.width
    igraph::E(g)$arrow.size <- arrow.size
    igraph::E(g)$label.color <- edge.label.color
    igraph::E(g)$label.cex <- edge.label.cex
    igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[, 
                                                                               1]], alpha.edge)
    if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
        igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                    1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                                   2] == edge.start[, 1]), 1]]
    }
    radian.rescale <- function(x, start = 0, direction = 1) {
        c.rotate <- function(x) (x + start)%%(2 * pi) * direction
        c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }
    label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                                 direction = -1, start = 0)
    
    set.seed(seed)
    label.dist <- vertex.weight/max(vertex.weight) + runif(length(igraph::V(g)),low_dist,high_dist)
    plot(g, edge.curved = edge.curved, vertex.shape = shape, 
         layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
         vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
         edge.label.family = "Helvetica")
    if (!is.null(title.name)) {
        text(0, 1.5, title.name, cex = 1.1)
    }
    gg <- recordPlot()
    return(gg)
}



netVisual_aggregate_modified <- function (object, signaling, signaling.name = NULL, color.use = NULL, 
                                          vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, 
                                          idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = 1, 
                                          vertex.weight.max = NULL, vertex.size.max = NULL, weight.scale = TRUE, 
                                          edge.weight.max = NULL, edge.width.max = 8, layout = c("circle", 
                                                                                                 "hierarchy", "chord"), thresh = 0.05, from = NULL, to = NULL, 
                                          bidirection = NULL, vertex.size = NULL, pt.title = 12, title.space = 6, 
                                          vertex.label.cex = 0.8, group = NULL, cell.order = NULL, 
                                          small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, 
                                          show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, 
                                          ...) 
{
    layout <- match.arg(layout)
    if (!is.null(vertex.size)) {
        warning("'vertex.size' is deprecated. Use `vertex.weight`")
    }
    if (is.null(vertex.weight)) {
        vertex.weight <- as.numeric(table(object@idents))
    }
    if (is.null(vertex.size.max)) {
        if (length(unique(vertex.weight)) == 1) {
            vertex.size.max <- 5
        }
        else {
            vertex.size.max <- 15
        }
    }
    pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, 
                         key = "pathway_name", matching.exact = T, pair.only = T)
    if (is.null(signaling.name)) {
        signaling.name <- signaling
    }
    net <- object@net
    pairLR.use.name <- dimnames(net$prob)[[3]]
    pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
    pairLR <- pairLR[pairLR.name, ]
    prob <- net$prob
    pval <- net$pval
    prob[pval > thresh] <- 0
    if (length(pairLR.name) > 1) {
        pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 
                                             3, sum) != 0]
    }
    else {
        pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 
                                           0]
    }
    if (length(pairLR.name.use) == 0) {
        stop(paste0("There is no significant communication of ", 
                    signaling.name))
    }
    else {
        pairLR <- pairLR[pairLR.name.use, ]
    }
    nRow <- length(pairLR.name.use)
    prob <- prob[, , pairLR.name.use]
    pval <- pval[, , pairLR.name.use]
    if (length(dim(prob)) == 2) {
        prob <- replicate(1, prob, simplify = "array")
        pval <- replicate(1, pval, simplify = "array")
    }
    if (layout == "hierarchy") {
        prob.sum <- apply(prob, c(1, 2), sum)
        if (is.null(edge.weight.max)) {
            edge.weight.max = max(prob.sum)
        }
        par(mfrow = c(1, 2), ps = pt.title)
        netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, 
                             sources.use = sources.use, targets.use = targets.use, 
                             remove.isolate = remove.isolate, top = top, color.use = color.use, 
                             vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, 
                             vertex.size.max = vertex.size.max, weight.scale = weight.scale, 
                             edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, 
                             title.name = NULL, vertex.label.cex = vertex.label.cex, 
                             ...)
        netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum), 
                                                                 vertex.receiver), sources.use = sources.use, targets.use = targets.use, 
                             remove.isolate = remove.isolate, top = top, color.use = color.use, 
                             vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, 
                             vertex.size.max = vertex.size.max, weight.scale = weight.scale, 
                             edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, 
                             title.name = NULL, vertex.label.cex = vertex.label.cex, 
                             ...)
        graphics::mtext(paste0(signaling.name, " signaling pathway network"), 
                        side = 3, outer = TRUE, cex = 1, line = -title.space)
        gg <- recordPlot()
    }
    else if (layout == "circle") {
        prob.sum <- apply(prob, c(1, 2), sum)
        gg <- netVisual_circle_modified(prob.sum, sources.use = sources.use, 
                                        targets.use = targets.use, idents.use = idents.use, 
                                        remove.isolate = remove.isolate, top = top, color.use = color.use, 
                                        vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, 
                                        vertex.size.max = vertex.size.max, weight.scale = weight.scale, 
                                        edge.weight.max = edge.weight.max, edge.width.max = edge.width.max, 
                                        title.name = paste0(signaling.name, " signaling pathway network"), 
                                        vertex.label.cex = vertex.label.cex, ...)
    }
    else if (layout == "chord") {
        prob.sum <- apply(prob, c(1, 2), sum)
        gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, 
                                            sources.use = sources.use, targets.use = targets.use, 
                                            remove.isolate = remove.isolate, group = group, cell.order = cell.order, 
                                            lab.cex = vertex.label.cex, small.gap = small.gap, 
                                            big.gap = big.gap, scale = scale, reduce = reduce, 
                                            title.name = paste0(signaling.name, " signaling pathway network"), 
                                            show.legend = show.legend, legend.pos.x = legend.pos.x, 
                                            legend.pos.y = legend.pos.y)
    }
    return(gg)
}






netVisual_bubble_modified <- 
    function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
              pairLR.use = NULL, color.heatmap = c("Spectral", "viridis"), 
              n.colors = 10, direction = -1, thresh = 0.05, comparison = NULL, 
              group = NULL, remove.isolate = FALSE, max.dataset = NULL, 
              min.dataset = NULL, min.quantile = 0, max.quantile = 1, line.on = TRUE, 
              line.size = 0.2, color.text.use = TRUE, color.text = NULL, 
              title.name = NULL, font.size = 10, font.size.title = 10, 
              show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
              angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE,
              point_scaler = 3, x.size = 10, y.size = 10) {
        color.heatmap <- match.arg(color.heatmap)
        if (is.list(object@net[[1]])) {
            message("Comparing communications on a merged object \n")
        }
        else {
            message("Comparing communications on a single object \n")
        }
        if (is.null(vjust.x) | is.null(hjust.x)) {
            angle = c(0, 45, 90)
            hjust = c(0, 1, 1)
            vjust = c(0, 1, 0.5)
            vjust.x = vjust[angle == angle.x]
            hjust.x = hjust[angle == angle.x]
        }
        if (length(color.heatmap) == 1) {
            color.use <- tryCatch({
                RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
            }, error = function(e) {
                (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
            })
        }
        else {
            color.use <- color.heatmap
        }
        if (direction == -1) {
            color.use <- rev(color.use)
        }
        if (is.null(comparison)) {
            cells.level <- levels(object@idents)
            if (is.numeric(sources.use)) {
                sources.use <- cells.level[sources.use]
            }
            if (is.numeric(targets.use)) {
                targets.use <- cells.level[targets.use]
            }
            df.net <- subsetCommunication(object, slot.name = "net", 
                                          sources.use = sources.use, targets.use = targets.use, 
                                          signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
            df.net$source.target <- paste(df.net$source, df.net$target, 
                                          sep = " -> ")
            source.target <- paste(rep(sources.use, each = length(targets.use)), 
                                   targets.use, sep = " -> ")
            source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
            if (length(source.target.isolate) > 0) {
                df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                                       ncol = ncol(df.net)))
                colnames(df.net.isolate) <- colnames(df.net)
                df.net.isolate$source.target <- source.target.isolate
                df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
                df.net.isolate$pval <- 1
                a <- stringr::str_split(df.net.isolate$source.target, 
                                        " -> ", simplify = T)
                df.net.isolate$source <- as.character(a[, 1])
                df.net.isolate$target <- as.character(a[, 2])
                df.net <- rbind(df.net, df.net.isolate)
            }
            df.net$pval[df.net$pval > 0.05] = 1
            df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
            df.net$pval[df.net$pval <= 0.01] = 3
            df.net$prob[df.net$prob == 0] <- NA
            df.net$prob.original <- df.net$prob
            df.net$prob <- -1/log(df.net$prob)
            idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
                              0)
            if (sum(idx1) > 0) {
                values.assign <- seq(max(df.net$prob, na.rm = T) * 
                                         1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
                position <- sort(prob.original[idx1], index.return = TRUE)$ix
                df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                                                         position)]
            }
            df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                            unique(df.net$source)])
            df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                            unique(df.net$target)])
            group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                                 levels(df.net$target), sep = " -> ")
            df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
            
            ##Modified here
            
            df.net$interaction_name_2 <- paste0(df.net$pathway_name, ": ",df.net$interaction_name_2)
            # End modification
            df.net <- with(df.net, df.net[order(interaction_name_2), 
            ])
            df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
                                                levels = unique(df.net$interaction_name_2))
            cells.order <- group.names
            df.net$source.target <- factor(df.net$source.target, 
                                           levels = cells.order)
            df <- df.net
        }
        else {
            dataset.name <- names(object@net)
            df.net.all <- subsetCommunication(object, slot.name = "net", 
                                              sources.use = sources.use, targets.use = targets.use, 
                                              signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
            df.all <- data.frame()
            for (ii in 1:length(comparison)) {
                cells.level <- levels(object@idents[[comparison[ii]]])
                if (is.numeric(sources.use)) {
                    sources.use <- cells.level[sources.use]
                }
                if (is.numeric(targets.use)) {
                    targets.use <- cells.level[targets.use]
                }
                df.net <- df.net.all[[comparison[ii]]]
                
                ##Modified here
                
                tryCatch({df.net$interaction_name_2 <- 
                    paste0(df.net$pathway_name, ": ",df.net$interaction_name_2)},
                    error = function(w) {
                        print("No early interactions")
                    })
                # End modification
                
                
                df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
                df.net$source.target <- paste(df.net$source, df.net$target, 
                                              sep = " -> ")
                source.target <- paste(rep(sources.use, each = length(targets.use)), 
                                       targets.use, sep = " -> ")
                source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
                if (length(source.target.isolate) > 0) {
                    df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                                           ncol = ncol(df.net)))
                    colnames(df.net.isolate) <- colnames(df.net)
                    df.net.isolate$source.target <- source.target.isolate
                    df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
                    df.net.isolate$pval <- 1
                    a <- stringr::str_split(df.net.isolate$source.target, 
                                            " -> ", simplify = T)
                    df.net.isolate$source <- as.character(a[, 1])
                    df.net.isolate$target <- as.character(a[, 2])
                    df.net <- rbind(df.net, df.net.isolate)
                }
                df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                                unique(df.net$source)])
                df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                                unique(df.net$target)])
                group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                                     levels(df.net$target), sep = " -> ")
                group.names0 <- group.names
                group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                                      ")")
                if (nrow(df.net) > 0) {
                    df.net$pval[df.net$pval > 0.05] = 1
                    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                                    0.05] = 2
                    df.net$pval[df.net$pval <= 0.01] = 3
                    df.net$prob[df.net$prob == 0] <- NA
                    df.net$prob.original <- df.net$prob
                    df.net$prob <- -1/log(df.net$prob)
                }
                else {
                    df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                                                   ncol = 5))
                    colnames(df.net) <- c("interaction_name_2", "source.target", 
                                          "prob", "pval", "prob.original")
                    df.net$source.target <- group.names0
                }
                df.net$group.names <- as.character(df.net$source.target)
                df.net$source.target <- paste0(df.net$source.target, 
                                               " (", dataset.name[comparison[ii]], ")")
                df.net$dataset <- dataset.name[comparison[ii]]
                df.all <- rbind(df.all, df.net)
            }
            if (nrow(df.all) == 0) {
                stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
            }
            idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
                              0)
            if (sum(idx1) > 0) {
                values.assign <- seq(max(df.all$prob, na.rm = T) * 
                                         1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
                position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
                df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                                                         position)]
            }
            df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
            df <- df.all
            df <- with(df, df[order(interaction_name_2), ])
            df$interaction_name_2 <- factor(df$interaction_name_2, 
                                            levels = unique(df$interaction_name_2))
            cells.order <- c()
            dataset.name.order <- c()
            for (i in 1:length(group.names0)) {
                for (j in 1:length(comparison)) {
                    cells.order <- c(cells.order, paste0(group.names0[i], 
                                                         " (", dataset.name[comparison[j]], ")"))
                    dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
                }
            }
            df$source.target <- factor(df$source.target, levels = cells.order)
        }
        min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
        max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
        df$prob[df$prob < min.cutoff] <- min.cutoff
        df$prob[df$prob > max.cutoff] <- max.cutoff
        if (remove.isolate) {
            df <- df[!is.na(df$prob), ]
            line.on <- FALSE
        }
        if (!is.null(max.dataset)) {
            signaling <- as.character(unique(df$interaction_name_2))
            for (i in signaling) {
                df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
                cell <- as.character(unique(df.i$group.names))
                for (j in cell) {
                    df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
                    values <- df.i.j$prob
                    idx.max <- which(values == max(values, na.rm = T))
                    idx.min <- which(values == min(values, na.rm = T))
                    dataset.na <- c(df.i.j$dataset[is.na(values)], 
                                    setdiff(dataset.name[comparison], df.i.j$dataset))
                    if (length(idx.max) > 0) {
                        if (!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset])) {
                            df.i.j$prob <- NA
                        }
                        else if ((idx.max != idx.min) & !is.null(min.dataset)) {
                            if (!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset])) {
                                df.i.j$prob <- NA
                            }
                            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                                                                    dataset.na)) > 0) {
                                df.i.j$prob <- NA
                            }
                        }
                    }
                    df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
                }
                df[df$interaction_name_2 == i, "prob"] <- df.i$prob
            }
        }
        if (remove.isolate)  {
            df <- df[!is.na(df$prob), ]
            line.on <- FALSE
        }
        if (nrow(df) == 0) {
            stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
        }
        df$interaction_name_2 <- factor(df$interaction_name_2, levels = unique(df$interaction_name_2))
        df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
                                                                          unique(df$source.target)))
        g <- ggplot(df, aes(x = source.target, y = interaction_name_2, 
                            color = prob, size = 2* pval)) + geom_point(pch = 16) + 
            theme_linedraw() + theme(panel.grid.major = element_blank()) + 
            theme(axis.text.x = element_blank(),#element_text(angle = angle.x, hjust = hjust.x, 
                  #            vjust = vjust.x, size = x.size),
                  axis.text.y = element_text(size = y.size),
                  axis.title.x = element_blank(), 
                  axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
        values <- c(1, 2, 3)
        names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
        g <- g + scale_radius(range = c(min(df$pval), point_scaler*max(df$pval)), 
                              breaks = sort(unique(df$pval)), labels = names(values)[values %in% 
                                                                                         sort(unique(df$pval))], name = "p-value")
        if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
            g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                            na.value = "white", limits = c(quantile(df$prob, 
                                                                                    0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
                                            breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 
                                                                                                 1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5, 
                                                                                                                                                                            title = "Commun. Prob."))
        }
        else {
            g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                            na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5, 
                                                                                                 title = "Commun. Prob."))
        }
        g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
            theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
        if (grid.on) {
            if (length(unique(df$source.target)) > 1) {
                g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                                                         0.5, 1), lwd = 0.1, colour = color.grid)
            }
            if (length(unique(df$interaction_name_2)) > 1) {
                g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                                         0.5, 1), lwd = 0.1, colour = color.grid)
            }
        }
        if (!is.null(title.name)) {
            g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
        }
        if (!is.null(comparison)) {
            if (line.on) {
                xintercept = seq(0.5 + length(dataset.name[comparison]), 
                                 length(group.names0) * length(dataset.name[comparison]), 
                                 by = length(dataset.name[comparison]))
                g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                                    color = "grey60", size = line.size)
            }
            if (color.text.use) {
                if (is.null(group)) {
                    group <- 1:length(comparison)
                    names(group) <- dataset.name[comparison]
                }
                if (is.null(color.text)) {
                    color <- ggPalette(length(unique(group)))
                }
                else {
                    color <- color.text
                }
                names(color) <- names(group[!duplicated(group)])
                color <- color[group]
                dataset.name.order <- levels(df$source.target)
                dataset.name.order <- stringr::str_match(dataset.name.order, 
                                                         "\\(.*\\)")
                dataset.name.order <- stringr::str_sub(dataset.name.order, 
                                                       2, stringr::str_length(dataset.name.order) - 
                                                           1)
                xtick.color <- color[dataset.name.order]
                g <- g + theme(axis.text.x = element_text(colour = xtick.color))
            }
        }
        if (!show.legend) {
            g <- g + theme(legend.position = "none")
        }
        if (return.data) {
            return(list(communication = df, gg.obj = g))
        }
        else {
            return(g)
        }
    }




netVisual_chord_gene_modified <- function (object, slot.name = "net", color.use = NULL, signaling = NULL, 
          pairLR.use = NULL, net = NULL, sources.use = NULL, targets.use = NULL, 
          lab.cex = 0.8, small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03), 
          link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, 
          reduce = -1, transparency = 0.4, link.border = NA, title.name = NULL, 
          legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE, 
          thresh = 0.05, ...) 
{
    if (!is.null(pairLR.use)) {
        if (!is.data.frame(pairLR.use)) {
            stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
        }
        else if ("pathway_name" %in% colnames(pairLR.use)) {
            message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
            slot.name = "netP"
        }
    }
    if (!is.null(pairLR.use) & !is.null(signaling)) {
        stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
    }
    if (is.null(net)) {
        prob <- slot(object, "net")$prob
        pval <- slot(object, "net")$pval
        prob[pval > thresh] <- 0
        net <- reshape2::melt(prob, value.name = "prob")
        colnames(net)[1:3] <- c("source", "target", "interaction_name")
        pairLR = dplyr::select(object@LR$LRsig, c("interaction_name_2", 
                                                  "pathway_name", "ligand", "receptor", "annotation", 
                                                  "evidence"))
        idx <- match(net$interaction_name, rownames(pairLR))
        temp <- pairLR[idx, ]
        net <- cbind(net, temp)
    }
    if (!is.null(signaling)) {
        pairLR.use <- data.frame()
        for (i in 1:length(signaling)) {
            pairLR.use.i <- searchPair(signaling = signaling[i], 
                                       pairLR.use = object@LR$LRsig, key = "pathway_name", 
                                       matching.exact = T, pair.only = T)
            pairLR.use <- rbind(pairLR.use, pairLR.use.i)
        }
    }
    if (!is.null(pairLR.use)) {
        if ("interaction_name" %in% colnames(pairLR.use)) {
            net <- subset(net, interaction_name %in% pairLR.use$interaction_name)
        }
        else if ("pathway_name" %in% colnames(pairLR.use)) {
            net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
        }
    }
    if (slot.name == "netP") {
        net <- dplyr::select(net, c("source", "target", "pathway_name", 
                                    "prob"))
        net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
        net <- net %>% dplyr::group_by(source_target, pathway_name) %>% 
            dplyr::summarize(prob = sum(prob))
        a <- stringr::str_split(net$source_target, "sourceTotarget", 
                                simplify = T)
        net$source <- as.character(a[, 1])
        net$target <- as.character(a[, 2])
        net$ligand <- net$pathway_name
        net$receptor <- " "
    }
    if (!is.null(sources.use)) {
        if (is.numeric(sources.use)) {
            sources.use <- levels(object@idents)[sources.use]
        }
        net <- subset(net, source %in% sources.use)
    }
    else {
        sources.use <- levels(object@idents)
    }
    if (!is.null(targets.use)) {
        if (is.numeric(targets.use)) {
            targets.use <- levels(object@idents)[targets.use]
        }
        net <- subset(net, target %in% targets.use)
    }
    else {
        targets.use <- levels(object@idents)
    }
    df <- subset(net, prob > 0)
    if (nrow(df) == 0) {
        stop("No signaling links are inferred! ")
    }
    if (length(unique(net$ligand)) == 1) {
        message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
    }
    df$id <- 1:nrow(df)
    ligand.uni <- unique(df$ligand)
    for (i in 1:length(ligand.uni)) {
        df.i <- df[df$ligand == ligand.uni[i], ]
        source.uni <- unique(df.i$source)
        for (j in 1:length(source.uni)) {
            df.i.j <- df.i[df.i$source == source.uni[j], ]
            df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(" ", 
                                                             j - 1), collapse = ""))
            df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
        }
    }
    receptor.uni <- unique(df$receptor)
    for (i in 1:length(receptor.uni)) {
        df.i <- df[df$receptor == receptor.uni[i], ]
        target.uni <- unique(df.i$target)
        for (j in 1:length(target.uni)) {
            df.i.j <- df.i[df.i$target == target.uni[j], ]
            df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(" ", 
                                                                 j - 1), collapse = ""))
            df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
        }
    }
    cell.order.sources <- levels(object@idents)[levels(object@idents) %in% 
                                                    sources.use]
    cell.order.targets <- levels(object@idents)[levels(object@idents) %in% 
                                                    targets.use]
    df$source <- factor(df$source, levels = cell.order.sources)
    df$target <- factor(df$target, levels = cell.order.targets)
    df.ordered.source <- df[with(df, order(source, -prob)), ]
    df.ordered.target <- df[with(df, order(target, -prob)), ]
    order.source <- unique(df.ordered.source[, c("ligand", "source")])
    order.target <- unique(df.ordered.target[, c("receptor", 
                                                 "target")])
    order.sector <- c(order.source$ligand, order.target$receptor)
    if (is.null(color.use)) {
        color.use = scPalette(nlevels(object@idents))
        names(color.use) <- levels(object@idents)
        color.use <- color.use[levels(object@idents) %in% as.character(union(df$source, 
                                                                             df$target))]
    }
    else if (is.null(names(color.use))) {
        names(color.use) <- levels(object@idents)
        color.use <- color.use[levels(object@idents) %in% as.character(union(df$source, 
                                                                             df$target))]
    }
    edge.color <- color.use[as.character(df.ordered.source$source)]
    names(edge.color) <- as.character(df.ordered.source$source)
    grid.col.ligand <- color.use[as.character(order.source$source)]
    names(grid.col.ligand) <- as.character(order.source$source)
    grid.col.receptor <- color.use[as.character(order.target$target)]
    names(grid.col.receptor) <- as.character(order.target$target)
    grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
    names(grid.col) <- order.sector
    df.plot <- df.ordered.source[, c("ligand", "receptor", "prob")]
    if (directional == 2) {
        link.arr.type = "triangle"
    }
    else {
        link.arr.type = "big.arrow"
    }
    circos.clear()
    chordDiagram(df.plot, order = order.sector, col = edge.color, 
                 grid.col = grid.col, transparency = transparency, link.border = link.border, 
                 directional = directional, direction.type = c("diffHeight", 
                                                               "arrows"), link.arr.type = link.arr.type, annotationTrack = "grid", 
                 annotationTrackHeight = annotationTrackHeight, preAllocateTracks = list(track.height = max(strwidth(order.sector))), 
                 small.gap = small.gap, big.gap = big.gap, link.visible = link.visible, 
                 scale = scale, link.target.prop = link.target.prop, reduce = reduce, 
                 ...)
    circos.track(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        xplot = get.cell.meta.data("xplot")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", 
                    niceFacing = TRUE, adj = c(0, 0.5), cex = lab.cex)
    }, bg.border = NA)
    if (show.legend) {
        lgd <- ComplexHeatmap::Legend(at = names(color.use), 
                                      type = "grid", legend_gp = grid::gpar(fill = color.use), 
                                      title = "Cell State")
        ComplexHeatmap::draw(lgd, x = unit(1, "npc") - unit(legend.pos.x, 
                                                            "mm"), y = unit(legend.pos.y, "mm"), just = c("right", 
                                                                                                          "bottom"))
    }
    circos.clear()
    if (!is.null(title.name)) {
        text(-0, 1.02, title.name, cex = 1)
    }
    gg <- recordPlot()
    return(list(plot = gg, df = df.plot))
}
