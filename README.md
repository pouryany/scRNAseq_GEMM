## Overview

This repository contains the codes for generating the figures and
results presented in **Tumor editing suppresses innate and adaptive
antitumor immunity and is reversed by inhibiting DNA methylation**.

please contact `pnaderiy [at] bidmc [dot] harvard [dot] edu` for any
questions you may have about the codes.

## Citations

Zhang, Y., Naderi Yeganeh, P. et al. **Tumor editing suppresses innate
and adaptive antitumor immunity and is reversed by inhibiting DNA
methylation** *Nature Immunology* (2024) doi:
10.1038/s41590-024-01932-8[link](https://www.nature.com/articles/s41590-024-01932-8).

## Code instructions

First, you need to download the filtered scRNA-Seq dataset from NCBI GEO
(GSE234917) and place it under a folder named
`filtered_feature_bc_matrix`.

All codes relating to scRNA-Seq and bulk RNA-Seq are placed in the
`codes` directory.

Raw RNA-Seq data can be accessed from NCBI GEO (GSE212029)

<table>
<colgroup>
<col style="width: 11%" />
<col style="width: 44%" />
<col style="width: 44%" />
</colgroup>
<thead>
<tr class="header">
<th>File name</th>
<th>Description</th>
<th>outputs/Figures</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>01_import_sc.R</code></td>
<td>Preprocessing raw scRNA-Seq reads from 10X assays.</td>
<td>processed RDS files</td>
</tr>
<tr class="even">
<td><code>02_01_clustering.R</code></td>
<td>Normalizing and clustering processed scRNA-Seq data.</td>
<td>processed RDS files</td>
</tr>
<tr class="odd">
<td><code>02_02_scPlots.R</code></td>
<td>Plotting clusters and cellular subclasses.</td>
<td>Figures 1B, 1C</td>
</tr>
<tr class="even">
<td><code>03_01_tumor_plots.R</code></td>
<td>Plotting tumor clusters, associated markers, and cellular
proportions.</td>
<td>Figures 1D, 1E, and 1G; Extended Data Fig 2</td>
</tr>
<tr class="odd">
<td><code>03_02_tumor_markers.R</code></td>
<td>Defining tumor markers and differentially expressed genes in
specific subclusters between early and late tumors.</td>
<td>Marker genes for epithelial cell types and DE genes.</td>
</tr>
<tr class="even">
<td><code>03_03_cytoTracing.R</code></td>
<td>CytoTrace algorithm generates cell differentiation status.</td>
<td>Figure 1G</td>
</tr>
<tr class="odd">
<td><code>03_04_dirchlet_tumor.R</code></td>
<td>Dirichlet regression determines changes in cellular proportions
between late and early tumors.</td>
<td>Cell proportion table and associated statistics.</td>
</tr>
<tr class="even">
<td><code>03_05_DEG_late_early_tumors.R</code></td>
<td>Differentially expressed genes in pooled tumor subgroups between
late and early samples.</td>
<td>Figures 2A-B</td>
</tr>
<tr class="odd">
<td><code>04_01_Tcell.R</code></td>
<td>Plots relating to T and NK cells, along with associated markers.
Cell proportion analysis of T and NK cells between early and late
tumors.</td>
<td>Figures 3A-B</td>
</tr>
<tr class="even">
<td><code>04_02_immune.R</code></td>
<td>Generating cell-type markers of immune cells. Plotting key markers
of immune cells</td>
<td>Extended Data Fig. 3A, marker genes for immune cells.</td>
</tr>
<tr class="odd">
<td><code>05_00_cell_cell_communication_early.R</code></td>
<td>CellChat algorithm defines significant communication events in early
tumors.</td>
<td></td>
</tr>
<tr class="even">
<td><code>05_00_cell_cell_communication_late.R</code></td>
<td>CellChat algorithm defines significant communication events in late
tumors.</td>
<td></td>
</tr>
<tr class="odd">
<td><code>05_01_cell_cell_communication_compare.R</code></td>
<td>Comparison of cell-cell communication events between late and early
tumors.</td>
<td>Figures 4A-B and Extended Data Fig. 5</td>
</tr>
<tr class="even">
<td><code>05_02_cell_cell_communication_secondary.R</code></td>
<td>Comparison of cell-cell communication events between late and early
tumors.</td>
<td>Figure 4</td>
</tr>
<tr class="odd">
<td><code>06_01_Bulk_comparison_tumor.R</code></td>
<td>Comparison of single-cell and bulk RNA-Seq datasets.</td>
<td>Figure 6 and Extended Data Fig. 7F</td>
</tr>
<tr class="even">
<td><code>07_updated_figures.R</code></td>
<td>Additional plots and analyses relating to immune cell
populations.</td>
<td>Extended Data Fig. 3C</td>
</tr>
<tr class="odd">
<td><code>08_vln_significant.R</code></td>
<td>Violin plots for functional markers across immune cell
clusters.</td>
<td>Figure 3 and 4D and Extended Data Fig. 3D</td>
</tr>
<tr class="even">
<td><code>09_Atp_genes.R</code></td>
<td>Heatmaps of ATP- and metabolism-associated genes and pathways in
tumors.</td>
<td>Extended Data Fig. 1E</td>
</tr>
<tr class="odd">
<td><code>10_ISG.R</code></td>
<td>Defining interferon stimulated genes among differentially expressed
genes in late versus early tumors.</td>
<td></td>
</tr>
<tr class="even">
<td><code>11_01_Tcell_specific_markers.R</code></td>
<td>Defining top markers of tumor infiltrating lymphocytes.</td>
<td>Extended Data Fig. 4</td>
</tr>
<tr class="odd">
<td><code>11_02_Checkpointmarkers.R</code></td>
<td>Expression of analysis of genes encoding for checkpoint ligand and
receptors.</td>
<td>Supplementary Figures</td>
</tr>
<tr class="even">
<td><code>11_03_exhausted.R</code></td>
<td>Defining the differentially expressed genes between CD8+ Ex vs
progenitor Ex.</td>
<td></td>
</tr>
<tr class="odd">
<td><code>11_04_Gasdermins.R</code></td>
<td>Differential expression analysis of Gasdermin gene family.</td>
<td></td>
</tr>
<tr class="even">
<td><code>11_05_TRM.R</code></td>
<td>Proportions of genes in interest expressed in TRM cell
population.</td>
<td></td>
</tr>
<tr class="odd">
<td><code>11_07_gdTcells.R</code></td>
<td>Expression of gdT cell related genes in the TRM population.</td>
<td></td>
</tr>
<tr class="even">
<td><code>11_09_gd_allT.R</code></td>
<td>Expression of gdT cell related genes in CD8+ population.</td>
<td></td>
</tr>
<tr class="odd">
<td><code>12_Saeki_data.R</code></td>
<td>Comparison of results with normal mammary epithelial scRNA-Seq.</td>
<td></td>
</tr>
<tr class="even">
<td><code>13_LiLab_data.R</code></td>
<td>Comparison of results with ILTCK scRNA-Seq.</td>
<td></td>
</tr>
<tr class="odd">
<td><code>14_Shiny_app.R</code></td>
<td>Building Shiny App.</td>
<td></td>
</tr>
<tr class="even">
<td><code>15_DESeq2_4T1.R</code></td>
<td>Bulk RNA-Seq Data Analysis.</td>
<td></td>
</tr>
<tr class="odd">
<td><code>utils.R</code></td>
<td>Predefined functions for plotting and analysis.</td>
<td></td>
</tr>
</tbody>
</table>

## Data instructions

Supplementary files that are necessary for reproducing the study are
provided in the `data/` directory

<table>
<colgroup>
<col style="width: 11%" />
<col style="width: 44%" />
<col style="width: 44%" />
</colgroup>
<thead>
<tr class="header">
<th>File name</th>
<th>Description</th>
<th>Reference/Resource</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>20220913071924_GeneSearchResults.txt</code></td>
<td>Interferon stimulated genes.</td>
<td>The Interferome database</td>
</tr>
<tr class="even">
<td><code>enrichment_neg_tumor_late_vs_early.csv</code></td>
<td>GO enrichment results for genes downregulated in late tumors
compared to early tumors</td>
<td>This study</td>
</tr>
<tr class="odd">
<td><code>Figure2B_table.csv</code></td>
<td>Top GO pathways, manually grouped.</td>
<td>This study</td>
</tr>
<tr class="even">
<td><code>RNA_Seq_B16_raw.csv</code></td>
<td>Gene expression profiles of decitabine-treated and control B16 cell
lines.</td>
<td>This study</td>
</tr>
<tr class="odd">
<td><code>RNA_Seq_4T1_raw.csv</code></td>
<td>Gene expression profiles of decitabine-treated and control 4T1 cell
lines.</td>
<td>This study</td>
</tr>
<tr class="even">
<td><code>Table5_DEGenes_4T1.csv</code></td>
<td>Differentially expressed genes between decitabine-treated and
control 4T1 cell lines.</td>
<td>This study</td>
</tr>
<tr class="odd">
<td><code>Tumor_late_early--clustered default  node.csv</code></td>
<td>GO clusters generated based on EnrichmentMap algorithm</td>
<td></td>
</tr>
<tr class="even">
<td><code>Cytoscape</code></td>
<td>Directory containing cytoscape sessions for analyzing the enrichment
terms</td>
<td></td>
</tr>
</tbody>
</table>

## SessionInfo

    R Under development (unstable) (2021-12-03 r81290)
    Platform: x86_64-apple-darwin17.0 (64-bit)
    Running under: macOS Monterey 12.1

    Matrix products: default
    LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    attached base packages:
    [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] CytoTRACE_0.3.3        pals_1.7               cartography_3.0.1      rcartocolor_2.0.0      palettetown_0.1.1      wesanderson_0.3.6     
     [7] ComplexHeatmap_2.11.0  CellChat_1.5.0         igraph_1.2.11          clusterProfiler_4.3.1  org.Mm.eg.db_3.14.0    AnnotationDbi_1.57.1  
    [13] IRanges_2.29.1         S4Vectors_0.33.9       Biobase_2.55.0         BiocGenerics_0.41.2    enrichplot_1.15.2      VennDiagram_1.7.1     
    [19] futile.logger_1.4.3    tidyr_1.1.4            tibble_3.1.7           DirichletReg_0.7-1     Formula_1.2-4          rlang_1.0.4           
    [25] viridis_0.6.2          viridisLite_0.4.0      dplyr_1.0.7            ggthemes_4.2.4         forcats_0.5.1          EnhancedVolcano_1.13.2
    [31] ggrepel_0.9.1          patchwork_1.1.1        RColorBrewer_1.1-2     ggplot2_3.3.5          stringr_1.4.0          harmony_0.1.0         
    [37] Rcpp_1.0.7             SeuratObject_4.0.4     Seurat_4.1.0          

    loaded via a namespace (and not attached):
      [1] scattermore_0.8        maxLik_1.5-2           coda_0.19-4            pkgmaker_0.32.2        bit64_4.0.5            irlba_2.3.5           
      [7] data.table_1.14.2      rpart_4.1-15           KEGGREST_1.35.0        RCurl_1.98-1.5         doParallel_1.0.16      generics_0.1.1        
     [13] cowplot_1.1.1          lambda.r_1.2.4         RSQLite_2.2.9          shadowtext_0.1.0       RANN_2.6.1             future_1.26.1         
     [19] bit_4.0.4              spatstat.data_2.1-2    httpuv_1.6.4           assertthat_0.2.1       promises_1.2.0.1       fansi_0.5.0           
     [25] DBI_1.1.2              htmlwidgets_1.5.4      spatstat.geom_2.3-2    purrr_0.3.4            ellipsis_0.3.2         RSpectra_0.16-0       
     [31] ggpubr_0.4.0           backports_1.4.1        gridBase_0.4-7         deldir_1.0-6           vctrs_0.4.1            ggalluvial_0.12.3     
     [37] here_1.0.1             ROCR_1.0-11            abind_1.4-5            cachem_1.0.6           withr_2.4.3            ggforce_0.3.3         
     [43] sctransform_0.3.3      sna_2.7                treeio_1.19.1          goftest_1.2-3          svglite_2.1.0          cluster_2.1.2         
     [49] DOSE_3.21.2            ape_5.6-1              lazyeval_0.2.2         crayon_1.4.2           pkgconfig_2.0.3        tweenr_1.0.2          
     [55] GenomeInfoDb_1.31.1    nlme_3.1-153           globals_0.15.1         lifecycle_1.0.1        miniUI_0.1.1.1         sandwich_3.0-1        
     [61] downloader_0.4         registry_0.5-1         dichromat_2.0-0.1      rprojroot_2.0.2        polyclip_1.10-0        matrixStats_0.61.0    
     [67] lmtest_0.9-39          rngtools_1.5.2         Matrix_1.4-0           aplot_0.1.1            carData_3.0-4          zoo_1.8-9             
     [73] ggridges_0.5.3         GlobalOptions_0.1.2    png_0.1-7              rjson_0.2.21           bitops_1.0-7           KernSmooth_2.23-20    
     [79] Biostrings_2.63.1      blob_1.2.2             shape_1.4.6            qvalue_2.27.0          parallelly_1.30.0      spatstat.random_2.1-0 
     [85] rstatix_0.7.0          gridGraphics_0.5-1     ggsignif_0.6.3         scales_1.1.1           memoise_2.0.1          magrittr_2.0.1        
     [91] plyr_1.8.6             ica_1.0-2              zlibbioc_1.41.0        compiler_4.2.0         scatterpie_0.1.7       miscTools_0.6-26      
     [97] clue_0.3-60            fitdistrplus_1.1-6     cli_3.3.0              XVector_0.35.0         listenv_0.8.0          pbapply_1.5-0         
    [103] formatR_1.11           MASS_7.3-54            mgcv_1.8-38            tidyselect_1.1.1       stringi_1.7.6          GOSemSim_2.21.0       
    [109] fastmatch_1.1-3        tools_4.2.0            future.apply_1.8.1     parallel_4.2.0         circlize_0.4.14        rstudioapi_0.13       
    [115] foreach_1.5.1          gridExtra_2.3          farver_2.1.0           Rtsne_0.15             ggraph_2.0.5           rgeos_0.5-9           
    [121] digest_0.6.29          FNN_1.1.3              shiny_1.7.1            car_3.0-12             broom_0.7.11           later_1.3.0           
    [127] RcppAnnoy_0.0.19       httr_1.4.2             colorspace_2.0-2       tensor_1.5             reticulate_1.24        splines_4.2.0         
    [133] uwot_0.1.14            yulab.utils_0.0.4      tidytree_0.3.6         spatstat.utils_2.3-0   graphlayouts_0.8.0     sp_1.4-6              
    [139] mapproj_1.2.9          ggplotify_0.1.0        systemfonts_1.0.4      plotly_4.10.0          xtable_1.8-4           jsonlite_1.7.2        
    [145] ggtree_3.3.1           futile.options_1.0.1   tidygraph_1.2.0        ggfun_0.0.4            R6_2.5.1               pillar_1.8.0          
    [151] htmltools_0.5.2        mime_0.12              NMF_0.24.0             glue_1.6.2             fastmap_1.1.0          BiocParallel_1.29.10  
    [157] codetools_0.2-18       maps_3.4.0             fgsea_1.21.0           utf8_1.2.2             lattice_0.20-45        spatstat.sparse_2.1-0 
    [163] network_1.17.2         leiden_0.3.9           GO.db_3.14.0           survival_3.2-13        statnet.common_4.6.0   munsell_0.5.0         
    [169] DO.db_2.9              GetoptLong_1.0.5       GenomeInfoDbData_1.2.7 iterators_1.0.13       reshape2_1.4.4         gtable_0.3.0          
    [175] spatstat.core_2.4-0   
