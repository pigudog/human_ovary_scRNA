####----scMetabolism_sinlge_cell_seq----####
options(
  showErrorCalls = TRUE,
  showWarnCalls = TRUE,
  warn = 1,
  warnPartialMatchAttr = TRUE,
  warnPartialMatchDollar = TRUE
)
install.packages(c("devtools", "data.table", "wesanderson", "Seurat", "devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
BiocManager::install("AUCell")
BiocManager::install("GSEABase",force = T)
BiocManager::install("GSVA",force = T)
devtools::install_github("YosefLab/VISION")
devtools::install_github("wu-yc/scMetabolism")

load(file = 'sce_immune_after_cluster.Rdata')
library(scMetabolism)
library(ggplot2)
library(rsvd)
countexp.Seurat <-sce
countexp.Seurat@meta.data$cell_type <- countexp.Seurat@active.ident
countexp.Seurat<-sc.metabolism.Seurat(obj = countexp.Seurat, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metabolism.matrix

p.dimplot<-DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Glycolysis / Gluconeogenesis", 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, size = 1)
ggsave(p.dimplot,filename=paste0("scMeatabolism",'/dimplot_immune_metabolism.png'),width = 12,height = 9)


input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:30]
p_dotplot <- DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "ovary_sample", 
                   norm = "y")
ggsave(p_dotplot,filename=paste0("scMeatabolism",'/dotplot_immune_metabolism.png'),width = 12,height = 9)


input.pathway <- c(rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[26],
                   rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[28])
p_boxplot <- BoxPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = "ovary_sample", ncol = 1)
ggsave(p_boxplot,filename=paste0("scMeatabolism",'/boxplot_immune_metabolism.png'),width = 12,height = 9)

sce_Metal_exp = countexp.Seurat
mscore_data = data.frame(t(sce_Metal_exp@assays[["METABOLISM"]][["score"]]),sce_Metal_exp$cell_type)
avg_sM=aggregate(mscore_data[,1:ncol(mscore_data)-1],list(mscore_data$sce_Metal_exp.cell_type),mean)
rownames(avg_sM) = avg_sM$Group.1
avg_sM=data.frame(t(avg_sM[,-1]))
avg_sM$KEGG = rownames(sce_Metal_exp@assays[["METABOLISM"]][["score"]])
rownames(avg_sM)=avg_sM$KEGG

c_k_l = c()
for(c in c(1:ncol(avg_sM))){
  c_k=avg_sM[order(avg_sM[,c]),]$KEGG[1:5]
  c_k_l=c(c_k_l,c_k)
}
c_k_l= unique(c_k_l)
c_k_d = avg_sM[avg_sM$KEGG %in%c_k_l,]

rownames(c_k_d) = c_k_d$KEGG
p_heatmap<-pheatmap::pheatmap(c_k_d[,-ncol(c_k_d)],show_colnames = T,scale='row')
ggsave(p_heatmap,filename=paste0("scMeatabolism",'/heatmap_immune_metabolism.png'),width = 12,height = 9)


# > sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.6 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] monocle3_1.3.1              sp_1.6-0                    magrittr_2.0.3              svglite_2.1.1              
# [5] ggalluvial_0.12.5           reticulate_1.28             clusterProfiler_4.6.2       wesanderson_0.3.6          
# [9] VISION_3.0.1                rsvd_1.0.5                  scMetabolism_0.2.1          Cairo_1.6-0                
# [13] ggsci_3.0.0                 ggsignif_0.6.4              ggpubr_0.6.0                harmony_0.1.1              
# [17] Rcpp_1.0.10                 DoubletFinder_2.0.3         hdf5r_1.3.8                 patchwork_1.1.2            
# [21] lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0               purrr_1.0.1                
# [25] readr_2.1.4                 tidyr_1.3.0                 tibble_3.2.1                tidyverse_2.0.0            
# [29] rhdf5_2.42.1                scran_1.26.2                scuttle_1.8.4               SingleCellExperiment_1.20.1
# [33] SummarizedExperiment_1.28.0 GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         MatrixGenerics_1.10.0      
# [37] matrixStats_0.63.0          genefilter_1.80.3           gplots_3.1.3                org.Hs.eg.db_3.16.0        
# [41] topGO_2.50.0                SparseM_1.81                GO.db_3.16.0                AnnotationDbi_1.60.2       
# [45] IRanges_2.32.0              S4Vectors_0.36.2            Biobase_2.58.0              graph_1.76.0               
# [49] BiocGenerics_0.44.0         Matrix_1.5-1                dplyr_1.1.2                 SeuratObject_4.1.3         
# [53] Seurat_4.3.0                plotly_4.10.1               ggplot2_3.4.2              
# 
# loaded via a namespace (and not attached):
#   [1] ica_1.0-3                 ps_1.7.5                  foreach_1.5.2             lmtest_0.9-40             rprojroot_2.0.3          
# [6] crayon_1.5.2              MASS_7.3-58.1             rhdf5filters_1.10.1       nlme_3.1-160              backports_1.4.1          
# [11] GOSemSim_2.24.0           rlang_1.1.1               argparse_2.2.2            HDO.db_0.99.1             XVector_0.38.0           
# [16] ROCR_1.0-11               loe_1.1                   irlba_2.3.5.1             nloptr_2.0.3              callr_3.7.3              
# [21] limma_3.54.2              BiocParallel_1.32.6       bit64_4.0.5               glue_1.6.2                pheatmap_1.0.12          
# [26] sctransform_0.3.5         parallel_4.2.2            processx_3.8.1            spatstat.sparse_3.0-1     DOSE_3.24.2              
# [31] spatstat.geom_3.2-1       tidyselect_1.2.0          usethis_2.1.6             fitdistrplus_1.1-11       XML_3.99-0.14            
# [36] zoo_1.8-12                wordspace_0.2-8           xtable_1.8-4              phyclust_0.1-33           cli_3.6.1                
# [41] zlibbioc_1.44.0           plumber_1.2.1             swagger_3.33.1            rstudioapi_0.14           miniUI_0.1.1.1           
# [46] rjags_4-14                parallelDist_0.2.6        fastmatch_1.1-3           pbmcapply_1.5.1           lambda.r_1.2.4           
# [51] treeio_1.22.0             shiny_1.7.4               BiocSingular_1.14.0       xfun_0.39                 gson_0.1.0               
# [56] pkgbuild_1.4.0            cluster_2.1.4             caTools_1.18.2            tidygraph_1.2.3           KEGGREST_1.38.0          
# [61] ggrepel_0.9.3             logging_0.10-108          ape_5.7-1                 listenv_0.9.0             Biostrings_2.66.0        
# [66] png_0.1-8                 permute_0.9-7             future_1.32.0-9105        withr_2.5.0               ggforce_0.4.1            
# [71] bitops_1.0-7              plyr_1.8.8                sparsesvd_0.2-2           dqrng_0.3.0               coda_0.19-4              
# [76] pillar_1.9.0              RcppParallel_5.1.7        cachem_1.0.8              multcomp_1.4-23           fs_1.6.2                 
# [81] DelayedMatrixStats_1.20.0 vctrs_0.6.2               ellipsis_0.3.2            generics_0.1.3            devtools_2.4.5           
# [86] tools_4.2.2               tweenr_2.0.2              munsell_0.5.0             fgsea_1.24.0              DelayedArray_0.24.0      
# [91] fastmap_1.1.1             compiler_4.2.2            pkgload_1.3.2             abind_1.4-5               httpuv_1.6.11            
# [96] sessioninfo_1.2.2         GenomeInfoDbData_1.2.9    gridExtra_2.3             edgeR_3.40.2              lattice_0.20-45          
# [101] deldir_1.0-9              utf8_1.2.3                later_1.3.1               jsonlite_1.8.4            scales_1.2.1             
# [106] ScaledMatrix_1.6.0        tidytree_0.4.2            pbapply_1.7-0             carData_3.0-5             sparseMatrixStats_1.10.0 
# [111] lazyeval_0.2.2            promises_1.2.0.1          car_3.1-2                 doParallel_1.0.17         goftest_1.2-3            
# [116] spatstat.utils_3.0-3      sandwich_3.0-2            cowplot_1.1.1             textshaping_0.3.6         statmod_1.5.0            
# [121] Rtsne_0.16                downloader_0.4            uwot_0.1.14               igraph_1.4.2              survival_3.4-0           
# [126] systemfonts_1.0.4         htmltools_0.5.5           memoise_2.0.1             profvis_0.3.8             modeltools_0.2-23        
# [131] locfit_1.5-9.7            graphlayouts_1.0.0        quadprog_1.5-8            viridisLite_0.4.2         digest_0.6.31            
# [136] mime_0.12                 futile.options_1.0.1      RSQLite_2.3.1             yulab.utils_0.0.6         future.apply_1.11.0      
# [141] remotes_2.4.2             data.table_1.14.8         urlchecker_1.0.1          blob_1.2.4                vegan_2.6-4              
# [146] futile.logger_1.4.3       ragg_1.2.5                labeling_0.4.2            fastICA_1.2-3             splines_4.2.2            
# [151] webutils_1.1              Rhdf5lib_1.20.0           iotools_0.3-2             RCurl_1.98-1.12           broom_1.0.4              
# [156] hms_1.1.3                 colorspace_2.1-0          BiocManager_1.30.20       aplot_0.1.10              libcoin_1.0-9            
# [161] coin_1.4-2                mclust_6.0.0              RANN_2.6.1                mvtnorm_1.1-3             enrichplot_1.18.4        
# [166] fansi_1.0.4               tzdb_0.3.0                parallelly_1.35.0         R6_2.5.1                  grid_4.2.2               
# [171] ggridges_0.5.4            lifecycle_1.0.3           formatR_1.14              bluster_1.8.0             curl_5.0.0               
# [176] minqa_1.2.5               leiden_0.4.3              phangorn_2.11.1           fastcluster_1.2.3         qvalue_2.30.0            
# [181] desc_1.4.2                RcppAnnoy_0.0.20          TH.data_1.1-2             RColorBrewer_1.1-3        iterators_1.0.14         
# [186] spatstat.explore_3.2-1    htmlwidgets_1.6.2         beachmat_2.14.2           polyclip_1.10-4           shadowtext_0.1.2         
# [191] terra_1.7-29              gridGraphics_0.5-1        timechange_0.2.0          mgcv_1.8-41               globals_0.16.2           
# [196] spatstat.random_3.1-5     progressr_0.13.0          codetools_0.2-18          metapod_1.6.0             gtools_3.9.4             
# [201] prettyunits_1.1.1         gtable_0.3.3              DBI_1.1.3                 ggfun_0.0.9               tensor_1.5               
# [206] httr_1.4.6                KernSmooth_2.23-20        stringi_1.7.12            farver_2.1.1              reshape2_1.4.4           
# [211] infercnv_1.14.2           viridis_0.6.3             annotate_1.76.0           ggtree_3.6.2              boot_1.3-28              
# [216] BiocNeighbors_1.16.0      lme4_1.1-33               ggplotify_0.1.0           scattermore_1.1           bit_4.0.5                
# [221] scatterpie_0.1.9          spatstat.data_3.0-1       ggraph_2.1.0              pkgconfig_2.0.3           rstatix_0.7.2            
# [226] knitr_1.42     