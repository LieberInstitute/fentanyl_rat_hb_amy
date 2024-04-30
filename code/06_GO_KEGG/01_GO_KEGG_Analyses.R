
library(here)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Rn.eg.db)
library(cowplot)
library(ggplot2)
library(sessioninfo)


#######################  Functional Enrichment Analysis  #######################

load(here('processed-data/05_DEA/de_genes_Substance_habenula.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/de_genes_Substance_amygdala.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/results_Substance_uncorr_vars_amygdala.Rdata'), verbose = TRUE)


## Groups of DEGs

########################
##    Habenula DEGs
########################
## Up and down habenula DEGs from model with uncorrelated sample variables
up_hab <- de_genes_habenula[which(de_genes_habenula$logFC>0), ]
down_hab <- de_genes_habenula[which(de_genes_habenula$logFC<0), ]

########################
##    Amygdala DEGs
########################
## Up and down amygdala DEGs from model with uncorrelated sample variables
up_amy <- de_genes_amygdala[which(de_genes_amygdala$logFC>0), ]
down_amy <- de_genes_amygdala[which(de_genes_amygdala$logFC<0), ]



## Function to find enriched GO and KEGG terms

GO_KEGG<- function(sigGeneList, geneUniverse, name){

    ## GO terms
    ## Obtain biological processes
    goBP_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = org.Rn.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    ## Save
    if (!is.null(goBP_Adj)){
        p1 <- dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes")
    }


    ## Obtain molecular functions
    goMF_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = org.Rn.eg.db,
        ont = "MF",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    if (!is.null(goMF_Adj)){
        p2 <- dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function")
    }


    ## Obtain cellular components
    goCC_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = org.Rn.eg.db,
        ont = "CC",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    if (!is.null(goCC_Adj)){
        p3 <- dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components")
    }


    ## KEGG terms
    kegg_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichKEGG",
        organism = 'rat',
        universe = geneUniverse,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    if (!is.null(kegg_Adj)){
        p4 <- dotplot(kegg_Adj, title="KEGG Enrichment Analysis")
    }

    ## Plots
    plot_grid(p1, p2, p3, p4, ncol=2, align = 'vh')
    ggsave(paste("plots/06_GO_KEGG/GO_KEGG_", name, ".pdf", sep=""), height = 10, width = 14)

    ## Save results
    goList <- list(
        BP = goBP_Adj,
        MF = goMF_Adj,
        CC = goCC_Adj,
        KEGG = kegg_Adj
    )

    return(goList)
}


## 1. Analysis for all DEGs from each brain region

######################
#      Habenula
######################

## List of DEGs
sigGeneList <- list("All"= de_genes_habenula[which(!is.na(de_genes_habenula$EntrezID) & !de_genes_habenula$EntrezID=='NULL'), 'EntrezID'])
## Background genes (all genes assessed for DGE)
geneUniverse <- as.character(results_Substance_uncorr_vars_amygdala[[1]]$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']

goList_habenula_all_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_all_DEGs')
save(goList_habenula_all_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_all_DEGs.Rdata")

######################
#      Amygdala
######################

sigGeneList <- list("All"= de_genes_amygdala[which(!is.na(de_genes_amygdala$EntrezID) & !de_genes_amygdala$EntrezID=='NULL'), 'EntrezID'])

goList_amygdala_all_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_all_DEGs')
save(goList_amygdala_all_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_all_DEGs.Rdata")





## 2. Analysis for up- and down-regulated DEGs from each brain region

######################
#      Habenula
######################

## List of DEG sets
sigGeneList <- list("Up"=up_hab[which(!is.na(up_hab$EntrezID) & !up_hab$EntrezID=='NULL'), 'EntrezID'],
                    "Down"=down_hab[which(!is.na(down_hab$EntrezID) & !down_hab$EntrezID=='NULL'), 'EntrezID'])

goList_habenula_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_up_down_DEGs')
save(goList_habenula_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_up_down_DEGs.Rdata")

######################
#      Amygdala
######################

sigGeneList <- list("Up"=up_amy[which(!is.na(up_amy$EntrezID) & !up_amy$EntrezID=='NULL'), 'EntrezID'],
                    "Down"=down_amy[which(!is.na(down_amy$EntrezID) & !down_amy$EntrezID=='NULL'), 'EntrezID'])

goList_amygdala_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_up_down_DEGs')
save(goList_amygdala_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_up_down_DEGs.Rdata")







## Reproducibility information

options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2024-04-23
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version    date (UTC) lib source
# abind                    1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
# AnnotationDbi          * 1.64.1     2023-11-02 [1] Bioconductor
# AnnotationHub            3.10.0     2023-10-26 [1] Bioconductor
# ape                      5.7-1      2023-03-13 [1] CRAN (R 4.3.0)
# aplot                    0.2.2      2023-10-06 [1] CRAN (R 4.3.1)
# Biobase                * 2.62.0     2023-10-26 [1] Bioconductor
# BiocFileCache            2.10.1     2023-10-26 [1] Bioconductor
# BiocGenerics           * 0.48.1     2023-11-02 [1] Bioconductor
# BiocManager              1.30.22    2023-08-08 [1] CRAN (R 4.3.0)
# BiocParallel             1.36.0     2023-10-26 [1] Bioconductor
# BiocVersion              3.18.1     2023-11-18 [1] Bioconductor 3.18 (R 4.3.2)
# Biostrings               2.70.2     2024-01-30 [1] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
# cachem                   1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
# cli                      3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# clusterProfiler        * 4.10.0     2023-11-06 [1] Bioconductor
# codetools                0.2-19     2023-02-01 [1] CRAN (R 4.3.2)
# colorspace               2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# cowplot                  1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon                   1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.2.1      2024-03-01 [1] CRAN (R 4.3.1)
# data.table               1.15.2     2024-02-29 [1] CRAN (R 4.3.1)
# DBI                      1.2.2      2024-02-16 [1] CRAN (R 4.3.2)
# dbplyr                   2.4.0      2023-10-26 [1] CRAN (R 4.3.1)
# DelayedArray             0.28.0     2023-11-06 [1] Bioconductor
# digest                   0.6.34     2024-01-11 [1] CRAN (R 4.3.1)
# DOSE                     3.28.2     2023-12-12 [1] Bioconductor 3.18 (R 4.3.2)
# dplyr                    1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# ellipsis                 0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
# enrichplot               1.22.0     2023-11-06 [1] Bioconductor
# evaluate                 0.23       2023-11-01 [1] CRAN (R 4.3.1)
# fansi                    1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                  1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# fastmatch                1.1-4      2023-08-18 [1] CRAN (R 4.3.0)
# fgsea                    1.28.0     2023-10-26 [1] Bioconductor
# filelock                 1.0.3      2023-12-11 [1] CRAN (R 4.3.1)
# fs                       1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.38.6     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11     2024-02-17 [1] Bioconductor
# GenomicRanges          * 1.54.1     2023-10-30 [1] Bioconductor
# ggforce                  0.4.2      2024-02-19 [1] CRAN (R 4.3.1)
# ggfun                    0.1.4      2024-01-19 [1] CRAN (R 4.3.1)
# ggplot2                  3.5.0      2024-02-23 [1] CRAN (R 4.3.1)
# ggplotify                0.1.2      2023-08-09 [1] CRAN (R 4.3.0)
# ggraph                   2.2.0      2024-02-27 [1] CRAN (R 4.3.1)
# ggrepel                  0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# ggtree                   3.10.1     2024-02-27 [1] Bioconductor 3.18 (R 4.3.2)
# glue                     1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# GO.db                    3.18.0     2024-02-17 [1] Bioconductor
# GOSemSim                 2.28.1     2024-01-20 [1] Bioconductor 3.18 (R 4.3.2)
# graphlayouts             1.1.0      2024-01-19 [1] CRAN (R 4.3.1)
# gridExtra                2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gridGraphics             0.5-1      2020-12-13 [1] CRAN (R 4.3.0)
# gson                     0.1.0      2023-03-07 [1] CRAN (R 4.3.0)
# gtable                   0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# HDO.db                   0.99.1     2023-05-28 [1] Bioconductor
# here                   * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# htmltools                0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
# httpuv                   1.6.14     2024-01-26 [1] CRAN (R 4.3.1)
# httr                     1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
# igraph                   2.0.2      2024-02-17 [1] CRAN (R 4.3.1)
# interactiveDisplayBase   1.40.0     2023-10-26 [1] Bioconductor
# IRanges                * 2.36.0     2023-10-26 [1] Bioconductor
# jsonlite                 1.8.8      2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.42.0     2023-10-26 [1] Bioconductor
# knitr                    1.45       2023-10-30 [1] CRAN (R 4.3.1)
# labeling                 0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# later                    1.3.2      2023-12-06 [1] CRAN (R 4.3.1)
# lattice                  0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
# lazyeval                 0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# magrittr                 2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                     7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                   1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.14.0     2023-10-26 [1] Bioconductor
# matrixStats            * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# memoise                  2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12       2021-09-28 [1] CRAN (R 4.3.0)
# munsell                  0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# org.Rn.eg.db           * 3.18.0     2024-02-17 [1] Bioconductor
# patchwork                1.2.0      2024-01-08 [1] CRAN (R 4.3.1)
# pillar                   1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plyr                     1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                      0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
# polyclip                 1.10-6     2023-09-27 [1] CRAN (R 4.3.1)
# promises                 1.2.1      2023-08-10 [1] CRAN (R 4.3.0)
# purrr                    1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
# qvalue                   2.34.0     2023-10-26 [1] Bioconductor
# R6                       2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                    1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# reshape2                 1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# rlang                    1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown                2.26       2024-03-05 [1] CRAN (R 4.3.1)
# rprojroot                2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
# RSQLite                  2.3.5      2024-01-21 [1] CRAN (R 4.3.1)
# rstudioapi               0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays                 1.2.0      2023-10-26 [1] Bioconductor
# S4Vectors              * 0.40.2     2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                   1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# scatterpie               0.2.1      2023-06-07 [1] CRAN (R 4.3.0)
# sessioninfo            * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# shadowtext               0.1.3      2024-01-19 [1] CRAN (R 4.3.1)
# shiny                    1.8.0      2023-11-17 [1] CRAN (R 4.3.1)
# SparseArray              1.2.4      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# stringi                  1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr                  1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.32.0     2023-11-06 [1] Bioconductor
# tibble                   3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidygraph                1.3.1      2024-01-30 [1] CRAN (R 4.3.1)
# tidyr                    1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect               1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tidytree                 0.4.6      2023-12-12 [1] CRAN (R 4.3.1)
# treeio                   1.26.0     2023-11-06 [1] Bioconductor
# tweenr                   2.0.3      2024-02-26 [1] CRAN (R 4.3.1)
# utf8                     1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# viridis                  0.6.5      2024-01-29 [1] CRAN (R 4.3.1)
# viridisLite              0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                    3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# xfun                     0.42       2024-02-08 [1] CRAN (R 4.3.1)
# xtable                   1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.42.0     2023-10-26 [1] Bioconductor
# yaml                     2.3.8      2023-12-11 [1] CRAN (R 4.3.1)
# yulab.utils              0.1.4      2024-01-28 [1] CRAN (R 4.3.1)
# zlibbioc                 1.48.0     2023-10-26 [1] Bioconductor
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
