
library(here)
library(rtracklayer)
library(sessioninfo)


################################################################################
##                     3. Preparation of gene location file
################################################################################

## GRCh37 human genome assembly
gtf = import("/dcs05/lieber/marmaypag/legacySingleCell_Tran_Maynard_Neuron_2021_LIBD001/refdata-cellranger-hg19-3.0.0/genes/genes.gtf")
## GRanges object with 2490452 range

## Only genes
gtf = gtf[gtf$type == "gene"]
length(gtf)
# [1] 32738

## Valid chr: 1-22, sexual and mitochondrial
length(which(! seqnames(gtf) %in% c(1:22,"X","Y","MT")))
# 29
gtf <- gtf[seqnames(gtf) %in% c(1:22,"X","Y","MT"),]
unique(seqnames(gtf))
# [1] 1  2  3  4  5  6  7  X  8  9  10 11 12 13 14 15 16 17 18 20 19 Y  22 21 MT

## Create file
geneloc_file <- data.frame(ensembl = gtf$gene_id,
                           chr = seqnames(gtf),
                           start = start(gtf),
                           end = end(gtf),
                           strand = strand(gtf),
                           symbol = gtf$gene_name)

head(geneloc_file)
#            ensembl chr  start    end strand       symbol
# 1 ENSG00000243485   1  29554  31109      +   MIR1302-10
# 2 ENSG00000237613   1  34554  36081      -      FAM138A
# 3 ENSG00000186092   1  69091  70008      +        OR4F5
# 4 ENSG00000238009   1  89295 133566      - RP11-34P13.7
# 5 ENSG00000239945   1  89551  91105      - RP11-34P13.8
# 6 ENSG00000237683   1 134901 139379      -   AL627309.1

dim(geneloc_file)
# [1] 32709     6

## Save
write.table(geneloc_file,
            file= here("processed-data/07_MAGMA/geneloc/human_GRCh37_Ensembl.gene.loc"),
            sep="\t",row.names=F, col.names=F, quote=F)







## Reproducibility information

options(width = 120)
session_info()
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# biomaRt              * 2.56.1    2023-06-09 [2] Bioconductor
# Biostrings             2.68.1    2023-05-16 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
# cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# curl                   5.0.2     2023-08-14 [2] CRAN (R 4.3.1)
# DBI                    1.1.3     2022-06-18 [2] CRAN (R 4.3.1)
# dbplyr                 2.3.3     2023-07-07 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# digest                 0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.1)
# filelock               1.0.2     2018-10-05 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicAlignments      1.36.0    2023-04-25 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.3.1)
# httr                   1.4.7     2023-08-15 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# KEGGREST               1.40.0    2023-04-25 [2] Bioconductor
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
# prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.3.1)
# progress               1.2.2     2019-05-16 [2] CRAN (R 4.3.1)
# purrr                  1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# Rsamtools              2.16.0    2023-04-25 [2] Bioconductor
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
# rtracklayer          * 1.60.1    2023-08-15 [2] Bioconductor
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# stringi                1.7.12    2023-01-11 [2] CRAN (R 4.3.1)
# stringr                1.5.0     2022-12-02 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XML                    3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
# xml2                   1.3.5     2023-07-06 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.1)
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
#
# [1] /users/dgonzale/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

