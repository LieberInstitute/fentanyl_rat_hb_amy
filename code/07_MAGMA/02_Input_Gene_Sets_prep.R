
library(here)
library(SummarizedExperiment)
library(biomaRt)
library(sessioninfo)


################################################################################
##                       2. Input DEG sets preparation
################################################################################

load(here('processed-data/05_DEA/de_genes_Substance_habenula.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/de_genes_Substance_amygdala.Rdata'), verbose = TRUE)


## Function to get orthologs genes of rat genes in human

human_orthologs <- function(de_genes_list){

    ## Mart
    mart <- useMart('ENSEMBL_MART_ENSEMBL', dataset='rnorvegicus_gene_ensembl')
    human_rat_ids <- getBM(values  = de_genes_list,
                           mart = mart,
                           attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
                           filters = "ensembl_gene_id")
    return(human_rat_ids)
}

load(here( here("dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data", "13_MAGMA", "gene_sets", "markerSets_broad_combo_ENSEMBL_FDR05.txt")), verbose=TRUE)


###################
##    Habenula
###################

## Get human orthologs of rat DEGs:

###################  All DEGs  ###################
DEGs_habenula_human_orthologs <- human_orthologs(de_genes_list = de_genes_habenula$ensemblID)[,'hsapiens_homolog_ensembl_gene']
DEGs_habenula_human_orthologs <- DEGs_habenula_human_orthologs[DEGs_habenula_human_orthologs!='']
which(duplicated(DEGs_habenula_human_orthologs))
# integer(0)

###############  Upregulated DEGs  ###############
up_de_genes_habenula <- subset(de_genes_habenula, logFC>0)
up_DEGs_habenula_human_orthologs <- human_orthologs(de_genes_list = up_de_genes_habenula$ensemblID)[,'hsapiens_homolog_ensembl_gene']
## Remove empty entries
up_DEGs_habenula_human_orthologs <- up_DEGs_habenula_human_orthologs[up_DEGs_habenula_human_orthologs!='']
## Unique genes only
up_DEGs_habenula_human_orthologs <- up_DEGs_habenula_human_orthologs[-which(duplicated(up_DEGs_habenula_human_orthologs))]
which(duplicated(up_DEGs_habenula_human_orthologs))
# integer(0)

##############  Downregulated DEGs  ##############
down_de_genes_habenula <- subset(de_genes_habenula, logFC<0)
down_DEGs_habenula_human_orthologs <- human_orthologs(de_genes_list = down_de_genes_habenula$ensemblID)[,'hsapiens_homolog_ensembl_gene']
down_DEGs_habenula_human_orthologs <- down_DEGs_habenula_human_orthologs[down_DEGs_habenula_human_orthologs!='']
which(duplicated(down_DEGs_habenula_human_orthologs))
# integer(0)



###################
##    Amygdala
###################

###################  All DEGs  ###################
DEGs_amygdala_human_orthologs <- human_orthologs(de_genes_list = de_genes_amygdala$ensemblID)[,'hsapiens_homolog_ensembl_gene']
DEGs_amygdala_human_orthologs <- DEGs_amygdala_human_orthologs[DEGs_amygdala_human_orthologs!='']
DEGs_amygdala_human_orthologs <- DEGs_amygdala_human_orthologs[-which(duplicated(DEGs_amygdala_human_orthologs))]
which(duplicated(DEGs_amygdala_human_orthologs))
# integer(0)

###############  Upregulated DEGs  ###############
up_de_genes_amygdala <- subset(de_genes_amygdala, logFC>0)
up_DEGs_amygdala_human_orthologs <- human_orthologs(de_genes_list = up_de_genes_amygdala$ensemblID)[,'hsapiens_homolog_ensembl_gene']
up_DEGs_amygdala_human_orthologs <- up_DEGs_amygdala_human_orthologs[up_DEGs_amygdala_human_orthologs!='']
up_DEGs_amygdala_human_orthologs <- up_DEGs_amygdala_human_orthologs[-which(duplicated(up_DEGs_amygdala_human_orthologs))]
which(duplicated(up_DEGs_amygdala_human_orthologs))

##############  Downregulated DEGs  ##############
down_de_genes_amygdala <- subset(de_genes_amygdala, logFC<0)
down_DEGs_amygdala_human_orthologs <- human_orthologs(de_genes_list = down_de_genes_amygdala$ensemblID)[,'hsapiens_homolog_ensembl_gene']
down_DEGs_amygdala_human_orthologs <- down_DEGs_amygdala_human_orthologs[down_DEGs_amygdala_human_orthologs!='']
down_DEGs_amygdala_human_orthologs <- down_DEGs_amygdala_human_orthologs[-which(duplicated(down_DEGs_amygdala_human_orthologs))]
which(duplicated(down_DEGs_amygdala_human_orthologs))


## Bind
gene_sets <- as.data.frame(rbind(cbind(DEGs_habenula_human_orthologs, 'All_DEGs_habenula'),
                                 cbind(up_DEGs_habenula_human_orthologs, 'Up_DEGs_habenula'),
                                 cbind(down_DEGs_habenula_human_orthologs, 'Down_DEGs_habenula'),
                                 cbind(DEGs_amygdala_human_orthologs, 'All_DEGs_amygdala'),
                                 cbind(up_DEGs_amygdala_human_orthologs, 'Up_DEGs_amygdala'),
                                 cbind(down_DEGs_amygdala_human_orthologs, 'Down_DEGs_amygdala')))
colnames(gene_sets) <- c('Gene', 'Set')
table(gene_sets$Set)
# All_DEGs_amygdala  All_DEGs_habenula Down_DEGs_amygdala Down_DEGs_habenula
#             2765                424                914                175
# Up_DEGs_amygdala   Up_DEGs_habenula
#             1853                250

save(gene_sets, file=here('processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt'))







## Reproducibility information

options(width = 120)
session_info()
# Biostrings             2.68.1    2023-05-16 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
# cachem                 1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
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
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# RSQLite                2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
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
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
#
# [1] /users/dgonzale/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

