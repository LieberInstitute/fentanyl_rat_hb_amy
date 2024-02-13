
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
library(edgeR)
library(ggplot2)
library(cowplot)
library(sessioninfo)

#######################   Data Preparation   #######################

## Load RSE objects
load(here('raw-data/count_objects/rse_gene_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('raw-data/count_objects/rse_exon_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('raw-data/count_objects/rse_jx_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('raw-data/count_objects/rpkmCounts_Jlab_experiment_n33.rda'), verbose=TRUE)
load(here('raw-data/count_objects/rawCounts_Jlab_experiment_n33.rda'), verbose=TRUE)
## Load sample data
sample_data <- as.data.frame(read_excel("raw-data/Sample_data.xlsx"))



## 1. Data exploration & correction

## Verify sample data is the same in gene, exon and jx datasets
identical(colData(rse_gene), colData(rse_exon))
# [1] TRUE
identical(colData(rse_gene), colData(rse_jx))
# [1] TRUE


## Add sample data to colData of gene RSE

## Revome last (empty) line and column
sample_data <- sample_data[-nrow(sample_data),-ncol(sample_data)]
## Correct colnames in sample data
colnames(sample_data) <- str_replace_all(colnames(sample_data), c(" "="_"))
## Change conflicting column names
colnames(sample_data)[1] <- 'Sample_Num'
colnames(sample_data)[4] <- 'Num_Fentanyl_Sessions_six_hrs'
colnames(sample_data)[5] <- 'Total_Num_Fentanyl_Sessions'
colnames(sample_data)[6] <- 'RNA_concentration'
colnames(sample_data)[7] <- 'RIN_LIBD'
colnames(sample_data)[9] <- 'Total_RNA_amount'

## Add correct sample IDs
sample_data$SAMPLE_ID <- str_replace_all(sample_data$Tissue_Punch_Label, c(" "="_", "-"="_"))
## Correct Sample_Num for the extra sample
sample_data[which(sample_data$SAMPLE_ID=="33_S_Amyg_20"), "Sample_Num"]="33.0"
## Numeric data
sample_data$Sample_Num <- as.numeric(sample_data$Sample_Num)

## Add information of batch for RNA extraction
sample_data$Batch_RNA_extraction <- apply(sample_data, 1, function(x){if(x['Sample_Num'] %in% c(10, 14, 16, 29, 33)){"3"} else if (x['Brain_Region']=='habenula'){"1"} else {"2"}})

## Add information of batch for library preparation
sample_data$Batch_lib_prep <- apply(sample_data, 1, function(x){if(x['Sample_Num'] %in% c(26, 30, 33)){"3"} else if (x['Brain_Region']=='habenula'){"1"} else {"2"}})

## Add information of batch for sequencing
sample_data$Batch_seq <- apply(sample_data, 1, function(x){if(x['Brain_Region']=='habenula'){"1"} else {"2"}})


## Verify that samples are the same in sample_data and rse
setdiff(sample_data$SAMPLE_ID, rse_gene$SAMPLE_ID)
# character(0)
setdiff(rse_gene$SAMPLE_ID, sample_data$SAMPLE_ID)
# character(0)

## Merge sample data in colData of rse_gene
colData(rse_gene) <- merge(colData(rse_gene), sample_data, by='SAMPLE_ID')

## Add full sample data to colData of exon and jxn RSE objects
colData(rse_exon) <- colData(rse_gene)
colData(rse_jx) <- colData(rse_gene)


## Save data
save(rse_gene, file="processed-data/03_Data_preparation/rse_gene_sample_info.Rdata")
save(rse_exon, file="processed-data/03_Data_preparation/rse_exon_sample_info.Rdata")
save(rse_jx, file="processed-data/03_Data_preparation/rse_jx_sample_info.Rdata")





## 2. Data normalization

## Percentage of zeros in each dataset
length(which(assay(rse_gene)==0))*100/(dim(rse_gene)[1]*dim(rse_gene)[2])
# 33.91736
length(which(assay(rse_exon)==0))*100/(dim(rse_exon)[1]*dim(rse_exon)[2])
# 21.2531
length(which(assay(rse_jx)==0))*100/(dim(rse_jx)[1]*dim(rse_jx)[2])
# 89.52425

## Transform read counts to log2(CPM + 0.5) -> Norm distribution
## Genes
assays(rse_gene, withDimnames=FALSE)$logcounts <- edgeR::cpm(calcNormFactors(rse_gene, method = "TMM"), log = TRUE, prior.count = 0.5)

## Exons
assays(rse_exon, withDimnames=FALSE)$logcounts<- edgeR::cpm(calcNormFactors(rse_exon, method = "TMM"), log = TRUE, prior.count = 0.5)

## Junctions
## TMMwsp method for >80% of zeros
assays(rse_jx, withDimnames=FALSE)$logcounts<- edgeR::cpm(calcNormFactors(rse_jx, method = "TMMwsp"), log = TRUE, prior.count = 0.5)

## Transcripts: TODO
## Scale TPM (Transcripts per million) to log2(TPM + 0.5)
## assays(rse_tx)$logcounts<-log2(assays(rse_tx)$tpm + 0.5)


## Save rse objects with the assays of normalized counts
save(rse_gene, file="processed-data/03_Data_preparation/rse_gene_logcounts.Rdata")
save(rse_exon, file="processed-data/03_Data_preparation/rse_exon_logcounts.Rdata")
save(rse_jx, file="processed-data/03_Data_preparation/rse_jx_logcounts.Rdata")
##save(rse_tx, file="processed-data/03_Data_preparation/rse_tx_logcounts.Rdata")





## 3. Feature filtering

## Feature filtering based on counts

## Keep genes with at least k CPM in n samples & with a minimum total number of counts across all samples
## Add design matrix to define n
rse_gene_filt<-rse_gene[which(filterByExpr(assay(rse_gene),
                                           design=with(colData(rse_gene), model.matrix(~ Brain_Region + Substance)))),]
dim(rse_gene_filt)
# 16708    33
## Percentage of genes that were kept
dim(rse_gene_filt)[1]*100/dim(rse_gene)[1]
# 54.86668


## Filter exons
rse_exon_filt<-rse_exon[which(filterByExpr(assay(rse_exon),
                                           design=with(colData(rse_exon), model.matrix(~ Brain_Region + Substance)))),]
dim(rse_exon_filt)
# 182291     33
## Percentage of exons that were kept
dim(rse_exon_filt)[1]*100/dim(rse_exon)[1]
#  66.76984


## Filter junctions
rse_jx_filt<-rse_jx[which(filterByExpr(assay(rse_jx),
                                       design=with(colData(rse_jx), model.matrix(~ Brain_Region + Substance)))),]
dim(rse_jx_filt)
# 150593     33
## Percentage of jxns that were kept
dim(rse_jx_filt)[1]*100/dim(rse_jx)[1]
# 4.338896


## Filter TPM TODO
## Identify potential cutoffs
# seed <- 2
# expression_cutoff(assays(rse_tx)$tpm, seed = seed, k=2)
# The suggested expression cutoff is 0.28:
# percent_features_cut  samples_nonzero_cut
#                 0.29                 0.27
# cutoff<-0.28
## Transcripts that pass cutoff
# rse_tx_filt<-rse_tx[rowMeans(assays(rse_tx)$tpm) > cutoff,]
# dim(rse_tx_filt)
#


## Add complete QC and sample filtering info to rse objects with filtered features and lognorm counts

## Load rse object with correct QC info
load(here('processed-data/04_EDA/01_QCA/rse_gene_complete.Rdata'), verbose=TRUE)
# Loading objects:
#     rse_gene

colData(rse_gene_filt) <- colData(rse_gene_complete)
colData(rse_exon_filt) <- colData(rse_gene_complete)
colData(rse_jx_filt) <- colData(rse_gene_complete)
#colData(rse_tx_filt) <- colData(rse_gene)

## Save
save(rse_gene_filt, file = 'processed-data/03_Data_preparation/rse_gene_filt.Rdata')
save(rse_exon_filt, file = 'processed-data/03_Data_preparation/rse_exon_filt.Rdata')
save(rse_jx_filt, file = 'processed-data/03_Data_preparation/rse_jx_filt.Rdata')
# save(rse_tx_filt, file = 'processed-data/03_Data_preparation/rse_tx_filt.Rdata')





## 4. Visualization

## Plots of the distribution of genes' counts before and after normalization and filtering

## Raw counts
counts_data <- data.frame(counts=as.vector(assays(rse_gene)$counts))
p1 <- ggplot(counts_data, aes(x=counts)) +
    geom_histogram(colour="black", fill="lightgray") +
    theme_classic() +
    labs(x='read counts', y='Frecuency')
ggsave(filename = 'plots/03_Data_preparation/Hist_counts.pdf', height = 3, width = 4)

## Normalized counts
logcounts_data <- data.frame(logcounts=as.vector(assays(rse_gene)$logcounts))
p2 <- ggplot(logcounts_data, aes(x=logcounts)) +
    geom_histogram(aes(y=..density..), colour="darkgray", fill="lightgray") +
    theme_classic() +
    geom_density(fill="#69b3a2", alpha=0.3) +
    labs(x='log(CPM+0.5)', y='Frecuency')
ggsave(filename = 'plots/03_Data_preparation/Hist_logcounts.pdf', height = 3, width = 4)

## Normalized counts and filtered genes
filt_logcounts_data <- data.frame(logcounts=as.vector(assays(rse_gene_filt)$logcounts))
p3 <- ggplot(filt_logcounts_data, aes(x=logcounts)) +
    geom_histogram(aes(y=..density..), colour="darkgray", fill="lightgray") +
    theme_classic() +
    geom_density(fill="#69b3a2", alpha=0.3) +
    labs(x='log(CPM+0.5)', y='Frecuency')
ggsave(filename = 'plots/03_Data_preparation/Hist_filt_logcounts.pdf', height = 3, width = 4)

plot_grid(p1, p2, p3, ncol=3)
ggsave(filename = 'plots/03_Data_preparation/Hist_count_dist.pdf', height = 3, width = 9)







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
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
# date     2024-02-13
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# Biobase              * 2.61.0    2023-06-02 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.26.6    2023-07-02 [1] Bioconductor
# digest                 0.6.34    2024-01-11 [1] CRAN (R 4.3.1)
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 3.43.7    2023-06-21 [1] Bioconductor
# evaluate               0.23      2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData       1.2.10    2023-05-28 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-30 [1] Bioconductor
# ggplot2              * 3.4.4     2023-10-12 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# knitr                  1.45      2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.21-9    2023-10-01 [1] CRAN (R 4.3.2)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.57.6    2023-06-21 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.0)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-1.1   2023-09-18 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.25      2023-09-18 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.41      2023-11-01 [1] CRAN (R 4.3.1)
# XVector                0.41.1    2023-06-02 [1] Bioconductor
# zlibbioc               1.47.0    2023-05-20 [1] Bioconductor
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
