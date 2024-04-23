
library(here)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(GGally)
library(sessioninfo)


################################################################################
##                                2. Comparisons
################################################################################

## Habenula results
load(here("processed-data/05_DEA/results_Substance_uncorr_vars_habenula.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_FirstHrIntakeSlope_habenula.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_TotalIntake_habenula.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_LastSessionIntake_habenula.Rdata"), verbose = TRUE)


## Amygdala results
load(here("processed-data/05_DEA/results_Substance_uncorr_vars_amygdala.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_FirstHrIntakeSlope_amygdala.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_TotalIntake_amygdala.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_LastSessionIntake_amygdala.Rdata"), verbose = TRUE)


## 2.1 Comparison of gene DE signal for substance and behavior in habenula

## DE t-stats of the genes for substance and behavioral covariates
t_stats_Substance_habenula <- results_Substance_uncorr_vars_habenula[[1]]
t_stats_FirstHrIntakeSlope_habenula <- results_FirstHrIntakeSlope_habenula[[1]]
t_stats_TotalIntake_habenula <- results_TotalIntake_habenula[[1]]
t_stats_LastSessionIntake_habenula <- results_LastSessionIntake_habenula[[1]]

## Df
df_habenula <- cbind(t_stats_Substance_habenula[,c("ensemblID", "Symbol", "logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_FirstHrIntakeSlope_habenula[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_TotalIntake_habenula[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_LastSessionIntake_habenula[,c("logFC", "t", "P.Value", "adj.P.Val")])

colnames(df_habenula)[3:18] <- paste(c("logFC", "t", "P.Value", "adj.P.Val"),
                               rep(c('Substance', 'FirstHrIntakeSlope', 'TotalIntake', 'LastSessionIntake'), c(4,4,4,4)),
                               sep='_')

## Pairwise correlations between t-stats
ggpairs(df_habenula, columns = c("t_Substance", "t_FirstHrIntakeSlope", "t_TotalIntake", "t_LastSessionIntake")) + theme_bw()
ggsave(here('plots/05_DEA/02_Comparisons/t_stats_pairs_habenula.pdf'))





## 2.2 Comparison of gene DE signal for substance and behavior in amygdala

t_stats_Substance_amygdala <- results_Substance_uncorr_vars_amygdala[[1]]
t_stats_FirstHrIntakeSlope_amygdala <- results_FirstHrIntakeSlope_amygdala[[1]]
t_stats_TotalIntake_amygdala <- results_TotalIntake_amygdala[[1]]
t_stats_LastSessionIntake_amygdala <- results_LastSessionIntake_amygdala[[1]]

df_amygdala <- cbind(t_stats_Substance_amygdala[,c("ensemblID", "Symbol", "logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_FirstHrIntakeSlope_amygdala[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_TotalIntake_amygdala[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_LastSessionIntake_amygdala[,c("logFC", "t", "P.Value", "adj.P.Val")])

colnames(df_amygdala)[3:18] <- paste(c("logFC", "t", "P.Value", "adj.P.Val"),
                                     rep(c('Substance', 'FirstHrIntakeSlope', 'TotalIntake', 'LastSessionIntake'), c(4,4,4,4)),
                                     sep='_')

## Pairwise correlations between t-stats
ggpairs(df_amygdala, columns = c("t_Substance", "t_FirstHrIntakeSlope", "t_TotalIntake", "t_LastSessionIntake")) + theme_bw()
ggsave(here('plots/05_DEA/02_Comparisons/t_stats_pairs_amygdala.pdf'))







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
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0    2023-10-26 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-11-06 [1] Bioconductor
# digest                 0.6.34    2024-01-11 [1] CRAN (R 4.3.1)
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# evaluate               0.23      2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.6    2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2024-02-17 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-30 [1] Bioconductor
# GGally               * 2.2.1     2024-02-14 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.0     2024-02-23 [1] CRAN (R 4.3.1)
# ggstats                0.6.0     2024-04-05 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# knitr                  1.45      2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1    2023-11-02 [1] Bioconductor
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-26 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.3.1)
# progress               1.2.3     2023-12-06 [1] CRAN (R 4.3.1)
# purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26      2024-03-05 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-26 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# SparseArray            1.2.4     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
# SummarizedExperiment * 1.32.0    2023-11-06 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.42      2024-02-08 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-26 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-26 [1] Bioconductor
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
