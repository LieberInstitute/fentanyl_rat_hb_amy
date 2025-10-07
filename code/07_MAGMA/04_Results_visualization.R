
library(here)
library(berryFunctions)
library(ComplexHeatmap)
library(sessioninfo)


################################################################################
##                       4. MAGMA results visualization
################################################################################

## Load outputs
SCZ_out <- read.table(here('processed-data/07_MAGMA/Output/SCZ/SCZ_MAGMA.gsa.out'))[-1,]
MDD_2019_out <- read.table(here('processed-data/07_MAGMA/Output/MDD_2019/MDD_2019_MAGMA.gsa.out'))[-1,]
Panic_out <- read.table(here('processed-data/07_MAGMA/Output/Panic/Panic_MAGMA.gsa.out'))[-1,]
SUD_out <- read.table(here('processed-data/07_MAGMA/Output/SUD/SUD_MAGMA.gsa.out'))[-1,]
MDD_out <- read.table(here('processed-data/07_MAGMA/Output/MDD/MDD_MAGMA.gsa.out'))[-1,]
OUD_out <- read.table(here('processed-data/07_MAGMA/Output/OUD/OUD_MAGMA.gsa.out'))[-1,]

## Colnames
colnames(SCZ_out) <- colnames(MDD_2019_out) <- colnames(Panic_out) <- colnames(SUD_out) <- colnames(MDD_out) <- colnames(OUD_out) <-
    c("VARIABLE", "TYPE", "NGENES" , "BETA",  "BETA_STD",  "SE",  "P")

## Missing info for OUD in Down_DEGs_habenula
OUD_out <- berryFunctions::insertRows(OUD_out, 4, c('Down_DEGs_habenula', rep(NA, 6)))

## Matrix with p-values
pvals <-  data.frame("SCZ"= as.numeric(SCZ_out$P), "MDD_2019"= as.numeric(MDD_2019_out$P), "Panic"= as.numeric(Panic_out$P),
                     "SUD"= as.numeric(SUD_out$P), "MDD"= as.numeric(MDD_out$P), "OUD"= as.numeric(OUD_out$P))
## -log10 of p-values
log_pvals <- -log10(pvals)
rownames(log_pvals) <- SCZ_out$VARIABLE
log_pvals <- as.matrix(log_pvals)


## Number of DEGs in each set
gene_sets <- read.table(here('processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt'))
colnames(gene_sets) <- gene_sets[1,]
gene_sets <- gene_sets[-1,]
num_DEGs <- as.data.frame(table(gene_sets$Set))


## Heatmap
## Row annotation
row_gene_anno <- ComplexHeatmap::rowAnnotation(
    'n genes' = ComplexHeatmap::anno_barplot(num_DEGs$Freq))

## With significance cutoff at 5%
pdf(here("plots/07_MAGMA/MAGMA_pval_heatmap_pval_05.pdf"), width = 5.5, height = 3)
Heatmap(log_pvals,
        name = "-log10(p-val)",
        col = colorRampPalette(c('azure2', 'dodgerblue4'))(50),
        na_col = 'gray90',
        rect_gp = grid::gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = 'GWAS',
        row_title = 'Groups of genes',
        right_annotation = row_gene_anno,
        ## Mark cells with p<0.05
        cell_fun = function(j, i, x, y,  w, h, col){
            if(! is.na(log_pvals[i,j])){
                if (log_pvals[i,j]>(-log10(0.05))){ grid.text('*', x, y, gp = gpar(fontsize = 15, col='magenta'))}
            }
            else {
                grid.text('NA', x, y, gp = gpar(fontsize = 10, col='gray40'))
            }
        })

dev.off()


## With significance cutoff at 10%
pdf(here("plots/07_MAGMA/MAGMA_pval_heatmap_pval_10.pdf"), width = 5.5, height = 3)
Heatmap(log_pvals,
        name = "-log10(p-val)",
        col = colorRampPalette(c('azure2', 'dodgerblue4'))(50),
        na_col = 'gray90',
        rect_gp = grid::gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = 'GWAS',
        row_title = 'Groups of genes',
        right_annotation = row_gene_anno,
        ## Mark cells with p<0.1
        cell_fun = function(j, i, x, y,  w, h, col){
            if(! is.na(log_pvals[i,j])){
                if (log_pvals[i,j]>(-log10(0.1))){ grid.text('*', x, y, gp = gpar(fontsize = 15, col='magenta'))}
            }
            else {
                grid.text('NA', x, y, gp = gpar(fontsize = 10, col='gray40'))
            }
        })

dev.off()







## Reproducibility information

options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 (2024-04-24)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Monterrey
# date     2024-07-30
# rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
# pandoc   3.1.11 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────
# package        * version date (UTC) lib source
# abind            1.4-5   2016-07-21 [1] CRAN (R 4.4.0)
# berryFunctions   1.22.5  2024-02-16 [1] CRAN (R 4.4.0)
# BiocGenerics     0.50.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# circlize         0.4.16  2024-02-20 [1] CRAN (R 4.4.0)
# cli              3.6.3   2024-06-21 [1] CRAN (R 4.4.0)
# clue             0.3-65  2023-09-23 [1] CRAN (R 4.4.0)
# cluster          2.1.6   2023-12-01 [1] CRAN (R 4.4.0)
# codetools        0.2-20  2024-03-31 [1] CRAN (R 4.4.0)
# colorspace       2.1-1   2024-07-26 [1] CRAN (R 4.4.0)
# ComplexHeatmap * 2.20.0  2024-04-30 [1] Bioconductor 3.19 (R 4.4.0)
# crayon           1.5.3   2024-06-20 [1] CRAN (R 4.4.0)
# digest           0.6.36  2024-06-23 [1] CRAN (R 4.4.0)
# doParallel       1.0.17  2022-02-07 [1] CRAN (R 4.4.0)
# evaluate         0.24.0  2024-06-10 [1] CRAN (R 4.4.0)
# fastmap          1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
# foreach          1.5.2   2022-02-02 [1] CRAN (R 4.4.0)
# GetoptLong       1.0.5   2020-12-15 [1] CRAN (R 4.4.0)
# GlobalOptions    0.1.2   2020-06-10 [1] CRAN (R 4.4.0)
# here           * 1.0.1   2020-12-13 [1] CRAN (R 4.4.0)
# htmltools        0.5.8.1 2024-04-04 [1] CRAN (R 4.4.0)
# IRanges          2.38.1  2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
# iterators        1.0.14  2022-02-05 [1] CRAN (R 4.4.0)
# knitr            1.48    2024-07-07 [1] CRAN (R 4.4.0)
# limma            3.60.4  2024-07-17 [1] Bioconductor 3.19 (R 4.4.1)
# matrixStats      1.3.0   2024-04-11 [1] CRAN (R 4.4.0)
# png              0.1-8   2022-11-29 [1] CRAN (R 4.4.0)
# RColorBrewer     1.1-3   2022-04-03 [1] CRAN (R 4.4.0)
# rjson            0.2.21  2022-01-09 [1] CRAN (R 4.4.0)
# rlang            1.1.4   2024-06-04 [1] CRAN (R 4.4.0)
# rmarkdown        2.27    2024-05-17 [1] CRAN (R 4.4.0)
# rprojroot        2.0.4   2023-11-05 [1] CRAN (R 4.4.0)
# S4Vectors        0.42.1  2024-07-03 [1] Bioconductor 3.19 (R 4.4.1)
# sessioninfo    * 1.2.2   2021-12-06 [1] CRAN (R 4.4.0)
# shape            1.4.6.1 2024-02-23 [1] CRAN (R 4.4.0)
# statmod          1.5.0   2023-01-06 [1] CRAN (R 4.4.0)
# xfun             0.46    2024-07-18 [1] CRAN (R 4.4.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────
