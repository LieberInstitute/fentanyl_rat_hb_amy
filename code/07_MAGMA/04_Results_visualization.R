
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

pdf(here("plots/07_MAGMA/MAGMA_pval_heatmap.pdf"), width = 6.4, height = 4)
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







## Reproducibility information

options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-05-07
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────
# package        * version date (UTC) lib source
# abind            1.4-5   2016-07-21 [2] CRAN (R 4.3.1)
# berryFunctions * 1.22.5  2024-02-16 [1] CRAN (R 4.3.1)
# BiocGenerics     0.46.0  2023-04-25 [2] Bioconductor
# Cairo            1.6-1   2023-08-18 [2] CRAN (R 4.3.1)
# circlize         0.4.15  2022-05-10 [2] CRAN (R 4.3.1)
# cli              3.6.1   2023-03-23 [2] CRAN (R 4.3.1)
# clue             0.3-64  2023-01-31 [2] CRAN (R 4.3.1)
# cluster          2.1.4   2022-08-22 [3] CRAN (R 4.3.1)
# codetools        0.2-19  2023-02-01 [3] CRAN (R 4.3.1)
# colorspace       2.1-0   2023-01-23 [2] CRAN (R 4.3.1)
# ComplexHeatmap * 2.16.0  2023-04-25 [2] Bioconductor
# crayon           1.5.2   2022-09-29 [2] CRAN (R 4.3.1)
# digest           0.6.33  2023-07-07 [2] CRAN (R 4.3.1)
# doParallel       1.0.17  2022-02-07 [2] CRAN (R 4.3.1)
# foreach          1.5.2   2022-02-02 [2] CRAN (R 4.3.1)
# GetoptLong       1.0.5   2020-12-15 [2] CRAN (R 4.3.1)
# GlobalOptions    0.1.2   2020-06-10 [2] CRAN (R 4.3.1)
# here           * 1.0.1   2020-12-13 [2] CRAN (R 4.3.1)
# IRanges          2.34.1  2023-06-22 [2] Bioconductor
# iterators        1.0.14  2022-02-05 [2] CRAN (R 4.3.1)
# magick           2.7.5   2023-08-07 [2] CRAN (R 4.3.1)
# magrittr         2.0.3   2022-03-30 [2] CRAN (R 4.3.1)
# matrixStats      1.0.0   2023-06-02 [2] C
