
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
library(ggplot2)
library(rlang)
library(ggstatsplot)
library(sessioninfo)


########################   Exploratory Data Analysis   #########################


################################################################################
##                        1. Quality Control Analysis
################################################################################
## (Note: genes were not filtered and counts were not normalized in these analyses)

## Load RSE object at gene level
load(here('processed-data/03_Data_preparation/rse_gene_sample_info.Rdata'), verbose=TRUE)



## 1.1 QC variable construction and correction

## Add library size for each sample
colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)

## Add detected number of genes (not zero-expressed genes) for each sample
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x){length(x[which(x>0)])})

## Create a single and complete RIN column (add LIBD RIN if Psomagen is not available)
colData(rse_gene)$RIN <- as.numeric(apply(colData(rse_gene), 1, function(x){if (x['Psomagen_RIN']!="error"){x['Psomagen_RIN']} else {x['RIN_LIBD']}}))

## Correct data for RNA concentration
colData(rse_gene)$RNA_concentration <- colData(rse_gene)$NanoDrop_or_BioA
colData(rse_gene)$RNA_concentration <- as.vector(sapply(rse_gene$NanoDrop_or_BioA, function(x){if (strsplit(x, '')[[1]][1]=='*'){strsplit(x, '\\*')[[1]][2]} else {x}}))
colData(rse_gene)$RNA_concentration <- as.numeric(colData(rse_gene)$RNA_concentration)



## 1.2 Evaluate QC metrics for groups of samples

## QC metrics of interest
qc_metrics <- c('mitoRate', 'overallMapRate', 'totalAssignedGene', 'concordMapRate', 'library_size', 'detected_num_genes', 'RIN', 'RNA_concentration', 'Psomagen_Total_RNA_amount')
## Sample variables of interest
sample_variables <- c("Brain_Region", "Substance", "Num_Fentanyl_Sessions_six_hrs", 'Total_Num_Fentanyl_Sessions')


## Function to create boxplots of QC metrics for groups of samples

QC_boxplots <- function(qc_metric, sample_var){
    plot <- ggstatsplot::ggbetweenstats(
        data = data.frame(colData(rse_gene)),
        x = !! rlang::sym(sample_var),
        y = !! rlang::sym(qc_metric),
        mean.plotting = FALSE,
        boxtype = "boxviolin",
        xlab = str_replace_all(sample_var, c("_"=" ")),
        ylab = str_replace_all(qc_metric, c("_"=" ")),
        ## turn off messages
        bf.message = FALSE,
        results.subtitle = FALSE,
        outlier.color = "red",
        ggtheme = ggplot2::theme_gray(),
        package = "yarrr", ## package for color palette
        palette = "info2", ## color palette
        point.args = list(alpha=0.7, size=2, position = ggplot2::position_jitterdodge(dodge.width = 0.6))
    )

    return(plot)
}


## Multiple plots
for (sample_var in sample_variables) {
    i=1
    plots = list()
    for (qc_metric in qc_metrics) {
        plots[[i]]<- QC_boxplots(qc_metric, sample_var)
        i=i+1
    }
    combine_plots(
        list(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]]),
        plotgrid.args = list(nrow = 3)
    )
    ggsave(paste("plots/01_EDA/01_QCA/QC_boxplots_", sample_var ,".pdf", sep=""), width=35, height=22, units = "cm")
}




## Check dark dots

## See https://mran.microsoft.com/snapshot/2018-05-29/web/packages/ggstatsplot/vignettes/ggbetweenstats.html








## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# setting  value
# version  R version 4.2.2 (2022-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-04-20
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA



