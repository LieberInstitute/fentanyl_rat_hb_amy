
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
library(ggplot2)
library(rlang)
library(ggstatsplot)
library(sessioninfo)


#######################   Exploratory Data Analysis   #######################


##  1. Quality Control Analysis
## Note: genes were not filtered and counts were not normalized in these analyses

## Load RSE object at gene level
load(here('raw-data/count_objects/rse_gene_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('raw-data/count_objects/rse_exon_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('raw-data/count_objects/rse_jx_Jlab_experiment_n33.Rdata'), verbose=TRUE)
## Load sample data
sample_data <- as.data.frame(read_excel("raw-data/FentanylvsSaline_SelfAdministration_RNAextraction.xlsx"))



## 1.1 Data exploration and preparation

## Verify sample data is the same in gene, exon and jx datasets
identical(colData(rse_gene), colData(rse_exon))
# [1] TRUE
identical(colData(rse_gene), colData(rse_jx))
# [1] TRUE


## Add sample data to colData of gene RSE

## Correct colnames in sample data
colnames(sample_data) <- str_replace_all(colnames(sample_data), c(" "="_"))
## Correct sample ID in sample data
sample_data$SAMPLE_ID <- str_replace_all(sample_data$Tissue_Punch_Label, c(" "="_", "-"="_"))
## Discard unused samples
sample_data <- sample_data[which(sample_data$SAMPLE_ID %in% rse_gene$SAMPLE_ID),]
## Merge data in colData
colData(rse_gene) <- merge(colData(rse_gene), sample_data, by='SAMPLE_ID')


## Add library size for each sample
colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)

## Add detected number of genes (not zero-expressed genes) for each sample
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x){length(x[which(x>0)])})




## 1.2 Evaluate QC metrics of samples groups

## Metrics of interest
qc_metrics <- c('mitoRate', 'overallMapRate', 'totalAssignedGene', 'concordMapRate', 'library_size', 'detected_num_genes')
## Sample variables of interest
sample_variables <- c("Brain_Region", "Substance")


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
        list(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]]),
        plotgrid.args = list(nrow = 2)
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



