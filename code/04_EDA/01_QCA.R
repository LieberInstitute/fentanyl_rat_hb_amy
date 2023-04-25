
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
library(ggplot2)
library(rlang)
library(smplot2)
library(Hmisc)
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
colData(rse_gene)$RNA_concentration <- as.vector(sapply(rse_gene$RNA_concentration, function(x){if (strsplit(x, '')[[1]][1]=='*'){strsplit(x, '\\*')[[1]][2]} else {x}}))
colData(rse_gene)$RNA_concentration <- as.numeric(colData(rse_gene)$RNA_concentration)

## Change numeric to char
colData(rse_gene)$Total_Num_Fentanyl_Sessions <- as.character(colData(rse_gene)$Total_Num_Fentanyl_Sessions)
colData(rse_gene)$Num_Fentanyl_Sessions_six_hrs <- as.character(colData(rse_gene)$Num_Fentanyl_Sessions_six_hrs)

## 1.2 Evaluate QC metrics for groups of samples

## QC metrics of interest
qc_metrics <- c('mitoRate', 'overallMapRate', 'totalAssignedGene', 'concordMapRate', 'library_size', 'detected_num_genes', 'RIN', 'RNA_concentration', 'Total_RNA_amount')

## Create variable of Brain Region + Substance
colData(rse_gene)$Brain_Region_and_Substance <- apply(colData(rse_gene), 1, function(x){ capitalize(paste(x['Brain_Region'], x['Substance'])) })
## Sample variables of interest
sample_variables <- c("Brain_Region", "Substance", "Brain_Region_and_Substance", "Num_Fentanyl_Sessions_six_hrs", 'Total_Num_Fentanyl_Sessions')


## Function to create boxplots of QC metrics for groups of samples

QC_boxplots <- function(qc_metric, sample_var){

    if (sample_var=="Brain_Region"){
        colors=c('amygdala'='palegreen3', 'habenula'='orchid1')
        violin_width=1
        jitter_width=0.1
    }
    else if (sample_var=="Substance"){
        colors=c('Fentanyl'='turquoise3', 'Saline'='yellow3')
        violin_width=0.7
        jitter_width=0.1
    }
    else if (sample_var=="Brain_Region_and_Substance"){
        colors=c('Amygdala Fentanyl'='springgreen3' , 'Amygdala Saline'='yellowgreen', 'Habenula Fentanyl'='hotpink1', 'Habenula Saline'='violet')
        violin_width=0.7
        jitter_width=0.09
    }
    else if (sample_var=='Total_Num_Fentanyl_Sessions'){
        colors=c('24'='salmon', '22'='pink2')
        violin_width=0.7
        jitter_width=0.1
    }
    else if (sample_var=='Num_Fentanyl_Sessions_six_hrs'){
        colors=c('18'='dodgerblue', '16'='lightskyblue')
        violin_width=0.7
        jitter_width=0.1
    }


    data <- data.frame(colData(rse_gene))
    plot <- ggplot(data = data, mapping = aes(x = !! rlang::sym(sample_var), y = !! rlang::sym(qc_metric), color = !! rlang::sym(sample_var))) +
                geom_violin(alpha = 0, size = 0.4, color='black', width=violin_width)+
                geom_jitter(width = jitter_width, alpha = 0.7, size = 2) +
                geom_boxplot(alpha = 0, size = 0.4, width=0.1, color='black') +
                scale_color_manual(values = colors) +
                sm_hgrid()

    return(plot)
}


## Multiple plots

for (sample_var in sample_variables){

    if (sample_var=="Brain_Region_and_Substance"){
        width=50
        height=40
    }
    else {
        width=25
        height=20
    }

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
    ggsave(paste("plots/01_EDA/01_QCA/QC_boxplots_", sample_var,".pdf", sep=""), width=width, height=height, units = "cm")
}









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



