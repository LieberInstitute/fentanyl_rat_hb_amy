
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)

#######################   Exploratory Data Analysis   #######################


##  1. Quality Control Analysis


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

## Correct sample ID in sample data
sample_data$SAMPLE_ID <- str_replace_all(sample_data$`Tissue Punch Label`, c(" "="_", "-"="_"))
## Discard unused samples
sample_data <- sample_data[which(sample_data$SAMPLE_ID %in% rse_gene$SAMPLE_ID),]
## Merge data in colData
colData(rse_gene) <- merge(colData(rse_gene), sample_data, by='SAMPLE_ID')


