
library(here)
library(SummarizedExperiment)
library(readxl)

#######################   Exploratory Data Analysis   #######################


##  1. Quality Control Analysis


## Load RSE object at gene level
load(here('raw-data/count_objects/rse_gene_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('raw-data/count_objects/rse_exon_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('raw-data/count_objects/rse_jx_Jlab_experiment_n33.Rdata'), verbose=TRUE)

## Load sample data
sample_data <- as.data.frame(read_excel("raw-data/FentanylvsSaline_SelfAdministration_RNAextraction.xlsx"))


## Verify sample data is the same in gene, exon and jx datasets
identical(colData(rse_gene), colData(rse_exon))
# [1] TRUE
identical(colData(rse_gene), colData(rse_jx))
# [1] TRUE




