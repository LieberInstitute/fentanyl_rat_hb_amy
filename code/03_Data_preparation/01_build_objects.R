
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
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
## Correct sample ID in sample data
sample_data$SAMPLE_ID <- str_replace_all(sample_data$Tissue_Punch_Label, c(" "="_", "-"="_"))
## Correct Sample_No. for the extra sample
sample_data[which(sample_data$SAMPLE_ID=="33_S_Amyg_20"), "Sample_No."]="33.0"

## Verify samples are the same in sample_data and rse
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
# date     2023-04-24
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA
