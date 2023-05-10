
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
library(edgeR)
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
## Add design matrix to account for group differences and define n
rse_gene_filt<-rse_gene[which(filterByExpr(assay(rse_gene),
                                           design=with(colData(rse_gene), model.matrix(~ Brain_Region + Substance)))),]
save(rse_gene_filt, file = 'processed-data/03_Data_preparation/rse_gene_filt.Rdata')
dim(rse_gene_filt)
# 16708    33
## Percentage of genes that were kept
dim(rse_gene_filt)[1]*100/dim(rse_gene)[1]
# 54.86668


## Filter exons
rse_exon_filt<-rse_exon[which(filterByExpr(assay(rse_exon),
                                           design=with(colData(rse_exon), model.matrix(~ Brain_Region + Substance)))),]
save(rse_exon_filt, file = 'processed-data/03_Data_preparation/rse_exon_filt.Rdata')
dim(rse_exon_filt)
# 182291     33
## Percentage of exons that were kept
dim(rse_exon_filt)[1]*100/dim(rse_exon)[1]
#  66.76984


## Filter junctions
rse_jx_filt<-rse_jx[which(filterByExpr(assay(rse_jx),
                                       design=with(colData(rse_jx), model.matrix(~ Brain_Region + Substance)))),]
save(rse_jx_filt, file = 'processed-data/03_Data_preparation/rse_jx_filt.Rdata')
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
# save(rse_tx_filt, file = 'processed-data/03_Data_preparation/rse_tx_filt.Rdata')
# dim(rse_tx_filt)
#






## 4. Visualization
## Plots of the distribution of genes' counts before and after normalization and filtering

## Raw counts
counts_data <- data.frame(counts=as.vector(assays(rse_gene)$counts))
ggplot(counts_data, aes(x=counts)) +
    geom_histogram(colour="black", fill="lightgray") +
    theme_classic() +
    labs(x='read counts', y='Frecuency')
ggsave(filename = 'plots/03_Data_preparation/Hist_counts.pdf', height = 3, width = 4)

## Normalized counts
logcounts_data <- data.frame(logcounts=as.vector(assays(rse_gene)$logcounts))
ggplot(logcounts_data, aes(x=logcounts)) +
    geom_histogram(aes(y=..density..), colour="darkgray", fill="lightgray") +
    theme_classic() +
    geom_density(fill="#69b3a2", alpha=0.3) +
    labs(x='log(CPM+0.5)', y='Frecuency')
ggsave(filename = 'plots/03_Data_preparation/Hist_logcounts.pdf', height = 3, width = 4)

## Normalized counts and filtered genes
filt_logcounts_data <- data.frame(logcounts=as.vector(assays(rse_gene_filt)$logcounts))
ggplot(filt_logcounts_data, aes(x=logcounts)) +
    geom_histogram(aes(y=..density..), colour="darkgray", fill="lightgray") +
    theme_classic() +
    geom_density(fill="#69b3a2", alpha=0.3) +
    labs(x='log(CPM+0.5)', y='Frecuency')
ggsave(filename = 'plots/03_Data_preparation/Hist_filt_logcounts.pdf', height = 3, width = 4)









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
