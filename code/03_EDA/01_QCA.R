
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
library(ggplot2)
library(rlang)
library(sessioninfo)


#######################   Exploratory Data Analysis   #######################


##  1. Quality Control Analysis
## Note: features are not filtered and counts are not normalized in this analysis

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
sample_data$SAMPLE_ID <- str_replace_all(sample_data$`Tissue Punch Label`, c(" "="_", "-"="_"))
## Discard unused samples
sample_data <- sample_data[which(sample_data$SAMPLE_ID %in% rse_gene$SAMPLE_ID),]
## Merge data in colData
colData(rse_gene) <- merge(colData(rse_gene), sample_data, by='SAMPLE_ID')
## Correct colnames


## Add library size for each sample
colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)

## Add detected number of genes (not zero-expressed genes) for each sample
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x){length(x[which(x>0)])})




## 1.2 Evaluate QC metrics of groups of samples

## Metrics of interest
qc_metrics <- c('mitoRate', 'rRNA_rate', 'overallMapRate', 'totalAssignedGene', 'concordMapRate', 'library_size', 'detected_num_genes')

plot <-  ggplot(data.frame(colData(rse_gene)), aes(x='`Brain Region`', y=eval(parse_expr('mitoRate')),
                          fill='`Brain Region`')) +
    geom_violin(width=0.98, aes(fill='Sex'), alpha=0.5, color=NA) +
    geom_boxplot(width=0.2, lwd=0.3, aes(col=Protocol, alpha=0.1), alpha=0, outlier.shape=NA, show.legend = F) +
    geom_point(data=subset(md, outlier), col='#707070', alpha=0.7, show.legend=F)














