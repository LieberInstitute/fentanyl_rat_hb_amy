# Exploratory Data Analysis (EDA)

This directory contains the following scripts with code for Exploratory Data Analyses prior to differential expression assessment. 

## 1. Quality Control Analysis (QCA)

* [`01_QCA.R`](01_QCA.R): 
    * **1\.1 QC variable construction and correction**: additional sample-level QC metrics were computed and variables created to assess the quality of specific sample categories of interest. 
    * **1\.2 Evaluate QC metrics for groups of samples**: compare `mitoRate`, `overallMapRate`, `totalAssignedGene`, `concordMapRate`, `library_size`, `detected_num_genes`, `RIN`, `RNA_concentration`, and `Total_RNA_amount` of samples from different brain regions, administered with different substances, with variable numbers of administration sessions and treated in different RNA extraction and library preparation batches. 
    * **1.\3 Compare QC metrics of samples**: assess the correlation between QC metrics and RNA concentration/RNA amount for all samples together and for amygdala and habenula samples separately. 
    * **1.\4  QC-based sample filtering**: identify poor-quality samples by detecting outlier samples based on their QC metrics, for all samples together and separately for habenula and amygdala samples; all metrics of these samples were examined to further support their removal as poor quality. 

See [TableSX](../../processed-data/Supplementary_Tables/TableSX_sample_metadata_and_QCmetrics.tsv) for sample variable meaning. 


## 2. Explore sample-level effects (PCA)



## 3. Explore gene-level effects (Variance Partition)
