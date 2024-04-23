# Script summary

## 01. Transfer data
FASTQ files were transferred from their original locations to [raw-data](../raw-data). 
    
See more details in [01_transfer_data](01_transfer_data/).

## 02. *SPEAQeasy* pipeline

*SPEAQeasy* pipeline was run to obtain the `RangedSummarizedExperiment` (RSE) objects.

See more details in [02_SPEAQeasy](02_SPEAQeasy/). 


## 03. Data processing

Sample metadata was cleaned and correctly formatted and added to RSE objects. Count normalization and gene filtering steps were carried out here. 

See more details in [03_Data_processing](03_Data_preparation/).


## 04. Exploratory Data Analysis (EDA)

Quality control analysis (QCA) was performed, as well as the exploration of sample-level and gene-level effects through dimensionality reduction  (Principal Component Analysis) and variance partition analyses, respectively. 

See more details in [04_EDA](04_EDA/).


## 05. Differential Expression Analysis (DEA)

Differential Expression Analyses were performed to assess the effect of substance, total fentanyl intake, last session fentanyl intake, and fentanyl infusion slope on gene expression in habenula and amygdala of rats. Results were visualized and compared within each brain region. 

See more details and specifications about each DGE analysis in [05_DEA](05_DEA/).


## 06. Functional Enrichment Analysis (GO & KEGG terms)



## 07. Multi-marker analysis of genomic annotation (MAGMA)


