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

Differential Expression Analyses were performed to assess the effect of substance, total fentanyl intake, last session fentanyl intake, and fentanyl infusion slope on gene expression in habenula and amygdala of rats. Results were visualized and compared within and between brain regions. 

See more details and specifications about each DGE analysis in [05_DEA](05_DEA/).


## 06. Functional Enrichment Analysis (GO & KEGG terms)

Analysis of overrepresentation of habenula and amygdala DEGs for substance among gene sets of [GO](https://geneontology.org) and [KEGG](https://www.genome.jp/kegg/) terms. 
 
See more details in [06_GO_KEGG](06_GO_KEGG/).


## 07. Multi-marker Analysis of GenoMic Annotation (MAGMA)

MAGMA analysis performed to find associations between the sets of DEGs with multiple substance use disorders (OUD and more broadly SUD) and mental disorders (MDD, panic, and schizophrenia) based on publicly available GWAS data. 

See more details in [07_MAGMA](07_MAGMA/).

## 08. Gene Set Enrichment Analysis (GSEA)

The enrichment of rat orthologs of cell type-specific marker genes in human epithalamus and amygdala, and mouse habenula, among habenula and amygdala DEGs, was assessed. 

See more details of the human and mouse marker sets used and the analysis in [08_GSEA](08_GSEA/).


