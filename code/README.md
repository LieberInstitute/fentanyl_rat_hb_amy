# Script summary

## 01. Transfer data




## 02. SPEAQeasy pipeline



## 03. Data processing

* [01_build_objects.R](03_Data_preparation/01_build_objects.R): 
    * 1. Sample metadata were added to `colData()` of the RSE objects in the correct format. 
    * 2. Raw gene expression counts were log-normalized to log(CPM+0.5).
    * 3. Lowly-expressed genes were filtered out. Additional QC metrics and sample filtering info from [04_EDA/01_QCA.R](04_EDA/01_QCA.R) analysis were included in these original datasets.
    * 4. Count distributions were visualized before and after normalization and feature filtering steps.

Note that exon, exon-exon junction, and transcript RSE objects were also processed but only gene-level data were used for downstream analyses.


## 04. Exploratory Data Analysis (EDA)



## 05. Differential Expression Analysis (DEA)



## 06. Functional Enrichment Analysis (GO & KEGG terms)



## 07. Multi-marker analysis of genomic annotation (MAGMA)



## 08. snRNA-data comparisons (*chage name in the future*)
