# Differential Expression Analysis (DEA)

This directory contains a single script in which all the following differential expression analyses were carried out. 

## 1. Modeling

* [`01_Modeling.R`](01_Modeling.R): 
    * **1\.1 DEA for *Fentanyl vs. Saline* with all samples from each brain region**: 
    gene expression was modeled using all sample covariates in both habenula and amygdala, and with the uncorrelated variables only (defined through variance partition analysis in [04_EDA/03_Explore_gene_level_effects.R](../04_EDA)):
        * Habenula:    ~ `Substance` + `` + `` + ``
        * Amygdala: 
            
    * **1\.2 DEA for *High vs Low fentanyl intake slope* in habenula and amygdala fentanyl samples**:
    
    * **1\.3 DEA for *Fentanyl vs. Saline* with behavioral covariates:**: 
        * *1\.3\.1 DEA with covariate 1st hour intake slope*:
        * *1\.3\.2 DEA with covariate total intake*:
        * *1\.3\.3 DEA with covariate last session intake*:
    
    Then subsetting to fentanyl samples only:
    * **1\.4 DEA for *1st Hour Intake Slope* in habenula and amygdala fentanyl samples**: 
        * *1\.4\.1 Analysis with all fentanyl samples*:   
        * *1\.4\.2 Analysis without samples from negative outlier fentanyl rat*:    
    
    * **1\.5 DEA for *Total Intake* in habenula and amygdala fentanyl samples**: 
        * *1\.5\.1 Analysis with all fentanyl samples*:   
        * *1\.5\.2 Analysis without samples from negative outlier fentanyl rat*:  
    
    * **1\.6 DEA for *Last Session Intake* in habenula and amygdala fentanyl samples**: 
        * *1\.6\.1 Analysis with all fentanyl samples*:   
        * *1\.6\.2 Analysis without samples from negative outlier fentanyl rat*:  
        
