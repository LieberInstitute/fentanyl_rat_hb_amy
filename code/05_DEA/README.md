# Differential Expression Analysis (DEA)

This directory contains the scripts in which all the following differential expression analyses were executed and their results visualized and compared. 

## 1. Modeling

* [`01_Modeling.R`](01_Modeling.R): 
    * **1\.1 DEA for *Fentanyl vs. Saline* with all samples from each brain region**: 
    gene expression was modeled using all sample covariates in both habenula and amygdala, and with the following uncorrelated and informative variables only (defined through variance partition analysis in [04_EDA/03_Explore_gene_level_effects.R](../04_EDA)):
        * Habenula: ~ **`Substance`** + `Batch_RNA_extraction` + `concordMapRate` + `RIN`
        * Amygdala: ~ **`Substance`** + `Batch_RNA_extraction` + `Batch_lib_prep` + `overallMapRate` + `RIN`
            
    * **1\.2 DEA for *High vs Low fentanyl intake slope* in habenula and amygdala fentanyl samples**: DGE analysis between fentanyl rats presenting high and low infusion slopes (computed with the old method based on 3 session averages) in habenula and amygdala. Gene expression was modeled as:
        * Habenula: ~ **`Intake_slope_binary`** + `concordMapRate` + `RIN`
        * Amygdala: ~ **`Intake_slope_binary`** + `overallMapRate` + `RIN`
    
    * **1\.3 DEA for *Fentanyl vs. Saline* with behavioral covariates in habenula and amygdala:**: DGE analysis for substance in habenula and amygdala adjusting gene expression for the final set of uncorrelated covariates as in **1.1** DEA but including also the rat behavioral covariates.  
        * *1\.3\.1 DEA with covariate 1st hour intake slope*:
            * Habenula: ~ `Substance` + `Batch_RNA_extraction` + `concordMapRate` + `RIN` + **`First_hr_infusion_slope`**
            * Amygdala: ~ `Substance` + `Batch_RNA_extraction` + `Batch_lib_prep` + `overallMapRate` + `RIN` + **`First_hr_infusion_slope`**
        * *1\.3\.2 DEA with covariate total intake*:
            * Habenula: ~ `Substance` + `Batch_RNA_extraction` + `concordMapRate` + `RIN` + **`Total_intake`**
            * Amygdala: ~ `Substance` + `Batch_RNA_extraction` + `Batch_lib_prep` + `overallMapRate` + `RIN` + **`Total_intake`**
        * *1\.3\.3 DEA with covariate last session intake*:
            * Habenula: ~ `Substance` + `Batch_RNA_extraction` + `concordMapRate` + `RIN` + **`Last_session_intake`**
            * Amygdala: ~ `Substance` + `Batch_RNA_extraction` + `Batch_lib_prep` + `overallMapRate` + `RIN` + **`Last_session_intake`**
    
    Then subsetting to fentanyl samples only we performed the following DGE analyses:
    * **1\.4 DEA for *1st Hour Intake Slope* in habenula and amygdala fentanyl samples**: DGE analysis for the first hour infusion slope among fentanyl rats in habenula and amygdala, adjusting for:
    
            * Habenula: ~ **`First_hr_infusion_slope`** + `concordMapRate` + `RIN` 
            * Amygdala: ~ **`First_hr_infusion_slope**` + `overallMapRate` + `RIN`    
            
        * *1\.4\.1 Analysis with all fentanyl samples*: all rats were considered for the analysis and adjusting for:

        * *1\.4\.2 Analysis without samples from negative outlier fentanyl rat*
    
    * **1\.5 DEA for *Total Intake* in habenula and amygdala fentanyl samples**: 
        * *1\.5\.1 Analysis with all fentanyl samples*:   
        * *1\.5\.2 Analysis without samples from negative outlier fentanyl rat*:  
    
    * **1\.6 DEA for *Last Session Intake* in habenula and amygdala fentanyl samples**: 
        * *1\.6\.1 Analysis with all fentanyl samples*:   
        * *1\.6\.2 Analysis without samples from negative outlier fentanyl rat*:  

See [Table S1](TableS1_sample_variables_dictionary.tsv) for the meaning of covariates. 


## 2. Comparisons
