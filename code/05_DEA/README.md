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
    
    * **1\.3 DEA for *Fentanyl vs. Saline* with behavioral covariates in habenula and amygdala**: DGE analysis for substance in habenula and amygdala adjusting gene expression for the final set of uncorrelated covariates as in **1.1** DEA but including also the rat behavioral covariates.  
        * *1\.3\.1 DEA with covariate 1st Hour Intake Slope*:
            * Habenula: ~ `Substance` + `Batch_RNA_extraction` + `concordMapRate` + `RIN` + **`First_Hour_Infusion_Slope`**
            * Amygdala: ~ `Substance` + `Batch_RNA_extraction` + `Batch_lib_prep` + `overallMapRate` + `RIN` + **`First_Hour_Infusion_Slope`**
        * *1\.3\.2 DEA with covariate Total Intake*:
            * Habenula: ~ `Substance` + `Batch_RNA_extraction` + `concordMapRate` + `RIN` + **`Total_Intake`**
            * Amygdala: ~ `Substance` + `Batch_RNA_extraction` + `Batch_lib_prep` + `overallMapRate` + `RIN` + **`Total_Intake`**
        * *1\.3\.3 DEA with covariate Last Session Intake*:
            * Habenula: ~ `Substance` + `Batch_RNA_extraction` + `concordMapRate` + `RIN` + **`Last_Session_Intake`**
            * Amygdala: ~ `Substance` + `Batch_RNA_extraction` + `Batch_lib_prep` + `overallMapRate` + `RIN` + **`Last_Session_Intake`**
            
    The following DGE analyses for rat behavior were performed on fentanyl samples only and were done with habenula and amygdala samples from all rats self-administered fentanyl and then excluding samples from a negative outlier fentanyl rat detected in **1.2**.
    * **1\.4 DEA for *1st Hour Intake Slope* in habenula and amygdala fentanyl samples**: DGE analysis for the first hour infusion slope among fentanyl rats in habenula and amygdala, adjusting for:

        * Habenula: ~ **`First_Hour_Infusion_Slope`** + `RIN` + `RNA_concentration` + `mitoRate`
        * Amygdala: ~ **`First_Hour_Infusion_Slope`** + `RIN` + `mitoRate`    
    
    * **1\.5 DEA for *Total Intake* in habenula and amygdala fentanyl samples**: DGE analysis for the total drug intake of the fentanyl rats in habenula and amygdala, adjusting gene expression for:
        * Habenula: ~ **`Total_Intake`** + `RIN` + `RNA_concentration` + `overallMapRate`
        * Amygdala: ~ **`Total_Intake`** + `RIN` + `mitoRate`    
    
    * **1\.6 DEA for *Last Session Intake* in habenula and amygdala fentanyl samples**: DGE analysis for the drug intake in the last session for each fentanyl rat in habenula and amygdala, adjusting for:
        * Habenula: ~ **`Last_Session_Intake`** + `RIN` + `RNA_concentration` + `mitoRate`
        * Amygdala: ~ **`Last_Session_Intake`** + `RIN` + `totalAssignedGene` + `concordMapRate`


See [Table S1](TableS1_sample_variables_dictionary.tsv) for the meaning of covariates. 


## 2. Comparisons
* [`02_Comparisons.R`](02_Comparisons.R): 
    * **2\.1 Comparison of gene DE signal for substance and behavior in habenula**:
    
    * **2\.2 Comparison of gene DE signal for substance and behavior in amygdala**:
        
