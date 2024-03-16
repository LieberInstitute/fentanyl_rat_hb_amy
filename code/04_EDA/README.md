# Exploratory Data Analysis (EDA)

This directory contains the following scripts with code for exploratory analyses prior to differential expression assessment. 

## 1. Quality Control Analysis (QCA)

* [`01_QCA.R`](01_QCA.R): 
    * **1\.1 QC variable construction and correction**: additional sample-level QC metrics were computed and variables created to assess the quality of specific sample categories of interest. 
    * **1\.2 Evaluate QC metrics for groups of samples**: compare `mitoRate`, `overallMapRate`, `totalAssignedGene`, `concordMapRate`, `library_size`, `detected_num_genes`, `RIN`, `RNA_concentration`, and `Total_RNA_amount` of samples from different brain regions, from rats administered with different substances and with variable numbers of self-administration sessions, and treated in different RNA extraction and library preparation batches. 
    * **1\.3 Compare QC metrics of samples**: assess the correlation between QC metrics and RNA concentration/RNA amount for all samples together and for amygdala and habenula samples separately. 
    * **1\.4  QC-based sample filtering**: identify bad-quality samples by detecting outlier samples based on their QC metrics, for all samples together and separately for habenula and amygdala samples; all metrics of these samples were examined to further support their removal as poor-quality. 

See [TableSA](../../processed-data/Supplementary_Tables/TableSX_sample_metadata_and_QCmetrics.tsv) for sample variable meaning. 


## 2. Explore sample-level effects (PCA)

* [`02_PCA.R`](02_PCA.R):
    * **2\.1 PCA data obtention**:
    PCA was performed using the following datasets with filtered genes and log-normalized counts: 
        - `RSE` with all samples
        - `RSE` with habenula samples only
        - `RSE` with amygdala samples only
    * **2\.2 PCA visualization**: the top 6 PCs were plotted to see how variable are the samples between them based on their gene expression.
    * **2\.3 Manual sample filtering**:
        * *2\.3\.1 Identify rare samples in PCA plots*: rare (segregated) and outlier samples from `01_QCA.R` were identified in PCx vs PCy plots.
        * *2\.3\.2 Explore QC metrics of rare samples*: QC metrics of the identified rare samples and outlier samples  were re-evaluated. 
        * *2\.3\.3 Remove rare/outlier samples*: (not performed; all samples were kept).
    * **2\.4 Explore differences within fentanyl and saline sample groups**: boxplots for PCs that separated the samples were created, coloring samples by substance to observe the differences within the groups of saline and fentanyl samples.


## 3. Explore gene-level effects (Variance Partition)

## INCLUDE ADDITIONAL COVARIATES!!!!

* [`03_Explore_gene_level_effects.R`](03_Explore_gene_level_effects.R):
    * **3\.1 Analysis of explanatory variables**: 
        * *3\.1\.1 Computation of gene-wise expression variance percentages*: the percentages of variance in the expression of each gene explained by each sample variable were computed.
        * *3\.1\.2 Expression exploration of most affected genes*: explore the expression levels of most affected genes by each sample variable. 
        
    * **3\.2 Variance Partition Analysis**: compute fraction of gene expression variation attributable to each variable after correcting for all other variables.
        * *3\.2\.1 Canonical Correlation Analysis (CCA)*: the pairwise correlations between sample variables were assessed and visualized.
        * *3\.2\.2 Model fit*: a linear mixed model (LMM) was fitted to the lognorm expression data of each gene to estimate the contribution in variance of each sample variable, taking all variables and subsetting to uncorrelated variables only.
        


