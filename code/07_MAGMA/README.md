
# Generalized Gene-Set Analysis of GWAS Data (MAGMA) ToDo

Describe here 

* [`01_Input_GWAS_data_prep.R`](01_Input_GWAS_data_prep.R): 
    GWAS results from 6 different studies previously formatted and prepared for MAGMA by [Louise Huuki-Myers](https://lahuuki.github.io) in the [Habenula Pilot LIBD project](https://github.com/LieberInstitute/Habenula_Pilot/tree/master) were explored and further processed to use as input for MAGMA analysis:

     * **GWAS SCZ data**:
     * **GWAS MDD data**:
     * **GWAS Panic data**:
     * **GWAS SUD data**:
     * **GWAS MDD data (2)**:
     * **GWAS OUD data**:
 
    For each dataset, the `.snploc` file with the position of the SNPs were generated, as well as the `.pval` file containing the *p*-values of the SNPs for their association with the phenotype. 

* [`02_Input_Gene_Sets_prep.R`](02_Input_Gene_Sets_prep.R): 
    the gene sets of human orthologues of rat habenula and amygdala DEGs (all, up-, and down-regulated) were obtained for MAGMA. 
    
* [`03_run_MAGMA.sh`](03_run_MAGMA.sh): script to run MAGMA on SNP *p*-values of the previous GWAS. The human gene data in the GRCh37 and GRCh38 human genome assemblies located in [geneloc/](../../processed-data/07_MAGMA/geneloc/) were used as the `.geneloc` files needed for the annotation step to map SNPs onto genes. As reference dataset in the Gene-level Analysis we used ...  The model SNP-wise mean was implemented and the ... for Gene Set Analysis. 
