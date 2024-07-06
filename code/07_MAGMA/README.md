
# Multi-marker Analysis of GenoMic Annotation (MAGMA)

The association between human orthologs of rat habenula and amygdala DEGs with SUD and psychiatric disorders was assessed based on the association significance of human genes according to their annotated SNPs. 

## 1. GWAS input data preparation

* [`01_Input_GWAS_data_prep.R`](01_Input_GWAS_data_prep.R): 
    * **1. GWAS input data preparation**: GWAS results from 6 different studies previously formatted and prepared for MAGMA by [Louise Huuki-Myers](https://lahuuki.github.io) in the [Habenula Pilot LIBD project](https://github.com/LieberInstitute/Habenula_Pilot/tree/master) were explored and further processed to use as input for MAGMA analysis:

         * **GWAS SCZ data**:
         * **GWAS MDD data**:
         * **GWAS Panic data**:
         * **GWAS SUD data**:
         * **GWAS MDD data (2)**:
         * **GWAS OUD data**:
         
        For each GWAS dataset, the `.snploc` file with the position of the SNPs was generated, as well as the `.pval` file containing the *p*-values of the same SNPs for their association with the phenotype.                                                      
        ***Note***: SNPs in sexual chromosomes were excluded from the analysis. 
       
       
       
## 2. Input DEG sets preparation

* [`02_Input_Gene_Sets_prep.R`](02_Input_Gene_Sets_prep.R): the gene sets of human orthologs of rat habenula and amygdala DEGs (all, up-, and down-regulated) were obtained using [Ensembl BioMart](http://www.ensembl.org/info/data/biomart/index.html) and their Ensembl IDs were used for subsequent MAGMA gene set analysis (`gene_sets.txt` file). 
    
    
    
## 3. Run MAGMA

* [`03_run_MAGMA.sh`](03_run_MAGMA.sh): script to run MAGMA on SNP *p*-values for each one of the previous GWASes:

    * **Step 1: Annotate SNPs onto genes**: the human gene data in the GRCh37/GRCh38* human genome assembly in the `.gene.loc` file were used for the annotation step in which SNPs from each GWAS are mapped onto genes according to their positions in the same genome build (in the `.snploc` file). The output is the `.annot` file containing the list of gene-wise SNPs.
    
        *\ GRCh37 (h19) human genome reference build was used for the SCZ (Trubetskoy et al., 2022), MDD (2019), Panic (2019), and SUD (Hatoum et al., 2023) SNPs, whereas GRCh38 (h38) build was used for the MDD (Levey et al., 2021) and OUD (Deak et al., 2022) SNPs.                       
    * **Step 2: Gene Analysis - SNP *p*-values (model SNP wise-mean)**: the SNP *p*-values (in the `.pval` file) of each gene (indicated in the `.annot` file) are transformed into $Z$-statistics and their mean (weighted sum) results in the gene test-statistic. The corresponding gene *p*-value for its association with the phenotype is computed through an approximation of the sampling distribution of the gene test-stats (**mean Z SNP statistic** default method) and then this *p*-value is converted into a *z*-score ($z_g$).                    
    As reference dataset the 1000 Genomes European panel data were used to obtain the SNP correlation matrix $R$ as most of the GWASes were conducted in Europeans. The sample size from each study is also provided (N). The output files `.genes.out` and `.genes.raw` have the gene-level association results.        
    
    * **Step 3: Gene Set Analysis**: the human orthologs of rat DEGs in each gene set (in `gene_sets.txt`) are assessed for their joint association with the phenotype based on the $z_g$ of each gene in the set (in `.genes.raw`). Specifically, the default **competitive analysis** was implemented to analyze if the genes in each set were more associated than the rest of genes with SNPs. The output file `.gsa.out` contains the gene set results. 



## 4. MAGMA results visualization

* [`04_Results_visualization.R`](04_Results_visualization.R): the resulting *p*-values for the Gene Set Analyses in the six GWAS were plotted in a heatmap. 
