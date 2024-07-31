
# Multi-marker Analysis of GenoMic Annotation (MAGMA)

The association between sets of human orthologs of rat habenula and amygdala DEGs with SUD and psychiatric disorders was assessed based on the association significance of human genes according to their annotated SNPs. 

## 1. GWAS input data preparation

* [`01_Input_GWAS_data_prep.R`](01_Input_GWAS_data_prep.R): GWAS results from 6 different studies previously formatted and prepared for MAGMA by [Louise Huuki-Myers](https://lahuuki.github.io) in the [Habenula Pilot LIBD project](https://github.com/LieberInstitute/Habenula_Pilot/tree/master) were explored and further processed to use as input for MAGMA analysis:

     * **Schizophrenia (SCZ) GWAS**: 
        * Trubetskoy, V. et al., *Nature*, 2022 (doi: [10.1038/s41586-022-04434-5](https://doi.org/10.1038/s41586-022-04434-5)).
        * This was a two-stage GWAS of up to 76,755 SCZ cases and 243,649 controls from multiple cohorts. 74.3% of the individuals in the primary GWAS (175,799: 74,776 cases and 101,023 controls) were EUR. Sample size (N) used was the one provided in the dataset.
        * In the primary GWAS a total of 7,585,078 SNPs were analyzed, of which 313 independent common SNPs were genome-wide significant (GWS) spanning 263 loci. In the data 7,659,767 total SNPs are present.
        * GRCh37 genome build used.                                                                                             
        
     * **Panic Disorder (PD) GWAS**:
        * Forstner, A. J. et al., *Molecular Psychiatry*, 2019 (doi: [10.1038/s41380-019-0590-2](https://doi.org/10.1038/s41380-019-0590-2)).
        * The GWAS discovery meta-analysis comprised 2,248 PD patients and 7,992 controls, all of European ancestry. 333 were excluded after QC, leading to N=9,907 total analyzed individuals: 2,147 PD patients and 7,760 controls.
        * In the data 10,151,300 total SNPs are provided but in the analysis only 8,757,275 SNPs that passed QC were included. No significant SNPs detected.
        * GRCh37 genome build used here.
        
    * **Opioid Use Disorder (OUD) GWAS**:
        * Deak, J. D. et al., *Molecular Psychiatry*, 2022 (doi: [10.1038/s41380-022-01709-1](https://doi.org/10.1038/s41380-022-01709-1)).
        * A total of N=639,063 individuals were analyzed: 20,686 cases and 618,377 controls. 554,186 (>86%) were EUR and the rest AFR. 
        * Data available for 1,322 SNPs. 3 GWS variants in EUR, none in AFR.
        * Genome Build GRCh38 used here.                                                                                       
     * **Substance Use Disorder (SUD) GWAS**:
        * Hatoum, A. S. et al., *Nature Mental Health*, 2023 (doi: [10.1038/s44220-023-00034-y](https://doi.org/10.1038/s44220-023-00034-y)).
        * Sample of 1,025,550 EUR and 92,630 AFR individuals. N used were the ones provided in the dataset.
        * Data available for 4,187,212 SNPs. 19 independent SNPs were GWS for general addition; other SNPs for specific substances.
        * GRCh37 genome build used here.
     
     * **Major Depressive Disorder (MDD) GWAS I**:
        * Howard, D. M. et al., *Nature neuroscience*, 2019 (doi: [10.1038/s41593-018-0326-7](https://doi.org/10.1038/s41593-018-0326-7)).
        * N=807,553 individuals (246,363 cases and 561,190 controls) from three previous studies of depression were used (one on UK Biobank individuals, the second from 23andMe, and the third using European-ancestry PGC cohorts).
        * Data provided for 8,481,694 SNPs but only 8,098,588 were tested. 102 independent genetic variants in 101 loci were GWS.
        * GRCh37 genome build used here.
     
    * **Major Depressive Disorder GWAS II**:
        * Levey, D. F. et al., *Nature neuroscience*, 2021 (doi: [10.1038/s41593-021-00860-2](https://doi.org/10.1038/s41593-021-00860-2)).
        * N=1,154,267 EUR individuals (340,591 cases; primary analysis) and 59,600 AFR individuals (25,843 cases).
        * Data for 11,700 SNPs. 178 genetic risk loci and 223 independently significant SNPs were identified.
        * GRCh38 genome build used here.
         

For each GWAS dataset, the `.snploc` file with the position of the SNPs was generated, as well as the `.pval` file containing the *p*-values of the same SNPs for their association with the phenotype.                                                      
        ***Note***: SNPs in sexual chromosomes were excluded from the analysis. 
       
       
       
## 2. Input DEG sets preparation

* [`02_Input_Gene_Sets_prep.R`](02_Input_Gene_Sets_prep.R): the gene sets of human orthologs of rat habenula and amygdala DEGs (all, up-, and down-regulated) were obtained using [Ensembl *biomaRt*](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) and their Ensembl IDs (`gene_sets.txt` file) were used for MAGMA gene set analysis. 
    
    
    
## 3. Run MAGMA

* [`03_run_MAGMA.sh`](03_run_MAGMA.sh): script to run MAGMA on SNP *p*-values for each one of the previous GWASes:

    * **Step 1: Annotate SNPs onto genes**: the human gene data in either the GRCh37/GRCh38** human genome assembly (in the `.gene.loc` file) were used for the annotation step in which SNPs from each GWAS are mapped onto genes according to their positions in the same genome build (in the `.snploc` file). The output is the `.annot` file containing the list of gene-wise SNPs.          
    
         <span style="color: gray;  font-size: 1"> ** GRCh37 (h19) human genome reference build was used for the SCZ (Trubetskoy et al., 2022), MDD (Howard et al., 2019), Panic (Forstner et al., 2019), and SUD (Hatoum et al., 2023) SNPs, whereas GRCh38 (h38) build was used for the MDD (Levey et al., 2021) and OUD (Deak et al., 2022) SNPs.  </span>                          
    * **Step 2: Gene Analysis - SNP *p*-values (model SNP wise-mean)**: the SNP *p*-values (in the `.pval` file) of each gene (indicated in the `.annot` file) are transformed into $Z$-statistics and their mean (weighted sum) results in the gene test-statistic. The corresponding gene *p*-value for its association with the phenotype is computed through an approximation of the sampling distribution of the gene test-stats (**mean Z SNP statistic** default method) and then this *p*-value is converted into a *z*-score ($z_g$).                    
    As reference dataset the 1000 Genomes European panel data were used to obtain the SNP correlation matrix $R$ as all/most individuals analyzed in the GWASes were Europeans. The sample size from each study is also provided (N). The output files `.genes.out` and `.genes.raw` have the gene-level association results.        
    
    * **Step 3: Gene Set Analysis**: the human orthologs of rat DEGs in each gene set (in `gene_sets.txt`) are assessed for their joint association with the phenotype based on the $z_g$ of each gene inside/outside the set (in `.genes.raw`). Specifically, the default **competitive analysis** was implemented to analyze if the genes in each set were more associated than the rest of genes (that also have mapped SNPs). The output file `.gsa.out` contains the gene set results. 



## 4. MAGMA results visualization

* [`04_Results_visualization.R`](04_Results_visualization.R): the resulting *p*-values from the Gene Set Analyses using the six GWASes were plotted in a heatmap. 
