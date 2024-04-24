
# Functional Enrichment Analysis (GO & KEGG terms)

Habenula and amygdala DEGs for substance obtained in [05_DEA/01_Modeling.R](../05_DEA) were assessed for their enrichment in [GO](https://geneontology.org) and [KEGG](https://www.genome.jp/kegg/) gene sets. 

* [`01_GO_KEGG_Analyses.R`](01_GO_KEGG_Analyses.R): GO biological processes (BP), molecular functions (MF), cellular components (CC), and KEGG pathways (KEGG) enriched in the clusters of DEGs were identified, specifically performing:  
    * **1\. Analysis for all DEGs from each brain region**: all habenula and amygdala DEGs were tested. 
    * **2\. Analysis for up- and down-regulated DEGs from each brain region**: habenula and amygdala up- and down-regulated DEGs were separately tested. 
