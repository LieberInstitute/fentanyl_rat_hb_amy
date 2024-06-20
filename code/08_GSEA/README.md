
# Gene Set Enrichment Analysis (GSEA)

The association between the differential expression of genes in rat habenula and amygdala and their cell type-specific expression was assessed based on the enrichment of such DEGs among rat orthologs of human marker genes for epithalamus (including the habenula) and amygdala cell types at broad and fine resolutions. 

Cell type marker genes in human epithalamus and amygdala were obtained applying the *MeanRatio* and *1vsALL* methods -described in [Huuki-Myers, L., et al. (2024) ](https://doi.org/10.1101/2024.02.09.579665)- using snRNA-seq data from [Yalcinbas, E., et al. (2024)]( https://doi.org/10.1101/2024.02.26.582081) and [Yu, B., et al. (2023)](https://doi.org/10.1038/s41421-022-00506-y), respectively. The following were the sets of markers used:

* *A) MeanRatio-based cell type marker genes*
    * i) Markers for cell types in human epithalamus:
    
        * Fine resolution cell type markers: the top 50 marker genes for fine cell types in human epithalamus. 
        
    * ii) Markers for cell types in human amygdala:
    
        * Broad resolution cell type markers: top 100 marker genes for broad cell types in human amygdala. 
        * Fine resolution cell type markers: top 100 marker genes for fine cell types in human amygdala. 
            
    ***Note***: only genes with mean expression ratios >1 were considered markers as these are, by definition, more expressed in the target cell type. 


* *B) 1vsALL-based cell type marker genes in human*
    * i) Markers for cell types in human epithalamus:
    
        * Broad resolution cell type markers: marker genes (DEGs with FDR<0.05) for broad cell types in human epithalamus.
        * Fine resolution cell type markers: marker genes (DEGs with FDR<0.05) for fine cell types in human epithalamus.
        
    * ii) Markers for cell types in human amygdala:
    
        * Broad resolution cell type markers: marker genes (DEGs with FDR<0.05) for broad cell types in human amygdala.
        * Fine resolution cell type markers: marker genes (DEGs with FDR<0.05) for fine cell types in human amygdala.


* [`01_enrich_DEGs_vs_cell_type_markers.R`](01_enrich_DEGs_vs_cell_type_markers.R): 

    * **1\. Obtain sets of orthologs of human marker genes in rat**: for the human marker genes in each of the sets previously described, their rat orthologs were obtained using ENSEMBL under the GRCh38 human genome version. 
        
    * **2\. Cell type enrichment analysis for rat DEGs**: all, up-, and down-regulated rat habenula and amygdala DEGs were separately assessed for their enrichment among rat genes that have at least one human ortholog that is a marker gene for a given cell type in epithalamus or amygdala, respectively; these human cell type markers genes are provided in the sets described above. Only rat genes that were expressed and assessed for DGE were considered in the analysis (n=16,708). Results are presented with heatmaps for each set of cell type markers and group of DEGs.  
    
