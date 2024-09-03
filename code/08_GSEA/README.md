
# Cell Type Enrichment Analysis

The association between the differential expression of genes in rat habenula and amygdala after fentanyl self-administration and their cell type-specific expression was assessed based on the enrichment among such rat DEGs, of rat marker genes for main and fine amygdala cell types, and rat orthologs of human marker genes for epithalamus (including the habenula) and amygdala cell types at broad and fine resolutions, as well as rat orthologs of mouse marker genes for all or only neuronal cell subpopulations of the habenula complex.

Cell type marker genes in human epithalamus and amygdala, mouse habenula, and rat amygdala were obtained applying the *MeanRatio* and *1vsALL* methods -described in [Huuki-Myers, L., et al. (2024) ](http://research.libd.org/DeconvoBuddies/index.html)- using snRNA-seq data from [Yalcinbas, E., et al. (2024)](https://doi.org/10.1101/2024.02.26.582081) for human epithalamus, [Yu, B., et al. (2023)](https://doi.org/10.1038/s41421-022-00506-y) for human amygdala, scRNA-seq data from [Hashikawa, Y., et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32272058/) for mouse habenula, and snRNA-seq data from [Zhou, J., et al. (2023)](https://www.nature.com/articles/s41593-023-01452-y) for rat amygdala. 

The following were the sets of markers used:

* *A) MeanRatio-based cell type marker genes*
    * i) Markers for cell types in human epithalamus:
    
        * Fine resolution cell type markers: the top 50 marker genes for each fine cell type in human epithalamus. 
        
    * ii) Markers for cell types in human amygdala:
    
        * Broad resolution cell type markers: the top 100 marker genes for each broad cell type in human amygdala. 
        * Fine resolution cell type markers: the top 100 marker genes for each fine cell type in human amygdala. 
        
    * iii) Markers for cell types in mouse habenula complex: 
        
        * All cell subpopulation markers: the top 100 marker genes for each of all cell subpopulations in mouse habenula complex.
        * Habenula neuronal subpopulation markers: the top 100 marker genes for each of the habenula neuronal subpopulations in mouse habenula complex.
        
    * iv) Markers for cell types in rat amygdala: 
        
        * Main cell type markers: the top 100 marker genes for each main cell type in rat amygdala. 
        * Fine cell type markers: the top 100 marker genes for each fine cell type in rat amygdala. 

            
    ***Note***: only genes with mean expression ratios >1 were considered markers as these are, by definition, more expressed in the target cell type than in any other cell type. 


* *B) 1vsALL-based cell type marker genes*
    * i) Markers for cell types in human epithalamus:
    
        * Broad resolution cell type markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each broad cell type in human epithalamus.
        * Fine resolution cell type markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each fine cell type in human epithalamus.
        
    * ii) Markers for cell types in human amygdala:
    
        * Broad resolution cell type markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each broad cell type in human amygdala.
        * Fine resolution cell type markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each fine cell type in human amygdala.
        
    * iii) Markers for cell types in mouse habenula complex:
    
        * All cell subpopulation markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each of all cell subpopulations in mouse habenula complex.
        * Habenula neuronal subpopulation markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each of the habenula neuronal subpopulations in mouse habenula complex.
        
    * iv) Markers for cell types in rat amygdala: 
        
        * Main cell type markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each main cell type in rat amygdala. 
        * Fine cell type markers: marker genes (DEGs with FDR<0.05 and logFC>0) for each fine cell type in rat amygdala. 


* [`01_enrich_DEGs_vs_cell_type_markers.R`](01_enrich_DEGs_vs_cell_type_markers.R): 
 
    * **1\. Obtain cell type marker genes in rodent habenula and amygdala**: 
        * *1\.1 Obtain cell type marker genes in mouse habenula*: scRNA-seq data from Hashikawa, Y., et al. (2020) were explored and log-normalized. Only data from control mice were used and cell types with <10 cells were discarded. These processed data were subsequently used to find mouse marker genes for (all/neuronal) habenula subpopulations with *MeanRatio* and *1vsALL*. 
        
        * *1\.2 Obtain cell type marker genes in rat amygdala*: snRNA-seq data from Zhou, J., et al. (2023) were explored and normalized counts of cells from control mice were used to find *MeanRatio* and *1vsALL* markers for main cell types and Inhib. neuronal subtypes.
    
        ***Note***: human habenula and amygdala cell type markers were generated and provided by [Louise Huuki-Myers](https://lahuuki.github.io) and [Michael Totty](https://github.com/MicTott), respectively. 
    
    * **2\. Obtain sets of orthologs of human/mouse marker genes in rat**: for the human and mouse marker genes in each of the sets previously described, their rat orthologs were obtained using ENSEMBL under the human GRCh38 and mouse GRCm39 genome versions, respectively. This was not necessary for rat amygdala markers. 
    
    * **3\. Compare *MeanRatio* vs *1vsALL* marker genes**: the ratio of *MeanRatio* marker genes for broad and fine human cell types, for all and neuronal mouse cell types, and for main and fine rat cell types, was compared against their *1vsALL* standard logFC (aka *t*-stats) for the same or equivalent cell types, identifying common markers with both methods.
        
    * **4\. Cell type enrichment analysis for rat DEGs**: rat amygdala markers and rat genes that have at least one human/mouse ortholog that is a marker gene for a given cell type in human epithalamus/amygdala, or mouse habenula, correspondingly, were assessed for their enrichment among all, up-, and down-regulated rat habenula and amygdala DEGs. These sets of human/mouse/rat cell type marker genes are described above. Only rat genes that were expressed and assessed for DGE were considered in the analysis (n=16,708). Results are presented in heatmaps for each set of cell type markers and group of DEGs.  
    
