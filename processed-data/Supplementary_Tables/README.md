# Supplementary Tables
TODO: add Tables S1,S2 to supp tables.
TODO: add number of figs where these tables were used. 


## [Supplementary Table 1](TableS1_sample_variables_dictionary.tsv)

**Dictionary of sample variables.**

Meaning of the sample (or rat) variables and sample-level quality control metrics used throughout this project.



## [Supplementary Table 2](TableS2_sample_metadata_and_QCmetrics.tsv)

**Sample metadata and QC metrics.**

Sample-level metadata and quality control data used for EDA and DGE, as well as data regarding the behavior of fentanyl and saline rats. See [Table S1](stabl_sample_variables_dictionary.tsv) for the description of these variables. 

Created [here](https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/e081ddd5a05334d9990002b9d441ba4bc8cd3651/code/04_EDA/01_QCA.R#L69).


## [Supplementary Table 3](TableS3_de_genes_Substance_habenula.tsv)

**DEGs for substance in habenula**

Metadata, *limma* DE statistics, and Ensembl associated phenotypes and descriptions of DEGs obtained for fentanyl vs. saline in habenula. See [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html) documentation for these statistics definitions. Volcano plot in **Figure XX** was created with the data provided in this table.

Created [here](https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/2ae73935abf555c6848a9a25c4ff3c322100923b/code/05_DEA/01_Modeling.R#L500).


## [Supplementary Table 4](TableS4_de_genes_Substance_amygdala.tsv)

**DEGs for substance in amygdala**

Metadata, *limma* DE statistics, and Ensembl associated phenotypes and descriptions of DEGs obtained for fentanyl vs. saline in amygdala. See [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html) documentation for these statistics definitions. Volcano plot in **Figure XX** was created with the data provided in this table.

Created [here](https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/2ae73935abf555c6848a9a25c4ff3c322100923b/code/05_DEA/01_Modeling.R#L513).


## [Supplementary Table 5](TableS5_de_genes_common_Substance_hab_amy.tsv)

**Common DEGs for substance in habenula and amygdala**

Metadata, region-specific *limma* DE statistics, and Ensembl associated phenotypes and descriptions of overlapping DEGs for fentanyl vs. saline in habenula and amygdala. See [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html) documentation for these statistics definitions. Volcano plot in **Figure XX** was created with the data provided in this table.

Created [here](https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/5f2104ac7dbfd4534db9f1bb88ab7e124b984cd5/code/05_DEA/01_Modeling.R#L570).


## [Supplementary Table 6](TableS6_DEAs_results_all_genes_hab_amyg.tsv)

**Results for all DGE analyses and genes in habenula and amygdala** 

Gene-level metadata and *limma* DE statistics of each gene for substance and rat behavior DGE analyses (fentanyl vs. saline, first hour infusion slope, total intake, and last session intake) in habenula and amygdala. See [_limma_](https://bioconductor.org/packages/release/bioc/html/limma.html) documentation for these statistics definitions. Scatter plots in **Figure XX** were created with the data provided in this table.

Created [here](https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/ccee0d3bb43a43490fe0d6da91253c3011a61640/code/05_DEA/02_Comparisons.R#L89).



## [Supplementary Table 7](TableS7_MeanRatio_markers_top100_hab_mouse.tsv)

**Top 100 *MeanRatio* marker genes per cell type in mouse habenula** 

For all cell types and habenula neuronal cell types in the habenula complex of control mice obtained in [Hashikawa et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32272058/), the top 100 most specific marker genes for each based on the *MeanRatio* method, are reported. See [_DeconvoBuddies_](http://research.libd.org/DeconvoBuddies/) documentation for column description. `Cell_type_resolution` column corresponds to the resolution of the cell type for which the gene is a marker. 

Created [here](https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/76baa6a0271249c89f57996f5aa6e9a58eeb40d8/code/08_GSEA/01_enrich_DEGs_vs_cell_type_markers.R#L269).



## [Supplementary Table 8](TableS8_MeanRatio_markers_top100_amy_rat.tsv)

**Top 100 *MeanRatio* marker genes per cell type in rat amygdala** 
For main cell types and inhibitory neuronal subtypes in the amygdala of control rats obtained in [Zhou et al., 2023](https://www.nature.com/articles/s41593-023-01452-y), the top 100 most specific marker genes for each based on the *MeanRatio* method, are reported. See [_DeconvoBuddies_](http://research.libd.org/DeconvoBuddies/) documentation for column description. `Cell_type_resolution` column corresponds to the resolution of the cell type for which the gene is a marker.

Created [here](https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/3d6ff603b91963e03458b2dd238e60be1a3465c6/code/08_GSEA/01_enrich_DEGs_vs_cell_type_markers.R#L572).






