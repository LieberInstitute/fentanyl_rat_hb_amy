
# Data processing

This directory contains the following script to process the original sample metadata and RSE objects:

* [`01_build_objects.R`](03_Data_preparation/01_build_objects.R): 
    * **1\.  Data exploration & correction**: Sample metadata were added to `colData()` of the RSE objects in the correct format.
    * **2\.  Data normalization**: Raw gene expression counts were log-normalized to $log_2(CPM+0.5)$.
    * **3\.  Feature filtering**: Lowly-expressed genes were filtered out. Additional QC metrics and sample filtering info from [`04_EDA/01_QCA.R`](04_EDA/01_QCA.R) analysis were included in these original datasets.
    * **4\.   Visualization**: Count distributions were visualized before and after normalization and feature filtering steps.

    Note that exon, exon-exon junction, and transcript RSE objects were also processed but only gene-level data were used for downstream analyses.
