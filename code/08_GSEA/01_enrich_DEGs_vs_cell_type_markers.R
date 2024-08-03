
library(here)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scater)
library(DeconvoBuddies)
library(biomaRt)
library(rlang)
library(readxl)
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
library(randomcoloR)
library(sessioninfo)


####################   Gene Set Enrichment Analysis   ######################

## All expressed genes as universe
all_genes <- eval(parse_expr(load(here('processed-data/05_DEA/results_Substance_uncorr_vars_amygdala.Rdata'), verbose = TRUE)))[[1]]$ensemblID
length(all_genes)
# [1] 16708

## Load lists of DEGs
load(here('processed-data/05_DEA/de_genes_Substance_habenula.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/de_genes_Substance_amygdala.Rdata'), verbose = TRUE)


############################################################################
##           1. Obtain cell type marker genes in mouse habenula*
############################################################################
# *From doi: 10.1016/j.neuron.2020.03.011

## Load single cell mouse data (internal LIBD paths)

## Count data for all cell subpopulations
all_counts <- read.csv("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data/count.csv", row.names = 1)
## Count data for habenula neuronal subpopulations
hab_counts <- read.csv("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data/count_neuron.csv", row.names = 1)

## colData for all subpopulations
all_meta <- read.csv("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data/meta.csv")
## colData for hab subpopulations
hab_meta <- read.csv("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data/meta_neuron.csv")


## Explore count data for all subpops
dim(all_counts)
# [1] 17726 11878

all_counts[1:5, 1:5]
#        AAACCTGAGGCCCTCAcntl AAACCTGAGTTTGCGTcntl AAACCTGTCGTAGGTTcntl
# Xkr4                      0                    0                    0
# Sox17                     0                    0                    0
# Mrpl15                    0                    0                    0
# Lypla1                    0                    0                    0
# Tcea1                     0                    0                    0
#        AAACGGGCAGCTCGCAcntl AAACGGGTCGCTGATAcntl
# Xkr4                      0                    0
# Sox17                     0                    0
# Mrpl15                    0                    0
# Lypla1                    0                    2
# Tcea1                     0                    0

dim(all_meta)
# [1] 11878    10

head(all_meta)
#                      X orig.ident nCount_RNA nFeature_RNA stim percent.mito
# 1 AAACCTGAGGCCCTCAcntl    10X_LHb       1138          609 cntl   0.10456942
# 2 AAACCTGAGTTTGCGTcntl    10X_LHb       1358          678 cntl   0.03681885
# 3 AAACCTGTCGTAGGTTcntl    10X_LHb       1292          669 cntl   0.07120743
# 4 AAACGGGCAGCTCGCAcntl    10X_LHb       1001          594 cntl   0.04995005
# 5 AAACGGGTCGCTGATAcntl    10X_LHb        712          400 cntl   0.03089888
# 6 AAAGATGAGACCGGATcntl    10X_LHb        840          493 cntl   0.04404762
#   nCount_integrated nFeature_integrated integrated_snn_res.0.8   celltype
# 1                NA                  NA                      0 Astrocyte1
# 2                NA                  NA                      0 Astrocyte1
# 3                NA                  NA                      0 Astrocyte1
# 4                NA                  NA                      0 Astrocyte1
# 5                NA                  NA                      0 Astrocyte1
# 6                NA                  NA                      0 Astrocyte1


## Explore count data for habenula subpops
dim(hab_counts)
# [1] 17726  5558

hab_counts[1:5, 1:5]
#        AACACGTGTGGCAAACcntl AACCGCGGTAGGCATGcntl AACCGCGTCGACCAGCcntl
# Xkr4                      0                    0                    0
# Sox17                     0                    0                    0
# Mrpl15                    0                    0                    0
# Lypla1                    0                    0                    0
# Tcea1                     0                    0                    1
#        AACTCCCCAAAGCAATcntl AACTCCCGTCTAACGTcntl
# Xkr4                      0                    0
# Sox17                     0                    0
# Mrpl15                    0                    0
# Lypla1                    0                    0
# Tcea1                     0                    0

dim(hab_meta)
# [1] 5558   10

head(hab_meta)
#                      X orig.ident nCount_RNA nFeature_RNA stim percent.mito
# 1 AACACGTGTGGCAAACcntl    10X_LHb       3809         1897 cntl   0.04988186
# 2 AACCGCGGTAGGCATGcntl    10X_LHb       4058         1801 cntl   0.04706752
# 3 AACCGCGTCGACCAGCcntl    10X_LHb       3066         1580 cntl   0.05120678
# 4 AACTCCCCAAAGCAATcntl    10X_LHb       3206         1690 cntl   0.03867748
# 5 AACTCCCGTCTAACGTcntl    10X_LHb       4712         2171 cntl   0.05369270
# 6 AACTGGTAGACTTTCGcntl    10X_LHb       4119         1884 cntl   0.03957271
#   nCount_integrated nFeature_integrated integrated_snn_res.0.8 celltype
# 1                NA                  NA                      0     MHb1
# 2                NA                  NA                      0     MHb1
# 3                NA                  NA                      0     MHb1
# 4                NA                  NA                      0     MHb1
# 5                NA                  NA                      0     MHb1
# 6                NA                  NA                      0     MHb1


## Create SingleCellExperiment objects
sce_mouse_all <- SingleCellExperiment(
                    rowData=DataFrame(gene_name=rownames(all_counts)),
                    colData=DataFrame(all_meta),
                    assays = list(counts = as(all_counts, "sparseMatrix")))
sce_mouse_hab <- SingleCellExperiment(
                    rowData=DataFrame(gene_name=rownames(hab_counts)),
                    colData=DataFrame(hab_meta),
                    assays = list(counts = as(hab_counts, "sparseMatrix")))


## Log-transform and normalize expression values

## Compute library size factors per cell
lib.sf_mouse_all <- librarySizeFactors(sce_mouse_all)
lib.sf_mouse_hab <- librarySizeFactors(sce_mouse_hab)
length(lib.sf_mouse_all)
# [1] 11878
length(lib.sf_mouse_hab)
# [1] 5558

## Mean of size factors is 1
summary(lib.sf_mouse_all)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.2592  0.4171  0.7504  1.0000  1.4278  5.5241
summary(lib.sf_mouse_hab)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1912  0.5623  1.0002  1.0000  1.2895  4.0759

## Log-normalize counts in a separate assay: log2(gene count/cell size factor)
## Provided size factors are the same computed internally by logNormCounts()
sce_mouse_all <- logNormCounts(sce_mouse_all, size.factors=lib.sf_mouse_all)
assays(sce_mouse_all)$logcounts[1:5, 1:5]
#        AAACCTGAGGCCCTCAcntl AAACCTGAGTTTGCGTcntl AAACCTGTCGTAGGTTcntl
# Xkr4                      .                    .                    .
# Sox17                     .                    .                    .
# Mrpl15                    .                    .                    .
# Lypla1                    .                    .                    .
# Tcea1                     .                    .                    .
#        AAACGGGCAGCTCGCAcntl AAACGGGTCGCTGATAcntl
# Xkr4                      .             .
# Sox17                     .             .
# Mrpl15                    .             .
# Lypla1                    .             3.103901
# Tcea1                     .             .

sce_mouse_hab <- logNormCounts(sce_mouse_hab, size.factors=lib.sf_mouse_hab)
assays(sce_mouse_hab)$logcounts[1:5, 1:5]
#        AACACGTGTGGCAAACcntl AACCGCGGTAGGCATGcntl AACCGCGTCGACCAGCcntl
# Xkr4                      .                    .             .
# Sox17                     .                    .             .
# Mrpl15                    .                    .             .
# Lypla1                    .                    .             .
# Tcea1                     .                    .             1.134619
# AACTCCCCAAAGCAATcntl AACTCCCGTCTAACGTcntl
# Xkr4                      .                    .
# Sox17                     .                    .
# Mrpl15                    .                    .
# Lypla1                    .                    .
# Tcea1                     .                    .


## Number of cells of each type from naive and exposed mice:
##  All cell types
table(colData(sce_mouse_all)[, c('stim', 'celltype')])
#        celltype
# stim   Astrocyte1 Astrocyte2 Endothelial Epen Microglia Mural Neuron1 Neuron2
# cntl        969         72          93    6       156    94     429     541
# stim        637         40          69   36       147   152     631     502
#
# stim   Neuron3 Neuron4 Neuron5 Neuron6 Neuron7 Neuron8 Oligo1 Oligo2 Oligo3
# cntl     395     461     372     344     181     225    782    304     77
# stim     620     326     399     377     130      61    685    210     50

# stim   OPC1 OPC2 OPC3
# cntl      9  359   73
# stim    581  226   57

sum(table(colData(sce_mouse_all)[, c('stim', 'celltype')])['cntl',])
# [1] 5942

##  Cells from habenula subtypes
table(colData(sce_mouse_hab)[, c('stim', 'celltype')])
#        celltype
# stim   LHb1 LHb2 LHb3 LHb4 LHb5 LHb6 MHb1 MHb2 MHb3 MHb4 MHb5 MHb6
# cntl    329  279  217  214  150  174  315  264  270  148  165  142
# stim    229  178  222  185  210  185  351  398  351  264  229   89

sum(table(colData(sce_mouse_hab)[, c('stim', 'celltype')])['cntl',])
# [1] 2667


## Subset to naive mice samples
sce_mouse_all_ctrl <- sce_mouse_all[, which(sce_mouse_all$stim=='cntl')]
sce_mouse_hab_ctrl <- sce_mouse_hab[, which(sce_mouse_hab$stim=='cntl')]

## Remove cell groups with <10 cells
sce_mouse_all_ctrl_filt <- sce_mouse_all_ctrl[, which(!sce_mouse_all_ctrl$celltype %in% names(which(table(sce_mouse_all_ctrl$celltype)<10)))]
## None in habenula subpops
sce_mouse_hab_ctrl_filt <- sce_mouse_hab_ctrl[, which(!sce_mouse_hab_ctrl$celltype %in% names(which(table(sce_mouse_hab_ctrl$celltype)<10)))]


## -----------------------------------------------------------------------------
##                    A) MeanRatio cell type marker genes
## -----------------------------------------------------------------------------

########################  All habenula subpopulations  #########################
MeanRatio_all_hab_mouse_genes <- as.data.frame(get_mean_ratio(
                                                sce_mouse_all_ctrl_filt,
                                                cellType_col = "celltype",
                                                assay_name = "logcounts"))

## Subset to top100 markers per cell type
MeanRatio_top100_all_hab_mouse_genes <- subset(MeanRatio_all_hab_mouse_genes, MeanRatio.rank<=100)


######################  Habenula neuronal subpopulations  ######################
MeanRatio_neu_hab_mouse_genes <- as.data.frame(get_mean_ratio(
                                                sce_mouse_hab_ctrl_filt,
                                                cellType_col = "celltype",
                                                assay_name = "logcounts"))
## Top100 only
MeanRatio_top100_neu_hab_mouse_genes <- subset(MeanRatio_neu_hab_mouse_genes, MeanRatio.rank<=100)


## -----------------------------------------------------------------------------
##                     B) 1vsALL cell type marker genes
## -----------------------------------------------------------------------------

########################  All habenula subpopulations  #########################
lvsALL_all_hab_mouse_genes <- as.data.frame(findMarkers_1vAll(
                                                sce_mouse_all_ctrl_filt,
                                                assay_name = "logcounts",
                                                cellType_col = "celltype",
                                                mod = NULL,
                                                verbose = TRUE))

## Subset to DEGs with logFC>0 and FDR<0.05
lvsALL_all_hab_mouse_genes <- subset(lvsALL_all_hab_mouse_genes, logFC>0 & log.FDR<log(0.05))

## Number of marker genes (DEGs) per cell subpopulation
table(lvsALL_all_hab_mouse_genes$cellType.target)
# Astrocyte1  Astrocyte2 Endothelial   Microglia       Mural     Neuron1
#       1086         286         746         713         356        4788
# Neuron2     Neuron3     Neuron4     Neuron5     Neuron6     Neuron7
#    5663        2376        1497        3182        3320         279
# Neuron8      Oligo1      Oligo2      Oligo3        OPC2        OPC3
#    6665        1592         605         802         646         615


######################  Habenula neuronal subpopulations  ######################
lvsALL_neu_hab_mouse_genes <- as.data.frame(findMarkers_1vAll(
                                                sce_mouse_hab_ctrl_filt,
                                                assay_name = "logcounts",
                                                cellType_col = "celltype",
                                                mod = NULL,
                                                verbose = TRUE))

lvsALL_neu_hab_mouse_genes <- subset(lvsALL_neu_hab_mouse_genes, logFC>0 & log.FDR<log(0.05))

## Number of marker genes (DEGs) per cell subpopulation
table(lvsALL_neu_hab_mouse_genes$cellType.target)
# LHb1 LHb2 LHb3 LHb4 LHb5 LHb6 MHb1 MHb2 MHb3 MHb4 MHb5 MHb6
# 4332  354 1728  231  412 2246  903 1381 2951  264  453   85





############################################################################
##     2. Obtain sets of orthologs of human/mouse marker genes in rat
############################################################################

## Obtain rat orthologs of human marker genes
obtain_rat_orthologs_human <- function(human_marker_genes){
    mart <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl",
                       GRCh = "GRCh38")

    human_rat_ids <- getBM(values = human_marker_genes,
                           mart = mart,
                           attributes = c("external_gene_name",
                                          "rnorvegicus_homolog_ensembl_gene",
                                          "rnorvegicus_homolog_associated_gene_name"),
                           filters = "external_gene_name")
    return(human_rat_ids)
}

## Obtain rat orthologs of mouse marker genes
obtain_rat_orthologs_in_mouse <- function(mouse_marker_genes){
    mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                       dataset = "mmusculus_gene_ensembl",
                       ## Use Ensembl release 112 (GRCm39)
                       version = "112")

    mouse_rat_ids <- getBM(values = mouse_marker_genes,
                           mart = mart,
                           attributes = c("external_gene_name",
                                          "rnorvegicus_homolog_ensembl_gene",
                                          "rnorvegicus_homolog_associated_gene_name"),
                           filters = "external_gene_name")
    return(mouse_rat_ids)
}


## -----------------------------------------------------------------------------
##                A) MeanRatio-based cell type marker genes
## -----------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#               i)  Markers for cell types in human epithalamus*
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# *From doi: 10.1101/2024.02.26.582081

####################  Fine resolution cell type markers  ######################

## All ranked marker genes
MeanRatio_genes <- as.data.frame(read_xlsx('processed-data/08_GSEA/Input_cell_type_markers/MeanRatio_Top50_MarkerGenes_hab.xlsx'))
MeanRatio_top50_fine_hab_human_genes <- MeanRatio_genes

## Cell types/clusters included
cell_types <- names(table(MeanRatio_genes$cellType.target))
cell_types
# [1] "Astrocyte"  "Endo"       "Excit.Thal" "Inhib.Thal" "LHb.1"      "LHb.2"      "LHb.3"      "LHb.4"      "LHb.5"
# [10] "LHb.6"      "LHb.7"      "MHb.1"      "MHb.2"      "MHb.3"      "Microglia"  "Oligo"      "OPC"

## Confirm there are top 50 markers per cell type
table(MeanRatio_genes$cellType.target)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3      LHb.4      LHb.5      LHb.6      LHb.7
#        50         50         50         50         50         50         50         50         50         50         50
# MHb.1      MHb.2      MHb.3  Microglia      Oligo        OPC
#    50         50         50         50         50         50

## Range of expression ratios of marker genes per cell type
for (cell_type in cell_types){
    ratios <- subset(MeanRatio_genes, cellType.target==cell_type)$ratio
    print(paste0('Range of ratios of marker genes for ', cell_type, ': ', signif(min(ratios), 2), ' - ', signif(max(ratios), 2)))
}

# [1] "Range of ratios of marker genes for Astrocyte: 2.8 - 80"
# [1] "Range of ratios of marker genes for Endo: 2.7 - 160"
# [1] "Range of ratios of marker genes for Excit.Thal: 1.1 - 14"
# [1] "Range of ratios of marker genes for Inhib.Thal: 2.1 - 41"
# [1] "Range of ratios of marker genes for LHb.1: 1.2 - 2.8"
# [1] "Range of ratios of marker genes for LHb.2: 1.1 - 14"
# [1] "Range of ratios of marker genes for LHb.3: 1.2 - 4"
# [1] "Range of ratios of marker genes for LHb.4: 1.1 - 1.7"
# [1] "Range of ratios of marker genes for LHb.5: 1.2 - 2.2"
# [1] "Range of ratios of marker genes for LHb.6: 1.4 - 8.8"
# [1] "Range of ratios of marker genes for LHb.7: 1.1 - 3"
# [1] "Range of ratios of marker genes for MHb.1: 1.8 - 17"
# [1] "Range of ratios of marker genes for MHb.2: 1.2 - 5.6"
# [1] "Range of ratios of marker genes for MHb.3: 1.6 - 14"
# [1] "Range of ratios of marker genes for Microglia: 6.4 - 9.8"
# [1] "Range of ratios of marker genes for Oligo: 6.1 - 24"
# [1] "Range of ratios of marker genes for OPC: 1.2 - 16"


## Divide marker genes per cell type and obtain rat IDs
MeanRatio_top50_fine_hab_ratIDs <- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs_human(markers$Symbol)
    ## Take unique rat ensembl IDs: rat genes with at least one human ortholog marker gene
    markers_rat_IDs <- unique(markers_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]
    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))

    MeanRatio_top50_fine_hab_ratIDs[[cell_type]] <- markers_rat_IDs

}
# [1] "Number of Astrocyte marker genes in rat: 40"
# [1] "Number of Endo marker genes in rat: 54"
# [1] "Number of Excit.Thal marker genes in rat: 36"
# [1] "Number of Inhib.Thal marker genes in rat: 35"
# [1] "Number of LHb.1 marker genes in rat: 32"
# [1] "Number of LHb.2 marker genes in rat: 43"
# [1] "Number of LHb.3 marker genes in rat: 28"
# [1] "Number of LHb.4 marker genes in rat: 43"
# [1] "Number of LHb.5 marker genes in rat: 38"
# [1] "Number of LHb.6 marker genes in rat: 37"
# [1] "Number of LHb.7 marker genes in rat: 39"
# [1] "Number of MHb.1 marker genes in rat: 33"
# [1] "Number of MHb.2 marker genes in rat: 37"
# [1] "Number of MHb.3 marker genes in rat: 38"
# [1] "Number of Microglia marker genes in rat: 37"
# [1] "Number of Oligo marker genes in rat: 40"
# [1] "Number of OPC marker genes in rat: 47"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                ii)  Markers for cell types in human amygdala*
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# *From https://doi.org/10.1038/s41421-022-00506-y

####################  Broad resolution cell type markers  ######################

## All ranked marker genes
MeanRatio_genes <- as.data.frame(read.csv('processed-data/08_GSEA/Input_cell_type_markers/MeanRatio_Top100_broadMarkerGenes_amyg.csv'))
MeanRatio_top100_broad_amy_human_genes <- MeanRatio_genes

## Broad cell types
cell_types <- names(table(MeanRatio_genes$cellType.target))
cell_types
# [1] "Astrocyte"       "Endothelial"     "ExN"             "InN"             "Microglia"       "Oligodendrocyte" "OPC"

## Top 100 markers per cell type
table(MeanRatio_genes$cellType.target)
# Astrocyte     Endothelial             ExN             InN       Microglia Oligodendrocyte             OPC
#      100             100             100             100             100             100             100

## Marker ratios
for (cell_type in cell_types){
    ratios <- subset(MeanRatio_genes, cellType.target==cell_type)$ratio
    print(paste0('Range of ratios of marker genes for ', cell_type, ': ', signif(min(ratios), 2), ' - ', signif(max(ratios), 2)))
}

# [1] "Range of ratios of marker genes for Astrocyte: 1.7 - 3.7"
# [1] "Range of ratios of marker genes for Endothelial: 1.3 - 3.5"
# [1] "Range of ratios of marker genes for ExN: 1.7 - 9.1"
# [1] "Range of ratios of marker genes for InN: 1.2 - 7.8"
# [1] "Range of ratios of marker genes for Microglia: 2.2 - 11"
# [1] "Range of ratios of marker genes for Oligodendrocyte: 2 - 11"
# [1] "Range of ratios of marker genes for OPC: 1.3 - 36"


## Marker genes per cell type and obtain rat IDs
MeanRatio_top100_broad_amyg_ratIDs <- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs_human(markers$gene)
    markers_rat_IDs <- unique(markers_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]
    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))

    MeanRatio_top100_broad_amyg_ratIDs[[cell_type]] <- markers_rat_IDs

}
# [1] "Number of Astrocyte marker genes in rat: 77"
# [1] "Number of Endothelial marker genes in rat: 145"
# [1] "Number of ExN marker genes in rat: 60"
# [1] "Number of InN marker genes in rat: 89"
# [1] "Number of Microglia marker genes in rat: 98"
# [1] "Number of Oligodendrocyte marker genes in rat: 86"
# [1] "Number of OPC marker genes in rat: 91"


####################  Fine resolution cell type markers  ######################

MeanRatio_genes <- as.data.frame(read.csv('processed-data/08_GSEA/Input_cell_type_markers/MeanRatio_Top100_fineMarkerGenes_amyg.csv'))
MeanRatio_top100_fine_amy_human_genes <- MeanRatio_genes

## Fine cell types
cell_types <- names(table(MeanRatio_genes$cellType.target))
cell_types
# [1] "Human_Astro_1 FGFR3"  "Human_Astro_2 FGFR3"  "Human_Astro_3 FGFR3"  "Human_Astro_4 FGFR3"  "Human_CALCR LHX8"     "Human_DRD2 ISL1"
# [7] "Human_DRD2 PAX6"      "Human_Endo NOSTRIN"   "Human_HGF C11orf87"   "Human_HGF ESR1"       "Human_HGF NPSR1"      "Human_HTR3A DRD2"
# [13] "Human_LAMP5 ABO"      "Human_LAMP5 BDNF"     "Human_LAMP5 COL14A1"  "Human_LAMP5 COL25A1"  "Human_LAMP5 NDNF"     "Human_Micro CTSS"
# [19] "Human_Oligo_1 OPALIN" "Human_Oligo_2 OPALIN" "Human_Oligo_3 OPALIN" "Human_Oligo_4 OPALIN" "Human_Oligo_5 OPALIN" "Human_Oligo_6 OPALIN"
# [25] "Human_OPC_1 PDGFRA"   "Human_OPC_2 PDGFRA"   "Human_OPC_3 PDGFRA"   "Human_OPC_4 PDGFRA"   "Human_PRKCD"          "Human_PVALB ADAMTS5"
# [31] "Human_RXFP2 RSPO2"    "Human_SATB2 CALCRL"   "Human_SATB2 IL15"     "Human_SATB2 ST8SIA2"  "Human_SOX11 EBF2"     "Human_SST EPYC"
# [37] "Human_SST HGF"        "Human_STRIP2"         "Human_TFAP2C"         "Human_TSHZ1 CALCRL"   "Human_TSHZ1 SEMA3C"   "Human_VGLL3 CNGB1"
# [43] "Human_VGLL3 MEPE"     "Human_VIP ABI3BP"     "Human_VIP NDNF"


## Top 100 markers per cell type
unique(table(MeanRatio_genes$cellType.target))
# [1] 100

## Ratios
for (cell_type in cell_types){
    ratios <- subset(MeanRatio_genes, cellType.target==cell_type)$ratio
    print(paste0('Range of ratios of marker genes for ', cell_type, ': ', signif(min(ratios), 2), ' - ', signif(max(ratios), 2)))
}

#- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#     We don't want markers with ratios =<1 as they are      |
#        more expressed in the non-target cell type !!!      |
#- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# [1] "Range of ratios of marker genes for Human_Astro_1 FGFR3: 1.1 - 2.4"
# [1] "Range of ratios of marker genes for Human_Astro_2 FGFR3: 1.1 - 1.9"
# [1] "Range of ratios of marker genes for Human_Astro_3 FGFR3: 0.99 - 2.5"
# [1] "Range of ratios of marker genes for Human_Astro_4 FGFR3: 0.92 - 1.1"
# [1] "Range of ratios of marker genes for Human_CALCR LHX8: 0.94 - 1"
# [1] "Range of ratios of marker genes for Human_DRD2 ISL1: 1.1 - 11"
# [1] "Range of ratios of marker genes for Human_DRD2 PAX6: 1 - 2.3"
# [1] "Range of ratios of marker genes for Human_Endo NOSTRIN: 1.1 - 1.8"
# [1] "Range of ratios of marker genes for Human_HGF C11orf87: 1.1 - 1.3"
# [1] "Range of ratios of marker genes for Human_HGF ESR1: 0.98 - 1.4"
# [1] "Range of ratios of marker genes for Human_HGF NPSR1: 0.97 - 1.2"
# [1] "Range of ratios of marker genes for Human_HTR3A DRD2: 0.99 - 2.5"
# [1] "Range of ratios of marker genes for Human_LAMP5 ABO: 1 - 2.1"
# [1] "Range of ratios of marker genes for Human_LAMP5 BDNF: 1 - 1.3"
# [1] "Range of ratios of marker genes for Human_LAMP5 COL14A1: 1 - 3.4"
# [1] "Range of ratios of marker genes for Human_LAMP5 COL25A1: 0.98 - 1.1"
# [1] "Range of ratios of marker genes for Human_LAMP5 NDNF: 1 - 2.1"
# [1] "Range of ratios of marker genes for Human_Micro CTSS: 1.2 - 1.6"
# [1] "Range of ratios of marker genes for Human_Oligo_1 OPALIN: 0.97 - 1.2"
# [1] "Range of ratios of marker genes for Human_Oligo_2 OPALIN: 1 - 1.6"
# [1] "Range of ratios of marker genes for Human_Oligo_3 OPALIN: 0.92 - 1"
# [1] "Range of ratios of marker genes for Human_Oligo_4 OPALIN: 1 - 1.4"
# [1] "Range of ratios of marker genes for Human_Oligo_5 OPALIN: 1 - 1.4"
# [1] "Range of ratios of marker genes for Human_Oligo_6 OPALIN: 1 - 1.1"
# [1] "Range of ratios of marker genes for Human_OPC_1 PDGFRA: 1.3 - 4.3"
# [1] "Range of ratios of marker genes for Human_OPC_2 PDGFRA: 0.94 - 1.1"
# [1] "Range of ratios of marker genes for Human_OPC_3 PDGFRA: 1 - 1.3"
# [1] "Range of ratios of marker genes for Human_OPC_4 PDGFRA: 0.96 - 2.9"
# [1] "Range of ratios of marker genes for Human_PRKCD: 1.1 - 3.1"
# [1] "Range of ratios of marker genes for Human_PVALB ADAMTS5: 1.1 - 4.3"
# [1] "Range of ratios of marker genes for Human_RXFP2 RSPO2: 0.98 - 1.5"
# [1] "Range of ratios of marker genes for Human_SATB2 CALCRL: 0.99 - 1.4"
# [1] "Range of ratios of marker genes for Human_SATB2 IL15: 1 - 2.3"
# [1] "Range of ratios of marker genes for Human_SATB2 ST8SIA2: 1.1 - 3"
# [1] "Range of ratios of marker genes for Human_SOX11 EBF2: 0.92 - 1.1"
# [1] "Range of ratios of marker genes for Human_SST EPYC: 1 - 23"
# [1] "Range of ratios of marker genes for Human_SST HGF: 0.97 - 1.5"
# [1] "Range of ratios of marker genes for Human_STRIP2: 1.1 - 4.1"
# [1] "Range of ratios of marker genes for Human_TFAP2C: 1.2 - 26"
# [1] "Range of ratios of marker genes for Human_TSHZ1 CALCRL: 0.98 - 1.3"
# [1] "Range of ratios of marker genes for Human_TSHZ1 SEMA3C: 1 - 1.8"
# [1] "Range of ratios of marker genes for Human_VGLL3 CNGB1: 1 - 1.3"
# [1] "Range of ratios of marker genes for Human_VGLL3 MEPE: 0.94 - 1.3"
# [1] "Range of ratios of marker genes for Human_VIP ABI3BP: 0.98 - 3.2"
# [1] "Range of ratios of marker genes for Human_VIP NDNF: 0.98 - 9.6"

## Subset to markers with ratio>1
MeanRatio_top100_fine_amy_human_genes <- subset(MeanRatio_top100_fine_amy_human_genes, ratio>1)

## Marker genes per cell type and obtain rat IDs
MeanRatio_top100_fine_amyg_ratIDs <- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Remove markers with ratio =<1
    markers <- markers[markers$ratio>1, ]

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs_human(markers$gene)
    markers_rat_IDs <- unique(markers_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]
    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))

    MeanRatio_top100_fine_amyg_ratIDs[[cell_type]] <- markers_rat_IDs

}
# [1] "Number of Human_Astro_1 FGFR3 marker genes in rat: 85"
# [1] "Number of Human_Astro_2 FGFR3 marker genes in rat: 85"
# [1] "Number of Human_Astro_3 FGFR3 marker genes in rat: 72"
# [1] "Number of Human_Astro_4 FGFR3 marker genes in rat: 11"
# [1] "Number of Human_CALCR LHX8 marker genes in rat: 2"
# [1] "Number of Human_DRD2 ISL1 marker genes in rat: 76"
# [1] "Number of Human_DRD2 PAX6 marker genes in rat: 92"
# [1] "Number of Human_Endo NOSTRIN marker genes in rat: 147"
# [1] "Number of Human_HGF C11orf87 marker genes in rat: 84"
# [1] "Number of Human_HGF ESR1 marker genes in rat: 50"
# [1] "Number of Human_HGF NPSR1 marker genes in rat: 23"
# [1] "Number of Human_HTR3A DRD2 marker genes in rat: 73"
# [1] "Number of Human_LAMP5 ABO marker genes in rat: 67"
# [1] "Number of Human_LAMP5 BDNF marker genes in rat: 69"
# [1] "Number of Human_LAMP5 COL14A1 marker genes in rat: 86"
# [1] "Number of Human_LAMP5 COL25A1 marker genes in rat: 47"
# [1] "Number of Human_LAMP5 NDNF marker genes in rat: 97"
# [1] "Number of Human_Micro CTSS marker genes in rat: 97"
# [1] "Number of Human_Oligo_1 OPALIN marker genes in rat: 59"
# [1] "Number of Human_Oligo_2 OPALIN marker genes in rat: 92"
# [1] "Number of Human_Oligo_3 OPALIN marker genes in rat: 2"
# [1] "Number of Human_Oligo_4 OPALIN marker genes in rat: 92"
# [1] "Number of Human_Oligo_5 OPALIN marker genes in rat: 75"
# [1] "Number of Human_Oligo_6 OPALIN marker genes in rat: 93"
# [1] "Number of Human_OPC_1 PDGFRA marker genes in rat: 102"
# [1] "Number of Human_OPC_2 PDGFRA marker genes in rat: 19"
# [1] "Number of Human_OPC_3 PDGFRA marker genes in rat: 90"
# [1] "Number of Human_OPC_4 PDGFRA marker genes in rat: 43"
# [1] "Number of Human_PRKCD marker genes in rat: 89"
# [1] "Number of Human_PVALB ADAMTS5 marker genes in rat: 79"
# [1] "Number of Human_RXFP2 RSPO2 marker genes in rat: 35"
# [1] "Number of Human_SATB2 CALCRL marker genes in rat: 65"
# [1] "Number of Human_SATB2 IL15 marker genes in rat: 91"
# [1] "Number of Human_SATB2 ST8SIA2 marker genes in rat: 86"
# [1] "Number of Human_SOX11 EBF2 marker genes in rat: 3"
# [1] "Number of Human_SST EPYC marker genes in rat: 84"
# [1] "Number of Human_SST HGF marker genes in rat: 41"
# [1] "Number of Human_STRIP2 marker genes in rat: 79"
# [1] "Number of Human_TFAP2C marker genes in rat: 73"
# [1] "Number of Human_TSHZ1 CALCRL marker genes in rat: 44"
# [1] "Number of Human_TSHZ1 SEMA3C marker genes in rat: 83"
# [1] "Number of Human_VGLL3 CNGB1 marker genes in rat: 81"
# [1] "Number of Human_VGLL3 MEPE marker genes in rat: 2"
# [1] "Number of Human_VIP ABI3BP marker genes in rat: 59"
# [1] "Number of Human_VIP NDNF marker genes in rat: 60"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#             iii)  Markers for cell types in mouse habenula*
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# *From doi: 10.1016/j.neuron.2020.03.011

####################  All cell subpopulations markers  ######################

MeanRatio_top100_all_hab_mouse_genes$ratio <- MeanRatio_top100_all_hab_mouse_genes$MeanRatio
MeanRatio_genes <- MeanRatio_top100_all_hab_mouse_genes

## Cell types/clusters included
cell_types <- names(table(MeanRatio_genes$cellType.target))
cell_types
# [1] "Astrocyte1"  "Astrocyte2"  "Endothelial" "Microglia"   "Mural"
# [6] "Neuron1"     "Neuron2"     "Neuron3"     "Neuron4"     "Neuron5"
# [11] "Neuron6"     "Neuron7"     "Neuron8"     "Oligo1"      "Oligo2"
# [16] "Oligo3"      "OPC2"        "OPC3"

## Number of top markers per cell type
table(MeanRatio_genes$cellType.target)
# Astrocyte1  Astrocyte2 Endothelial   Microglia       Mural     Neuron1
#         78         100         100         100         100         100
# Neuron2     Neuron3     Neuron4     Neuron5     Neuron6     Neuron7
#     100         100         100         100         100         100
# Neuron8      Oligo1      Oligo2      Oligo3        OPC2        OPC3
#     100         100         100         100         100         100

## Range of expression ratios of marker genes per cell type
for (cell_type in cell_types){
    ratios <- subset(MeanRatio_genes, cellType.target==cell_type)$ratio
    print(paste0('Range of ratios of marker genes for ', cell_type, ': ', signif(min(ratios), 2), ' - ', signif(max(ratios), 2)))
}

# [1] "Range of ratios of marker genes for Astrocyte1: 0.19 - 1.5"
# [1] "Range of ratios of marker genes for Astrocyte2: 0.85 - 21"
# [1] "Range of ratios of marker genes for Endothelial: 0.88 - 430"
# [1] "Range of ratios of marker genes for Microglia: 1.1 - 390"
# [1] "Range of ratios of marker genes for Mural: 0.66 - 16"
# [1] "Range of ratios of marker genes for Neuron1: 1.3 - 4.1"
# [1] "Range of ratios of marker genes for Neuron2: 0.97 - 1.5"
# [1] "Range of ratios of marker genes for Neuron3: 0.97 - 2.1"
# [1] "Range of ratios of marker genes for Neuron4: 0.64 - 1.6"
# [1] "Range of ratios of marker genes for Neuron5: 0.91 - 2.6"
# [1] "Range of ratios of marker genes for Neuron6: 1 - 1.7"
# [1] "Range of ratios of marker genes for Neuron7: 0.35 - 1.3"
# [1] "Range of ratios of marker genes for Neuron8: 1.4 - 1.9"
# [1] "Range of ratios of marker genes for Oligo1: 0.92 - 1.6"
# [1] "Range of ratios of marker genes for Oligo2: 0.7 - 2"
# [1] "Range of ratios of marker genes for Oligo3: 1.2 - 3.7"
# [1] "Range of ratios of marker genes for OPC2: 0.45 - 5"
# [1] "Range of ratios of marker genes for OPC3: 1.3 - 18"

## Subset to markers with ratio>1
MeanRatio_top100_all_hab_mouse_genes <- subset(MeanRatio_top100_all_hab_mouse_genes, ratio>1)

## Number of markers per cell type with ratio>1
table(MeanRatio_top100_all_hab_mouse_genes$cellType.target)
# Astrocyte1  Astrocyte2 Endothelial   Microglia       Mural     Neuron1
#          7          80          76         100          27         100
# Neuron2     Neuron3     Neuron4     Neuron5     Neuron6     Neuron7
#      61          81          44          41         100          24
# Neuron8      Oligo1      Oligo2      Oligo3        OPC2        OPC3
#     100          63           3         100          22         100

## Obtain rat IDs per cell type
MeanRatio_top100_all_hab_ratIDs<- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Remove markers with ratio =<1
    markers <- markers[markers$ratio>1, ]

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs_in_mouse(markers$gene)
    ## Take unique rat ensembl IDs: rat genes with at least one mouse ortholog marker gene
    markers_rat_IDs <- unique(markers_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]
    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))

    MeanRatio_top100_all_hab_ratIDs[[cell_type]] <- markers_rat_IDs
}
# [1] "Number of Astrocyte1 marker genes in rat: 4"
# [1] "Number of Astrocyte2 marker genes in rat: 85"
# [1] "Number of Endothelial marker genes in rat: 78"
# [1] "Number of Microglia marker genes in rat: 123"
# [1] "Number of Mural marker genes in rat: 24"
# [1] "Number of Neuron1 marker genes in rat: 89"
# [1] "Number of Neuron2 marker genes in rat: 60"
# [1] "Number of Neuron3 marker genes in rat: 89"
# [1] "Number of Neuron4 marker genes in rat: 46"
# [1] "Number of Neuron5 marker genes in rat: 38"
# [1] "Number of Neuron6 marker genes in rat: 107"
# [1] "Number of Neuron7 marker genes in rat: 24"
# [1] "Number of Neuron8 marker genes in rat: 99"
# [1] "Number of Oligo1 marker genes in rat: 65"
# [1] "Number of Oligo2 marker genes in rat: 3"
# [1] "Number of Oligo3 marker genes in rat: 90"
# [1] "Number of OPC2 marker genes in rat: 22"
# [1] "Number of OPC3 marker genes in rat: 94"


####################  Neuronal cell subpopulations markers  ######################

MeanRatio_top100_neu_hab_mouse_genes$ratio <- MeanRatio_top100_neu_hab_mouse_genes$MeanRatio
MeanRatio_genes <- MeanRatio_top100_neu_hab_mouse_genes

## Cell types/clusters included
cell_types <- names(table(MeanRatio_genes$cellType.target))
cell_types
# [1] "LHb1" "LHb2" "LHb3" "LHb4" "LHb5" "LHb6" "MHb1" "MHb2" "MHb3" "MHb4"
# [11] "MHb5" "MHb6"

## Number of top markers per cell type
table(MeanRatio_genes$cellType.target)
# LHb1 LHb2 LHb3 LHb4 LHb5 LHb6 MHb1 MHb2 MHb3 MHb4 MHb5 MHb6
# 100  100  100  100  100  100  100  100  100  100  100  100

## Range of expression ratios of marker genes per cell type
for (cell_type in cell_types){
    ratios <- subset(MeanRatio_genes, cellType.target==cell_type)$ratio
    print(paste0('Range of ratios of marker genes for ', cell_type, ': ', signif(min(ratios), 2), ' - ', signif(max(ratios), 2)))
}
# [1] "Range of ratios of marker genes for LHb1: 1.2 - 2.3"
# [1] "Range of ratios of marker genes for LHb2: 0.5 - 1.9"
# [1] "Range of ratios of marker genes for LHb3: 1 - 2"
# [1] "Range of ratios of marker genes for LHb4: 0.67 - 1.3"
# [1] "Range of ratios of marker genes for LHb5: 0.98 - 1.7"
# [1] "Range of ratios of marker genes for LHb6: 1.2 - 2.3"
# [1] "Range of ratios of marker genes for MHb1: 1 - 1.9"
# [1] "Range of ratios of marker genes for MHb2: 0.99 - 1.4"
# [1] "Range of ratios of marker genes for MHb3: 1.3 - 4.2"
# [1] "Range of ratios of marker genes for MHb4: 0.9 - 3"
# [1] "Range of ratios of marker genes for MHb5: 1 - 5.3"
# [1] "Range of ratios of marker genes for MHb6: 0.39 - 1.5"

## Subset to markers with ratio>1
MeanRatio_top100_neu_hab_mouse_genes <- subset(MeanRatio_top100_neu_hab_mouse_genes, ratio>1)

## Number of markers per cell type with ratio>1
table(MeanRatio_top100_neu_hab_mouse_genes$cellType.target)
# LHb1 LHb2 LHb3 LHb4 LHb5 LHb6 MHb1 MHb2 MHb3 MHb4 MHb5 MHb6
#  100   34  100   22   65  100  100   74  100   49  100   28

## Obtain rat IDs per cell type
MeanRatio_top100_neu_hab_ratIDs<- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Remove markers with ratio =<1
    markers <- markers[markers$ratio>1, ]

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs_in_mouse(markers$gene)
    ## Take unique rat ensembl IDs: rat genes with at least one mouse ortholog marker gene
    markers_rat_IDs <- unique(markers_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]
    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))

    MeanRatio_top100_neu_hab_ratIDs[[cell_type]] <- markers_rat_IDs
}
# [1] "Number of LHb1 marker genes in rat: 96"
# [1] "Number of LHb2 marker genes in rat: 30"
# [1] "Number of LHb3 marker genes in rat: 99"
# [1] "Number of LHb4 marker genes in rat: 27"
# [1] "Number of LHb5 marker genes in rat: 63"
# [1] "Number of LHb6 marker genes in rat: 97"
# [1] "Number of MHb1 marker genes in rat: 86"
# [1] "Number of MHb2 marker genes in rat: 90"
# [1] "Number of MHb3 marker genes in rat: 95"
# [1] "Number of MHb4 marker genes in rat: 52"
# [1] "Number of MHb5 marker genes in rat: 121"
# [1] "Number of MHb6 marker genes in rat: 26"



## -----------------------------------------------------------------------------
##             B) 1vsALL-based cell type marker genes in human
## -----------------------------------------------------------------------------

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                i)  Markers for cell types in human epithalamus
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

####################  Broad resolution cell type markers  ######################
## (DEGs (FDR<0.05 and logFC>0) based on the enrichment model for one cell type vs the rest were taken as markers)

lvsALL_broad_genes <- eval(parse_expr(load(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broad_MarkerGenes_hab.Rdata'))))
lvsALL_broad_hab_human_genes <- lvsALL_broad_genes_enrich_stats <- lvsALL_broad_genes$enrichment

## Cell types
cell_types <- gsub('fdr_', '', colnames(lvsALL_broad_genes_enrich_stats)[grep('fdr', colnames(lvsALL_broad_genes_enrich_stats))])
cell_types
# [1] "Astrocyte"   "Endo"   "Excit.Thal"  "Inhib.Thal"  "LHb"   "MHb"   "Microglia"  "Oligo"   "OPC"

## Cell type-specific DEGs and corresponding orthologs in rat
lvsALL_broad_hab_ratIDs <- list()
for (cell_type in cell_types){

    ## Cell type data
    cell_type_stats <- lvsALL_broad_genes_enrich_stats[,  c(paste0(c('t_stat_', 'p_value_', 'fdr_', 'logFC_'), cell_type), 'gene')]
    ## DEGs symbols
    cell_type_DEGs <- subset(cell_type_stats, eval(parse_expr(paste0('fdr_', cell_type)))<0.05  &
                                              eval(parse_expr(paste0('logFC_', cell_type)))>0)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs_human(cell_type_DEGs)
    markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
    lvsALL_broad_hab_ratIDs[[cell_type]] <- markers_rat_IDs
}

# [1] "Number of Astrocyte human DEGs: 669"
# [1] "Number of Astrocyte marker genes in rat: 409"
# [1] "Number of Endo human DEGs: 427"
# [1] "Number of Endo marker genes in rat: 364"
# [1] "Number of Excit.Thal human DEGs: 304"
# [1] "Number of Excit.Thal marker genes in rat: 156"
# [1] "Number of Inhib.Thal human DEGs: 443"
# [1] "Number of Inhib.Thal marker genes in rat: 251"
# [1] "Number of LHb human DEGs: 22"
# [1] "Number of LHb marker genes in rat: 8"
# [1] "Number of MHb human DEGs: 77"
# [1] "Number of MHb marker genes in rat: 30"
# [1] "Number of Microglia human DEGs: 994"
# [1] "Number of Microglia marker genes in rat: 864"
# [1] "Number of Oligo human DEGs: 252"
# [1] "Number of Oligo marker genes in rat: 164"
# [1] "Number of OPC human DEGs: 196"
# [1] "Number of OPC marker genes in rat: 118"


####################  Fine resolution cell type markers  #######################

lvsALL_fine_genes <- eval(parse_expr(load(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_fine_MarkerGenes_hab.Rdata'))))
lvsALL_fine_hab_human_genes <- lvsALL_fine_genes_enrich_stats <- lvsALL_fine_genes$enrichment

## Cell types
cell_types <- gsub('fdr_', '', colnames(lvsALL_fine_genes_enrich_stats)[grep('fdr', colnames(lvsALL_fine_genes_enrich_stats))])
cell_types
# [1] "Astrocyte"  "Endo"       "Excit.Thal" "Inhib.Thal" "LHb.1"      "LHb.2"      "LHb.3"      "LHb.4"      "LHb.5"
# [10] "LHb.6"      "LHb.7"      "MHb.1"      "MHb.2"      "MHb.3"      "Microglia"  "Oligo"      "OPC"

lvsALL_fine_hab_ratIDs <- list()
for (cell_type in cell_types){

    ## Cell type data
    cell_type_stats <- lvsALL_fine_genes_enrich_stats[,  c(paste0(c('t_stat_', 'p_value_', 'fdr_', 'logFC_'), cell_type), 'gene')]
    ## DEGs symbols
    cell_type_DEGs <- subset(cell_type_stats, eval(parse_expr(paste0('fdr_', cell_type)))<0.05  &
                                              eval(parse_expr(paste0('logFC_', cell_type)))>0)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs (if DEGs exist)
    if (length(cell_type_DEGs)>0){
        cell_type_DEGs_rat_IDs <- obtain_rat_orthologs_human(cell_type_DEGs)
        markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
        markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

        print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
        lvsALL_fine_hab_ratIDs[[cell_type]] <- markers_rat_IDs
    }
}

# [1] "Number of Astrocyte human DEGs: 1193"
# [1] "Number of Astrocyte marker genes in rat: 801"
# [1] "Number of Endo human DEGs: 788"
# [1] "Number of Endo marker genes in rat: 621"
# [1] "Number of Excit.Thal human DEGs: 523"
# [1] "Number of Excit.Thal marker genes in rat: 296"
# [1] "Number of Inhib.Thal human DEGs: 1066"
# [1] "Number of Inhib.Thal marker genes in rat: 693"
# [1] "Number of LHb.1 human DEGs: 22"
# [1] "Number of LHb.1 marker genes in rat: 7"
# [1] "Number of LHb.2 human DEGs: 20"
# [1] "Number of LHb.2 marker genes in rat: 11"
# [1] "Number of LHb.3 human DEGs: 8"
# [1] "Number of LHb.3 marker genes in rat: 1"
# [1] "Number of LHb.4 human DEGs: 0"
# [1] "Number of LHb.5 human DEGs: 1"
# [1] "Number of LHb.5 marker genes in rat: 1"
# [1] "Number of LHb.6 human DEGs: 15"
# [1] "Number of LHb.6 marker genes in rat: 6"
# [1] "Number of LHb.7 human DEGs: 5"
# [1] "Number of LHb.7 marker genes in rat: 2"
# [1] "Number of MHb.1 human DEGs: 54"
# [1] "Number of MHb.1 marker genes in rat: 20"
# [1] "Number of MHb.2 human DEGs: 42"
# [1] "Number of MHb.2 marker genes in rat: 10"
# [1] "Number of MHb.3 human DEGs: 31"
# [1] "Number of MHb.3 marker genes in rat: 6"
# [1] "Number of Microglia human DEGs: 1757"
# [1] "Number of Microglia marker genes in rat: 1499"
# [1] "Number of Oligo human DEGs: 578"
# [1] "Number of Oligo marker genes in rat: 387"
# [1] "Number of OPC human DEGs: 440"
# [1] "Number of OPC marker genes in rat: 295"


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                  ii)  Markers for cell types in human amygdala
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

####################  Broad resolution cell type markers  ######################

## Cell type-specific enrichment stats for all genes
lvsALL_broad_genes <- read.csv('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broad_MarkerGenes_amyg.csv')
## (DEGs (FDR<0.05) based on the enrichment model as marker genes)

lvsALL_broad_amy_human_genes <- lvsALL_broad_genes

## Cell types
cell_types <- names(table(lvsALL_broad_genes$cellType.target))
cell_types
# [1] "Astrocyte"       "Endothelial"     "ExN"             "InN"             "Microglia"       "Oligodendrocyte" "OPC"

## Cell type-specific DEGs and corresponding orthologs in rat
lvsALL_broad_amy_ratIDs <- list()
for (cell_type in cell_types){

    ## Cell type data
    cell_type_stats <- subset(lvsALL_broad_genes, cellType.target==cell_type)
    ## DEGs
    cell_type_DEGs <- subset(cell_type_stats, log.FDR<log(0.05) & logFC>0)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs_human(cell_type_DEGs)
    markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
    lvsALL_broad_amy_ratIDs[[cell_type]] <- markers_rat_IDs
}

# [1] "Number of Astrocyte human DEGs: 6440"
# [1] "Number of Astrocyte marker genes in rat: 4105"
# [1] "Number of Endothelial human DEGs: 4784"
# [1] "Number of Endothelial marker genes in rat: 4349"
# [1] "Number of ExN human DEGs: 19235"
# [1] "Number of ExN marker genes in rat: 11378"
# [1] "Number of InN human DEGs: 11411"
# [1] "Number of InN marker genes in rat: 7889"
# [1] "Number of Microglia human DEGs: 4522"
# [1] "Number of Microglia marker genes in rat: 3615"
# [1] "Number of Oligodendrocyte human DEGs: 4306"
# [1] "Number of Oligodendrocyte marker genes in rat: 2945"
# [1] "Number of OPC human DEGs: 3125"
# [1] "Number of OPC marker genes in rat: 2223"


####################  Fine resolution cell type markers  ######################

lvsALL_fine_genes <- read.csv('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_fine_MarkerGenes_amyg.csv')
lvsALL_fine_amy_human_genes <- lvsALL_fine_genes

## Cell types
cell_types <- names(table(lvsALL_fine_genes$cellType.target))
cell_types
# [1] "Human_Astro_1 FGFR3"  "Human_Astro_2 FGFR3"  "Human_Astro_3 FGFR3"  "Human_Astro_4 FGFR3"  "Human_CALCR LHX8"     "Human_DRD2 ISL1"
# [7] "Human_DRD2 PAX6"      "Human_Endo NOSTRIN"   "Human_HGF C11orf87"   "Human_HGF ESR1"       "Human_HGF NPSR1"      "Human_HTR3A DRD2"
# [13] "Human_LAMP5 ABO"      "Human_LAMP5 BDNF"     "Human_LAMP5 COL14A1"  "Human_LAMP5 COL25A1"  "Human_LAMP5 NDNF"     "Human_Micro CTSS"
# [19] "Human_Oligo_1 OPALIN" "Human_Oligo_2 OPALIN" "Human_Oligo_3 OPALIN" "Human_Oligo_4 OPALIN" "Human_Oligo_5 OPALIN" "Human_Oligo_6 OPALIN"
# [25] "Human_OPC_1 PDGFRA"   "Human_OPC_2 PDGFRA"   "Human_OPC_3 PDGFRA"   "Human_OPC_4 PDGFRA"   "Human_PRKCD"          "Human_PVALB ADAMTS5"
# [31] "Human_RXFP2 RSPO2"    "Human_SATB2 CALCRL"   "Human_SATB2 IL15"     "Human_SATB2 ST8SIA2"  "Human_SOX11 EBF2"     "Human_SST EPYC"
# [37] "Human_SST HGF"        "Human_STRIP2"         "Human_TFAP2C"         "Human_TSHZ1 CALCRL"   "Human_TSHZ1 SEMA3C"   "Human_VGLL3 CNGB1"
# [43] "Human_VGLL3 MEPE"     "Human_VIP ABI3BP"     "Human_VIP NDNF"

## Cell type-specific DEGs and corresponding orthologs in rat
lvsALL_fine_amy_ratIDs <- list()
for (cell_type in cell_types){

    ## Cell type data
    cell_type_stats <- subset(lvsALL_fine_genes, cellType.target==cell_type)
    ## DEGs
    cell_type_DEGs <- subset(cell_type_stats, log.FDR<log(0.05) & logFC>0)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs_human(cell_type_DEGs)
    markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
    lvsALL_fine_amy_ratIDs[[cell_type]] <- markers_rat_IDs
}

# [1] "Number of Human_Astro_1 FGFR3 human DEGs: 4601"
# [1] "Number of Human_Astro_1 FGFR3 marker genes in rat: 3117"
# [1] "Number of Human_Astro_2 FGFR3 human DEGs: 5844"
# [1] "Number of Human_Astro_2 FGFR3 marker genes in rat: 3813"
# [1] "Number of Human_Astro_3 FGFR3 human DEGs: 4902"
# [1] "Number of Human_Astro_3 FGFR3 marker genes in rat: 3243"
# [1] "Number of Human_Astro_4 FGFR3 human DEGs: 5022"
# [1] "Number of Human_Astro_4 FGFR3 marker genes in rat: 3545"
# [1] "Number of Human_CALCR LHX8 human DEGs: 10327"
# [1] "Number of Human_CALCR LHX8 marker genes in rat: 7857"
# [1] "Number of Human_DRD2 ISL1 human DEGs: 7457"
# [1] "Number of Human_DRD2 ISL1 marker genes in rat: 5044"
# [1] "Number of Human_DRD2 PAX6 human DEGs: 6078"
# [1] "Number of Human_DRD2 PAX6 marker genes in rat: 4321"
# [1] "Number of Human_Endo NOSTRIN human DEGs: 4784"
# [1] "Number of Human_Endo NOSTRIN marker genes in rat: 4349"
# [1] "Number of Human_HGF C11orf87 human DEGs: 15423"
# [1] "Number of Human_HGF C11orf87 marker genes in rat: 10238"
# [1] "Number of Human_HGF ESR1 human DEGs: 7994"
# [1] "Number of Human_HGF ESR1 marker genes in rat: 5775"
# [1] "Number of Human_HGF NPSR1 human DEGs: 8327"
# [1] "Number of Human_HGF NPSR1 marker genes in rat: 5973"
# [1] "Number of Human_HTR3A DRD2 human DEGs: 5155"
# [1] "Number of Human_HTR3A DRD2 marker genes in rat: 3912"
# [1] "Number of Human_LAMP5 ABO human DEGs: 14766"
# [1] "Number of Human_LAMP5 ABO marker genes in rat: 9507"
# [1] "Number of Human_LAMP5 BDNF human DEGs: 15129"
# [1] "Number of Human_LAMP5 BDNF marker genes in rat: 9683"
# [1] "Number of Human_LAMP5 COL14A1 human DEGs: 5830"
# [1] "Number of Human_LAMP5 COL14A1 marker genes in rat: 4259"
# [1] "Number of Human_LAMP5 COL25A1 human DEGs: 7947"
# [1] "Number of Human_LAMP5 COL25A1 marker genes in rat: 5300"
# [1] "Number of Human_LAMP5 NDNF human DEGs: 5558"
# [1] "Number of Human_LAMP5 NDNF marker genes in rat: 4223"
# [1] "Number of Human_Micro CTSS human DEGs: 4522"
# [1] "Number of Human_Micro CTSS marker genes in rat: 3615"
# [1] "Number of Human_Oligo_1 OPALIN human DEGs: 6097"
# [1] "Number of Human_Oligo_1 OPALIN marker genes in rat: 4572"
# [1] "Number of Human_Oligo_2 OPALIN human DEGs: 1783"
# [1] "Number of Human_Oligo_2 OPALIN marker genes in rat: 1381"
# [1] "Number of Human_Oligo_3 OPALIN human DEGs: 3365"
# [1] "Number of Human_Oligo_3 OPALIN marker genes in rat: 2923"
# [1] "Number of Human_Oligo_4 OPALIN human DEGs: 2863"
# [1] "Number of Human_Oligo_4 OPALIN marker genes in rat: 1977"
# [1] "Number of Human_Oligo_5 OPALIN human DEGs: 4268"
# [1] "Number of Human_Oligo_5 OPALIN marker genes in rat: 2882"
# [1] "Number of Human_Oligo_6 OPALIN human DEGs: 2325"
# [1] "Number of Human_Oligo_6 OPALIN marker genes in rat: 1673"
# [1] "Number of Human_OPC_1 PDGFRA human DEGs: 918"
# [1] "Number of Human_OPC_1 PDGFRA marker genes in rat: 706"
# [1] "Number of Human_OPC_2 PDGFRA human DEGs: 3359"
# [1] "Number of Human_OPC_2 PDGFRA marker genes in rat: 2575"
# [1] "Number of Human_OPC_3 PDGFRA human DEGs: 2935"
# [1] "Number of Human_OPC_3 PDGFRA marker genes in rat: 2086"
# [1] "Number of Human_OPC_4 PDGFRA human DEGs: 2797"
# [1] "Number of Human_OPC_4 PDGFRA marker genes in rat: 2080"
# [1] "Number of Human_PRKCD human DEGs: 4780"
# [1] "Number of Human_PRKCD marker genes in rat: 3655"
# [1] "Number of Human_PVALB ADAMTS5 human DEGs: 5547"
# [1] "Number of Human_PVALB ADAMTS5 marker genes in rat: 4180"
# [1] "Number of Human_RXFP2 RSPO2 human DEGs: 14945"
# [1] "Number of Human_RXFP2 RSPO2 marker genes in rat: 9629"
# [1] "Number of Human_SATB2 CALCRL human DEGs: 12705"
# [1] "Number of Human_SATB2 CALCRL marker genes in rat: 8449"
# [1] "Number of Human_SATB2 IL15 human DEGs: 11643"
# [1] "Number of Human_SATB2 IL15 marker genes in rat: 8140"
# [1] "Number of Human_SATB2 ST8SIA2 human DEGs: 9954"
# [1] "Number of Human_SATB2 ST8SIA2 marker genes in rat: 7817"
# [1] "Number of Human_SOX11 EBF2 human DEGs: 6645"
# [1] "Number of Human_SOX11 EBF2 marker genes in rat: 4775"
# [1] "Number of Human_SST EPYC human DEGs: 3874"
# [1] "Number of Human_SST EPYC marker genes in rat: 2949"
# [1] "Number of Human_SST HGF human DEGs: 9730"
# [1] "Number of Human_SST HGF marker genes in rat: 7217"
# [1] "Number of Human_STRIP2 human DEGs: 9528"
# [1] "Number of Human_STRIP2 marker genes in rat: 7457"
# [1] "Number of Human_TFAP2C human DEGs: 7106"
# [1] "Number of Human_TFAP2C marker genes in rat: 5104"
# [1] "Number of Human_TSHZ1 CALCRL human DEGs: 6388"
# [1] "Number of Human_TSHZ1 CALCRL marker genes in rat: 4471"
# [1] "Number of Human_TSHZ1 SEMA3C human DEGs: 7266"
# [1] "Number of Human_TSHZ1 SEMA3C marker genes in rat: 5278"
# [1] "Number of Human_VGLL3 CNGB1 human DEGs: 16668"
# [1] "Number of Human_VGLL3 CNGB1 marker genes in rat: 10485"
# [1] "Number of Human_VGLL3 MEPE human DEGs: 6832"
# [1] "Number of Human_VGLL3 MEPE marker genes in rat: 4378"
# [1] "Number of Human_VIP ABI3BP human DEGs: 5165"
# [1] "Number of Human_VIP ABI3BP marker genes in rat: 3955"
# [1] "Number of Human_VIP NDNF human DEGs: 6878"
# [1] "Number of Human_VIP NDNF marker genes in rat: 5179"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#             iii)  Markers for cell types in mouse habenula
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

####################  Broad resolution cell type markers  ######################
## (DEGs (FDR<0.05 and logFC>0) based on the enrichment model for one cell type vs the rest were taken as markers)

lvsALL_broad_genes <- eval(parse_expr(load(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broad_MarkerGenes_hab.Rdata'))))
lvsALL_broad_hab_human_genes <- lvsALL_broad_genes_enrich_stats <- lvsALL_broad_genes$enrichment

## Cell types
cell_types <- gsub('fdr_', '', colnames(lvsALL_broad_genes_enrich_stats)[grep('fdr', colnames(lvsALL_broad_genes_enrich_stats))])
cell_types
# [1] "Astrocyte"   "Endo"   "Excit.Thal"  "Inhib.Thal"  "LHb"   "MHb"   "Microglia"  "Oligo"   "OPC"

## Cell type-specific DEGs and corresponding orthologs in rat
lvsALL_broad_hab_ratIDs <- list()
for (cell_type in cell_types){

    ## Cell type data
    cell_type_stats <- lvsALL_broad_genes_enrich_stats[,  c(paste0(c('t_stat_', 'p_value_', 'fdr_', 'logFC_'), cell_type), 'gene')]
    ## DEGs symbols
    cell_type_DEGs <- subset(cell_type_stats, eval(parse_expr(paste0('fdr_', cell_type)))<0.05  &
                                 eval(parse_expr(paste0('logFC_', cell_type)))>0)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs_human(cell_type_DEGs)
    markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
    lvsALL_broad_hab_ratIDs[[cell_type]] <- markers_rat_IDs
}







# Pasas 5
############################################################################
##             3. Compare MeanRatio vs 1vsALL marker genes
############################################################################

compare_markers <- function(region, species, resolution_MR, resolution_lvsALL){

    if(region== 'habenula' & species=='human'){
        top_n <- 'top50'
        region_name <- 'hab'
    }
    else if(region== 'amygdala' & species=='human'){
        top_n <- 'top100'
        region_name <- 'amy'
    }
    else if(region== 'habenula' & species=='mouse'){
        top_n <- 'top100'
        region_name <- 'hab'
    }

    ## Markers
    MR_markers <- eval(parse_expr(paste('MeanRatio', top_n, resolution_MR, region_name, species, 'genes', sep='_')))
    lvsALL_markers <- eval(parse_expr(paste('lvsALL', resolution_lvsALL, region_name, species, 'genes', sep='_')))

    ## Confirm all MeanRatios are >1
    stopifnot(length(which(MR_markers$ratio<=1))==0)

    ## Define cell types per set
    if(region == 'habenula' & species=='human'){
        cell_types_MR <- unique(MR_markers$cellType.target)
        cell_types_lvsALL <- gsub('fdr_', '', colnames(lvsALL_markers)[grep('fdr', colnames(lvsALL_markers))])
    }
    else if (region== 'amygdala' & species=='human'){
        cell_types_MR <- unique(MR_markers$cellType.target)
        cell_types_lvsALL <- unique(lvsALL_markers$cellType.target)
    }
    else if(region == 'habenula' & species=='mouse'){

    }

    ## Colors and alphas for plots
    colors <- randomColor(50, luminosity="dark")
    alphas <- c('TRUE'=1, 'FALSE'=0.3)

    ##_____________________________________________________________________________
    ##       1. Compare ratio vs standard logFC at same cell type resolution
    ##_____________________________________________________________________________
    if(resolution_MR == resolution_lvsALL){
        ## Verify same cell types at same resolution
        stopifnot(length(setdiff(cell_types_MR, cell_types_lvsALL))==0)

        ## Define needed column for plot
        if(resolution_MR=='fine'){
            MR_markers$cellType.target_fine <- MR_markers$cellType.target
        }
        else {
            MR_markers$cellType.target_fine <- NA
        }

        cell_type_colors_res <- ifelse(resolution_MR=='fine', 'cellType.target_fine', 'cellType.target')


        ## --------------------------------------------------
        ##  A) MeanRatio vs 1vsALL markers in human habenula
        ## --------------------------------------------------
        if(region=='habenula' & species=='human'){

            ########### i) Fine MeanRatio vs Fine 1vsALL cell types ###########
            if(resolution_MR=='fine' & resolution_lvsALL=='fine'){
                ## Use fine cell types for both methods
                cell_types <- cell_types_MR
                cell_type_MR <- "cell_type"
                cell_type_1vsALL <- "cell_type"

                ## Broad cell types: set to NA
                broad_cell_types <- NA
            }
        }

        ## --------------------------------------------------
        ##  B) MeanRatio vs 1vsALL markers in human amygdala
        ## --------------------------------------------------
        else if(region=='amygdala' & species=='human'){

            ########### i) Fine MeanRatio vs Fine 1vsALL cell types ###########
            if(resolution_MR=='fine' & resolution_lvsALL=='fine'){
                ## Use fine cell types for both methods
                cell_types <- cell_types_MR
                cell_type_MR <- "cell_type"
                cell_type_1vsALL <- "cell_type"

                ## Broad cell types: set to NA
                broad_cell_types <- NA
            }

            ########### ii) Broad MeanRatio vs Broad 1vsALL cell types ###########
            if(resolution_MR=='broad' & resolution_lvsALL=='broad'){
                ## Use broad cell types for both methods
                cell_types <- cell_types_MR
                cell_type_MR <- "cell_type"
                cell_type_1vsALL <- "cell_type"

                ## Broad cell types: set to NA (as we are not comparing against fine)
                broad_cell_types <- NA
            }
        }

        ## --------------------------------------------------
        ##  C) MeanRatio vs 1vsALL markers in mouse habenula
        ## --------------------------------------------------
        else if(region=='habenula' & species=='mouse'){

        }
    }


    ##_____________________________________________________________________________
    ##    2. Compare ratio vs standard logFC at different cell type resolutions
    ##_____________________________________________________________________________
    else{

        ## Use fine cell types for plotting
        cell_type_colors_res <- 'cellType.target_fine'

        ## --------------------------------------------------
        ##  A) MeanRatio vs 1vsALL markers in human habenula
        ## --------------------------------------------------
        if(region=='habenula' & species=='human'){

            ########### i) Fine MeanRatio vs Broad 1vsALL cell types ###########

            if(resolution_MR=='fine' & resolution_lvsALL=='broad'){

                MR_markers$cellType.target_fine <- MR_markers$cellType.target

                ## Broad cell type for each fine cell type
                broad_cell_types <- replace(replace(unique(MR_markers$cellType.target_fine),
                                                    grep('LHb.', unique(MR_markers$cellType.target_fine)), 'LHb'),
                                            grep('MHb.', unique(MR_markers$cellType.target_fine)), 'MHb')
                names(broad_cell_types) <- unique(MR_markers$cellType.target_fine)

                ## Verify fine cell types have a broad type in the broad data
                stopifnot(broad_cell_types %in% cell_types_lvsALL)

                ## Use fine cell types for plots
                cell_types <- cell_types_MR
                cell_type_MR <- "cell_type"
                ## Broad cell types for 1vsALL
                cell_type_1vsALL <- "broad_cell_type"

            }
        }

        ## --------------------------------------------------
        ##  B) MeanRatio vs 1vsALL markers in human amygdala
        ## --------------------------------------------------
        else if(region=='amygdala' & species=='human'){

            ########### i) Fine MeanRatio vs Broad 1vsALL cell types ###########

            if(resolution_MR=='fine' & resolution_lvsALL=='broad'){

                ## Fine MeanRatio markers
                MR_markers$cellType.target_fine <- MR_markers$cellType.target
                ## Broad cell types of the fine ones
                fine_and_broad_cell_types <- unique(MR_markers[,c('cellType.target_fine', 'broadCellType')])
                broad_cell_types <- fine_and_broad_cell_types$broadCellType
                names(broad_cell_types) <- fine_and_broad_cell_types$cellType.target_fine

                ## Verify fine cell types have a broad type in the broad data
                stopifnot(broad_cell_types %in% cell_types_lvsALL)

                ## Use fine cell types for plots
                cell_types <- cell_types_MR
                cell_type_MR <- "cell_type"
                ## Broad cell types for 1vsALL
                cell_type_1vsALL <- "broad_cell_type"
            }

            ########### ii) Broad MeanRatio vs Fine 1vsALL cell types ###########

            if(resolution_MR=='broad' & resolution_lvsALL=='fine'){

                ## Fine 1vsALL markers
                lvsALL_markers$cellType.target_fine <- lvsALL_markers$cellType.target

                ## Broad MeanRatio markers
                MR_markers$cellType.target_fine <- NA
                cell_type_colors_res <- 'cellType.target'

                ## Broad cell types of the fine ones
                fine_and_broad_cell_types <- unique(lvsALL_markers[,c('cellType.target_fine', 'broadCellType')])
                broad_cell_types <- fine_and_broad_cell_types$broadCellType
                names(broad_cell_types) <- fine_and_broad_cell_types$cellType.target_fine

                ## Verify fine cell types have a broad type in the broad data
                stopifnot(broad_cell_types %in% cell_types_MR)

                ## Use fine cell types for plots
                cell_types <- cell_types_lvsALL
                cell_type_1vsALL <- "cell_type"
                ## Broad cell types for MeanRatio
                cell_type_MR <- "broad_cell_type"
            }

        }

        ## --------------------------------------------------
        ##  C) MeanRatio vs 1vsALL markers in mouse habenula
        ## --------------------------------------------------
        else if(region=='habenula' & species=='mouse'){

        }

    }

    ## Plot per cell type
    plots <- list()
    for (i in 1:length(cell_types)){
        cell_type <- cell_types[i]
        ## Only to compare fine vs broad
        broad_cell_type <- broad_cell_types[cell_type]

        ## Genes' MeanRatio for the cell type
        MR_markers_cell_type <- subset(MR_markers, get(cell_type_colors_res)==get(cell_type_MR))[,c("gene", "ratio",
                                                                                                    "cellType.target",
                                                                                                    "cellType.target_fine")]

        ## Genes' std logFC and p-val for the cell type
        if(region=='habenula' & species=='human'){
            data <- cbind(MR_markers_cell_type, lvsALL_markers[match(MR_markers_cell_type$gene, lvsALL_markers$gene),
                                                               c(paste0('logFC_', get(cell_type_1vsALL)),
                                                                 paste0('fdr_', get(cell_type_1vsALL)),
                                                                 paste0('t_stat_', get(cell_type_1vsALL)))])
            data$logFC <- data[,paste0('logFC_', get(cell_type_1vsALL))]
            data[,paste0('logFC_', get(cell_type_1vsALL))] <- NULL
            data$FDR <- data[, paste0('fdr_', get(cell_type_1vsALL))]
            data[, paste0('fdr_', get(cell_type_1vsALL))] <- NULL
            data$t_stat <- data[, paste0('t_stat_', get(cell_type_1vsALL))]
            data[, paste0('t_stat_', get(cell_type_1vsALL))] <- NULL

            ## Discard NAs (MeanRatio markers without 1vsALL metrics available)
            if(length(which(is.na(data$t_stat)))>0){
                data <- data[-which(is.na(data$t_stat)), ]
            }

            ## Add if gene is MeanRatio marker and also lvsALL marker
            data$lvsALL_marker <- apply(data, 1, function(x){if(as.numeric(x['FDR'])<0.05 &
                                                                as.numeric(x['logFC'])>0 &
                                                                as.numeric(x['ratio'])>1){TRUE}
                                                             else{FALSE}})
        }

        else if(region=='amygdala' & species=='human'){
            lvsALL_data <- subset(lvsALL_markers, cellType.target==get(cell_type_1vsALL))
            data <- cbind(MR_markers_cell_type, lvsALL_data[match(MR_markers_cell_type$gene, lvsALL_data$gene),
                                                               c('logFC', 'log.FDR', 'std.logFC')])
            data$t_stat <- data$std.logFC

            ## Discard NAs (MeanRatio markers without 1vsALL metrics available)
            if(length(which(is.na(data$t_stat)))>0){
                data <- data[-which(is.na(data$t_stat)), ]
            }

            ## Add same info of common markers
            data$lvsALL_marker <- apply(data, 1, function(x){if(as.numeric(x['log.FDR'])<log(0.05) &
                                                                as.numeric(x['logFC'])>0 &
                                                                as.numeric(x['ratio'])>1){TRUE}
                                                             else{FALSE}})
        }

        ## % of MeanRatio markers that are also 1vsALL markers
        num_common <- ifelse(!is.na(table(data$lvsALL_marker)['TRUE']), table(data$lvsALL_marker)['TRUE'], 0)
        percent <- signif(num_common / dim(data)[1] *100, 4)

        ## Plot ratio of MeanRatio markers per cell type vs their std logFC for same/equivalent cell type
        plots[[i]] <- ggplot(data, aes(x=ratio, y=t_stat,
                                       color=get(cell_type_colors_res),
                                       alpha=lvsALL_marker)) +
            geom_point(size = 2, color=colors[i]) +
            scale_alpha_manual(values = alphas) +
            theme_bw() +
            labs(x=paste0("Ratio of MeanRatio marker genes for ", get(cell_type_MR)),
                 y=paste0("1vsALL standard logFC for ", get(cell_type_1vsALL)),
                 color= "Cell type",
                 alpha="MeanRatio & 1vsALL marker",
                 subtitle = paste0(percent , '% of ', dim(data)[1],
                                   ' MeanRatio markers are 1vsALL markers')) +
            theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"),
                  axis.title = element_text(size = 11),
                  axis.text = element_text(size = (10)),
                  legend.text = element_text(size=8.4),
                  legend.title = element_text(size=10))
    }

    return(plots)
}


## -----------------------------------------------------------------------------
##              A) MeanRatio vs 1vsALL markers in human habenula
## -----------------------------------------------------------------------------

####################  Fine resolution cell type markers  #######################

## Compare MeanRatio vs 1vsALL std logFC for fine cell types
region <- 'habenula'
species <- 'human'
resolution_MR <- 'fine'
resolution_lvsALL <- 'fine'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=5)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 24, height = 20)

###############  Fine and broad resolution cell type markers  ##################

## Compare MeanRatio for fine cell types vs 1vsALL std logFC for the respective broad cell types
region <- 'habenula'
species <- 'human'
resolution_MR <- 'fine'
resolution_lvsALL <- 'broad'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=5)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 20, height = 15)


## -----------------------------------------------------------------------------
##              B) MeanRatio vs 1vsALL markers in human amygdala
## -----------------------------------------------------------------------------

####################  Fine resolution cell type markers  #######################

## Compare MeanRatio vs 1vsALL std logFC for fine cell types
region <- 'amygdala'
species <- 'human'
resolution_MR <- 'fine'
resolution_lvsALL <- 'fine'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=5)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 65, height = 26, limitsize = FALSE)

####################  Broad resolution cell type markers  ######################

## Compare MeanRatio vs 1vsALL std logFC for broad cell types
region <- 'amygdala'
species <- 'human'
resolution_MR <- 'broad'
resolution_lvsALL <- 'broad'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=2)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 24, height = 8)

###############  Fine and broad resolution cell type markers  ##################

## Compare MeanRatio for fine cell types vs 1vsALL std logFC for the respective broad cell types
region <- 'amygdala'
species <- 'human'
resolution_MR <- 'fine'
resolution_lvsALL <- 'broad'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=5)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 65, height = 26, limitsize = FALSE)

## Compare 1vsALL std logFC for fine cell types vs MeanRatio for the respective broad cell types
region <- 'amygdala'
species <- 'human'
resolution_MR <- 'broad'
resolution_lvsALL <- 'fine'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=5)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 65, height = 26, limitsize = FALSE)


## -----------------------------------------------------------------------------
##              C) MeanRatio vs 1vsALL markers in mouse habenula
## -----------------------------------------------------------------------------

####################  All cell type markers  #######################

## Compare MeanRatio vs 1vsALL std logFC for all cell types
region <- 'habenula'
species <- 'mouse'
resolution_MR <- 'all'
resolution_lvsALL <- 'all'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=5)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 65, height = 26, limitsize = FALSE)

####################  Habenula cell type markers  ######################

## Compare MeanRatio vs 1vsALL std logFC for habenula cell types
region <- 'habenula'
species <- 'mouse'
resolution_MR <- 'hab'
resolution_lvsALL <- 'hab'
p <- compare_markers(region, species, resolution_MR, resolution_lvsALL)
plot_grid(plotlist = p, nrow=2)
ggsave(filename = paste0('plots/08_GSEA/MeanRatio_', resolution_MR, '_vs_lvsALL_',
                         resolution_lvsALL, '_', region, '_', species, '.pdf'), width = 24, height = 8)





############################################################################
##             4. Cell type enrichment analysis for rat DEGs
############################################################################

## Enrichment analysis
enrichment_analysis<- function(region, method, top_n, resolution, DEGs_region){

    ## Define marker genes
    if(method == "lvsALL"){
        markers <- eval(parse_expr(paste(method, resolution, region, 'ratIDs', sep='_')))
    }
    else{
        markers <- eval(parse_expr(paste(method, top_n, resolution, region, 'ratIDs', sep='_')))
    }

    ## Define groups of rat DEGs
    de_genes <- eval(parse_expr(paste0('de_genes_', DEGs_region)))
    all_DEGs <- de_genes$ensemblID
    up_DEGs <- subset(de_genes, logFC>0)$ensemblID
    down_DEGs <- subset(de_genes, logFC<0)$ensemblID

    DEGs_groups <- list("all_DEGs"=all_DEGs, "up_DEGs"=up_DEGs, "down_DEGs"=down_DEGs)

    ## Create contingency table and perform Fisher test for cell type and DEGs group
    p_values <- data.frame(matrix(nrow=length(DEGs_groups), ncol=length(names(markers))))
    colnames(p_values) <- names(markers)
    rownames(p_values) <-  names(DEGs_groups)

    ## Individual matrices
    ms <- list()

    for (group in names(DEGs_groups)){
        for (cell_type in names(markers)){

            ## Universe = DEGs + non-DEGs
            ## * DEGs
            de_genes <- DEGs_groups[[group]]
            ## * Non-DE genes
            non_de_genes <- all_genes[! all_genes %in% de_genes]

            ## Marker genes for the cell type
            cell_type_markers <-  markers[[cell_type]]
            ## Subset to marker genes in universe
            cell_type_markers <- cell_type_markers[cell_type_markers %in% all_genes]

            ## Intersections
            ## * Cell type markers and DEGs
            DEGs_and_cellTypeMarkers <- intersect(de_genes, cell_type_markers)
            ## * Cell type markers among non-DEGs
            nonDEGs_and_cellTypeMarkers <- intersect(non_de_genes, cell_type_markers)
            ## * Rest of DEGs that are not markers
            DEGs_and_no_cellTypeMarkers <- de_genes[! de_genes %in% DEGs_and_cellTypeMarkers]
            ## * Rest of non-DEGs that are not markers
            nonDEGs_and_no_cellTypeMarkers <- non_de_genes[! non_de_genes %in% nonDEGs_and_cellTypeMarkers]


            ## Contingency table
            m <- matrix(c(length(DEGs_and_cellTypeMarkers), length(nonDEGs_and_cellTypeMarkers),
                          length(DEGs_and_no_cellTypeMarkers), length(nonDEGs_and_no_cellTypeMarkers)), ncol=2)
            rownames(m) <- c('DEGs', 'non-DEGs')
            colnames(m) <- c('Markers', 'non-Markers')

            ms[[paste(group, cell_type, sep='-')]] <- m

            ## Confirm number of DEGs, cell type markers and universe size
            if (sum(m[1,])==length(de_genes) & sum(m[,1])==length(cell_type_markers) & sum(m)==length(all_genes)){
                ## Assess enrichment of markers among DEGs
                p <- fisher.test(m, alternative = "greater")$p.value
                p_values[group, cell_type] <- p
            }
                else {
                print('error')
            }

        }

    }

    return(list(p_values, ms))
}


## Create heatmaps for results
heatmap_pvals <- function(p_values, DEGs_region, marker_set_name, filename, width){

    ## Heatmap for -log(p-values)
    log_p_values <- -log10(p_values)

    ## Num of DEGs in each group
    de_genes <- eval(parse_expr(paste0('de_genes_', DEGs_region)))
    all_DEGs <- de_genes$ensemblID
    up_DEGs <- subset(de_genes, logFC>0)$ensemblID
    down_DEGs <- subset(de_genes, logFC<0)$ensemblID
    num_DEGs_group <- data.frame('Freq'=c('all_DEGs'=length(all_DEGs), 'up_DEGs'=length(up_DEGs), 'down_DEGs'=length(down_DEGs)))

    ## Num of marker genes for each cell type
    markers <- eval(parse_expr(paste0(filename, '_ratIDs')))
    num_markers_cell_type <- data.frame('Freq'=unlist(lapply(markers, length)))

    ## Row and col bars
    row_gene_anno <- ComplexHeatmap::rowAnnotation('n genes' = ComplexHeatmap::anno_barplot(num_DEGs_group$Freq))
    col_gene_anno <- ComplexHeatmap::columnAnnotation('n markers' = ComplexHeatmap::anno_barplot(num_markers_cell_type$Freq))

    h <- Heatmap(log_p_values,
                 name='-log10(p-value)',
                 col= colorRampPalette(c('aliceblue', 'dodgerblue3'))(50),
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 border_gp = gpar(col = "gray20", lty = 1),
                 rect_gp = gpar(col = "gray20", lwd = 1),
                 column_title = marker_set_name,
                 row_title = paste0('Groups of ', DEGs_region, ' DEGs'),
                 column_title_gp = gpar(fontsize = 11),
                 row_title_gp = gpar(fontsize = 11),
                 column_names_gp = gpar(fontsize = 9),
                 row_names_gp = gpar(fontsize = 9),
                 right_annotation = row_gene_anno,
                 top_annotation = col_gene_anno,
                 ## Add '*' if p<0.05
                 cell_fun = function(j, i, x, y,  w, h, col)
                 { if(log_p_values[i,j]>(-log10(0.05))){ grid.text('*', x, y, gp = gpar(fontsize = 17, col='yellow1'))} }
    )

    pdf(file=paste0('plots/08_GSEA/', filename, '_vs_', DEGs_region, 'DEGs.pdf'), height = 5, width = width)
    print(h)
    dev.off()
}


####  Compare amygdala DEGs vs cell type marker genes in amygdala:

#   * MeanRatio-based cell type marker genes at fine resolution
results_MeanRatio_100_fine_amyg <- enrichment_analysis('amyg', 'MeanRatio', 'top100', 'fine', 'amygdala')
p_values_MeanRatio_100_fine_amyg <- results_MeanRatio_100_fine_amyg[[1]]
heatmap_pvals(p_values_MeanRatio_100_fine_amyg, 'amygdala', 'Top100 MeanRatio-based amygdala fine cell type markers',
              'MeanRatio_top100_fine_amyg', 11)
ms_MeanRatio_100_fine_amyg <- results_MeanRatio_100_fine_amyg[[2]]

#   * MeanRatio-based cell type marker genes at broad resolution
results_MeanRatio_100_broad_amyg <- enrichment_analysis('amyg', 'MeanRatio', 'top100', 'broad', 'amygdala')
p_values_MeanRatio_100_broad_amyg <- results_MeanRatio_100_broad_amyg[[1]]
heatmap_pvals(p_values_MeanRatio_100_broad_amyg, 'amygdala', 'Top100 MeanRatio-based amygdala broad cell type markers',
              'MeanRatio_top100_broad_amyg', 6.7)
ms_MeanRatio_100_broad_amyg <- results_MeanRatio_100_broad_amyg[[2]]

#   * 1vsALL-based cell type marker genes at fine resolution
results_lvsALL_fine_amyg <- enrichment_analysis('amy', 'lvsALL', NULL, 'fine', 'amygdala')
p_values_lvsALL_fine_amyg <- results_lvsALL_fine_amyg[[1]]
heatmap_pvals(p_values_lvsALL_fine_amyg, 'amygdala', '1vsALL-based amygdala fine cell type markers',
              'lvsALL_fine_amy', 11)
ms_lvsALL_fine_amyg <- results_lvsALL_fine_amyg[[2]]

#   * 1vsALL-based cell type marker genes at broad resolution
results_lvsALL_broad_amyg <- enrichment_analysis('amy', 'lvsALL', NULL, 'broad', 'amygdala')
p_values_lvsALL_broad_amyg <- results_lvsALL_broad_amyg[[1]]
heatmap_pvals(p_values_lvsALL_broad_amyg, 'amygdala', '1vsALL-based amygdala broad cell type markers',
              'lvsALL_broad_amy', 6.7)
ms_lvsALL_broad_amyg <- results_lvsALL_broad_amyg[[2]]



####  Compare habenula DEGs vs cell type marker genes in epithalamus:

#   * MeanRatio-based cell type marker genes at fine resolution
results_MeanRatio_50_fine_hab <- enrichment_analysis('hab', 'MeanRatio', 'top50', 'fine', 'habenula')
p_values_MeanRatio_50_fine_hab <- results_MeanRatio_50_fine_hab[[1]]
heatmap_pvals(p_values_MeanRatio_50_fine_hab, 'habenula', 'Top50 MeanRatio-based habenula fine cell type markers',
              'MeanRatio_top50_fine_hab', 8)
ms_MeanRatio_50_fine_hab <- results_MeanRatio_50_fine_hab[[2]]

#   * 1vsALL-based cell type marker genes at fine resolution
results_lvsALL_fine_hab <- enrichment_analysis('hab', 'lvsALL', NULL, 'fine', 'habenula')
p_values_lvsALL_fine_hab <- results_lvsALL_fine_hab[[1]]
heatmap_pvals(p_values_lvsALL_fine_hab, 'habenula', '1vsALL-based habenula fine cell type markers',
              'lvsALL_fine_hab', 8)
ms_lvsALL_fine_hab <- results_lvsALL_fine_hab[[2]]

#   * 1vsALL-based cell type marker genes at broad resolution
results_lvsALL_broad_hab <- enrichment_analysis('hab', 'lvsALL', NULL, 'broad', 'habenula')
p_values_lvsALL_broad_hab <- results_lvsALL_broad_hab[[1]]
heatmap_pvals(p_values_lvsALL_broad_hab, 'habenula', '1vsALL-based habenula broad cell type markers',
              'lvsALL_broad_hab', 6)
ms_lvsALL_broad_hab <- results_lvsALL_broad_hab[[2]]







## Reproducibility information

options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 (2024-04-24)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2024-06-20
# rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
# pandoc   3.1.11 @ /private/var/folders/r6/94s0dsks4m3298b0d1mrp5tw0000gn/T/AppTranslocation/73D2BAF6-C47F-4B84-B5C9-D80EC60F9834/d/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version date (UTC) lib source
# abind                  1.4-5   2016-07-21 [1] CRAN (R 4.4.0)
# AnnotationDbi          1.67.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# Biobase              * 2.65.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# BiocFileCache          2.13.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# BiocGenerics         * 0.51.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# biomaRt              * 2.61.1  2024-06-12 [1] Bioconductor 3.20 (R 4.4.1)
# Biostrings             2.73.1  2024-06-02 [1] Bioconductor 3.20 (R 4.4.0)
# bit                    4.0.5   2022-11-15 [1] CRAN (R 4.4.0)
# bit64                  4.0.5   2020-08-30 [1] CRAN (R 4.4.0)
# blob                   1.2.4   2023-03-17 [1] CRAN (R 4.4.0)
# cachem                 1.1.0   2024-05-16 [1] CRAN (R 4.4.0)
# cellranger             1.1.0   2016-07-27 [1] CRAN (R 4.4.0)
# circlize               0.4.16  2024-02-20 [1] CRAN (R 4.4.0)
# cli                    3.6.2   2023-12-11 [1] CRAN (R 4.4.0)
# clue                   0.3-65  2023-09-23 [1] CRAN (R 4.4.0)
# cluster                2.1.6   2023-12-01 [1] CRAN (R 4.4.0)
# codetools              0.2-20  2024-03-31 [1] CRAN (R 4.4.0)
# colorspace             2.1-0   2023-01-23 [1] CRAN (R 4.4.0)
# ComplexHeatmap       * 2.21.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# crayon                 1.5.2   2022-09-29 [1] CRAN (R 4.4.0)
# curl                   5.2.1   2024-03-01 [1] CRAN (R 4.4.0)
# DBI                    1.2.3   2024-06-02 [1] CRAN (R 4.4.0)
# dbplyr                 2.5.0   2024-03-19 [1] CRAN (R 4.4.0)
# DelayedArray           0.31.1  2024-05-07 [1] Bioconductor 3.20 (R 4.4.0)
# digest                 0.6.35  2024-03-11 [1] CRAN (R 4.4.0)
# doParallel             1.0.17  2022-02-07 [1] CRAN (R 4.4.0)
# dplyr                  1.1.4   2023-11-17 [1] CRAN (R 4.4.0)
# evaluate               0.24.0  2024-06-10 [1] CRAN (R 4.4.0)
# fansi                  1.0.6   2023-12-08 [1] CRAN (R 4.4.0)
# fastmap                1.2.0   2024-05-15 [1] CRAN (R 4.4.0)
# filelock               1.0.3   2023-12-11 [1] CRAN (R 4.4.0)
# foreach                1.5.2   2022-02-02 [1] CRAN (R 4.4.0)
# generics               0.1.3   2022-07-05 [1] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.41.1  2024-05-24 [1] Bioconductor 3.20 (R 4.4.0)
# GenomeInfoDbData       1.2.12  2024-06-12 [1] Bioconductor
# GenomicRanges        * 1.57.1  2024-06-12 [1] Bioconductor 3.20 (R 4.4.1)
# GetoptLong             1.0.5   2020-12-15 [1] CRAN (R 4.4.0)
# GlobalOptions          0.1.2   2020-06-10 [1] CRAN (R 4.4.0)
# glue                   1.7.0   2024-01-09 [1] CRAN (R 4.4.0)
# here                 * 1.0.1   2020-12-13 [1] CRAN (R 4.4.0)
# hms                    1.1.3   2023-03-21 [1] CRAN (R 4.4.0)
# htmltools              0.5.8.1 2024-04-04 [1] CRAN (R 4.4.0)
# httr                   1.4.7   2023-08-15 [1] CRAN (R 4.4.0)
# httr2                  1.0.1   2024-04-01 [1] CRAN (R 4.4.0)
# IRanges              * 2.39.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# iterators              1.0.14  2022-02-05 [1] CRAN (R 4.4.0)
# jsonlite               1.8.8   2023-12-04 [1] CRAN (R 4.4.0)
# KEGGREST               1.45.1  2024-06-18 [1] Bioconductor 3.20 (R 4.4.1)
# knitr                  1.47    2024-05-29 [1] CRAN (R 4.4.0)
# lattice                0.22-6  2024-03-20 [1] CRAN (R 4.4.0)
# lifecycle              1.0.4   2023-11-07 [1] CRAN (R 4.4.0)
# magrittr               2.0.3   2022-03-30 [1] CRAN (R 4.4.0)
# Matrix                 1.7-0   2024-03-22 [1] CRAN (R 4.4.0)
# MatrixGenerics       * 1.17.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# matrixStats          * 1.3.0   2024-04-11 [1] CRAN (R 4.4.0)
# memoise                2.0.1   2021-11-26 [1] CRAN (R 4.4.0)
# pillar                 1.9.0   2023-03-22 [1] CRAN (R 4.4.0)
# pkgconfig              2.0.3   2019-09-22 [1] CRAN (R 4.4.0)
# png                    0.1-8   2022-11-29 [1] CRAN (R 4.4.0)
# prettyunits            1.2.0   2023-09-24 [1] CRAN (R 4.4.0)
# progress               1.2.3   2023-12-06 [1] CRAN (R 4.4.0)
# purrr                  1.0.2   2023-08-10 [1] CRAN (R 4.4.0)
# R6                     2.5.1   2021-08-19 [1] CRAN (R 4.4.0)
# rappdirs               0.3.3   2021-01-31 [1] CRAN (R 4.4.0)
# RColorBrewer           1.1-3   2022-04-03 [1] CRAN (R 4.4.0)
# readxl               * 1.4.3   2023-07-06 [1] CRAN (R 4.4.0)
# rjson                  0.2.21  2022-01-09 [1] CRAN (R 4.4.0)
# rlang                * 1.1.4   2024-06-04 [1] CRAN (R 4.4.0)
# rmarkdown              2.27    2024-05-17 [1] CRAN (R 4.4.0)
# rprojroot              2.0.4   2023-11-05 [1] CRAN (R 4.4.0)
# RSQLite                2.3.7   2024-05-27 [1] CRAN (R 4.4.0)
# S4Arrays               1.5.1   2024-05-29 [1] Bioconductor 3.20 (R 4.4.0)
# S4Vectors            * 0.43.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# sessioninfo          * 1.2.2   2021-12-06 [1] CRAN (R 4.4.0)
# shape                  1.4.6.1 2024-02-23 [1] CRAN (R 4.4.0)
# SparseArray            1.5.8   2024-06-14 [1] Bioconductor 3.20 (R 4.4.1)
# stringi                1.8.4   2024-05-06 [1] CRAN (R 4.4.0)
# stringr                1.5.1   2023-11-14 [1] CRAN (R 4.4.0)
# SummarizedExperiment * 1.35.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# tibble                 3.2.1   2023-03-20 [1] CRAN (R 4.4.0)
# tidyselect             1.2.1   2024-03-11 [1] CRAN (R 4.4.0)
# UCSC.utils             1.1.0   2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# utf8                   1.2.4   2023-10-22 [1] CRAN (R 4.4.0)
# vctrs                  0.6.5   2023-12-01 [1] CRAN (R 4.4.0)
# withr                  3.0.0   2024-01-16 [1] CRAN (R 4.4.0)
# xfun                   0.45    2024-06-16 [1] CRAN (R 4.4.0)
# xml2                   1.3.6   2023-12-04 [1] CRAN (R 4.4.0)
# XVector                0.45.0  2024-05-04 [1] Bioconductor 3.20 (R 4.4.0)
# zlibbioc               1.51.1  2024-06-05 [1] Bioconductor 3.20 (R 4.4.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
