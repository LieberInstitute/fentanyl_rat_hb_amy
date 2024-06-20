
library(here)
library(SummarizedExperiment)
library(biomaRt)
library(rlang)
library(readxl)
library(ComplexHeatmap)
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
##         1. Obtain sets of orthologs of human marker genes in rat
############################################################################

## Obtain rat orthologs of human marker genes
obtain_rat_orthologs <- function(human_marker_genes){
    mart <- useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl",
                       GRCh = "GRCh38")

    human_rat_ids <- getBM(values  = human_marker_genes,
                           mart  = mart,
                           attributes = c("external_gene_name", "rnorvegicus_homolog_ensembl_gene", "rnorvegicus_homolog_associated_gene_name"),
                           filters    = "external_gene_name")
    return(human_rat_ids)
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

## Divide marker genes per cell type and obtain rat IDs
MeanRatio_top50_fine_hab_ratIDs <- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs(markers$Symbol)
    ## Take unique rat ensembl IDs: human marker genes with at least one ortholog in rat
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

## Broad cell types
cell_types <- names(table(MeanRatio_genes$cellType.target))
cell_types
# [1] "Astrocyte"       "Endothelial"     "ExN"             "InN"             "Microglia"       "Oligodendrocyte" "OPC"

## Top 100 markers per cell type
table(MeanRatio_genes$cellType.target)
# Astrocyte     Endothelial             ExN             InN       Microglia Oligodendrocyte             OPC
#      100             100             100             100             100             100             100

## Marker genes per cell type and obtain rat IDs
MeanRatio_top100_broad_amyg_ratIDs <- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs(markers$gene)
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

## Marker genes per cell type and obtain rat IDs
MeanRatio_top100_fine_amyg_ratIDs <- list()
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs(markers$gene)
    markers_rat_IDs <- unique(markers_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]
    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))

    MeanRatio_top100_fine_amyg_ratIDs[[cell_type]] <- markers_rat_IDs

}
# [1] "Number of Human_Astro_1 FGFR3 marker genes in rat: 85"
# [1] "Number of Human_Astro_2 FGFR3 marker genes in rat: 85"
# [1] "Number of Human_Astro_3 FGFR3 marker genes in rat: 79"
# [1] "Number of Human_Astro_4 FGFR3 marker genes in rat: 74"
# [1] "Number of Human_CALCR LHX8 marker genes in rat: 96"
# [1] "Number of Human_DRD2 ISL1 marker genes in rat: 76"
# [1] "Number of Human_DRD2 PAX6 marker genes in rat: 93"
# [1] "Number of Human_Endo NOSTRIN marker genes in rat: 147"
# [1] "Number of Human_HGF C11orf87 marker genes in rat: 84"
# [1] "Number of Human_HGF ESR1 marker genes in rat: 97"
# [1] "Number of Human_HGF NPSR1 marker genes in rat: 88"
# [1] "Number of Human_HTR3A DRD2 marker genes in rat: 89"
# [1] "Number of Human_LAMP5 ABO marker genes in rat: 67"
# [1] "Number of Human_LAMP5 BDNF marker genes in rat: 69"
# [1] "Number of Human_LAMP5 COL14A1 marker genes in rat: 88"
# [1] "Number of Human_LAMP5 COL25A1 marker genes in rat: 96"
# [1] "Number of Human_LAMP5 NDNF marker genes in rat: 97"
# [1] "Number of Human_Micro CTSS marker genes in rat: 97"
# [1] "Number of Human_Oligo_1 OPALIN marker genes in rat: 98"
# [1] "Number of Human_Oligo_2 OPALIN marker genes in rat: 92"
# [1] "Number of Human_Oligo_3 OPALIN marker genes in rat: 104"
# [1] "Number of Human_Oligo_4 OPALIN marker genes in rat: 92"
# [1] "Number of Human_Oligo_5 OPALIN marker genes in rat: 75"
# [1] "Number of Human_Oligo_6 OPALIN marker genes in rat: 93"
# [1] "Number of Human_OPC_1 PDGFRA marker genes in rat: 102"
# [1] "Number of Human_OPC_2 PDGFRA marker genes in rat: 103"
# [1] "Number of Human_OPC_3 PDGFRA marker genes in rat: 90"
# [1] "Number of Human_OPC_4 PDGFRA marker genes in rat: 92"
# [1] "Number of Human_PRKCD marker genes in rat: 89"
# [1] "Number of Human_PVALB ADAMTS5 marker genes in rat: 79"
# [1] "Number of Human_RXFP2 RSPO2 marker genes in rat: 77"
# [1] "Number of Human_SATB2 CALCRL marker genes in rat: 84"
# [1] "Number of Human_SATB2 IL15 marker genes in rat: 91"
# [1] "Number of Human_SATB2 ST8SIA2 marker genes in rat: 86"
# [1] "Number of Human_SOX11 EBF2 marker genes in rat: 85"
# [1] "Number of Human_SST EPYC marker genes in rat: 84"
# [1] "Number of Human_SST HGF marker genes in rat: 89"
# [1] "Number of Human_STRIP2 marker genes in rat: 79"
# [1] "Number of Human_TFAP2C marker genes in rat: 73"
# [1] "Number of Human_TSHZ1 CALCRL marker genes in rat: 87"
# [1] "Number of Human_TSHZ1 SEMA3C marker genes in rat: 94"
# [1] "Number of Human_VGLL3 CNGB1 marker genes in rat: 81"
# [1] "Number of Human_VGLL3 MEPE marker genes in rat: 92"
# [1] "Number of Human_VIP ABI3BP marker genes in rat: 94"
# [1] "Number of Human_VIP NDNF marker genes in rat: 90"





## -----------------------------------------------------------------------------
##             B) 1vsALL-based cell type marker genes in human
## -----------------------------------------------------------------------------

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                i)  Markers for cell types in human epithalamus
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

####################  Broad resolution cell type markers  ######################
## (DEGs (FDR<0.05) based in the enrichment model for one cell type vs the rest were taken as markers)

lvsALL_broad_genes <- eval(parse_expr(load(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broad_MarkerGenes_hab.Rdata'))))
lvsALL_broad_genes_enrich_stats <- lvsALL_broad_genes$enrichment

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
    cell_type_DEGs <- subset(cell_type_stats, eval(parse_expr(paste0('fdr_', cell_type)))<0.05)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs(cell_type_DEGs)
    markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
    lvsALL_broad_hab_ratIDs[[cell_type]] <- markers_rat_IDs
}

# [1] "Number of Astrocyte human DEGs: 763"
# [1] "Number of Astrocyte marker genes in rat: 485"
# [1] "Number of Endo human DEGs: 5125"
# [1] "Number of Endo marker genes in rat: 3884"
# [1] "Number of Excit.Thal human DEGs: 322"
# [1] "Number of Excit.Thal marker genes in rat: 173"
# [1] "Number of Inhib.Thal human DEGs: 512"
# [1] "Number of Inhib.Thal marker genes in rat: 316"
# [1] "Number of LHb human DEGs: 25"
# [1] "Number of LHb marker genes in rat: 11"
# [1] "Number of MHb human DEGs: 127"
# [1] "Number of MHb marker genes in rat: 71"
# [1] "Number of Microglia human DEGs: 5959"
# [1] "Number of Microglia marker genes in rat: 4717"
# [1] "Number of Oligo human DEGs: 269"
# [1] "Number of Oligo marker genes in rat: 181"
# [1] "Number of OPC human DEGs: 196"
# [1] "Number of OPC marker genes in rat: 118"


####################  Fine resolution cell type markers  #######################

lvsALL_fine_genes <- eval(parse_expr(load(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_fine_MarkerGenes_hab.Rdata'))))
lvsALL_fine_genes_enrich_stats <- lvsALL_fine_genes$enrichment

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
    cell_type_DEGs <- subset(cell_type_stats, eval(parse_expr(paste0('fdr_', cell_type)))<0.05)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs (if DEGs exist)
    if (length(cell_type_DEGs)>0){
        cell_type_DEGs_rat_IDs <- obtain_rat_orthologs(cell_type_DEGs)
        markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
        markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

        print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
        lvsALL_fine_hab_ratIDs[[cell_type]] <- markers_rat_IDs
    }
}

# [1] "Number of Astrocyte human DEGs: 1448"
# [1] "Number of Astrocyte marker genes in rat: 1008"
# [1] "Number of Endo human DEGs: 6602"
# [1] "Number of Endo marker genes in rat: 4961"
# [1] "Number of Excit.Thal human DEGs: 552"
# [1] "Number of Excit.Thal marker genes in rat: 322"
# [1] "Number of Inhib.Thal human DEGs: 1223"
# [1] "Number of Inhib.Thal marker genes in rat: 821"
# [1] "Number of LHb.1 human DEGs: 25"
# [1] "Number of LHb.1 marker genes in rat: 10"
# [1] "Number of LHb.2 human DEGs: 25"
# [1] "Number of LHb.2 marker genes in rat: 15"
# [1] "Number of LHb.3 human DEGs: 11"
# [1] "Number of LHb.3 marker genes in rat: 4"
# [1] "Number of LHb.4 human DEGs: 0"
# [1] "Number of LHb.5 human DEGs: 2"
# [1] "Number of LHb.5 marker genes in rat: 1"
# [1] "Number of LHb.6 human DEGs: 17"
# [1] "Number of LHb.6 marker genes in rat: 7"
# [1] "Number of LHb.7 human DEGs: 5"
# [1] "Number of LHb.7 marker genes in rat: 2"
# [1] "Number of MHb.1 human DEGs: 64"
# [1] "Number of MHb.1 marker genes in rat: 29"
# [1] "Number of MHb.2 human DEGs: 131"
# [1] "Number of MHb.2 marker genes in rat: 94"
# [1] "Number of MHb.3 human DEGs: 308"
# [1] "Number of MHb.3 marker genes in rat: 271"
# [1] "Number of Microglia human DEGs: 7440"
# [1] "Number of Microglia marker genes in rat: 5888"
# [1] "Number of Oligo human DEGs: 668"
# [1] "Number of Oligo marker genes in rat: 460"
# [1] "Number of OPC human DEGs: 440"
# [1] "Number of OPC marker genes in rat: 295"



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                  ii)  Markers for cell types in human amygdala
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

####################  Broad resolution cell type markers  ######################

## Cell type-specific enrichment stats for all genes
lvsALL_broad_genes <- read.csv('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broad_MarkerGenes_amyg.csv')
## (DEGs (FDR<0.05) based on the enrichment model as marker genes)

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
    cell_type_DEGs <- subset(cell_type_stats, log.FDR<log(0.05))$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs(cell_type_DEGs)
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
    cell_type_DEGs <- subset(cell_type_stats, log.FDR<log(0.05))$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs(cell_type_DEGs)
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





############################################################################
##             2. Cell type enrichment analysis for rat DEGs
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

            ## Universe:
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

#   * Top 100 MeanRatio-based cell type marker genes at fine resolution
results_MeanRatio_100_fine_amyg <- enrichment_analysis('amyg', 'MeanRatio', 'top100', 'fine', 'amygdala')
p_values_MeanRatio_100_fine_amyg <- results_MeanRatio_100_fine_amyg[[1]]
heatmap_pvals(p_values_MeanRatio_100_fine_amyg, 'amygdala', 'Top100 MeanRatio-based amygdala fine cell type markers',
              'MeanRatio_top100_fine_amyg', 11)
ms_MeanRatio_100_fine_amyg <- results_MeanRatio_100_fine_amyg[[2]]

#   * Top 100 MeanRatio-based cell type marker genes at broad resolution
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



####  Compare habenula DEGs vs cell type marker genes in epithalamus

#   * Top 50 MeanRatio-based cell type marker genes at fine resolution
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






