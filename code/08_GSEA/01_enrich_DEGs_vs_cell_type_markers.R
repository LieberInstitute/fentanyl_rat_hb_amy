
library(here)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(biomaRt)
library(sessioninfo)


####################   Gene Set Enrichment Analysis   ######################

## All genes as universe
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
                       dataset="hsapiens_gene_ensembl")

    human_rat_ids <- getBM(values  = human_marker_genes,
                           mart  = mart,
                           attributes = c("external_gene_name", "rnorvegicus_homolog_ensembl_gene", "rnorvegicus_homolog_associated_gene_name"),
                           filters    = "external_gene_name")
    return(human_rat_ids)
}


## -----------------------------------------------------------------------------
##                 MeanRatio-based cell type marker genes
## -----------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                 Markers for cell types in human epithalamus*
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# *From https://github.com/LieberInstitute/Habenula_Pilot/tree/master

## All ranked marker genes
MeanRatio_genes <- as.data.frame(read_xlsx('processed-data/08_GSEA/Input_cell_type_markers/MeanRatio_Top50MarkerGenes.xlsx'))

## Cell types/clusters included
cell_types <- names(table(MeanRatio_genes$cellType.target))
# [1] "Astrocyte"  "Endo"       "Excit.Thal" "Inhib.Thal" "LHb.1"      "LHb.2"      "LHb.3"      "LHb.4"      "LHb.5"
# [10] "LHb.6"      "LHb.7"      "MHb.1"      "MHb.2"      "MHb.3"      "Microglia"  "Oligo"      "OPC"

## Confirm there are top 50 markers per cell type
table(MeanRatio_genes$cellType.target)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2      LHb.3      LHb.4      LHb.5      LHb.6      LHb.7
#        50         50         50         50         50         50         50         50         50         50         50
# MHb.1      MHb.2      MHb.3  Microglia      Oligo        OPC
#    50         50         50         50         50         50

## Divide marker genes per cell type
for (cell_type in cell_types){
    markers <- subset(MeanRatio_genes, cellType.target==cell_type)

    ## Find rat orthologs
    markers_rat_IDs <- obtain_rat_orthologs(markers$Symbol)
    ## Take unique rat ensembl IDs: human marker genes with at least one ortholog in rat
    markers_rat_IDs <- unique(markers_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]
    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))

    assign(paste0(cell_type, "_top50_marker_genes"), markers)
    assign(paste0(cell_type, "_top50_markers_ratIDs"), markers_rat_IDs)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                 Markers for cell types in human amygdala*
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# *From https://doi.org/10.1038/s41421-022-00506-y












## -----------------------------------------------------------------------------
##                 1vsALL-based cell type marker genes in human
## -----------------------------------------------------------------------------

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#               Markers for cell types in human epithalamus
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

####################  Broad resolution cell type markers  ######################
## (DEGs (FDR<0.05) based in the enrichment model for one cell type vs the rest were taken as markers)

lvsALL_broad_genes <- eval(parse_expr(load(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broadMarkerGenes.Rdata'))))
lvsALL_broad_genes_enrich_stats <- lvsALL_broad_genes$enrichment

## Cell types
cell_types <- gsub('fdr_', '', colnames(lvsALL_broad_genes_enrich_stats)[grep('fdr', colnames(lvsALL_broad_genes_enrich_stats))])
cell_types
# [1] "Astrocyte"   "Endo"   "Excit.Thal"  "Inhib.Thal"  "LHb"   "MHb"   "Microglia"  "Oligo"   "OPC"

## Cell type-specific DEGs and corresponding orthologs in rat
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
    assign(paste0(cell_type, "_1vsALL_broad_markers_ratIDs"), markers_rat_IDs)
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

lvsALL_fine_genes <- eval(parse_expr(load(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_finalMarkerGenes.Rdata'))))
lvsALL_fine_genes_enrich_stats <- lvsALL_fine_genes$enrichment

## Cell types
cell_types <- gsub('fdr_', '', colnames(lvsALL_fine_genes_enrich_stats)[grep('fdr', colnames(lvsALL_fine_genes_enrich_stats))])
cell_types
# [1] "Astrocyte"  "Endo"       "Excit.Thal" "Inhib.Thal" "LHb.1"      "LHb.2"      "LHb.3"      "LHb.4"      "LHb.5"
# [10] "LHb.6"      "LHb.7"      "MHb.1"      "MHb.2"      "MHb.3"      "Microglia"  "Oligo"      "OPC"

## Cell type-specific DEGs and corresponding orthologs in rat
for (cell_type in cell_types){

    ## Cell type data
    cell_type_stats <- lvsALL_fine_genes_enrich_stats[,  c(paste0(c('t_stat_', 'p_value_', 'fdr_', 'logFC_'), cell_type), 'gene')]
    ## DEGs symbols
    cell_type_DEGs <- subset(cell_type_stats, eval(parse_expr(paste0('fdr_', cell_type)))<0.05)$gene
    print(paste0('Number of ', cell_type, ' human DEGs: ', length(cell_type_DEGs)))

    ## Rat orthologs
    cell_type_DEGs_rat_IDs <- obtain_rat_orthologs(cell_type_DEGs)
    markers_rat_IDs <- unique(cell_type_DEGs_rat_IDs$rnorvegicus_homolog_ensembl_gene)
    markers_rat_IDs <- markers_rat_IDs[markers_rat_IDs!=""]

    print(paste0('Number of ', cell_type, ' marker genes in rat: ', length(markers_rat_IDs)))
    assign(paste0(cell_type, "_1vsALL_fine_markers_ratIDs"), markers_rat_IDs)
}



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                   Markers for cell types in human amygdala
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

lvsALL_broad_genes <- readRDS(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broadMarkerGenes.rds'))


####################  Fine resolution cell type markers  ######################

lvsALL_broad_genes <- readRDS(here('processed-data/08_GSEA/Input_cell_type_markers/lvsALL_broadMarkerGenes.rds'))












############################################################################
##                2. Enrichment analysis of rat DEGs
############################################################################





