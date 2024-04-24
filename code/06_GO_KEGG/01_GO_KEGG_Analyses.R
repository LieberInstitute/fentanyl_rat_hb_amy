
library(here)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Rn.eg.db)
library(sessioninfo)


#######################  Functional Enrichment Analysis  #######################

load(here('processed-data/05_DEA/de_genes_Substance_habenula.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/de_genes_Substance_amygdala.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/results_Substance_uncorr_vars_amygdala.Rdata'), verbose = TRUE)


## Groups of DEGs

########################
##    Habenula DEGs
########################
## Up and down habenula DEGs from model with uncorrelated sample variables
up_hab <- de_genes_habenula[which(de_genes_habenula$logFC>0), ]
down_hab <- de_genes_habenula[which(de_genes_habenula$logFC<0), ]

########################
##    Amygdala DEGs
########################
## Up and down amygdala DEGs from model with uncorrelated sample variables
up_amy <- de_genes_amygdala[which(de_genes_amygdala$logFC>0), ]
down_amy <- de_genes_amygdala[which(de_genes_amygdala$logFC<0), ]



## Function to find enriched GO and KEGG terms

GO_KEGG<- function(sigGeneList, geneUniverse, name){

    height=5
    width=7

    ## GO terms
    ## Obtain biological processes
    goBP_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = org.Rn.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    ## Save
    if (!is.null(goBP_Adj)){
        pdf(paste("plots/06_GO_KEGG/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
        print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes"))
        dev.off()
    }


    ## Obtain molecular functions
    goMF_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = org.Rn.eg.db,
        ont = "MF",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    if (!is.null(goMF_Adj)){
        pdf(paste("plots/06_GO_KEGG/GO_MF_", name, ".pdf", sep=""), height = height, width = width)
        print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function"))
        dev.off()
    }


    ## Obtain cellular components
    goCC_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = org.Rn.eg.db,
        ont = "CC",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    if (!is.null(goCC_Adj)){
        pdf(paste("plots/06_GO_KEGG/GO_CC_", name, ".pdf", sep=""), height = height, width = width)
        print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components"))
        dev.off()
    }


    ## KEGG terms
    kegg_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichKEGG",
        organism = 'rat',
        universe = geneUniverse,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    if (!is.null(kegg_Adj)){
        pdf(paste("plots/06_GO_KEGG/KEGG_", name, ".pdf", sep=""), height = height, width = width)
        print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis"))
        dev.off()
    }


    goList <- list(
        BP = goBP_Adj,
        MF = goMF_Adj,
        CC = goCC_Adj,
        KEGG = kegg_Adj
    )

    return(goList)
}


## 1. Analysis for all DEGs from each brain region

######################
#      Habenula
######################

## List of DEGs
sigGeneList <- list("All"= de_genes_habenula[which(!is.na(de_genes_habenula$EntrezID) & !de_genes_habenula$EntrezID=='NULL'), 'EntrezID'])
## Background genes (all genes assessed for DGE)
geneUniverse <- as.character(results_Substance_uncorr_vars_amygdala[[1]]$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']

goList_habenula_all_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_all_DEGs')
save(goList_habenula_all_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_all_DEGs.Rdata")

######################
#      Amygdala
######################

sigGeneList <- list("All"= de_genes_amygdala[which(!is.na(de_genes_amygdala$EntrezID) & !de_genes_amygdala$EntrezID=='NULL'), 'EntrezID'])

goList_amygdala_all_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_all_DEGs')
save(goList_amygdala_all_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_all_DEGs.Rdata")





## 2. Analysis for up- and down-regulated DEGs from each brain region

######################
#      Habenula
######################

## List of DEG sets
sigGeneList <- list("Up"=up_hab[which(!is.na(up_hab$EntrezID) & !up_hab$EntrezID=='NULL'), 'EntrezID'],
                    "Down"=down_hab[which(!is.na(down_hab$EntrezID) & !down_hab$EntrezID=='NULL'), 'EntrezID'])

goList_habenula_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_up_down_DEGs')
save(goList_habenula_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_up_down_DEGs.Rdata")

######################
#      Amygdala
######################

sigGeneList <- list("Up"=up_amy[which(!is.na(up_amy$EntrezID) & !up_amy$EntrezID=='NULL'), 'EntrezID'],
                    "Down"=down_amy[which(!is.na(down_amy$EntrezID) & !down_amy$EntrezID=='NULL'), 'EntrezID'])

goList_amygdala_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_up_down_DEGs')
save(goList_amygdala_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_up_down_DEGs.Rdata")







## Reproducibility information

options(width = 120)
session_info()

