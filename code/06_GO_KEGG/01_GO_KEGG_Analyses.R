
library(here)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Rn.eg.db)
library(sessioninfo)


#######################   GO & KEGG Analyses   #######################

load(here('processed-data/05_DEA/results_all_vars_habenula.Rdata'))
load(here('processed-data/05_DEA/results_uncorr_vars_habenula.Rdata'))
load(here('processed-data/05_DEA/results_all_vars_amygdala.Rdata'))
load(here('processed-data/05_DEA/results_uncorr_vars_amygdala.Rdata'))


## Groups of DEG

## Habenula DEGs from model with uncorrelated sample variables
top_genes_hab <- results_uncorr_vars_habenula[[1]]
## All DEGs
DEGs_hab <- top_genes_hab[which(top_genes_hab$adj.P.Val<0.05), ]
## Up and down DEGs
up_hab <- top_genes_hab[which(top_genes_hab$adj.P.Val<0.05 & top_genes_hab$logFC>0), ]
down_hab <- top_genes_hab[which(top_genes_hab$adj.P.Val<0.05 & top_genes_hab$logFC<0), ]

## Amygdala DEGs from model with uncorrelated sample variables
top_genes_amy <- results_uncorr_vars_amygdala[[1]]
DEGs_amy <- top_genes_amy[which(top_genes_amy$adj.P.Val<0.05), ]
up_amy <- top_genes_amy[which(top_genes_amy$adj.P.Val<0.05 & top_genes_amy$logFC>0), ]
down_amy <- top_genes_amy[which(top_genes_amy$adj.P.Val<0.05 & top_genes_amy$logFC<0), ]


## Function to do GO and KEGG analyses

GO_KEGG<- function(sigGeneList, geneUniverse, name){

    height=10
    width=9

    ## Do GO
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

    ## Save
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

    ## Save
    if (!is.null(goCC_Adj)){
        pdf(paste("plots/06_GO_KEGG/GO_CC_", name, ".pdf", sep=""), height = height, width = width)
        print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components"))
        dev.off()
    }


    ## Do KEGG
    kegg_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichKEGG",
        organism = 'rat',
        universe = geneUniverse,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    ## Save
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


##################################  All DEGs  ##################################

######################
#      Habenula
######################

## List of DEGs
sigGeneList <- list("All"= DEGs_hab[which(!is.na(DEGs_hab$EntrezID) & !DEGs_hab$EntrezID=='NULL'), 'EntrezID'])
## Background genes
geneUniverse <- as.character(top_genes_hab$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']

goList_habenula_all_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_all_DEGs')
save(goList_habenula_all_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_all_DEGs.Rdata")

######################
#      Amygdala
######################

sigGeneList <- list("All"= DEGs_amy[which(!is.na(DEGs_amy$EntrezID) & !DEGs_amy$EntrezID=='NULL'), 'EntrezID'])
geneUniverse <- as.character(top_genes_amy$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']

goList_amygdala_all_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_all_DEGs')
save(goList_amygdala_all_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_all_DEGs.Rdata")



################################  Up/Down DEGs  ################################

######################
#      Habenula
######################

## List of DEG sets
sigGeneList <- list("Up"=up_hab[which(!is.na(up_hab$EntrezID) & !up_hab$EntrezID=='NULL'), 'EntrezID'],
                    "Down"=down_hab[which(!is.na(down_hab$EntrezID) & !down_hab$EntrezID=='NULL'), 'EntrezID'])
## Background genes
geneUniverse <- as.character(top_genes_hab$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']

goList_habenula_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_up_down_DEGs')
save(goList_habenula_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_up_down_DEGs.Rdata")


######################
#      Amygdala
######################

sigGeneList <- list("Up"=up_amy[which(!is.na(up_amy$EntrezID) & !up_amy$EntrezID=='NULL'), 'EntrezID'],
                    "Down"=down_amy[which(!is.na(down_amy$EntrezID) & !down_amy$EntrezID=='NULL'), 'EntrezID'])

geneUniverse <- as.character(top_genes_amy$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']

goList_amygdala_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_up_down_DEGs')
save(goList_amygdala_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_up_down_DEGs.Rdata")


