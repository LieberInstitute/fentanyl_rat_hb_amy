
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

## Habenula DEGs from model with all sample variables
top_genes_hab <- results_all_vars_habenula[[1]]
up_hab <- top_genes_hab[which(top_genes_hab$adj.P.Val<0.1 & top_genes_hab$logFC>0), ]
down_hab <- top_genes_hab[which(top_genes_hab$adj.P.Val<0.1 & top_genes_hab$logFC<0), ]


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
        OrgDb =
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    ## Save
    pdf(paste("plots/06_GO_KEGG/GO_BP_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes"))
    dev.off()


    ## Obtain molecular functions
    goMF_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = ,
        ont = "MF",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    ## Save
    pdf(paste("plots/06_GO_KEGG/GO_MF_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function"))
    dev.off()


    ## Obtain cellular components
    goCC_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichGO",
        universe = geneUniverse,
        OrgDb = ,
        ont = "CC",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    ## Save
    pdf(paste("plots/06_GO_KEGG/GO_CC_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components"))
    dev.off()


    ## Do KEGG
    kegg_Adj <- compareCluster(
        sigGeneList,
        fun = "enrichKEGG",
        organism = "",
        universe = geneUniverse,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    ## Save
    pdf(paste("plots/06_GO_KEGG/KEGG_", name, ".pdf", sep=""), height = height, width = width)
    print(dotplot(kegg_Adj, title="KEGG Enrichment Analysis"))
    dev.off()


    goList <- list(
        BP = goBP_Adj,
        MF = goMF_Adj,
        CC = goCC_Adj,
        KEGG = kegg_Adj
    )

    return(goList)
}




######################
#    Up/Down DEGs
######################

## List of DEG sets
sigGeneList <- list("Up"=up_hab[which(!is.na(up_hab$EntrezID) & !up_hab$EntrezID=='NULL'), 'EntrezID'],
                    "Down"=down_hab[which(!is.na(down_hab$EntrezID) & !down_hab$EntrezID=='NULL'), 'EntrezID'])
## Background genes
geneUniverse <- as.character(top_genes_hab$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse) & !geneUniverse=='NULL']

goList_habenula_all_vars<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_all_vars')
save(goList_habenula_all_vars, file="processed-data/06_GO_KEGG/goList_habenula_all_vars.Rdata")


