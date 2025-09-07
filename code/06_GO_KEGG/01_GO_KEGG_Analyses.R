
library(dplyr)
library(tidyr)
library(tibble)
library(here)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Rn.eg.db)
library(cowplot)
library(ggplot2)
library(biomaRt)
library(sessioninfo)


#######################  Functional Enrichment Analysis  #######################

load(here('processed-data/05_DEA/de_genes_Substance_habenula.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/de_genes_Substance_amygdala.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/results_Substance_uncorr_vars_habenula.Rdata'), verbose = TRUE)
load(here('processed-data/05_DEA/results_Substance_uncorr_vars_amygdala.Rdata'), verbose = TRUE)


## Groups of DEGs

########################
##    Habenula DEGs
########################
## Up and down habenula DEGs from model with uncorrelated sample variables
up_hab <- de_genes_habenula[which(de_genes_habenula$logFC>0),]
down_hab <- de_genes_habenula[which(de_genes_habenula$logFC<0),]

########################
##    Amygdala DEGs
########################
## Up and down amygdala DEGs from model with uncorrelated sample variables
up_amy <- de_genes_amygdala[which(de_genes_amygdala$logFC>0),]
down_amy <- de_genes_amygdala[which(de_genes_amygdala$logFC<0),]

#########################################
##   Up/Down unique/shared in Hb/Amyg
#########################################
only_up_hab <- de_genes_habenula[which(!de_genes_habenula$ensemblID %in% de_genes_amygdala$ensemblID & de_genes_habenula$logFC>0),]
only_down_hab <- de_genes_habenula[which(!de_genes_habenula$ensemblID %in% de_genes_amygdala$ensemblID & de_genes_habenula$logFC<0),]

only_up_amy <- de_genes_amygdala[which(!de_genes_amygdala$ensemblID %in% de_genes_habenula$ensemblID & de_genes_amygdala$logFC>0),]
only_down_amy <- de_genes_amygdala[which(!de_genes_amygdala$ensemblID %in% de_genes_habenula$ensemblID & de_genes_amygdala$logFC<0),]

shared_hab_amy <- inner_join(de_genes_habenula, de_genes_amygdala, by = colnames(de_genes_amygdala)[1:9], suffix = c(".hb", ".amyg"))

shared_up_hab_up_amy <- shared_hab_amy[shared_hab_amy$logFC.hb>0 & shared_hab_amy$logFC.amyg>0, ]
shared_up_hab_down_amy <- shared_hab_amy[shared_hab_amy$logFC.hb>0 & shared_hab_amy$logFC.amyg<0, ]
shared_down_hab_up_amy <- shared_hab_amy[shared_hab_amy$logFC.hb<0 & shared_hab_amy$logFC.amyg>0, ]
shared_down_hab_down_amy <- shared_hab_amy[shared_hab_amy$logFC.hb<0 & shared_hab_amy$logFC.amyg<0, ]


## Retrieve valid Entrez IDs
only_up_hab_genes <- only_up_hab %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()
only_down_hab_genes <- only_down_hab %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()
only_up_amy_genes <- only_up_amy %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()
only_down_amy_genes <- only_down_amy %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()
shared_up_hab_up_amy_genes <- shared_up_hab_up_amy %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()
shared_up_hab_down_amy_genes <- shared_up_hab_down_amy %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()
shared_down_hab_up_amy_genes <- shared_down_hab_up_amy %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()
shared_down_hab_down_amy_genes <- shared_down_hab_down_amy %>% dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()


## Background genes (all genes assessed for DGE) -- same genes in Hb and Amyg DGE
geneUniverse <- results_Substance_uncorr_vars_amygdala[[1]] %>%
    dplyr::filter(!is.na(EntrezID) & !is.null(EntrezID) & EntrezID != "NULL" & EntrezID != "") %>% pull(EntrezID) %>% unique()


## Function to find enriched GO and KEGG terms
GO_KEGG<- function(sigGeneList, geneUniverse, name){

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
        p1 <- dotplot(goBP_Adj, title="GO Enrichment Analysis: Biological processes")
        goBP_Adj <- as.data.frame(goBP_Adj)
        goBP_Adj$geneID <- sapply(goBP_Adj$geneID, function(row){gsub("/", ", ", row)})
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
        p2 <- dotplot(goMF_Adj, title="GO Enrichment Analysis: Molecular function")
        goMF_Adj <- as.data.frame(goMF_Adj)
        goMF_Adj$geneID <- sapply(goMF_Adj$geneID, function(row){gsub("/", ", ", row)})
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
        p3 <- dotplot(goCC_Adj, title="GO Enrichment Analysis: Cellular components")
        goCC_Adj <- as.data.frame(goCC_Adj)
        goCC_Adj$geneID <- sapply(goCC_Adj$geneID, function(row){gsub("/", ", ", row)})
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
        p4 <- dotplot(kegg_Adj, title="KEGG Enrichment Analysis")
        ## Add symbols
        kegg_Adj <- as.data.frame(kegg_Adj)
        genes <-sapply(kegg_Adj$geneID, function(term_genes){unlist(strsplit(term_genes, "/"))})
        names(genes) <- kegg_Adj$ID
        mart = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
        term_symbols <- lapply(genes, function(term_genes){

            symbols <- getBM(attributes = c("entrezgene_id", "external_gene_name", "ensembl_gene_id"),
                             filters = "entrezgene_id",
                             values = term_genes,
                             mart = mart)
            symbols <- apply(symbols, 1, function(gene){if(!is.na(gene["external_gene_name"])){gene["external_gene_name"]}
                else if(!is.na(gene["ensembl_gene_id"])){gene["ensembl_gene_id"]}
                else{gene["entrezgene_id"]}})

            paste(symbols, collapse = ", ")

        })
        kegg_Adj$geneID <- do.call(rbind, term_symbols)
    }

    ## Plots
    if(name != "Hb_and_Amyg_Up_and_Down_unique_and_shared_DEGs"){
        h = 10
        w = 14
    } else{
        h = 35
        w = 49
    }
    plot_grid(p1, p2, p3, p4, ncol=2, align = 'vh')
    ggsave(paste("plots/06_GO_KEGG/GO_KEGG_", name, ".pdf", sep=""), height = h, width = w)

    ## Save results
    goList <- list(
        BP = goBP_Adj,
        MF = goMF_Adj,
        CC = goCC_Adj,
        KEGG = kegg_Adj
    )

    return(goList)
}

#-------------------------------------------------------------------------------
## A. Analysis for all DEGs from each brain region

######################
#      Habenula
######################
sigGeneList <- list("All"= unique(de_genes_habenula[which(!is.na(de_genes_habenula$EntrezID) & !de_genes_habenula$EntrezID=='NULL' & !de_genes_habenula$EntrezID==''), 'EntrezID']))
goList_habenula_all_DEGs <- GO_KEGG(sigGeneList, geneUniverse, 'habenula_all_DEGs')
save(goList_habenula_all_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_all_DEGs.Rdata")

######################
#      Amygdala
######################
sigGeneList <- list("All"= unique(de_genes_amygdala[which(!is.na(de_genes_amygdala$EntrezID) & !de_genes_amygdala$EntrezID=='NULL' & !de_genes_amygdala$EntrezID==''), 'EntrezID']))
goList_amygdala_all_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_all_DEGs')
save(goList_amygdala_all_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_all_DEGs.Rdata")

#-------------------------------------------------------------------------------
## B. Analysis for up- and down-regulated DEGs from each brain region

######################
#      Habenula
######################

## List of DEG sets
sigGeneList <- list("Up"=up_hab[which(!is.na(up_hab$EntrezID) & !up_hab$EntrezID=='NULL' & !up_hab$EntrezID==''), 'EntrezID'],
                    "Down"=down_hab[which(!is.na(down_hab$EntrezID) & !down_hab$EntrezID=='NULL' & !down_hab$EntrezID==''), 'EntrezID'])

goList_habenula_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'habenula_up_down_DEGs')
save(goList_habenula_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_habenula_up_down_DEGs.Rdata")

## Merge
go_kegg_results_hab <- rbind(cbind(goList_habenula_up_down_DEGs$BP, Ontology = "BP"),
                             cbind(goList_habenula_up_down_DEGs$MF, Ontology = "MF"),
                             cbind(goList_habenula_up_down_DEGs$CC, Ontology = "CC"),
                             cbind(goList_habenula_up_down_DEGs$KEGG[colnames(goList_habenula_up_down_DEGs$BP)], Ontology = "KEGG"))
go_kegg_results_hab$DEGs_set <- go_kegg_results_hab$Cluster
go_kegg_results_hab$Cluster <- NULL

go_kegg_results_hab <- go_kegg_results_hab[, c("Ontology", "DEGs_set", "ID", "Description", "Count", "GeneRatio",
                                               "BgRatio", "FoldEnrichment", "pvalue", "p.adjust", "geneID")]
go_kegg_results_hab <- go_kegg_results_hab[order(go_kegg_results_hab$Ontology, go_kegg_results_hab$DEGs_set, go_kegg_results_hab$p.adjust), ]


write.table(go_kegg_results_hab, "processed-data/Supplementary_Tables/TableS8_GO_KEGG_results_hab.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')

######################
#      Amygdala
######################

sigGeneList <- list("Up"=up_amy[which(!is.na(up_amy$EntrezID) & !up_amy$EntrezID=='NULL' & !up_amy$EntrezID==''), 'EntrezID'],
                    "Down"=down_amy[which(!is.na(down_amy$EntrezID) & !down_amy$EntrezID=='NULL' & !down_amy$EntrezID==''), 'EntrezID'])

goList_amygdala_up_down_DEGs<-GO_KEGG(sigGeneList, geneUniverse, 'amygdala_up_down_DEGs')
save(goList_amygdala_up_down_DEGs, file="processed-data/06_GO_KEGG/goList_amygdala_up_down_DEGs.Rdata")

go_kegg_results_amy <- rbind(cbind(goList_amygdala_up_down_DEGs$BP, Ontology = "BP"),
                             cbind(goList_amygdala_up_down_DEGs$MF, Ontology = "MF"),
                             cbind(goList_amygdala_up_down_DEGs$CC, Ontology = "CC"),
                             cbind(goList_amygdala_up_down_DEGs$KEGG[colnames(goList_amygdala_up_down_DEGs$BP)], Ontology = "KEGG"))
go_kegg_results_amy$DEGs_set <- go_kegg_results_amy$Cluster
go_kegg_results_amy$Cluster <- NULL

go_kegg_results_amy <- go_kegg_results_amy[, c("Ontology", "DEGs_set", "ID", "Description", "Count", "GeneRatio",
                                               "BgRatio", "FoldEnrichment", "pvalue", "p.adjust", "geneID")]
go_kegg_results_amy <- go_kegg_results_amy[order(go_kegg_results_amy$Ontology, go_kegg_results_amy$DEGs_set, go_kegg_results_amy$p.adjust), ]
write.table(go_kegg_results_amy, "processed-data/Supplementary_Tables/TableS9_GO_KEGG_results_amy.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')

#-------------------------------------------------------------------------------
## C. Analysis for up/down DEGs unique/shared in Hb and Amyg

sigGeneList <- list("Unique in Hb - Up" = only_up_hab_genes,
                    "Unique in Hb - Down" = only_down_hab_genes,
                    "Unique in Amyg - Up" = only_up_amy_genes,
                    "Unique in Amyg - Down" = only_down_amy_genes,
                    "Shared: Up in Hb, Up in Amyg" = shared_up_hab_up_amy_genes,
                    "Shared: Up in Hb, Down in Amyg" = shared_up_hab_down_amy_genes,
                    "Shared: Down in Hb, Up in Amyg" = shared_down_hab_up_amy_genes,
                    "Shared: Down in Hb, Down in Amyg" = shared_down_hab_down_amy_genes)

goList_hb_and_amyg_DEGs <- GO_KEGG(sigGeneList, geneUniverse,
                                   'Hb_and_Amyg_Up_and_Down_unique_and_shared_DEGs')
save(goList_hb_and_amyg_DEGs, file="processed-data/06_GO_KEGG/goList_Hb_and_Amyg_Up_and_Down_unique_and_shared_DEGs.Rdata")

## Merge
go_kegg_results <- rbind(cbind(goList_hb_and_amyg_DEGs$BP, Ontology = "BP"),
                         cbind(goList_hb_and_amyg_DEGs$MF, Ontology = "MF"),
                         cbind(goList_hb_and_amyg_DEGs$CC, Ontology = "CC"),
                         cbind(goList_hb_and_amyg_DEGs$KEGG[colnames(goList_hb_and_amyg_DEGs$BP)], Ontology = "KEGG"))

go_kegg_results$DEGs_set <- go_kegg_results$Cluster
go_kegg_results$Cluster <- NULL

go_kegg_results <- go_kegg_results[, c("Ontology", "DEGs_set", "ID", "Description", "Count", "GeneRatio",
                                               "BgRatio", "FoldEnrichment", "pvalue", "p.adjust", "geneID")]
go_kegg_results <- go_kegg_results[order(go_kegg_results$Ontology, go_kegg_results$DEGs_set, go_kegg_results$p.adjust), ]

write.table(go_kegg_results, "processed-data/Supplementary_Tables/TableS10_GO_KEGG_results_hab_vs_amyg.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')

## ------
## Heatmap with GO & KEGG enrichment results for Hb vs Amyg

## Find enriched terms of interest (list provided by Kristen and Robin)

## For Hb:
go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "regulation of synapse structure or activity") %>% .[,1:4]
#   Ontology            DEGs_set         ID                                 Description
# 1       BP   Unique in Hb - Up GO:0050803 regulation of synapse structure or activity
# 2       BP Unique in Amyg - Up GO:0050803 regulation of synapse structure or activity

go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "regulation of presynaptic membrane potential") %>% .[,1:4]
#   Ontology                     DEGs_set         ID                                  Description
# 1       BP            Unique in Hb - Up GO:0099505 regulation of presynaptic membrane potential
# 2       BP Shared: Up in Hb, Up in Amyg GO:0099505 regulation of presynaptic membrane potential

go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "regulation of postsynaptic membrane potential") %>% .[,1:4]
#   Ontology          DEGs_set         ID                                   Description
# 1       BP Unique in Hb - Up GO:0060078 regulation of postsynaptic membrane potential

go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "extracellular matrix organization") %>% .[,1:4]
#   Ontology                         DEGs_set         ID                       Description
# 1       BP              Unique in Hb - Down GO:0030198 extracellular matrix organization
# 2       BP            Unique in Amyg - Down GO:0030198 extracellular matrix organization
# 3       BP Shared: Down in Hb, Down in Amyg GO:0030198 extracellular matrix organization

go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "oligodendrocyte differentiation") %>% .[,1:4]
#   Ontology                         DEGs_set         ID                     Description
# 1       BP            Unique in Amyg - Down GO:0048709 oligodendrocyte differentiation
# 2       BP Shared: Down in Hb, Down in Amyg GO:0048709 oligodendrocyte differentiation


go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "presynaptic membrane") %>% .[,1:4]
#   Ontology                     DEGs_set         ID          Description
# 1       CC            Unique in Hb - Up GO:0042734 presynaptic membrane
# 2       CC Shared: Up in Hb, Up in Amyg GO:0042734 presynaptic membrane

go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "postsynaptic membrane") %>% .[,1:4]
#   Ontology                     DEGs_set         ID           Description
# 1       CC            Unique in Hb - Up GO:0045211 postsynaptic membrane
# 2       CC Shared: Up in Hb, Up in Amyg GO:0045211 postsynaptic membrane

go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "collagen-containing extracellular matrix") %>% .[,1:4]
#   Ontology                         DEGs_set         ID                              Description
# 1       CC              Unique in Hb - Down GO:0062023 collagen-containing extracellular matrix
# 2       CC            Unique in Amyg - Down GO:0062023 collagen-containing extracellular matrix
# 3       CC Shared: Down in Hb, Down in Amyg GO:0062023 collagen-containing extracellular matrix

go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "myelin sheath") %>% .[,1:4]
# Ontology                         DEGs_set         ID   Description
# 1       CC            Unique in Amyg - Down GO:0043209 myelin sheath
# 2       CC Shared: Down in Hb, Down in Amyg GO:0043209 myelin sheath
#

go_kegg_results %>% dplyr::filter(Ontology == "MF", Description == "voltage-gated monoatomic cation channel activity") %>% .[,1:4]
#   Ontology                     DEGs_set         ID                                      Description
# 1       MF            Unique in Hb - Up GO:0022843 voltage-gated monoatomic cation channel activity
# 2       MF Shared: Up in Hb, Up in Amyg GO:0022843 voltage-gated monoatomic cation channel activity

go_kegg_results %>% dplyr::filter(Ontology == "MF", Description == "metal ion transmembrane transporter activity") %>% .[,1:4]
# Ontology          DEGs_set         ID                                  Description
# 1       MF Unique in Hb - Up GO:0046873 metal ion transmembrane transporter activity


go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "Morphine addiction") %>% .[,1:4]
#   Ontology          DEGs_set       ID        Description
# 1     KEGG Unique in Hb - Up rno05032 Morphine addiction

go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "Glutamatergic synapse") %>% .[,1:4]
#   Ontology                       DEGs_set       ID           Description
# 1     KEGG              Unique in Hb - Up rno04724 Glutamatergic synapse
# 2     KEGG            Unique in Amyg - Up rno04724 Glutamatergic synapse
# 3     KEGG Shared: Up in Hb, Down in Amyg rno04724 Glutamatergic synapse

go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "Gastric acid secretion") %>% .[,1:4]
# Ontology                       DEGs_set       ID            Description
# 1     KEGG              Unique in Hb - Up rno04971 Gastric acid secretion
# 2     KEGG            Unique in Hb - Down rno04971 Gastric acid secretion
# 3     KEGG Shared: Up in Hb, Down in Amyg rno04971 Gastric acid secretion

go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "ECM-receptor interaction") %>% .[,1:4]
#   Ontology                       DEGs_set       ID              Description
# 1     KEGG            Unique in Hb - Down rno04512 ECM-receptor interaction
# 2     KEGG          Unique in Amyg - Down rno04512 ECM-receptor interaction
# 3     KEGG Shared: Up in Hb, Down in Amyg rno04512 ECM-receptor interaction


## For Amyg:
go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "aerobic respiration") %>% .[,1:4]
#   Ontology            DEGs_set         ID         Description
# 1       BP Unique in Amyg - Up GO:0009060 aerobic respiration

go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "vesicle-mediated transport in synapse") %>% .[,1:4]
#   Ontology            DEGs_set         ID                           Description
# 1       BP Unique in Amyg - Up GO:0099003 vesicle-mediated transport in synapse

go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "glial cell differentiation") %>% .[,1:4]
#   Ontology                         DEGs_set         ID                Description
# 1       BP            Unique in Amyg - Down GO:0010001 glial cell differentiation
# 2       BP Shared: Down in Hb, Down in Amyg GO:0010001 glial cell differentiation

go_kegg_results %>% dplyr::filter(Ontology == "BP", Description == "extracellular matrix organization") %>% .[,1:4]
#   Ontology                         DEGs_set         ID                       Description
# 1       BP              Unique in Hb - Down GO:0030198 extracellular matrix organization
# 2       BP            Unique in Amyg - Down GO:0030198 extracellular matrix organization
# 3       BP Shared: Down in Hb, Down in Amyg GO:0030198 extracellular matrix organization


go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "mitochondrial inner membrane") %>% .[,1:4]
#   Ontology            DEGs_set         ID                  Description
# 1       CC Unique in Amyg - Up GO:0005743 mitochondrial inner membrane

go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "exocytic vesicle") %>% .[,1:4]
#   Ontology            DEGs_set         ID      Description
# 1       CC Unique in Amyg - Up GO:0070382 exocytic vesicle

go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "extracellular matrix") %>% .[,1:4]
#   Ontology              DEGs_set         ID          Description
# 1       CC   Unique in Hb - Down GO:0031012 extracellular matrix
# 2       CC Unique in Amyg - Down GO:0031012 extracellular matrix

go_kegg_results %>% dplyr::filter(Ontology == "CC", Description == "myelin sheath") %>% .[,1:4]
#   Ontology                         DEGs_set         ID   Description
# 1       CC            Unique in Amyg - Down GO:0043209 myelin sheath
# 2       CC Shared: Down in Hb, Down in Amyg GO:0043209 myelin sheath


go_kegg_results %>% dplyr::filter(Ontology == "MF", Description == "structural constituent of ribosome") %>% .[,1:4]
#   Ontology            DEGs_set         ID                        Description
# 1       MF Unique in Amyg - Up GO:0003735 structural constituent of ribosome

go_kegg_results %>% dplyr::filter(Ontology == "MF", Description == "proton-transporting ATP synthase activity, rotational mechanism") %>% .[,1:4]
#   Ontology            DEGs_set         ID                                                     Description
# 1       MF Unique in Amyg - Up GO:0046933 proton-transporting ATP synthase activity, rotational mechanism

go_kegg_results %>% dplyr::filter(Ontology == "MF", Description == "cytoskeletal motor activity") %>% .[,1:4]
#   Ontology              DEGs_set         ID                 Description
# 1       MF Unique in Amyg - Down GO:0003774 cytoskeletal motor activity

go_kegg_results %>% dplyr::filter(Ontology == "MF", Description == "extracellular matrix structural constituent") %>% .[,1:4]
#   Ontology                       DEGs_set         ID                                 Description
# 1       MF          Unique in Amyg - Down GO:0005201 extracellular matrix structural constituent
# 2       MF Shared: Up in Hb, Down in Amyg GO:0005201 extracellular matrix structural constituent


go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "Oxidative phosphorylation") %>% .[,1:4]
#   Ontology            DEGs_set       ID               Description
# 1     KEGG Unique in Amyg - Up rno00190 Oxidative phosphorylation

go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "Parkinson disease") %>% .[,1:4]
#   Ontology            DEGs_set       ID       Description
# 1     KEGG Unique in Amyg - Up rno05012 Parkinson disease

go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "ECM-receptor interaction") %>% .[,1:4]
#   Ontology                       DEGs_set       ID              Description
# 1     KEGG            Unique in Hb - Down rno04512 ECM-receptor interaction
# 2     KEGG          Unique in Amyg - Down rno04512 ECM-receptor interaction
# 3     KEGG Shared: Up in Hb, Down in Amyg rno04512 ECM-receptor interaction

go_kegg_results %>% dplyr::filter(Ontology == "KEGG", Description == "Regulation of actin cytoskeleton") %>% .[,1:4]
#   Ontology              DEGs_set       ID                      Description
# 1     KEGG   Unique in Amyg - Up rno04810 Regulation of actin cytoskeleton
# 2     KEGG Unique in Amyg - Down rno04810 Regulation of actin cytoskeleton


enriched_terms <- list("only_up_hab" = c("regulation of synapse structure or activity", "regulation of presynaptic membrane potential",
                                         "regulation of postsynaptic membrane potential",
                                            "presynaptic membrane", "postsynaptic membrane", "voltage-gated monoatomic cation channel activity",
                                            "metal ion transmembrane transporter activity", "Morphine addiction", "Glutamatergic synapse",
                                            "Gastric acid secretion"),
                           "only_down_hab" = c("extracellular matrix organization", "collagen-containing extracellular matrix", "Gastric acid secretion",
                                               "ECM-receptor interaction", "extracellular matrix"),
                           "only_up_amy" = c("regulation of synapse structure or activity", "Glutamatergic synapse", "aerobic respiration",
                                             "vesicle-mediated transport in synapse", "mitochondrial inner membrane", "exocytic vesicle",
                                             "structural constituent of ribosome", "proton-transporting ATP synthase activity, rotational mechanism",
                                             "Oxidative phosphorylation", "Parkinson disease", "Regulation of actin cytoskeleton"),
                           "only_down_amy" = c("extracellular matrix organization", "oligodendrocyte differentiation",
                                               "collagen-containing extracellular matrix", "myelin sheath", "ECM-receptor interaction",
                                               "glial cell differentiation", "extracellular matrix", "cytoskeletal motor activity",
                                               "extracellular matrix structural constituent", "Regulation of actin cytoskeleton"),
                           "shared_up_hab_up_amy" = c("regulation of presynaptic membrane potential", "presynaptic membrane",
                                                      "postsynaptic membrane", "voltage-gated monoatomic cation channel activity"),
                           "shared_up_hab_down_amy" = c("Glutamatergic synapse", "Gastric acid secretion", "ECM-receptor interaction",
                                                        "extracellular matrix structural constituent"),
                           "shared_down_hab_up_amy" = c(),
                           "shared_down_hab_down_amy" = c("extracellular matrix organization", "oligodendrocyte differentiation",
                                                          "collagen-containing extracellular matrix", "myelin sheath", "glial cell differentiation")
                       )

## Annotate DEGs of interest in each term of interest (Kristen and Robin list) in specific group(s) (up/down unique/shared in Hb/Amyg)

only_up_hab_genes_2_show = list()
only_down_hab_genes_2_show = list()
only_up_amy_genes_2_show = list()
only_down_amy_genes_2_show = list()

shared_up_hab_up_amy_genes_2_show = list("regulation of presynaptic membrane potential" = c("Kcnj3", "Kcnc2", "Kcnj9", "Kctd16"),
                                         "presynaptic membrane" = c("Kcnj3", "Kcnc2", "Kcnj9", "Kctd16"),
                                         "postsynaptic membrane" = c("Kcnc2", "LRRTM1", "Epha4", "Lrrtm2"),
                                         "voltage-gated monoatomic cation channel activity" = c("Kcnj3", "Kcnc2", "Kcnj9"))

shared_up_hab_down_amy_genes_2_show = list("Glutamatergic synapse",
                                           "Gastric acid secretion",
                                           "extracellular matrix structural constituent" = c("Col4a3"),
                                           "ECM-receptor interaction" = c("Col4a3"))

shared_down_hab_down_amy_genes_2_show  = list("extracellular matrix organization" = c("Tgfbi", "Antxr1", "Col9a3", "Loxl4", "Sox9"),
                                              "oligodendrocyte differentiation" = c("Cnp", "Sox8", "Gsn", "Sox9", "Opalin"),
                                              "collagen-containing extracellular matrix" = c("Tgfbi", "Col9a3", "Loxl4", "Fgfr2"),
                                              "myelin sheath" = c("Cnp", "Tubb4a", "Gsn", "Tspan2"),
                                              "glial cell differentiation" = c("Cnp", "Gsn", "Sox9", "Opalin", "Tspan2"))



## Show additional top 5 most signif DEGs in each group
only_up_hab_genes_additional_top = list("regulation of synapse structure or activity" = c(),
                                "regulation of presynaptic membrane potential" = c(),
                                "regulation of postsynaptic membrane potential" = c(),
                                "presynaptic membrane" = c(),
                                "postsynaptic membrane" = c(),
                                "voltage-gated monoatomic cation channel activity" = c(),
                                "metal ion transmembrane transporter activity" = c(),
                                "Morphine addiction" = c(),
                                "Glutamatergic synapse" = c(),
                                "Gastric acid secretion" = c())


only_down_hab_genes_additional_top = list("extracellular matrix organization" = c(),
                                  "collagen-containing extracellular matrix" = c(),
                                  "Gastric acid secretion" = c(),
                                  "ECM-receptor interaction" = c(),
                                  "extracellular matrix" = c(),
                                  "ECM-receptor interaction" = c())

only_up_amy_genes_additional_top = list("regulation of synapse structure or activity" = c(),
                                "Glutamatergic synapse" = c(),
                                "aerobic respiration" = c(),
                                "vesicle-mediated transport in synapse" = c(),
                                "mitochondrial inner membrane" = c(),
                                "exocytic vesicle" = c(),
                                "structural constituent of ribosome" = c(),
                                "proton-transporting ATP synthase activity, rotational mechanism" = c(),
                                "Oxidative phosphorylation" = c(),
                                "Parkinson disease" = c(),
                                "Regulation of actin cytoskeleton" = c())

only_down_amy_genes_2_show = list("extracellular matrix organization",
                                  "oligodendrocyte differentiation",
                                  "collagen-containing extracellular matrix",
                                  "myelin sheath",
                                  "ECM-receptor interaction",
                                  "glial cell differentiation",
                                  "extracellular matrix",
                                  "cytoskeletal motor activity",
                                  "extracellular matrix structural constituent",
                                  "ECM-receptor interaction",
                                  "Regulation of actin cytoskeleton")


## Extract DEGs of interest terms across DEG groups (include interest DEGs from Kristen and Robin list)
genes_2_show_x_term_x_group <- list(list(), list(), list(), list(), list(), list(), list(), list())
names(genes_2_show_x_term_x_group) <- unique(go_kegg_results$DEGs_set)

## Specific DEGs to highlight
interest_genes <- unique(c("Ptpn3", "Kcnc2", "Kcnj9", "Kcnj3", "Kctd16",
                    "LRRTM1", "Lrrtm1", "Epha4", "Lrrn3", "Fam107a", "LRRTM2", "Lrrtm2",
                    "Loxl4", "Antxr1", "Tgfbi", "Col9a3", "Sox9",
                    "Cnp", "Gsn", "Sox8", "Sox9", "Opalin",
                    "Kcnc2", "Npy1r", "Kcnj9", "Epha4", "Kcnj3",
                    "Kcnc2", "LRRTM1", "Epha4", "Kctd16", "Cdh10",
                    "Loxl4", "Col9a3", "Tgfbi", "Fgfr2",
                    "Cnp", "Gsn", "Tspan2", "Tubb4a",
                    "Kcnc2", "Kcnj9", "Kcnj3",
                    "Slc13a5", "Kcnc2", "Kcnj9", "Kcnj3",
                    "Col9a3",
                    "Ap3s1",
                    "Col4a3"))

## Terms of interest
terms <- c("regulation of presynaptic membrane potential",  "regulation of postsynaptic membrane potential",
           "regulation of synapse structure or activity", "extracellular matrix organization", "oligodendrocyte differentiation",
           "presynaptic membrane", "postsynaptic membrane", "collagen-containing extracellular matrix", "myelin sheath",
           "voltage-gated monoatomic cation channel activity", "metal ion transmembrane transporter activity",
           "Morphine addiction", "Glutamatergic synapse", "Gastric acid secretion", "ECM-receptor interaction",
           "aerobic respiration", "vesicle-mediated transport in synapse", "glial cell differentiation", "extracellular matrix organization",
           "mitochondrial inner membrane", "exocytic vesicle", "extracellular matrix", "structural constituent of ribosome",
           "proton-transporting ATP synthase activity, rotational mechanism", "cytoskeletal motor activity",
           "extracellular matrix structural constituent", "Oxidative phosphorylation", "Parkinson disease",
           "Regulation of actin cytoskeleton")

for(term in terms){

    ## Extract groups of DEGs where term is enriched
    DEG_groups_with_term <- go_kegg_results %>% dplyr::filter(Description == term) %>% pull(DEGs_set)

    for(group in DEG_groups_with_term){

        ## Genes in term and group
        intersection_genes <- strsplit(go_kegg_results %>% dplyr::filter(Description == term & DEGs_set == group) %>% pull(geneID), ", ") %>% unlist
        ## Subset to the specific genes of interest
        genes_of_interest_2_show <- intersect(intersection_genes, interest_genes)

        if(length(genes_of_interest_2_show) >= 3){
            genes_2_show_x_term_x_group[[group]][[term]] = genes_of_interest_2_show
        }

        else{
            ## Add top n most signif genes in group and term
            n = 3 - length(genes_of_interest_2_show)

            if(group == "Unique in Hb - Up"){
                top3_intersection_genes <- only_up_hab %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>% arrange(adj.P.Val) %>% slice(1:n) %>% pull(Symbol)
            }
            else if(group == "Unique in Hb - Down"){
                top3_intersection_genes <- only_down_hab %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>% arrange(adj.P.Val) %>% slice(1:n) %>% pull(Symbol)
            }
            else if(group == "Unique in Amyg - Up"){
                top3_intersection_genes <- only_up_amy %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>% arrange(adj.P.Val) %>% slice(1:n) %>% pull(Symbol)
            }
            else if(group == "Unique in Amyg - Down"){
                top3_intersection_genes <- only_down_amy %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>% arrange(adj.P.Val) %>% slice(1:n) %>% pull(Symbol)
            }
            else if(group == "Shared: Up in Hb, Up in Amyg"){
                top3_intersection_genes <- shared_up_hab_up_amy %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>%
                    mutate(min_p = pmin(adj.P.Val.hb, adj.P.Val.amyg)) %>% arrange(min_p) %>% slice(1:n) %>% pull(Symbol)
            }
            else if(group == "Shared: Up in Hb, Down in Amyg"){
                top3_intersection_genes <- shared_up_hab_down_amy %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>%
                    mutate(min_p = pmin(adj.P.Val.hb, adj.P.Val.amyg)) %>% arrange(min_p) %>% slice(1:n) %>% pull(Symbol)
            }
            else if(group == "Shared: Down in Hb, Up in Amyg"){
                top3_intersection_genes <- shared_down_hab_up_amy %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>%
                    mutate(min_p = pmin(adj.P.Val.hb, adj.P.Val.amyg)) %>% arrange(min_p) %>% slice(1:n) %>% pull(Symbol)
            }
            else if(group == "Shared: Down in Hb, Down in Amyg"){
                top3_intersection_genes <- shared_down_hab_down_amy %>% dplyr::filter(Symbol %in% intersection_genes, !Symbol %in% interest_genes) %>%
                    mutate(min_p = pmin(adj.P.Val.hb, adj.P.Val.amyg)) %>% arrange(min_p) %>% slice(1:n) %>% pull(Symbol)
            }

            genes_2_show_x_term_x_group[[group]][[term]] = c(genes_of_interest_2_show, top3_intersection_genes)
        }

    }
}

l <- sapply(names(genes_2_show_x_term_x_group), function(set){
    imap_dfr(genes_2_show_x_term_x_group[[set]], ~ enframe(.x, value = "Symbol") %>% mutate(term = .y, DEGs_set = set))
})

l <- do.call(rbind, l)


## Add logFC of selected genes in Hb and Amyg
df <- l %>% left_join(results_Substance_uncorr_vars_habenula[[1]][, c("Symbol", "logFC", "t",  "P.Value", "adj.P.Val")], by = "Symbol", multiple = "any") %>% left_join(results_Substance_uncorr_vars_amygdala[[1]][, c("Symbol", "logFC", "t",  "P.Value", "adj.P.Val")], by = "Symbol", multiple = "any", suffix = c(".Hb", ".Amy"))

df$DEGs_set <- factor(df$DEGs_set, levels =  unique(df$DEGs_set))

## Order terms by alp order x group
df <- df %>% arrange(DEGs_set, term) %>% as.data.frame()

num_genes_x_term_x_group <- df %>% group_by(DEGs_set, term) %>% summarise(count = length(unique(Symbol)))
term_indices <- map2(.x = c(1, head(cumsum(num_genes_x_term_x_group$count), -1)+1),
                     .y = cumsum(num_genes_x_term_x_group$count),
                     ~(1:dim(df)[1])[.x:.y])
names(term_indices) <- num_genes_x_term_x_group$term

la = rowAnnotation(Group = df$DEGs_set)
ra = rowAnnotation(term = anno_block(align_to = term_indices,
                                     panel_fun = function(index, nm){
                                         grid.text(nm, rot = 0, just = "left", name = "term",
                                                   gp = gpar(fontsize = 7), x = 3.2)}))

h <- Heatmap(as.matrix(df[, c("logFC.Hb", "logFC.Amy")]),
        name = "logFC",
        border = T,
        row_labels = df$Symbol,
        row_names_side = "right",
        column_labels = c("logFC in Hb", "logFC in Amyg"),
        column_names_rot = 90,
        column_names_centered = F,
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 5, fontface = "italic"),
        row_split = df[, c("DEGs_set", "term")],
        gap = unit(0.75, "mm"),
        left_annotation = la,
        right_annotation = ra,
        row_title_gp = gpar(fontsize = 0),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        height = unit(35, "cm"),
        heatmap_width = unit(5, "cm"))


pdf(file = paste0("plots/06_GO_KEGG/GO_KEGG_heatmap_Hb_vs_Amyg.pdf"), height = 15, width = 10)
draw(h, heatmap_legend_side = "left")
dev.off()




df_wide <- l %>% select(Symbol, DEGs_set, term) %>%
    mutate(present = 1) %>%
    pivot_wider(
        names_from = term,
        values_from = present,
        values_fill = 0
    )

df_longer <- df_wide %>% pivot_longer(cols = setdiff(colnames(df_wide), c("Symbol", "DEGs_set")))

ggplot(df_longer, aes(x = Symbol, y = term, fill = value)) +
    geom_tile() +
    scale_fill_continuous(low="blue", high="red", name = "fold change")



## Reproducibility information

options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2024-04-23
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version    date (UTC) lib source
# abind                    1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
# AnnotationDbi          * 1.64.1     2023-11-02 [1] Bioconductor
# AnnotationHub            3.10.0     2023-10-26 [1] Bioconductor
# ape                      5.7-1      2023-03-13 [1] CRAN (R 4.3.0)
# aplot                    0.2.2      2023-10-06 [1] CRAN (R 4.3.1)
# Biobase                * 2.62.0     2023-10-26 [1] Bioconductor
# BiocFileCache            2.10.1     2023-10-26 [1] Bioconductor
# BiocGenerics           * 0.48.1     2023-11-02 [1] Bioconductor
# BiocManager              1.30.22    2023-08-08 [1] CRAN (R 4.3.0)
# BiocParallel             1.36.0     2023-10-26 [1] Bioconductor
# BiocVersion              3.18.1     2023-11-18 [1] Bioconductor 3.18 (R 4.3.2)
# Biostrings               2.70.2     2024-01-30 [1] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
# cachem                   1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
# cli                      3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# clusterProfiler        * 4.10.0     2023-11-06 [1] Bioconductor
# codetools                0.2-19     2023-02-01 [1] CRAN (R 4.3.2)
# colorspace               2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# cowplot                  1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon                   1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.2.1      2024-03-01 [1] CRAN (R 4.3.1)
# data.table               1.15.2     2024-02-29 [1] CRAN (R 4.3.1)
# DBI                      1.2.2      2024-02-16 [1] CRAN (R 4.3.2)
# dbplyr                   2.4.0      2023-10-26 [1] CRAN (R 4.3.1)
# DelayedArray             0.28.0     2023-11-06 [1] Bioconductor
# digest                   0.6.34     2024-01-11 [1] CRAN (R 4.3.1)
# DOSE                     3.28.2     2023-12-12 [1] Bioconductor 3.18 (R 4.3.2)
# dplyr                    1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# ellipsis                 0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
# enrichplot               1.22.0     2023-11-06 [1] Bioconductor
# evaluate                 0.23       2023-11-01 [1] CRAN (R 4.3.1)
# fansi                    1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# farver                   2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                  1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# fastmatch                1.1-4      2023-08-18 [1] CRAN (R 4.3.0)
# fgsea                    1.28.0     2023-10-26 [1] Bioconductor
# filelock                 1.0.3      2023-12-11 [1] CRAN (R 4.3.1)
# fs                       1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.38.6     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11     2024-02-17 [1] Bioconductor
# GenomicRanges          * 1.54.1     2023-10-30 [1] Bioconductor
# ggforce                  0.4.2      2024-02-19 [1] CRAN (R 4.3.1)
# ggfun                    0.1.4      2024-01-19 [1] CRAN (R 4.3.1)
# ggplot2                  3.5.0      2024-02-23 [1] CRAN (R 4.3.1)
# ggplotify                0.1.2      2023-08-09 [1] CRAN (R 4.3.0)
# ggraph                   2.2.0      2024-02-27 [1] CRAN (R 4.3.1)
# ggrepel                  0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# ggtree                   3.10.1     2024-02-27 [1] Bioconductor 3.18 (R 4.3.2)
# glue                     1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# GO.db                    3.18.0     2024-02-17 [1] Bioconductor
# GOSemSim                 2.28.1     2024-01-20 [1] Bioconductor 3.18 (R 4.3.2)
# graphlayouts             1.1.0      2024-01-19 [1] CRAN (R 4.3.1)
# gridExtra                2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gridGraphics             0.5-1      2020-12-13 [1] CRAN (R 4.3.0)
# gson                     0.1.0      2023-03-07 [1] CRAN (R 4.3.0)
# gtable                   0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# HDO.db                   0.99.1     2023-05-28 [1] Bioconductor
# here                   * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# htmltools                0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
# httpuv                   1.6.14     2024-01-26 [1] CRAN (R 4.3.1)
# httr                     1.4.7      2023-08-15 [1] CRAN (R 4.3.0)
# igraph                   2.0.2      2024-02-17 [1] CRAN (R 4.3.1)
# interactiveDisplayBase   1.40.0     2023-10-26 [1] Bioconductor
# IRanges                * 2.36.0     2023-10-26 [1] Bioconductor
# jsonlite                 1.8.8      2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.42.0     2023-10-26 [1] Bioconductor
# knitr                    1.45       2023-10-30 [1] CRAN (R 4.3.1)
# labeling                 0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# later                    1.3.2      2023-12-06 [1] CRAN (R 4.3.1)
# lattice                  0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
# lazyeval                 0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
# lifecycle                1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# magrittr                 2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                     7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                   1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.14.0     2023-10-26 [1] Bioconductor
# matrixStats            * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# memoise                  2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12       2021-09-28 [1] CRAN (R 4.3.0)
# munsell                  0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# org.Rn.eg.db           * 3.18.0     2024-02-17 [1] Bioconductor
# patchwork                1.2.0      2024-01-08 [1] CRAN (R 4.3.1)
# pillar                   1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plyr                     1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                      0.1-8      2022-11-29 [1] CRAN (R 4.3.0)
# polyclip                 1.10-6     2023-09-27 [1] CRAN (R 4.3.1)
# promises                 1.2.1      2023-08-10 [1] CRAN (R 4.3.0)
# purrr                    1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
# qvalue                   2.34.0     2023-10-26 [1] Bioconductor
# R6                       2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                    1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# reshape2                 1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# rlang                    1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown                2.26       2024-03-05 [1] CRAN (R 4.3.1)
# rprojroot                2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
# RSQLite                  2.3.5      2024-01-21 [1] CRAN (R 4.3.1)
# rstudioapi               0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays                 1.2.0      2023-10-26 [1] Bioconductor
# S4Vectors              * 0.40.2     2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                   1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# scatterpie               0.2.1      2023-06-07 [1] CRAN (R 4.3.0)
# sessioninfo            * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# shadowtext               0.1.3      2024-01-19 [1] CRAN (R 4.3.1)
# shiny                    1.8.0      2023-11-17 [1] CRAN (R 4.3.1)
# SparseArray              1.2.4      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# stringi                  1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr                  1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.32.0     2023-11-06 [1] Bioconductor
# tibble                   3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidygraph                1.3.1      2024-01-30 [1] CRAN (R 4.3.1)
# tidyr                    1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect               1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tidytree                 0.4.6      2023-12-12 [1] CRAN (R 4.3.1)
# treeio                   1.26.0     2023-11-06 [1] Bioconductor
# tweenr                   2.0.3      2024-02-26 [1] CRAN (R 4.3.1)
# utf8                     1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# viridis                  0.6.5      2024-01-29 [1] CRAN (R 4.3.1)
# viridisLite              0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                    3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# xfun                     0.42       2024-02-08 [1] CRAN (R 4.3.1)
# xtable                   1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.42.0     2023-10-26 [1] Bioconductor
# yaml                     2.3.8      2023-12-11 [1] CRAN (R 4.3.1)
# yulab.utils              0.1.4      2024-01-28 [1] CRAN (R 4.3.1)
# zlibbioc                 1.48.0     2023-10-26 [1] Bioconductor
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
