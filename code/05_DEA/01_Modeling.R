
library(here)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(rlang)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(ggrepel)
library(smplot2)
library(httr)
library(jsonlite)
library(xml2)
library(readxl)
library(sessioninfo)


####################   Differential Expression Analysis   ######################


################################################################################
##                                1. Modeling
################################################################################

load(here('raw-data/count_objects/rse_gene_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('processed-data/04_EDA/03_Explore_gene_level_effects/rse_gene_habenula_filt.Rdata'), verbose = TRUE)
load(here('processed-data/04_EDA/03_Explore_gene_level_effects/rse_gene_amygdala_filt.Rdata'), verbose = TRUE)
## Load sample data for additional covariates in the model
rat_behavioral_data <- as.data.frame(read.csv("raw-data/rat_behavioral_data.csv"))


## Extract previous output from calcNormFactors for count normalization (with all samples and genes)
norm_factors<-calcNormFactors(rse_gene, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)


## -----------------------------------------------------------------------------
##   1.1  DEA for Fentanyl vs. Saline with all samples from each brain region
## -----------------------------------------------------------------------------

DEA<- function(RSE, brain_region, formula, name, coef){

    ## Set "Saline" as the reference group for model.matrix()
    RSE$Substance <- factor(RSE$Substance, levels = c("Saline", "Fentanyl"))

    ## Previous lib sizes of each sample
    match_samples <- match(RSE$SAMPLE_ID, samples_factors$SAMPLE_ID)
    factors<-samples_factors[match_samples, ]

    pdf(file = paste("plots/05_DEA/01_Modeling/DEA_plots_", brain_region, "_", name, ".pdf", sep="" ))
    par(mfrow=c(2,3), mar = c(8,3.9,8,3.9))

    ## Model matrix
    model=model.matrix(formula, data=colData(RSE))
    ## Use previous norm factors to scale the raw library sizes
    RSE_scaled = calcNormFactors(RSE)
    RSE_scaled$samples$lib.size<-factors$lib.size
    RSE_scaled$samples$norm.factors<-factors$norm.factors

    ## Transform counts to log2(CPM)
    ## Estimate mean-variance relationship for each gene
    vGene = voom(RSE_scaled, design=model, plot=TRUE)

    ## Fit linear model for each gene
    fitGene = lmFit(vGene)

    ## Empirical Bayesian calculation to obtain the significant genes:
    ## compute moderated F and t-statistics, and log-odds of DE
    eBGene = eBayes(fitGene)

    ## Plot average log expression vs logFC
    limma::plotMA(eBGene, coef = coef, xlab = "Mean of normalized counts",
                  ylab="logFC")
    ## Plot -log(p-value) vs logFC
    volcanoplot(eBGene, coef = coef)

    ## Top-ranked genes for Substance (cases vs ctrls)
    top_genes = topTable(eBGene, coef=coef, p.value = 1, number=nrow(RSE), sort.by="none")
    ## Histogram of adjusted p-values
    hist(top_genes$adj.P.Val, xlab="FDR", main="")
    ## Histogram of p-values
    hist(top_genes$P.Value, xlab="p-value", main="")

    dev.off()

    return(list(top_genes, vGene, eBGene))

}

#######################
##  Habenula samples
#######################

## Model with all variables (without behavioral ones)
formula<- ~ Substance + Batch_RNA_extraction + Total_Num_Fentanyl_Sessions + mitoRate + totalAssignedGene +
            overallMapRate + concordMapRate + detected_num_genes + library_size + RIN + Total_RNA_amount +
            RNA_concentration
name<-"Substance_all_variables"
coef<-"SubstanceFentanyl"
results_all_vars_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_all_vars_habenula, file = 'processed-data/05_DEA/results_all_vars_habenula.Rdata')
## Number of DEGs (FDR<0.05)
length(which(results_all_vars_habenula[[1]]$adj.P.Val<0.05))
#  0

## Model with uncorrelated variables only
## formula <-  ~ Substance + Batch_RNA_extraction + overallMapRate + RIN + mitoRate (previous)
formula <-  ~ Substance + Batch_RNA_extraction + concordMapRate + RIN
name<-"Substance_uncorr_variables"
coef<-"SubstanceFentanyl"
results_uncorr_vars_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_uncorr_vars_habenula, file = 'processed-data/05_DEA/results_uncorr_vars_habenula.Rdata')
length(which(results_uncorr_vars_habenula[[1]]$adj.P.Val<0.05))
#  453


#######################
##  Amygdala samples
#######################

## Model with all variables
formula<- ~ Substance + Batch_RNA_extraction + Batch_lib_prep + Total_Num_Fentanyl_Sessions + mitoRate +
            totalAssignedGene + overallMapRate + concordMapRate + detected_num_genes + library_size + RIN +
            Total_RNA_amount + RNA_concentration
name<-"Substance_all_variables"
coef<-"SubstanceFentanyl"
results_all_vars_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_all_vars_amygdala, file = 'processed-data/05_DEA/results_all_vars_amygdala.Rdata')
## Number of DEGs (FDR<0.10)
length(which(results_all_vars_amygdala[[1]]$adj.P.Val<0.10))
#  0

## Model with uncorrelated variables only
## formula<- ~ Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate  + RIN
name<-"Substance_uncorr_variables"
coef<-"SubstanceFentanyl"
results_uncorr_vars_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_uncorr_vars_amygdala, file = 'processed-data/05_DEA/results_uncorr_vars_amygdala.Rdata')
length(which(results_uncorr_vars_amygdala[[1]]$adj.P.Val<0.05))
#  3041



## Plots for DEGs
plots_DEGs<-function(brain_region, top_genes, vGene, FDR, name) {

    if(name=='First_hr_infusion_slope'){
        FClab = 'log2FC(1st Hour Infusion Slope)'
    }
    else if (name=='Total_intake'){
        FClab = 'log2FC(Total Drug Intake)'
    }
    else if(name=='Last_session_intake'){
        FClab = 'log2FC(Last Session Intake)'
    }
    else{
        FClab='log2FC(Fentanyl vs. Saline)'
    }

    ## n.s./Down/Upregulated genes
    DE<-vector()
    for (i in 1:dim(top_genes)[1]) {
        if (top_genes$adj.P.Val[i]>FDR) {
            DE<-append(DE, "n.s.")
        }
        else {
            if (top_genes$logFC[i]>0) {
                DE<-append(DE, "Up")
            }
            else {
                DE<-append(DE, "Down")
            }
        }
    }
    top_genes$DE<- DE


    ## Gene symbols (or ensemblID if missing) for top 10 most up and downregulated DEGs
    top_genes$symbol_or_ensemblID <- unlist(apply(top_genes, 1, function(x){if(is.na(x['Symbol'])){x['ensemblID']} else{x['Symbol']}}))
    top_upDEGs <- subset(top_genes, logFC>0)[order(subset(top_genes, logFC>0)$adj.P.Val, decreasing = FALSE), 'symbol_or_ensemblID'][1:5]
    top_downDEGs <- subset(top_genes, logFC<0)[order(subset(top_genes, logFC<0)$adj.P.Val, decreasing = FALSE), 'symbol_or_ensemblID'][1:5]
    top_DEGs <- c(top_upDEGs, top_downDEGs)

    ## Position of caption in plot
    caption_x_units <- 0.55
    caption_y_units1 <- 0.15
    caption_y_units2 <- 0.08

    if (brain_region=='amygdala' & !name=='Total_intake'){
        caption_x_units <- 0.55
        caption_y_units1 <- 0.1
        caption_y_units2 <- -0.09
    }

    # ## Gene symbols for all 6 DEGs in habenula DEA for 1st hr infusion slope
    # else if(name=='First_hr_infusion_slope'){
    #     DEGs <- de_genes_FirstHrIntakeSlopeDEA_habenula$Symbol
    #     caption_x_units <- 1.2
    #     caption_y_units1 <- 0.4
    #     caption_y_units2 <- 0.31
    # }
    #
    # else if(name=='Total_intake'){
    #     caption_x_units <- 1.7
    #     caption_y_units1 <- 0.13
    #     caption_y_units2 <- 0.08
    # }
    top_genes$DEG_symbol<- sapply(top_genes$symbol_or_ensemblID, function(x){ if(x %in% top_DEGs){x} else {NA}})


    ## Plots
    cols <- c("Up" = "indianred2", "Down" = "steelblue2", "n.s." = "grey")
    sizes <- c("Up" = 1.3, "Down" = 1.3, "n.s." = 0.8)
    alphas <- c("Up" = 0.4, "Down" = 0.6, "n.s." = 0.5)

    ## MA plot for DE genes
    top_genes$mean_log_expr<-apply(vGene$E, 1, mean)
    p1<-ggplot(data = top_genes,
               aes(x = mean_log_expr,y = logFC,
                   fill = DE,
                   color = DE,
                   size = DE,
                   alpha = DE)) +
        sm_hgrid(legends = TRUE) +
        geom_point(shape = 21) +
        scale_color_manual(values = cols, name=NULL) +
        scale_fill_manual(values = cols, name=NULL) +
        scale_size_manual(values = sizes, name=NULL) +
        scale_alpha_manual(values = alphas, name=NULL) +
        labs(x="Mean of normalized counts", y=FClab) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"),
              legend.position = c(0.82, 0.15),
              legend.background = element_rect(fill=NA),
              legend.key.height = unit(0.15,"cm"),
              axis.title = element_text(size = (10)),
              legend.text = element_text(size=10))


    ## Volcano plot for DE genes
    p2<-ggplot(data = top_genes,
               aes(x = logFC,y = -log10(adj.P.Val),
                   color = DE,
                   fill = DE,
                   size = DE,
                   alpha = DE,
                   label= DEG_symbol)) +
        sm_hgrid(legends = TRUE) +
        geom_point(shape = 21) +
        geom_hline(yintercept = -log10(FDR),
                   linetype = "dashed", color = 'gray35', linewidth=0.5) +
        geom_vline(xintercept = c(-1,1),
                   linetype = "dashed", color = 'gray35', linewidth=0.5) +
        geom_text_repel(data = subset(top_genes, symbol_or_ensemblID %in% top_downDEGs), aes(fontface = 'bold'),
                        size=2.8,
                        color='black',
                        alpha = 1,
                        max.overlaps = Inf,
                        box.padding = 0.15, nudge_y = -0.1, nudge_x = -0.4,
                        segment.size = unit(0.4, 'mm'),
                        segment.alpha = 0.5,
                        show.legend=FALSE) +
        geom_text_repel(data = subset(top_genes, symbol_or_ensemblID %in% top_upDEGs), aes(fontface = 'bold'),
                        size=2.8,
                        color='black',
                        alpha = 1,
                        max.overlaps = Inf,
                        box.padding = 0.15, nudge_y = 0.1, nudge_x = 0.4,
                        segment.size = unit(0.4, 'mm'),
                        segment.alpha = 0.5,
                        show.legend=FALSE) +
        labs(y="-log10(FDR)", x=FClab)+
        scale_color_manual(values = cols, name=NULL) +
        scale_fill_manual(values = cols, name=NULL) +
        scale_size_manual(values = sizes, name=NULL) +
        scale_alpha_manual(values = alphas, name=NULL) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"),
              legend.position = c(0.13, 0.15),
              legend.background = element_rect(fill=NA),
              legend.key.height = unit(0.15,"cm"),
              axis.title = element_text(size = (10)),
              legend.text = element_text(size=10)) +
        ## Caption: number of DEGs
        annotate("text", x=max(top_genes$logFC)-caption_x_units, y=caption_y_units1, label= paste0(length(which(top_genes$adj.P.Val<FDR)), ' DEGs'),
                 color='gray40', size=2.8, fontface = 'bold') +
        ## Caption: FDR threshold
        annotate("text", x=max(top_genes$logFC)-caption_x_units, y=caption_y_units2, label= paste0("(FDR<", FDR, ")"),
                 color='gray40', size=2.5)

    plot_grid(p1, p2, ncol=2)
    ggsave(paste("plots/05_DEA/01_Modeling/DEG_plots_", brain_region, '_', name, ".pdf", sep=""),
           width = 22, height = 11, units = "cm")
}


## Plots for habenula DEGs from the model without correlated variables
plots_DEGs('habenula', top_genes = results_uncorr_vars_habenula[[1]], vGene = results_uncorr_vars_habenula[[2]], FDR = 0.05,
           name='Substance')
## Extract DEGs
de_genes_habenula <- results_uncorr_vars_habenula[[1]][which(results_uncorr_vars_habenula[[1]]$adj.P.Val<0.05),]


## Plots for amygdala DEGs from the model without correlated variables
plots_DEGs('amygdala', top_genes = results_uncorr_vars_amygdala[[1]], vGene = results_uncorr_vars_amygdala[[2]], FDR = 0.05,
           name='Substance')
de_genes_amygdala <- results_uncorr_vars_amygdala[[1]][which(results_uncorr_vars_amygdala[[1]]$adj.P.Val<0.05),]



## Add Ensembl info of DEGs

add_phenotypes <- function(de_genes){

    de_genes$associated_phenotypes <- NA

    for (i in 1:dim(de_genes)[1]){

        ## URL for the gene
        server <- "https://rest.ensembl.org"
        ext1 <- paste0("/phenotype/gene/Rattus_norvegicus/", de_genes$ensemblID[i], "?include_associated=1")
        r1 <- GET(paste(server, ext1, sep = ""), content_type("application/json"))

        ## Search and save associated phenotypes associated with variants reporting this gene
        associated_phenotypes <- unlist(fromJSON(toJSON(content(r1)))$description)
        de_genes$associated_phenotypes[i] <- paste(as.list(associated_phenotypes), collapse = '; ')
    }

    return(de_genes)
}

add_description <- function(de_genes){

    de_genes$gene_description <- NA

    for (i in 1:dim(de_genes)[1]){

        ## URL for the gene
        server <- "https://rest.ensembl.org"
        ext2 <- paste0("/lookup/id/", de_genes$ensemblID[i], "?expand=1;content-type=application/json")

        r2 <- GET(paste(server, ext2, sep = ""), content_type("application/json"))

        ## Search and save gene description
        description <- fromJSON(toJSON(content(r2)))$description
        if (length(description)==0){
            de_genes$gene_description[i] <- NA
        }
        else{
            de_genes$gene_description[i] <- description
        }

    }

    return(de_genes)
}


## Habenula
## Add associated phenotypes
de_genes_habenula <- add_phenotypes(de_genes_habenula)
## Add descriptions
de_genes_habenula <- add_description(de_genes_habenula)
de_genes_habenula$EntrezID <- as.character(de_genes_habenula$EntrezID)
## Order first by increasing FDR and secondly by decreasing |logFC|
de_genes_habenula <- de_genes_habenula[order(de_genes_habenula$adj.P.Val, -abs(de_genes_habenula$logFC)),]
save(de_genes_habenula, file = 'processed-data/05_DEA/de_genes_habenula.Rdata')
write.csv(de_genes_habenula, "generated_data/de_genes_habenula.csv")

## Amygdala
de_genes_amygdala <- add_phenotypes(de_genes_amygdala)
de_genes_amygdala <- add_description(de_genes_amygdala)
de_genes_amygdala$EntrezID <- as.character(de_genes_amygdala$EntrezID)
de_genes_amygdala <- de_genes_amygdala[order(de_genes_amygdala$adj.P.Val, -abs(de_genes_amygdala$logFC)),]
save(de_genes_amygdala, file = 'processed-data/05_DEA/de_genes_amygdala.Rdata')
write.csv(de_genes_amygdala, "generated_data/de_genes_amygdala.csv")





## -----------------------------------------------------------------------------
##      1.2  DEA for high/low fentanyl intake slope in fentanyl samples
## -----------------------------------------------------------------------------

##############################
## Habenula fentanyl samples
##############################

rse_gene_habenula_fent <- rse_gene_habenula_filt[,which(rse_gene_habenula_filt$Substance=='Fentanyl')]

## High/low intake samples
rse_gene_habenula_fent$Intake_slope_binary <- sapply(rse_gene_habenula_fent$Sample_Num, function(x){if(x %in% c(1,4,7)){'High'}
                                                                                             else if(x %in% c(2,5,6,8)){'Low'}
                                                                                             else {NA}})
## Remove sample from the outlier rat (in intake slope)
## Identify it
outlier_fent_sample <- colData(rse_gene_habenula_fent)[ which(is.na(rse_gene_habenula_fent$Intake_slope_binary)), 'Rat_ID']
# [1] "LgA 09"
rse_gene_habenula_fent <- rse_gene_habenula_fent[,-which(is.na(rse_gene_habenula_fent$Intake_slope_binary))]

## Same variables for habenula DEA
## Substance and batches are the same for all these samples
formula<- ~ Intake_slope + overallMapRate + RIN + mitoRate
name<-"Intake_slope_Habenula"
coef<-"Intake_slopeLow"
results_intake_slope_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
## DEGs (FDR<0.05)
de_genes_intake_slope_habenula <- results_intake_slope_habenula[[1]][which(results_intake_slope_habenula[[1]]$adj.P.Val<0.05), ]
save(results_intake_slope_habenula, file = 'processed-data/05_DEA/results_intake_slope_habenula.Rdata')
## Number of DEGs
dim(de_genes_intake_slope_habenula)[1]
##  65

## Add associated phenotypes and descriptions of DEGs
de_genes_intake_slope_habenula <- add_phenotypes(de_genes_intake_slope_habenula)
de_genes_intake_slope_habenula <- add_description(de_genes_intake_slope_habenula)
de_genes_intake_slope_habenula$EntrezID <- as.character(de_genes_intake_slope_habenula$EntrezID)
## Order by FDR and |logFC|
de_genes_intake_slope_habenula <- de_genes_intake_slope_habenula[order(de_genes_intake_slope_habenula$adj.P.Val, -abs(de_genes_intake_slope_habenula$logFC)),]
save(de_genes_intake_slope_habenula, file = 'processed-data/05_DEA/de_genes_intake_slope_habenula.Rdata')
write.csv(de_genes_intake_slope_habenula, "generated_data/de_genes_intake_slope_habenula.csv")


##############################
## Amygdala fentanyl samples
##############################

rse_gene_amygdala_fent <- rse_gene_amygdala_filt[,which(rse_gene_amygdala_filt$Substance=='Fentanyl')]

## High/low intake samples
rse_gene_amygdala_fent$Intake_slope <- sapply(rse_gene_amygdala_fent$Sample_Num, function(x){if(x %in% c(17,20,23)){'High'}
                                                                                             else if(x %in% c(18,21,22,24)){'Low'}
                                                                                             else {NA}})
## Remove sample from the outlier rat (in intake slope)
rse_gene_amygdala_fent <- rse_gene_amygdala_fent[,-which(is.na(rse_gene_amygdala_fent$Intake_slope))]

## Same variables for amygdala DEA
## Substance and batches are the same
formula<- ~ Intake_slope + overallMapRate + RIN + mitoRate
name<-"Intake_slope_Amygdala"
coef<-"Intake_slopeLow"
results_intake_slope_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
## DEGs (FDR<0.05)
de_genes_intake_slope_amygdala <- results_intake_slope_amygdala[[1]][which(results_intake_slope_amygdala[[1]]$adj.P.Val<0.05), ]
save(results_intake_slope_amygdala, file = 'processed-data/05_DEA/results_intake_slope_amygdala.Rdata')
## Number of DEGs
dim(de_genes_intake_slope_amygdala)[1]
##  2

## Add associated phenotypes and descriptions of DEGs
de_genes_intake_slope_amygdala <- add_phenotypes(de_genes_intake_slope_amygdala)
de_genes_intake_slope_amygdala <- add_description(de_genes_intake_slope_amygdala)
de_genes_intake_slope_amygdala$EntrezID <- as.character(de_genes_intake_slope_amygdala$EntrezID)
## Order by FDR and |logFC|
de_genes_intake_slope_amygdala <- de_genes_intake_slope_amygdala[order(de_genes_intake_slope_amygdala$adj.P.Val, -abs(de_genes_intake_slope_amygdala$logFC)),]
save(de_genes_intake_slope_amygdala, file = 'processed-data/05_DEA/de_genes_intake_slope_amygdala.Rdata')
write.csv(de_genes_intake_slope_amygdala, "generated_data/de_genes_intake_slope_amygdala.csv")





## -----------------------------------------------------------------------------
##                 1.3  DEA for fentanyl vs. saline with covariates:
##       '1st hour intake slope', 'Total intake' and 'Last session intake'
## -----------------------------------------------------------------------------

colnames(covariate_data) <- gsub(' ', '_', colnames(covariate_data))


############### 1.3.1 DEA with covariate '1st hour intake slope' ###############

#####################
## Habenula samples
#####################

## Add info of 1st hour intake slope for each sample (CAPITALIZE behav variables)
rse_gene_habenula_filt$First_hr_infusion_slope <- sapply(rse_gene_habenula_filt$Rat_ID,
                                                         function(x){covariate_data[which(covariate_data$Rat_ID==x), '1st_Hour_Infusion_Slope']})

## Same uncorrelated variables as before
formula <-  ~ Substance + First_hr_infusion_slope + Batch_RNA_extraction + overallMapRate + RIN + mitoRate
name <-"with_First_hr_infusion_slope"
coef <-"SubstanceFentanyl"
results_SubstanceDEA_FirstHrIntakeSlope_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_SubstanceDEA_FirstHrIntakeSlope_habenula, file = 'processed-data/05_DEA/results_SubstanceDEA_FirstHrIntakeSlope_habenula.Rdata')
length(which(results_SubstanceDEA_FirstHrIntakeSlope_habenula[[1]]$adj.P.Val<0.05))
#  0

#####################
## Amygdala samples
#####################

rse_gene_amygdala_filt$First_hr_infusion_slope <- sapply(rse_gene_amygdala_filt$Rat_ID,
                                                         function(x){covariate_data[which(covariate_data$Rat_ID==x), '1st_Hour_Infusion_Slope']})

## Same uncorrelated variables as before
formula <-  ~ Substance + First_hr_infusion_slope + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate
name <-"with_First_hr_infusion_slope"
coef <-"SubstanceFentanyl"
results_SubstanceDEA_FirstHrIntakeSlope_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_SubstanceDEA_FirstHrIntakeSlope_amygdala, file = 'processed-data/05_DEA/results_SubstanceDEA_FirstHrIntakeSlope_amygdala.Rdata')
length(which(results_SubstanceDEA_FirstHrIntakeSlope_amygdala[[1]]$adj.P.Val<0.05))
#  0



##################   1.3.2 DEA with covariate 'total intake'  ##################

#####################
## Habenula samples
#####################

## Add info of total intake for each sample
rse_gene_habenula_filt$Total_intake <- sapply(rse_gene_habenula_filt$Rat_ID,
                                                         function(x){covariate_data[which(covariate_data$Rat_ID==x), 'Total_Intake_(mg)']})

## Same uncorrelated variables as before
formula <-  ~ Substance + Total_intake + Batch_RNA_extraction + overallMapRate + RIN + mitoRate
name <-"with_Total_intake"
coef <-"SubstanceFentanyl"
results_SubstanceDEA_TotalIntake_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_SubstanceDEA_TotalIntake_habenula, file = 'processed-data/05_DEA/results_SubstanceDEA_TotalIntake_habenula.Rdata')
length(which(results_SubstanceDEA_TotalIntake_habenula[[1]]$adj.P.Val<0.05))
#  0

#####################
## Amygdala samples
#####################

rse_gene_amygdala_filt$Total_intake <- sapply(rse_gene_amygdala_filt$Rat_ID,
                                              function(x){covariate_data[which(covariate_data$Rat_ID==x), 'Total_Intake_(mg)']})

formula <-  ~ Substance + Total_intake + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate
name <-"with_Total_intake"
coef <-"SubstanceFentanyl"
results_SubstanceDEA_TotalIntake_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_SubstanceDEA_TotalIntake_amygdala, file = 'processed-data/05_DEA/results_SubstanceDEA_TotalIntake_amygdala.Rdata')
length(which(results_SubstanceDEA_TotalIntake_amygdala[[1]]$adj.P.Val<0.05))
#  0



##############   1.3.3 DEA with covariate 'last session intake'  ###############

#####################
## Habenula samples
#####################

## Add info of last session intake for each sample
rse_gene_habenula_filt$Last_session_intake<- sapply(rse_gene_habenula_filt$Rat_ID,
                                              function(x){covariate_data[which(covariate_data$Rat_ID==x), 'Last_Session_Intake_(mg)']})

## Same uncorrelated variables as before
formula <-  ~ Substance + Last_session_intake + Batch_RNA_extraction + overallMapRate + RIN + mitoRate
name <-"with_Last_session_intake"
coef <-"SubstanceFentanyl"
results_SubstanceDEA_LastSessionIntake_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_SubstanceDEA_LastSessionIntake_habenula, file = 'processed-data/05_DEA/results_SubstanceDEA_LastSessionIntake_habenula.Rdata')
length(which(results_SubstanceDEA_LastSessionIntake_habenula[[1]]$adj.P.Val<0.05))
#  0

#####################
## Amygdala samples
#####################

rse_gene_amygdala_filt$Last_session_intake <- sapply(rse_gene_amygdala_filt$Rat_ID,
                                              function(x){covariate_data[which(covariate_data$Rat_ID==x), 'Last_Session_Intake_(mg)']})

formula <-  ~ Substance + Last_session_intake + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate
name <-"with_Last_session_intake"
coef <-"SubstanceFentanyl"
results_SubstanceDEA_LastSessionIntake_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_SubstanceDEA_LastSessionIntake_amygdala, file = 'processed-data/05_DEA/results_SubstanceDEA_LastSessionIntake_amygdala.Rdata')
length(which(results_SubstanceDEA_LastSessionIntake_amygdala[[1]]$adj.P.Val<0.05))
#  0



## Plots of correlation between the additional covariates and the rest of the sample variables
covariates_CCA<- function(brain_region){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## Define variables
    if (brain_region == 'habenula'){
        formula = ~ Substance + Batch_RNA_extraction + Total_Num_Fentanyl_Sessions +
            mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
            library_size + detected_num_genes + RNA_concentration + Total_RNA_amount +
            First_hr_infusion_slope + Total_intake + Last_session_intake

    }
    else {
        formula = ~ Substance + Batch_RNA_extraction + Batch_lib_prep + Total_Num_Fentanyl_Sessions +
            mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
            library_size + detected_num_genes + RNA_concentration + Total_RNA_amount +
            First_hr_infusion_slope + Total_intake + Last_session_intake
    }

    ## Correlations
    C=canCorPairs(formula, colData(RSE))

    ## Heatmap
    pheatmap(
        C,
        color = hcl.colors(50, "YlOrRd", rev = TRUE),
        fontsize=8,
        border_color = "black",
        height = 6,
        width = 6.5,
        cluster_rows=TRUE,
        cluster_cols=TRUE,
        filename=paste("plots/04_EDA/03_Explore_gene_level_effects/02_Var_Partition/Covariates_corr_plot_", brain_region, ".pdf", sep="")
    )

    return(C)
}

Covariates_corr_habenula<- covariates_CCA('habenula')
Covariates_corr_amygdala <- covariates_CCA('amygdala')





## ---------------------------------------------------------------------------------
##   1.4  DEA for '1st hour intake slope' in habenula and amygdala fentanyl samples
## ---------------------------------------------------------------------------------

####################  1.4.1 Analysis with all fentanyl samples  ####################

##############################
## Habenula fentanyl samples
##############################

## Fentanyl samples only
rse_gene_habenula_fent <- rse_gene_habenula_filt[,which(rse_gene_habenula_filt$Substance=='Fentanyl')]

## formula <- ~ First_hr_infusion_slope + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  First_Hour_Infusion_Slope + RIN + RNA_concentration + overallMapRate + totalAssignedGene
name <-"for_First_Hour_Infusion_Slope"
## New contrast of interest
coef <-"First_Hour_Infusion_Slope"
results_FirstHrIntakeSlopeDEA_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
save(results_FirstHrIntakeSlopeDEA_habenula, file = 'processed-data/05_DEA/results_FirstHrIntakeSlopeDEA_habenula.Rdata')

## DEGs (FDR<0.05)
length(which(results_FirstHrIntakeSlopeDEA_habenula[[1]]$adj.P.Val<0.05))
#  6 (previous)
#  0

# de_genes_FirstHrIntakeSlopeDEA_habenula<- results_FirstHrIntakeSlopeDEA_habenula[[1]][which(results_FirstHrIntakeSlopeDEA_habenula[[1]]$adj.P.Val<0.05), ]
#
# ## Add associated phenotypes and descriptions of DEGs
# de_genes_FirstHrIntakeSlopeDEA_habenula <- add_phenotypes(de_genes_FirstHrIntakeSlopeDEA_habenula)
# de_genes_FirstHrIntakeSlopeDEA_habenula <- add_description(de_genes_FirstHrIntakeSlopeDEA_habenula)
# de_genes_FirstHrIntakeSlopeDEA_habenula$EntrezID <- as.character(de_genes_FirstHrIntakeSlopeDEA_habenula$EntrezID)
# ## Order by FDR and |logFC|
# de_genes_FirstHrIntakeSlopeDEA_habenula <- de_genes_FirstHrIntakeSlopeDEA_habenula[order(de_genes_FirstHrIntakeSlopeDEA_habenula$adj.P.Val, -abs(de_genes_FirstHrIntakeSlopeDEA_habenula$logFC)),]
# save(de_genes_FirstHrIntakeSlopeDEA_habenula, file = 'processed-data/05_DEA/de_genes_FirstHrIntakeSlopeDEA_habenula.Rdata')
# write.csv(de_genes_FirstHrIntakeSlopeDEA_habenula, "generated_data/de_genes_FirstHrIntakeSlopeDEA_habenula.csv")
#
# ## Plots for DEGs
# plots_DEGs('habenula', top_genes = results_FirstHrIntakeSlopeDEA_habenula[[1]], vGene = results_FirstHrIntakeSlopeDEA_habenula[[2]], FDR = 0.05,
#            name='habenula_for_First_hr_infusion_slope')


##############################
## Amygdala fentanyl samples
##############################

rse_gene_amygdala_fent <- rse_gene_amygdala_filt[,which(rse_gene_amygdala_filt$Substance=='Fentanyl')]

## DEA
## formula <-  ~ First_hr_infusion_slope + overallMapRate + RIN + mitoRate
formula <-  ~  First_Hour_Infusion_Slope + RIN + mitoRate
name <-"for_First_Hour_Infusion_Slope"
coef <-"First_Hour_Infusion_Slope"
results_FirstHrIntakeSlopeDEA_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
save(results_FirstHrIntakeSlopeDEA_amygdala, file = 'processed-data/05_DEA/results_FirstHrIntakeSlopeDEA_amygdala.Rdata')
length(which(results_FirstHrIntakeSlopeDEA_amygdala[[1]]$adj.P.Val<0.1))
#  0



#########  1.4.2 Analysis without samples from negative outlier fentanyl rat ########

##############################
## Habenula fentanyl samples
##############################

## Remove outlier rat sample in habenula
rse_gene_habenula_fent_withoutOutlier <- rse_gene_habenula_fent[,-which(rse_gene_habenula_fent$Rat_ID==outlier_fent_sample)]

## DEA
formula <- ~ First_hr_infusion_slope + overallMapRate + RIN + mitoRate
name <-"for_First_hr_infusion_slope_withoutOutlier"
coef <-"First_hr_infusion_slope"
results_FirstHrIntakeSlopeDEA_habenula_withoutOutlier<-DEA(rse_gene_habenula_fent_withoutOutlier, 'habenula', formula, name, coef)
save(results_FirstHrIntakeSlopeDEA_habenula_withoutOutlier, file = 'processed-data/05_DEA/results_FirstHrIntakeSlopeDEA_habenula_withoutOutlier.Rdata')

## DEGs (FDR<0.1)
length(which(results_FirstHrIntakeSlopeDEA_habenula_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0

##############################
## Amygdala fentanyl samples
##############################

rse_gene_amygdala_fent_withoutOutlier <- rse_gene_amygdala_fent[,-which(rse_gene_amygdala_fent$Rat_ID==outlier_fent_sample)]

formula <- ~ First_hr_infusion_slope + overallMapRate + RIN + mitoRate
name <-"for_First_hr_infusion_slope_withoutOutlier"
coef <-"First_hr_infusion_slope"
results_FirstHrIntakeSlopeDEA_amygdala_withoutOutlier<-DEA(rse_gene_amygdala_fent_withoutOutlier, 'amygdala', formula, name, coef)
save(results_FirstHrIntakeSlopeDEA_amygdala_withoutOutlier, file = 'processed-data/05_DEA/results_FirstHrIntakeSlopeDEA_amygdala_withoutOutlier.Rdata')
length(which(results_FirstHrIntakeSlopeDEA_amygdala_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0





## -----------------------------------------------------------------------------
##    1.5  DEA for 'Total intake' in habenula and amygdala fentanyl samples
## -----------------------------------------------------------------------------

##################  1.5.1 Analysis with all fentanyl samples  ##################

#####################
## Habenula samples
#####################

##formula <-  ~ Total_Intake + overallMapRate + RIN + mitoRate (previous)
formula <- ~  Total_Intake + RIN + RNA_concentration + overallMapRate
name <-"for_Total_intake"
coef <-"Total_Intake"
results_TotalIntakeDEA_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
save(results_TotalIntakeDEA_habenula, file = 'processed-data/05_DEA/results_TotalIntakeDEA_habenula.Rdata')
length(which(results_TotalIntakeDEA_habenula[[1]]$adj.P.Val<0.1))
#  0

#####################
## Amygdala samples
#####################

## formula <-  ~ Total_intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Total_Intake + RIN + mitoRate
name <-"for_Total_intake"
coef <-"Total_Intake"
results_TotalIntakeDEA_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
save(results_TotalIntakeDEA_amygdala, file = 'processed-data/05_DEA/results_TotalIntakeDEA_amygdala.Rdata')
length(which(results_TotalIntakeDEA_amygdala[[1]]$adj.P.Val<0.1))
#  26
#. 0
length(which(results_TotalIntakeDEA_amygdala[[1]]$adj.P.Val<0.05))
#  0

## DEGs (FDR<0.1)
# de_genes_TotalIntakeDEA_amygdala<- results_TotalIntakeDEA_amygdala[[1]][which(results_TotalIntakeDEA_amygdala[[1]]$adj.P.Val<0.1), ]
# ## Add associated phenotypes and descriptions of DEGs
# de_genes_TotalIntakeDEA_amygdala<- add_phenotypes(de_genes_TotalIntakeDEA_amygdala)
# de_genes_TotalIntakeDEA_amygdala <- add_description(de_genes_TotalIntakeDEA_amygdala)
# de_genes_TotalIntakeDEA_amygdala$EntrezID <- as.character(de_genes_TotalIntakeDEA_amygdala$EntrezID)
# de_genes_TotalIntakeDEA_amygdala<- de_genes_TotalIntakeDEA_amygdala[order(de_genes_TotalIntakeDEA_amygdala$adj.P.Val, -abs(de_genes_TotalIntakeDEA_amygdala$logFC)),]
# save(de_genes_TotalIntakeDEA_amygdala, file = 'processed-data/05_DEA/de_genes_TotalIntakeDEA_amygdala.Rdata')
# write.csv(de_genes_TotalIntakeDEA_amygdala, "generated_data/de_genes_TotalIntakeDEA_amygdala.csv")
#
# ## Plots for DEGs
# plots_DEGs('amygdala', top_genes = results_TotalIntakeDEA_amygdala[[1]], vGene = results_TotalIntakeDEA_amygdala[[2]], FDR = 0.1,
#            name='amygdala_for_Total_intake')



########  1.5.2 Analysis without samples from negative outlier fentanyl rat #######

##############################
## Habenula fentanyl samples
##############################

formula <-  ~ Total_intake + overallMapRate + RIN + mitoRate
name <-"for_Total_intake_withoutOutlier"
coef <-"Total_intake"
results_TotalIntakeDEA_habenula_withoutOutlier<-DEA(rse_gene_habenula_fent_withoutOutlier, 'habenula', formula, name, coef)
save(results_TotalIntakeDEA_habenula_withoutOutlier, file = 'processed-data/05_DEA/results_TotalIntakeDEA_habenula_withoutOutlier.Rdata')
length(which(results_TotalIntakeDEA_habenula_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0


##############################
## Amygdala fentanyl samples
##############################

formula <-  ~ Total_intake + overallMapRate + RIN + mitoRate
name <-"for_Total_intake_withoutOutlier"
coef <-"Total_intake"
results_TotalIntakeDEA_amygdala_withoutOutlier<-DEA(rse_gene_amygdala_fent_withoutOutlier, 'amygdala', formula, name, coef)
save(results_TotalIntakeDEA_amygdala_withoutOutlier, file = 'processed-data/05_DEA/results_TotalIntakeDEA_amygdala_withoutOutlier.Rdata')
length(which(results_TotalIntakeDEA_amygdala_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0





## -----------------------------------------------------------------------------
#  1.6  DEA for 'Last Session Intake' in habenula and amygdala fentanyl samples
## -----------------------------------------------------------------------------

##################  1.6.1 Analysis with all fentanyl samples  ##################

#####################
## Habenula samples
#####################

# formula <-  ~ Last_session_intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Last_Session_Intake + RIN + RNA_concentration + overallMapRate + totalAssignedGene
name <-"for_Last_Session_Intake"
coef <-"Last_Session_Intake"
results_LastSessionIntakeDEA_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
save(results_LastSessionIntakeDEA_habenula, file = 'processed-data/05_DEA/results_LastSessionIntakeDEA_habenula.Rdata')
length(which(results_LastSessionIntakeDEA_habenula[[1]]$adj.P.Val<0.1))
#  0
# 257

#####################
## Amygdala samples
#####################

# formula <-  ~ Last_session_intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Last_Session_Intake + RIN + totalAssignedGene + concordMapRate
name <-"for_Last_Session_Intake"
coef <-"Last_Session_Intake"
results_LastSessionIntakeDEA_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
save(results_LastSessionIntakeDEA_amygdala, file = 'processed-data/05_DEA/results_LastSessionIntakeDEA_amygdala.Rdata')
length(which(results_LastSessionIntakeDEA_amygdala[[1]]$adj.P.Val<0.1))
#  0



#########  1.6.2 Analysis without samples from negative outlier fentanyl rat ########

##############################
## Habenula fentanyl samples
##############################

formula <-  ~ Last_session_intake + overallMapRate + RIN + mitoRate
name <-"for_Last_session_intake_withoutOutlier"
coef <-"Last_session_intake"
results_LastSessionIntakeDEA_habenula_withoutOutlier<-DEA(rse_gene_habenula_fent_withoutOutlier, 'habenula', formula, name, coef)
save(results_LastSessionIntakeDEA_habenula_withoutOutlier, file = 'processed-data/05_DEA/results_LastSessionIntakeDEA_habenula_withoutOutlier.Rdata')
length(which(results_LastSessionIntakeDEA_habenula_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0

#####################
## Amygdala samples
#####################

formula <-  ~ Last_session_intake + overallMapRate + RIN + mitoRate
name <-"for_Last_session_intake_withoutOutlier"
coef <-"Last_session_intake"
results_LastSessionIntakeDEA_amygdala_withoutOutlier<-DEA(rse_gene_amygdala_fent_withoutOutlier, 'amygdala', formula, name, coef)
save(results_LastSessionIntakeDEA_amygdala_withoutOutlier, file = 'processed-data/05_DEA/results_LastSessionIntakeDEA_amygdala_withoutOutlier.Rdata')
length(which(results_LastSessionIntakeDEA_amygdala_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0




## Plot gene expression vs total/last/slope of drug intake

plot_gene_expr_vs_intake <- function(brain_region, sample_var, gene_id){

    rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'fent', sep="_")))

    if(sample_var=='First_hr_infusion_slope'){
        x_label="First hr infusion slope"
    }

    else if(sample_var=='Total_intake'){
        x_label="Total drug intake"
    }

    else if(sample_var=='Last_session_intake'){
        x_label="Last session drug intake"
    }

    ## Lognorm counts of the gene across samples
    data <- colData(rse_gene)
    data$gene_expr <- assays(rse_gene)$logcounts[gene_id,]

    ## Gene symbol
    gene_symbol <- rowData(rse_gene)[which(rowData(rse_gene)$ensemblID==gene_id), 'Symbol']
    if (is.na(gene_symbol)){
        gene_ids <- gene_id
    }
    else{
        gene_ids <- paste(gene_symbol, gene_id, sep='-')
    }

    ## Color for the outlier sample and the rest
    data$is_outlier <- sapply(data$Rat_ID, function(x){if(x==outlier_fent_sample){'coral2'} else {'gray50'}})

    plot <- ggplot(as.data.frame(data), aes(x=eval(parse_expr(sample_var)), y=gene_expr)) +
        geom_point(aes(color=is_outlier)) +
        stat_smooth (geom="line", alpha=0.4, size=0.6, span=0.3, method = lm, color='orangered4') +
        theme_bw() +
        guides(color="none") +
        labs(title = gene_ids,
             y= 'lognorm counts', x = x_label) +
        theme(plot.margin=unit (c (1,1,1,1), 'cm'),
              axis.title = element_text(size = (7)),
              axis.text = element_text(size = (6)),
              plot.title = element_text(hjust=0.5, size=7.5, face="bold"),
              legend.text = element_text(size=6),
              legend.title = element_text(size=7))


    return(plot)
}


##############################
## Habenula fentanyl samples
##############################

## Top most significant gene from DEA for 1st hr infusion slope
dea_results <- results_FirstHrIntakeSlopeDEA_habenula[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p1 <- plot_gene_expr_vs_intake('habenula', 'First_hr_infusion_slope', gene_id)

## Top most significant gene from DEA for total intake
dea_results <- results_TotalIntakeDEA_habenula[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p2 <- plot_gene_expr_vs_intake('habenula', 'Total_intake', gene_id)

## Top most significant gene from DEA for last intake
dea_results <- results_LastSessionIntakeDEA_habenula[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p3 <- plot_gene_expr_vs_intake('habenula', 'Last_session_intake', gene_id)

plot_grid(p1, p2, p3, nrow=1)
ggsave('plots/05_DEA/01_Modeling/geneExpr_VS_drugIntake_habenula_withOutlier.pdf', width = 20, height = 7, units = "cm")



##############################
## Amygdala fentanyl samples
##############################

## Top most significant gene from DEA for 1st hr infusion slope
dea_results <- results_FirstHrIntakeSlopeDEA_amygdala[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p1 <- plot_gene_expr_vs_intake('amygdala', 'First_hr_infusion_slope', gene_id)

## Top most significant gene from DEA for total intake
dea_results <- results_TotalIntakeDEA_amygdala[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p2 <- plot_gene_expr_vs_intake('amygdala', 'Total_intake', gene_id)

## Top most significant gene from DEA for last intake
dea_results <- results_LastSessionIntakeDEA_amygdala[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p3 <- plot_gene_expr_vs_intake('amygdala', 'Last_session_intake', gene_id)

plot_grid(p1, p2, p3, nrow=1)
ggsave('plots/05_DEA/01_Modeling/geneExpr_VS_drugIntake_amygdala_withOutlier.pdf', width = 20, height = 7, units = "cm")







## Reproducibility information

options(width = 120)
session_info()

# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Tijuana
# date     2023-06-05
# rstudio  2023.03.1+446 Cherry Blossom (desktop)
# pandoc   NA

