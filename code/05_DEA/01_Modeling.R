
library(here)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(rlang)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(smplot2)
library(sessioninfo)


########################   Differential Expression Analysis   #########################

load(here('raw-data/count_objects/rse_gene_Jlab_experiment_n33.Rdata'), verbose=TRUE)
load(here('processed-data/04_EDA/02_PCA/rse_gene_habenula_filt.Rdata'), verbose = TRUE)
load(here('processed-data/04_EDA/02_PCA/rse_gene_amygdala_filt.Rdata'), verbose = TRUE)


################################################################################
##                                 1. Modeling
################################################################################

## Extract previous output from calcNormFactors for count normalization (with all samples and genes)
norm_factors<-calcNormFactors(rse_gene, method = "TMM")
samples_factors<-data.frame(SAMPLE_ID=norm_factors$samples$SAMPLE_ID,
                            norm.factors=norm_factors$samples$norm.factors,
                            lib.size=norm_factors$samples$lib.size)


## DEA for habenula and amygdala samples
DEA<- function(brain_region, formula, name, coef){

    RSE <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## Previous lib sizes of each sample
    match_samples <- match(RSE$SAMPLE_ID, samples_factors$SAMPLE_ID)
    factors<-samples_factors[match_samples, ]

    pdf(file = paste("plots/05_DEA/01_Modeling/DEA_plots_", brain_region, "_", name, ".pdf", sep="" ))
    par(mfrow=c(2,2))

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
    ## Histogram of adjusted p values
    hist(top_genes$adj.P.Val, xlab="FDR", main="")

    dev.off()

    return(list(top_genes, vGene, eBGene))

}


## Habenula samples

## Model with all variables
formula<- ~ Substance + Batch_RNA_extraction + Total_Num_Fentanyl_Sessions + mitoRate + totalAssignedGene + overallMapRate + concordMapRate + detected_num_genes + library_size + RIN + Total_RNA_amount + RNA_concentration
name<-"all_variables"
coef<-"SubstanceSaline"
results_all_vars_habenula<-DEA('habenula', formula, name, coef)
save(results_all_vars_habenula, file = 'processed-data/05_DEA/results_all_vars_habenula.Rdata')
## Number of DEG (FDR<0.10)
length(which(results_all_vars_habenula[[1]]$adj.P.Val<0.05))
#  0

## Model with uncorrelated variables only
formula <-  ~ Substance + Batch_RNA_extraction + overallMapRate + RIN + mitoRate
name<-"uncorr_variables"
coef<-"SubstanceSaline"
results_uncorr_vars_habenula<-DEA('habenula', formula, name, coef)
save(results_uncorr_vars_habenula, file = 'processed-data/05_DEA/results_uncorr_vars_habenula.Rdata')
length(which(results_uncorr_vars_habenula[[1]]$adj.P.Val<0.05))
#  88



## Amygdala samples

## Model with all variables
formula<- ~ Substance + Batch_RNA_extraction + Batch_lib_prep + Total_Num_Fentanyl_Sessions + mitoRate + totalAssignedGene + overallMapRate + concordMapRate + detected_num_genes + library_size + RIN + Total_RNA_amount + RNA_concentration
name<-"all_variables"
coef<-"SubstanceSaline"
results_all_vars_amygdala<-DEA('amygdala', formula, name, coef)
save(results_all_vars_amygdala, file = 'processed-data/05_DEA/results_all_vars_amygdala.Rdata')
## Number of DEG (FDR<0.10)
length(which(results_all_vars_amygdala[[1]]$adj.P.Val<0.10))
#  0

## Model with uncorrelated variables only
formula<- ~ ~ Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate
name<-"uncorr_variables"
coef<-"SubstanceSaline"
results_uncorr_vars_amygdala<-DEA('amygdala', formula, name, coef)
save(results_uncorr_vars_amygdala, file = 'processed-data/05_DEA/results_uncorr_vars_amygdala.Rdata')
length(which(results_uncorr_vars_amygdala[[1]]$adj.P.Val<0.05))
#  2728





## Plots for DEGs
plots_DEGs<-function(top_genes, vGene, FDR, name) {

    ## |logFC| for gene labels in volcanoPlot
    if (length(which(top_genes$adj.P.Val<0.05))<100){
        logFC_abs = 1
    }
    else{
        logFC_abs = 1.5
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

    ## Gene symbols for DEGs with |logFC|>logFC_abs
    DEG_symbol<-vector()
    for (i in 1:dim(top_genes)[1]) {
        if (top_genes$DE[i]!="n.s." & abs(top_genes$logFC[i])>logFC_abs) {
            DEG_symbol<-append(DEG_symbol, top_genes$Symbol[i])
        }
        else {
            DEG_symbol<-append(DEG_symbol, NA)
        }
    }
    top_genes$DEG_symbol<- DEG_symbol

    ## Plots
    cols <- c("Up" = "red3", "Down" = "steelblue2", "n.s." = "grey")
    sizes <- c("Up" = 2, "Down" = 2, "n.s." = 1)
    alphas <- c("Up" = 1, "Down" = 1, "n.s." = 0.5)

    ## MA plot for DE genes
    top_genes$mean_log_expr<-apply(vGene$E, 1, mean)
    p1<-ggplot(data = top_genes,
               aes(x = mean_log_expr,y = logFC,
                   fill = DE,
                   size = DE,
                   alpha = DE)) +
        sm_hgrid(legends = TRUE) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) +
        geom_point(shape = 21) +
        scale_fill_manual(values = cols) +
        scale_size_manual(values = sizes) +
        scale_alpha_manual(values = alphas) +
        labs(x="Mean of normalized counts")


    ## Volcano plot for DE genes
    p2<-ggplot(data = top_genes,
               aes(x = logFC,y = -log10(adj.P.Val),
                   fill = DE,
                   size = DE,
                   alpha = DE,
                   label= DEG_symbol)) +
        sm_hgrid(legends = TRUE) +
        geom_point(shape =21) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) +
        geom_hline(yintercept = -log10(FDR),
                   linetype = "dashed") +
        geom_vline(xintercept = c(-logFC_abs,logFC_abs),
                   linetype = "dashed") +
        geom_label_repel(fill="white", size=2, max.overlaps = Inf,
                         box.padding = 0.2,
                         show.legend=FALSE) +
        labs(y="-log10(FDR)")+
        scale_fill_manual(values = cols) +
        scale_size_manual(values = sizes) +
        scale_alpha_manual(values = alphas)

    plot_grid(p1, p2, ncol=2)
    ggsave(paste("plots/05_DEA/01_Modeling/DEG_plots_", name, ".pdf", sep=""),
           width = 35, height = 15, units = "cm")
}


## Plots for habenula DEGs from the model without correlated variables
plots_DEGs(top_genes = results_uncorr_vars_habenula[[1]], vGene = results_uncorr_vars_habenula[[2]], FDR = 0.05,
           name='habenula_uncorr_variables')

## Plots for amygdala DEGs from the model without correlated variables
plots_DEGs(top_genes = results_uncorr_vars_amygdala[[1]], vGene = results_uncorr_vars_amygdala[[2]], FDR = 0.05,
           name='amygdala_uncorr_variables')






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

