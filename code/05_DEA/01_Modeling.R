
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
results_Substance_all_vars_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_Substance_all_vars_habenula, file = 'processed-data/05_DEA/results_Substance_all_vars_habenula.Rdata')
## Number of DEGs (FDR<0.05)
length(which(results_Substance_all_vars_habenula[[1]]$adj.P.Val<0.05))
#  0

## Model with uncorrelated variables only
## formula <-  ~ Substance + Batch_RNA_extraction + overallMapRate + RIN + mitoRate (previous)
formula <-  ~ Substance + Batch_RNA_extraction + concordMapRate + RIN
name<-"Substance_uncorr_variables"
coef<-"SubstanceFentanyl"
results_Substance_uncorr_vars_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_Substance_uncorr_vars_habenula, file = 'processed-data/05_DEA/results_Substance_uncorr_vars_habenula.Rdata')
length(which(results_Substance_uncorr_vars_habenula[[1]]$adj.P.Val<0.05))
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
results_Substance_all_vars_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_Substance_all_vars_amygdala, file = 'processed-data/05_DEA/results_Substance_all_vars_amygdala.Rdata')
## Number of DEGs (FDR<0.10)
length(which(results_Substance_all_vars_amygdala[[1]]$adj.P.Val<0.10))
#  0

## Model with uncorrelated variables only
## formula<- ~ Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate  + RIN
name<-"Substance_uncorr_variables"
coef<-"SubstanceFentanyl"
results_Substance_uncorr_vars_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_Substance_uncorr_vars_amygdala, file = 'processed-data/05_DEA/results_Substance_uncorr_vars_amygdala.Rdata')
length(which(results_Substance_uncorr_vars_amygdala[[1]]$adj.P.Val<0.05))
#  3041


#-----------------
## Plots for DEGs
top_genes_habenula <- results_Substance_uncorr_vars_habenula[[1]]
top_genes_habenula$symbol_or_ensemblID <- unlist(apply(top_genes_habenula, 1, function(x){if(is.na(x['Symbol'])){x['ensemblID']} else{x['Symbol']}}))
de_genes_habenula <- top_genes_habenula[which(top_genes_habenula$adj.P.Val<0.05),]

top_genes_amygdala <- results_Substance_uncorr_vars_amygdala[[1]]
top_genes_amygdala$symbol_or_ensemblID <- unlist(apply(top_genes_amygdala, 1, function(x){if(is.na(x['Symbol'])){x['ensemblID']} else{x['Symbol']}}))
de_genes_amygdala <- top_genes_amygdala[which(top_genes_amygdala$adj.P.Val<0.05),]

## Common DEGs
common_DEGs <- intersect(de_genes_habenula$ensemblID, de_genes_amygdala$ensemblID)


plots_DEGs<-function(brain_region, top_genes, vGene, name, DEGs_list) {

    if(name=='First_Hour_Infusion_Slope'){
        FClab = 'log2FC(1st Hour Infusion Slope)'
    }
    else if (name=='Total_Intake'){
        FClab = 'log2FC(Total Drug Intake)'
    }
    else if(name=='Last_Session_Intake'){
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


    ## Gene symbols (or ensemblID if missing) for specific DEGs
    top_genes$symbol_or_ensemblID <- unlist(apply(top_genes, 1, function(x){if(is.na(x['Symbol'])){x['ensemblID']} else{x['Symbol']}}))
    up_genes <- subset(top_genes, logFC>0)$symbol_or_ensemblID
    down_genes <- subset(top_genes, logFC<0)$symbol_or_ensemblID

    ## DEGs
    de_genes <- subset(top_genes, adj.P.Val < 0.05)
    up_de_genes <- subset(de_genes, logFC>0)
    top_up_de_genes <- up_de_genes[order(up_de_genes$adj.P.Val), ][1:5, "symbol_or_ensemblID"]
    down_de_genes <- subset(de_genes, logFC<0)
    top_down_de_genes <- down_de_genes[order(down_de_genes$adj.P.Val), ][1:5, "symbol_or_ensemblID"]
    genes_to_label <- c(top_up_de_genes, top_down_de_genes)

    ## Position of caption in plot
    caption_x_units <- 0.55
    caption_y_units1 <- 0.15
    caption_y_units2 <- 0.08

    if (brain_region=='amygdala' & !name=='Total_Intake'){
        caption_x_units <- 0.55
        caption_y_units1 <- 0.1
        caption_y_units2 <- -0.09
    }

    ## Label specific DEGs in plot
    # top_genes$DEG_symbol<- sapply(top_genes$symbol_or_ensemblID, function(x){ if(x %in% DEGs_list){x} else {NA}})

    ## Plots
    cols <- c("Up" = "indianred2", "Down" = "steelblue2", "n.s." = "grey")
    sizes <- c("Up" = 1.3, "Down" = 1.3, "n.s." = 0.8)
    alphas <- c("Up" = 0.4, "Down" = 0.6, "n.s." = 0.5)

    ## Max pval of DEGs
    maxP <- max(de_genes$P.Value)

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
              legend.position.inside = c(0.82, 0.15),
              legend.background = element_rect(fill=NA),
              legend.key.height = unit(0.15,"cm"),
              axis.title = element_text(size = (10)),
              legend.text = element_text(size=10))


    ## Volcano plot for DE genes
    p2<-ggplot(data = top_genes,
               aes(x = logFC,y = -log10(P.Value),
                   color = DE,
                   fill = DE,
                   size = DE,
                   alpha = DE)) +
        sm_hgrid(legends = TRUE) +
        geom_point(shape = 21) +
        geom_hline(yintercept = -log10(maxP),
                   linetype = "dashed", color = 'gray35', linewidth=0.5) +
        geom_vline(xintercept = c(-1,1),
                   linetype = "dashed", color = 'gray35', linewidth=0.5) +
        ## Label top DEGs
        geom_text_repel(data = subset(top_genes, symbol_or_ensemblID %in% genes_to_label),
                        aes(label = symbol_or_ensemblID),
                        size=2.7,
                        color='black',
                        alpha = 1,
                        max.overlaps = Inf,
                        box.padding = 0.15,
                        segment.size = unit(0.35, 'mm'),
                        segment.alpha = 0.4,
                        show.legend=FALSE) +
        ## Label down DEGs (not common)
        # geom_text_repel(data = subset(top_genes, symbol_or_ensemblID %in% down_genes &
        #                                   ! symbol_or_ensemblID %in% common_DEGs),
        #                 aes(fontface = 'bold'),
        #                 size=2.3,
        #                 color='black',
        #                 alpha = 1,
        #                 max.overlaps = Inf,
        #                 box.padding = 0.15, nudge_y = -0.1, nudge_x = -0.4,
        #                 segment.size = unit(0.35, 'mm'),
        #                 segment.alpha = 0.4,
        #                 show.legend=FALSE) +
        ## Label up DEGs (not common)
        # geom_text_repel(data = subset(top_genes, symbol_or_ensemblID %in% up_genes &
        #                                   ! symbol_or_ensemblID %in% common_DEGs),
        #                 aes(fontface = 'bold'),
        #                 size=2.3,
        #                 color='black',
        #                 alpha = 1,
        #                 max.overlaps = Inf,
        #                 box.padding = 0.15, nudge_y = 0.1, nudge_x = 0.4,
        #                 segment.size = unit(0.35, 'mm'),
        #                 segment.alpha = 0.4,
        #                 show.legend=FALSE) +
        ## Label common DEGs (down)
        # geom_text_repel(data = subset(top_genes, symbol_or_ensemblID %in% common_DEGs &
        #                                          symbol_or_ensemblID %in% DEGs_list &
        #                                          symbol_or_ensemblID %in% down_genes),
        #                  aes(fontface = 'bold'),
        #                  size=2.3,
        #                  color='darkorange3',
        #                  alpha = 1,
        #                  max.overlaps = Inf,
        #                  box.padding = 0.15, nudge_y = -0.3, nudge_x = -0.7,
        #                  segment.size = unit(0.4, 'mm'),
        #                  segment.alpha = 0.48,
        #                  show.legend=FALSE)+
        ## Label common DEGs (up)
        # geom_text_repel(data = subset(top_genes, symbol_or_ensemblID %in% common_DEGs &
        #                                   symbol_or_ensemblID %in% DEGs_list &
        #                                   symbol_or_ensemblID %in% up_genes),
        #                 aes(fontface = 'bold'),
        #                 size=2.3,
        #                 color='darkorange3',
        #                 alpha = 1,
        #                 max.overlaps = Inf,
        #                 box.padding = 0.15, nudge_y = 0.3, nudge_x = 0.7,
        #                 segment.size = unit(0.4, 'mm'),
        #                 segment.alpha = 0.48,
        #                 show.legend=FALSE)+

        labs(y="-log10(P)", x=FClab)+
        scale_color_manual(values = cols, name=NULL) +
        scale_fill_manual(values = cols, name=NULL) +
        scale_size_manual(values = sizes, name=NULL) +
        scale_alpha_manual(values = alphas, name=NULL) +
        theme(plot.margin = unit(c(1,1,1,1), "cm"),
              legend.position = c(0.13, 0.15),
              legend.background = element_rect(fill=NA),
              legend.key.height = unit(0.15,"cm"),
              axis.title = element_text(size = (10)),
              legend.text = element_text(size=10))
        ## Caption: number of DEGs
        # annotate("text", x=max(top_genes$logFC)-caption_x_units, y=caption_y_units1, label= paste0(length(which(top_genes$adj.P.Val<FDR)), ' DEGs'),
         #        color='gray40', size=3.1, fontface = 'bold')
        ## Caption: FDR threshold
        # annotate("text", x=max(top_genes$logFC)-caption_x_units, y=caption_y_units2, label= paste0("(FDR<", FDR, ")"),
        #          color='gray40', size=2.5)

    plot_grid(p1, p2, ncol=2)
    ggsave(paste("plots/05_DEA/01_Modeling/DEG_plots_", brain_region, '_', name, ".pdf", sep=""),
           width = 24, height = 12, units = "cm")
}


## Plots for habenula DEGs from the model without correlated variables
plots_DEGs('habenula', top_genes = results_Substance_uncorr_vars_habenula[[1]],
           vGene = results_Substance_uncorr_vars_habenula[[2]],
           name='Substance', DEGs_list = hab_DEGs)

## Plots for amygdala DEGs from the model without correlated variables
plots_DEGs('amygdala', top_genes = results_Substance_uncorr_vars_amygdala[[1]],
           vGene = results_Substance_uncorr_vars_amygdala[[2]],
           name='Substance', DEGs_list = amy_DEGs)


## Regress out covariates from gene expr
regress_covs <- function(ensemblID, variable, brain_region){

    if(variable == "Substance"){
        if(brain_region == "habenula"){
            RSE <- rse_gene_habenula_filt
            sample_data <- matrix(as.numeric(as.matrix(colData(RSE)[, c("Batch_RNA_extraction", "concordMapRate", "RIN")])), ncol = 3)
            vGene <- results_Substance_uncorr_vars_habenula[[2]]
            eGene <- results_Substance_uncorr_vars_habenula[[3]]
            coeffs <- eGene$coefficients[ensemblID,  c("Batch_RNA_extraction3", "concordMapRate", "RIN")]
        }
        else{
            RSE <- rse_gene_amygdala_filt
            sample_data <- matrix(as.numeric(as.matrix(colData(RSE)[, c("Batch_RNA_extraction", "Batch_lib_prep",
                                                                         "overallMapRate", "RIN")])), ncol = 4)
            vGene <- results_Substance_uncorr_vars_amygdala[[2]]
            eGene <- results_Substance_uncorr_vars_amygdala[[3]]
            coeffs <- eGene$coefficients[ensemblID,  c("Batch_RNA_extraction3", "Batch_lib_prep3", "overallMapRate", "RIN")]

        }
    }

    y = vGene$E[ensemblID, ]
    return(y - (sample_data %*% coeffs))
}

## Boxplot/Scatterplot for a DEG of interest
DEG_expression_plot <- function (variable, brain_region, gene){

    if (variable=='Substance'){
        rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
    }

    else {
        rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
        rse_gene <- rse_gene[which(rse_gene$Substance=='Fentanyl'), ]
    }

    top_genes <- eval(parse_expr(paste0("top_genes_", brain_region)))
    vGene <- eval(parse_expr(paste0("results_", variable, "_uncorr_vars_", brain_region, "[[2]]")))

    ## Sample colors
    colors = list("Substance"= c('Fentanyl'="#851C84", 'Saline'="#505050"),
                  "Total_Intake" = "yellow2",
                  "Last_Session_Intake" = "pink",
                  "First_Hour_Infusion_Slope" = "lightblue2")
    x_labs = c("Substance"= "Substance",
               "Total_Intake" = "Total intake",
               "Last_Session_Intake" = "Last session intake",
               "First_Hour_Infusion_Slope" = "First hour infusion slope")

    ## q-value for the gene
    q_value <-signif(top_genes[which(top_genes$symbol_or_ensemblID==gene)[1], "adj.P.Val"], digits = 2)

    ## FC
    FC <-signif(2**(top_genes[which(top_genes$symbol_or_ensemblID==gene)[1], "logFC"]), digits=2)

    # Gene symbol + ensemblID
    ensemblID <- top_genes[which(top_genes$symbol_or_ensemblID==gene)[1], "ensemblID"]
    if (length(ensemblID)!=0 && gene!=ensemblID){
        gene_title <- paste(gene, ensemblID, sep="-")
    }else {
        gene_title <- gene
    }

    ## Merge lognorm counts of DEG with sample data
    lognorm_DE <- regress_covs(ensemblID, variable, brain_region)
    lognorm_DE <- as.vector(t(lognorm_DE))
    lognorm_DE <- data.frame(deg = lognorm_DE, variable=colData(rse_gene)[variable])
    lognorm_DE$Substance <- factor(lognorm_DE$Substance, levels = c("Saline", "Fentanyl"))

    ## Boxplot for Substance DGE
    if (variable == 'Substance'){
        plot <- ggplot(data = lognorm_DE,
                   aes(x=eval(parse_expr(variable)),
                       y=deg)) +
                    ## Hide outliers
                    geom_boxplot(outlier.color = "#FFFFFFFF", width=0.3, color = "gray30", alpha = 0.6) +
                    ## Samples colored by variable of interest
                    geom_jitter(aes(fill=eval(parse_expr(variable))), color = "black", alpha = 0.85,
                                shape=21,
                                position=position_jitter(0.1),
                                size=2.2) +
                    theme_classic() +
                    scale_fill_manual(values=colors[[variable]]) +
                    labs(x = x_labs[variable], y = "log(cpm) - covariates",
                         title = gene,
                         subtitle = paste("FDR-adjusted p:", q_value, '    ', 'FC:', FC)) +
                    theme(legend.position = "none",
                          plot.title = element_text(hjust=0.5, size=10, face="bold"),
                          plot.subtitle = element_text(size = 9),
                          axis.title = element_text(size = (10)),
                          axis.text = element_text(size = 8.5))
    }

    return(plot)
}

## Plot expression of multiple DEGs of interest
DEG_expression_plots <- function(DEGs, variable, brain_region, name, w = 15, h = 11){

    plots<-list()
    for (i in 1:length(DEGs)){
        p<-DEG_expression_plot(variable, brain_region, DEGs[i])
        plots[[i]] <- p
    }
    plot_grid(plotlist = plots, ncol = 3, align = 'hv')
    ggsave(here(paste("plots/05_DEA/01_Modeling/", name, "DEGs_expression_", brain_region, "_", variable, ".pdf", sep="")),
           width = w, height = h, units = "cm")
}


## Top 5 most downregulated DEGs for Substance in Hb
down_DEGs <- de_genes_habenula[which(de_genes_habenula$logFC<0), ]
down_DEGs <- down_DEGs[order(down_DEGs$adj.P.Val, decreasing = FALSE), 'symbol_or_ensemblID'][1:5]
DEG_expression_plots(down_DEGs, 'Substance', 'habenula', 'down_', w = 21, h = 14)
## Top 5 most upregulated DEGs for Substance in Hb
up_DEGs <- de_genes_habenula[which(de_genes_habenula$logFC>0), ]
up_DEGs <- up_DEGs[order(up_DEGs$adj.P.Val, decreasing = FALSE), 'symbol_or_ensemblID'][1:5]
DEG_expression_plots(up_DEGs, 'Substance', 'habenula', 'up_', w = 21, h = 14)

## Top 5 most downregulated DEGs for Substance in Amyg
down_DEGs <- de_genes_amygdala[which(de_genes_amygdala$logFC<0), ]
down_DEGs <- down_DEGs[order(down_DEGs$adj.P.Val, decreasing = FALSE), 'symbol_or_ensemblID'][1:5]
DEG_expression_plots(down_DEGs, 'Substance', 'amygdala', 'down_', w = 21, h = 14)
## Top 5 most upregulated DEGs for Substance in Amyg
up_DEGs <- de_genes_amygdala[which(de_genes_amygdala$logFC>0), ]
up_DEGs <- up_DEGs[order(up_DEGs$adj.P.Val, decreasing = FALSE), 'symbol_or_ensemblID'][1:5]
DEG_expression_plots(up_DEGs, 'Substance', 'amygdala', 'up_', w = 21, h = 14)

## Unique/shared up/down DEGs in Hb and Amyg (relevant from GO & KEGG terms)

unique_up_Hb <- c("Grik1", "Scn1a", "Gria4", "Grm1", "Rgs7", "Cttnbp2", "Zfp804a", "Cacna1i",
                  "Kcnj10", "Kcnq3", "Adora1", "Gabra4", "Adcy8", "Chrm3")
unique_down_Hb <- c("Sulf1", "Pbxip1", "Spint2", "Lama3", "Col8a2", "Sod3", "Ezr", "Slc4a2",
                    "Calml4", "Itgb6", "Itga2")
unique_up_Amyg <- c("Snca", "Itpka", "Lrfn5", "Gria2", "Ppp3cb", "Adcy1", "Cox6b1", "Uqcrq", "Calm1",
                    "Arpc3", "Mpv17l2", "Mal2", "Slc2a3", "Mrpl45", "Mrps18a", "Mrpl21", "Atp5mg", "Atp5me",
                    "Atp5mf", "Atp6v1e1", "Atp5mc2", "Atp6v0e2", "Slc39a10", "Nckap1", "Cfl1")
unique_down_Amyg <- c("Phldb1", "Col15a1", "Col14a1", "Abca2", "Sox13", "Notch1", "Lama4", "Myo1d", "Scrib", "Llgl1",
                      "Tnc", "Lama5", "Arsg", "Abcc1", "Stard9", "Kif13a", "Myo9b", "Col5a3", "Vcl", "C6", "Wasf2")
shared_up_Hb_up_Amyg <- c("Kcnj3", "Kcnc2", "Kcnj9", "Kctd16", "Cdh10", "Npy1r", "Epha4",
                          "Lrrtm2", "LRRTM1")
shared_up_Hb_down_Amyg <- c("Plcb4", "Col4a3")
shared_down_Hb_down_Amyg <- c("Tgfbi", "Antxr1", "Col9a3", "Loxl4", "Sox9", "Cnp", "Sox8",
                              "Gsn", "Opalin", "Tspan2", "Fgfr2", "Tubb4a")

## Unique up DEGs in Hb
DEG_expression_plots(unique_up_Hb, 'Substance', 'habenula', 'unique_up_Hb_', w = 21, h = 35)
DEG_expression_plots(unique_up_Hb, 'Substance', 'amygdala', 'unique_up_Hb_', w = 21, h = 35)
## Unique down DEGs in Hb
DEG_expression_plots(unique_down_Hb, 'Substance', 'habenula', 'unique_down_Hb_', w = 21, h = 28)
DEG_expression_plots(unique_down_Hb, 'Substance', 'amygdala', 'unique_down_Hb_', w = 21, h = 28)
## Unique up DEGs in Amyg
DEG_expression_plots(unique_up_Amyg, 'Substance', 'habenula', 'unique_up_Amyg_', w = 21, h = 63)
DEG_expression_plots(unique_up_Amyg, 'Substance', 'amygdala', 'unique_up_Amyg_', w = 21, h = 63)
## Unique down DEGs in Amyg
DEG_expression_plots(unique_down_Amyg, 'Substance', 'habenula', 'unique_down_Amyg_', w = 21, h = 49)
DEG_expression_plots(unique_down_Amyg, 'Substance', 'amygdala', 'unique_down_Amyg_', w = 21, h = 49)
## Shared up DEGs in Hb and Amyg
DEG_expression_plots(shared_up_Hb_up_Amyg, 'Substance', 'habenula', 'shared_up_Hb_up_Amyg_', w = 21, h = 21)
DEG_expression_plots(shared_up_Hb_up_Amyg, 'Substance', 'amygdala', 'shared_up_Hb_up_Amyg_', w = 21, h = 21)
## Shared DEGs: up in Hb and down in Amyg
DEG_expression_plots(shared_up_Hb_down_Amyg, 'Substance', 'habenula', 'shared_up_Hb_down_Amyg_', w = 21, h = 7)
DEG_expression_plots(shared_up_Hb_down_Amyg, 'Substance', 'amygdala', 'shared_up_Hb_down_Amyg_', w = 21, h = 7)
## Shared down DEGs in Hb and Amyg
DEG_expression_plots(shared_down_Hb_down_Amyg, 'Substance', 'habenula', 'shared_down_Hb_down_Amyg_', w = 21, h = 28)
DEG_expression_plots(shared_down_Hb_down_Amyg, 'Substance', 'amygdala', 'shared_down_Hb_down_Amyg_', w = 21, h = 28)
#-----------------


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
## Create supp table with results
de_genes_habenula$chr <- de_genes_habenula$seqnames
de_genes_habenula <- de_genes_habenula[,c("chr", "start", "end", "width", "strand", "Length", "ensemblID",
                                          "EntrezID", "Symbol", "meanExprs", "logFC", "t", "P.Value", "adj.P.Val",
                                          "gene_description", "associated_phenotypes")]
save(de_genes_habenula, file = 'processed-data/05_DEA/de_genes_Substance_habenula.Rdata')
write.table(de_genes_habenula, "processed-data/Supplementary_Tables/TableS5_de_genes_Substance_habenula.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')

## Amygdala
de_genes_amygdala <- add_phenotypes(de_genes_amygdala)
de_genes_amygdala <- add_description(de_genes_amygdala)
de_genes_amygdala$EntrezID <- as.character(de_genes_amygdala$EntrezID)
de_genes_amygdala <- de_genes_amygdala[order(de_genes_amygdala$adj.P.Val, -abs(de_genes_amygdala$logFC)),]
## Create supp table with results
de_genes_amygdala$chr <- de_genes_amygdala$seqnames
de_genes_amygdala <- de_genes_amygdala[,c("chr", "start", "end", "width", "strand", "Length", "ensemblID",
                                          "EntrezID", "Symbol", "meanExprs", "logFC", "t", "P.Value", "adj.P.Val",
                                          "gene_description", "associated_phenotypes")]
save(de_genes_amygdala, file = 'processed-data/05_DEA/de_genes_Substance_amygdala.Rdata')
write.table(de_genes_amygdala, "processed-data/Supplementary_Tables/TableS6_de_genes_Substance_amygdala.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')


## Supp table for common DEGs
de_genes_common <- de_genes_habenula[c(which(de_genes_habenula$Symbol %in% common_DEGs),
                                       which(de_genes_habenula$ensemblID %in% common_DEGs)), ]
colnames(de_genes_common)[11:14] <- paste0(colnames(de_genes_common)[11:14], "_habenula")
## Bind DGE metrics in amygdala
de_genes_common_in_amy <- de_genes_amygdala[c(which(de_genes_amygdala$Symbol %in% common_DEGs),
                                       which(de_genes_amygdala$ensemblID %in% common_DEGs)), ]
de_genes_common <- merge(de_genes_common, de_genes_common_in_amy[, c("ensemblID", "logFC", "t", "P.Value", "adj.P.Val")], by = "ensemblID")
colnames(de_genes_common)[17:20] <- paste0(colnames(de_genes_common)[17:20], "_amygdala")
de_genes_common <- de_genes_common[,c("chr", "start", "end", "width", "strand", "Length", "ensemblID",
                                      "EntrezID", "Symbol", "meanExprs",
                                      paste0(c("logFC", "t", "P.Value", "adj.P.Val"), "_habenula"),
                                      paste0(c("logFC", "t", "P.Value", "adj.P.Val"), "_amygdala"),
                                      "gene_description", "associated_phenotypes")]
save(de_genes_common, file = 'processed-data/05_DEA/de_genes_common_Substance_hab_amy.Rdata')
write.table(de_genes_common, "processed-data/Supplementary_Tables/TableS7_de_genes_common_Substance_hab_amy.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')





## --------------------------------------------------------------------------------------------------
##      1.2 DEA for High vs Low fentanyl intake slope in habenula and amygdala fentanyl samples
## --------------------------------------------------------------------------------------------------
## Note: slopes used to categorize rats in this analysis were computed with an old method based on 3 session averages per rat

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

## Same variables for habenula Substance DEA
## (Substance and batches are the same for all these samples)
## formula<- ~ Intake_slope_binary + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Intake_slope_binary + concordMapRate + RIN
name<-"Intake_slope_binary"
coef<-"Intake_slope_binaryLow"
results_Intake_Slope_binary_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
## DEGs (FDR<0.05)
de_genes_Intake_Slope_binary_habenula <- results_Intake_Slope_binary_habenula[[1]][which(results_Intake_Slope_binary_habenula[[1]]$adj.P.Val<0.05), ]
save(results_Intake_Slope_binary_habenula, file = 'processed-data/05_DEA/results_Intake_Slope_binary_habenula.Rdata')
## No DEGs
dim(de_genes_Intake_Slope_binary_habenula)[1]
##  0


##############################
## Amygdala fentanyl samples
##############################

rse_gene_amygdala_fent <- rse_gene_amygdala_filt[,which(rse_gene_amygdala_filt$Substance=='Fentanyl')]

## High/low intake samples
rse_gene_amygdala_fent$Intake_slope_binary <- sapply(rse_gene_amygdala_fent$Sample_Num, function(x){if(x %in% c(17,20,23)){'High'}
                                                                                             else if(x %in% c(18,21,22,24)){'Low'}
                                                                                             else {NA}})
## Remove sample from the outlier rat (in intake slope)
rse_gene_amygdala_fent <- rse_gene_amygdala_fent[,-which(is.na(rse_gene_amygdala_fent$Intake_slope_binary))]

## Same variables for amygdala DEA
## formula<- ~ Intake_slope_binary + overallMapRate + RIN + mitoRate  (previous)
formula <-  ~  Intake_slope_binary + overallMapRate  + RIN
name<-"Intake_slope_binary"
coef<-"Intake_slope_binaryLow"
results_Intake_Slope_binary_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
de_genes_Intake_Slope_binary_amygdala <- results_Intake_Slope_binary_amygdala[[1]][which(results_Intake_Slope_binary_amygdala[[1]]$adj.P.Val<0.05), ]
save(results_Intake_Slope_binary_amygdala, file = 'processed-data/05_DEA/results_Intake_Slope_binary_amygdala.Rdata')
## No DEGs
dim(de_genes_Intake_Slope_binary_amygdala)[1]
##  0





## ----------------------------------------------------------------------------------------------
##      1.3  DEA for Fentanyl vs. Saline with behavioral covariates in habenula and amygdala
##                (1st hour intake/infusion slope, Total intake, and Last session intake)
## ----------------------------------------------------------------------------------------------

############### 1.3.1 DEA with covariate '1st Hour Intake Slope' ###############

#####################
## Habenula samples
#####################

## Same uncorrelated variables as before
## formula <-  ~ Substance + First_Hour_Infusion_Slope + Batch_RNA_extraction + overallMapRate + RIN + mitoRate (previous)
formula <-  ~ Substance + Batch_RNA_extraction + concordMapRate + RIN + First_Hour_Infusion_Slope
name <-"Substance_with_First_Hour_Infusion_Slope"
coef <-"SubstanceFentanyl"
results_Substance_with_FirstHrIntakeSlope_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_Substance_with_FirstHrIntakeSlope_habenula, file = 'processed-data/05_DEA/results_Substance_with_FirstHrIntakeSlope_habenula.Rdata')
length(which(results_Substance_with_FirstHrIntakeSlope_habenula[[1]]$adj.P.Val<0.05))
#  0

#####################
## Amygdala samples
#####################

## formula <-  ~ Substance + First_Hour_Infusion_Slope + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate  + RIN + First_Hour_Infusion_Slope
name <-"Substance_with_First_Hour_Infusion_Slope"
coef <-"SubstanceFentanyl"
results_Substance_with_FirstHrIntakeSlope_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_Substance_with_FirstHrIntakeSlope_amygdala, file = 'processed-data/05_DEA/results_Substance_with_FirstHrIntakeSlope_amygdala.Rdata')
length(which(results_Substance_with_FirstHrIntakeSlope_amygdala[[1]]$adj.P.Val<0.05))
#  0



##################   1.3.2 DEA with covariate 'Total Intake'  ##################

#####################
## Habenula samples
#####################

## Same uncorrelated variables as before
## formula <-  ~ Substance + Total_Intake + Batch_RNA_extraction + overallMapRate + RIN + mitoRate (previous)
formula <-  ~ Substance + Batch_RNA_extraction + concordMapRate + RIN + Total_Intake
name <-"Substance_with_Total_Intake"
coef <-"SubstanceFentanyl"
results_Substance_with_TotalIntake_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_Substance_with_TotalIntake_habenula, file = 'processed-data/05_DEA/results_Substance_with_TotalIntake_habenula.Rdata')
length(which(results_Substance_with_TotalIntake_habenula[[1]]$adj.P.Val<0.05))
#  0

#####################
## Amygdala samples
#####################

## formula <-  ~ Substance + Total_Intake + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate  + RIN + Total_Intake
name <-"Substance_with_Total_Intake"
coef <-"SubstanceFentanyl"
results_Substance_with_TotalIntake_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_Substance_with_TotalIntake_amygdala, file = 'processed-data/05_DEA/results_Substance_with_TotalIntake_amygdala.Rdata')
length(which(results_Substance_with_TotalIntake_amygdala[[1]]$adj.P.Val<0.05))
#  0



##############   1.3.3 DEA with covariate 'Last Session Intake'  ###############

#####################
## Habenula samples
#####################

## Same uncorrelated variables as before
## formula <-  ~ Substance + Last_Session_Intake + Batch_RNA_extraction + overallMapRate + RIN + mitoRate
formula <-  ~ Substance + Batch_RNA_extraction + concordMapRate + RIN + Last_Session_Intake
name <-"Substance_with_Last_Session_Intake"
coef <-"SubstanceFentanyl"
results_Substance_with_LastSessionIntake_habenula<-DEA(rse_gene_habenula_filt, 'habenula', formula, name, coef)
save(results_Substance_with_LastSessionIntake_habenula, file = 'processed-data/05_DEA/results_Substance_with_LastSessionIntake_habenula.Rdata')
length(which(results_Substance_with_LastSessionIntake_habenula[[1]]$adj.P.Val<0.05))
#  0

#####################
## Amygdala samples
#####################

## formula <-  ~ Substance + Last_Session_Intake + Batch_RNA_extraction + Batch_lib_prep + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Substance + Batch_RNA_extraction + Batch_lib_prep + overallMapRate  + RIN + Last_Session_Intake
name <-"Substance_with_Last_Session_Intake"
coef <-"SubstanceFentanyl"
results_Substance_with_LastSessionIntake_amygdala<-DEA(rse_gene_amygdala_filt, 'amygdala', formula, name, coef)
save(results_Substance_with_LastSessionIntake_amygdala, file = 'processed-data/05_DEA/results_Substance_with_LastSessionIntake_amygdala.Rdata')
length(which(results_Substance_with_LastSessionIntake_amygdala[[1]]$adj.P.Val<0.05))
#  0





## ---------------------------------------------------------------------------------
##   1.4  DEA for 1st Hour Intake Slope in habenula and amygdala fentanyl samples
## ---------------------------------------------------------------------------------
## Note: infusion slopes used for this analysis were calculated fitting a regression line

####################  1.4.1 Analysis with all fentanyl samples  ####################

#####################
## Habenula samples
#####################

## Fentanyl samples only
rse_gene_habenula_fent <- rse_gene_habenula_filt[,which(rse_gene_habenula_filt$Substance=='Fentanyl')]

## formula <- ~ First_Hour_Infusion_Slope + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  First_Hour_Infusion_Slope + RIN + RNA_concentration + mitoRate
name <-"for_First_Hour_Infusion_Slope"
## New contrast of interest
coef <-"First_Hour_Infusion_Slope"
results_FirstHrIntakeSlope_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
save(results_FirstHrIntakeSlope_habenula, file = 'processed-data/05_DEA/results_FirstHrIntakeSlope_habenula.Rdata')
## No DEGs (FDR<0.05)
length(which(results_FirstHrIntakeSlope_habenula[[1]]$adj.P.Val<0.05))
#  0

#####################
## Amygdala samples
#####################

rse_gene_amygdala_fent <- rse_gene_amygdala_filt[,which(rse_gene_amygdala_filt$Substance=='Fentanyl')]

## formula <-  ~ First_Hour_Infusion_Slope + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  First_Hour_Infusion_Slope + RIN + mitoRate
name <-"for_First_Hour_Infusion_Slope"
coef <-"First_Hour_Infusion_Slope"
results_FirstHrIntakeSlope_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
save(results_FirstHrIntakeSlope_amygdala, file = 'processed-data/05_DEA/results_FirstHrIntakeSlope_amygdala.Rdata')
length(which(results_FirstHrIntakeSlope_amygdala[[1]]$adj.P.Val<0.1))
#  0



#########  1.4.2 Analysis without samples from negative outlier fentanyl rat ########

#####################
## Habenula samples
#####################

## Remove outlier rat sample in habenula
rse_gene_habenula_fent_withoutOutlier <- rse_gene_habenula_fent[,-which(rse_gene_habenula_fent$Rat_ID==outlier_fent_sample)]

formula <-  ~  First_Hour_Infusion_Slope + RIN + RNA_concentration + mitoRate
name <-"for_First_Hour_Infusion_Slope_withoutOutlier"
coef <-"First_Hour_Infusion_Slope"
results_FirstHrIntakeSlope_habenula_withoutOutlier<-DEA(rse_gene_habenula_fent_withoutOutlier, 'habenula', formula, name, coef)
save(results_FirstHrIntakeSlope_habenula_withoutOutlier, file = 'processed-data/05_DEA/results_FirstHrIntakeSlope_habenula_withoutOutlier.Rdata')
length(which(results_FirstHrIntakeSlope_habenula_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0

#####################
## Amygdala samples
#####################

rse_gene_amygdala_fent_withoutOutlier <- rse_gene_amygdala_fent[,-which(rse_gene_amygdala_fent$Rat_ID==outlier_fent_sample)]

formula <- ~ First_Hour_Infusion_Slope + RIN + mitoRate
name <-"for_First_Hour_Infusion_Slope_withoutOutlier"
coef <-"First_Hour_Infusion_Slope"
results_FirstHrIntakeSlope_amygdala_withoutOutlier<-DEA(rse_gene_amygdala_fent_withoutOutlier, 'amygdala', formula, name, coef)
save(results_FirstHrIntakeSlope_amygdala_withoutOutlier, file = 'processed-data/05_DEA/results_FirstHrIntakeSlope_amygdala_withoutOutlier.Rdata')
length(which(results_FirstHrIntakeSlope_amygdala_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0





## -----------------------------------------------------------------------------
##     1.5  DEA for Total Intake in habenula and amygdala fentanyl samples
## -----------------------------------------------------------------------------

##################  1.5.1 Analysis with all fentanyl samples  ##################

#####################
## Habenula samples
#####################

##formula <-  ~ Total_Intake + overallMapRate + RIN + mitoRate (previous)
formula <- ~  Total_Intake + RIN + RNA_concentration + overallMapRate
name <-"for_Total_Intake"
coef <-"Total_Intake"
results_TotalIntake_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
save(results_TotalIntake_habenula, file = 'processed-data/05_DEA/results_TotalIntake_habenula.Rdata')
length(which(results_TotalIntake_habenula[[1]]$adj.P.Val<0.1))
#  0

#####################
## Amygdala samples
#####################

## formula <-  ~ Total_Intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Total_Intake + RIN + mitoRate
name <-"for_Total_Intake"
coef <-"Total_Intake"
results_TotalIntake_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
save(results_TotalIntake_amygdala, file = 'processed-data/05_DEA/results_TotalIntake_amygdala.Rdata')
length(which(results_TotalIntake_amygdala[[1]]$adj.P.Val<0.1))
#  0



########  1.5.2 Analysis without samples from negative outlier fentanyl rat #######

#####################
## Habenula samples
#####################

## formula <-  ~ Total_Intake + overallMapRate + RIN + mitoRate (previous)
formula <- ~  Total_Intake + RIN + RNA_concentration + overallMapRate
name <-"for_Total_Intake_withoutOutlier"
coef <-"Total_Intake"
results_TotalIntake_habenula_withoutOutlier<-DEA(rse_gene_habenula_fent_withoutOutlier, 'habenula', formula, name, coef)
save(results_TotalIntake_habenula_withoutOutlier, file = 'processed-data/05_DEA/results_TotalIntake_habenula_withoutOutlier.Rdata')
length(which(results_TotalIntake_habenula_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0

#####################
## Amygdala samples
#####################

## formula <-  ~ Total_intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Total_Intake + RIN + mitoRate
name <-"for_Total_Intake_withoutOutlier"
coef <-"Total_Intake"
results_TotalIntake_amygdala_withoutOutlier<-DEA(rse_gene_amygdala_fent_withoutOutlier, 'amygdala', formula, name, coef)
save(results_TotalIntake_amygdala_withoutOutlier, file = 'processed-data/05_DEA/results_TotalIntake_amygdala_withoutOutlier.Rdata')
length(which(results_TotalIntake_amygdala_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0





## -----------------------------------------------------------------------------
#   1.6  DEA for Last Session Intake in habenula and amygdala fentanyl samples
## -----------------------------------------------------------------------------

##################  1.6.1 Analysis with all fentanyl samples  ##################

#####################
## Habenula samples
#####################

# formula <-  ~ Last_Session_Intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Last_Session_Intake + RIN + RNA_concentration + mitoRate
name <-"for_Last_Session_Intake"
coef <-"Last_Session_Intake"
results_LastSessionIntake_habenula<-DEA(rse_gene_habenula_fent, 'habenula', formula, name, coef)
save(results_LastSessionIntake_habenula, file = 'processed-data/05_DEA/results_LastSessionIntake_habenula.Rdata')
length(which(results_LastSessionIntake_habenula[[1]]$adj.P.Val<0.1))
#  0

#####################
## Amygdala samples
#####################

# formula <-  ~ Last_Session_Intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Last_Session_Intake + RIN + totalAssignedGene + concordMapRate
name <-"for_Last_Session_Intake"
coef <-"Last_Session_Intake"
results_LastSessionIntake_amygdala<-DEA(rse_gene_amygdala_fent, 'amygdala', formula, name, coef)
save(results_LastSessionIntake_amygdala, file = 'processed-data/05_DEA/results_LastSessionIntake_amygdala.Rdata')
length(which(results_LastSessionIntake_amygdala[[1]]$adj.P.Val<0.1))
#  0



#########  1.6.2 Analysis without samples from negative outlier fentanyl rat ########

#####################
## Habenula samples
#####################

## formula <-  ~ Last_Session_Intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Last_Session_Intake + RIN + RNA_concentration + mitoRate
name <-"for_Last_Session_Intake_withoutOutlier"
coef <-"Last_Session_Intake"
results_LastSessionIntake_habenula_withoutOutlier<-DEA(rse_gene_habenula_fent_withoutOutlier, 'habenula', formula, name, coef)
save(results_LastSessionIntake_habenula_withoutOutlier, file = 'processed-data/05_DEA/results_LastSessionIntake_habenula_withoutOutlier.Rdata')
length(which(results_LastSessionIntake_habenula_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0

#####################
## Amygdala samples
#####################

## formula <-  ~ Last_Session_Intake + overallMapRate + RIN + mitoRate (previous)
formula <-  ~  Last_Session_Intake + RIN + totalAssignedGene + concordMapRate
name <-"for_Last_Session_Intake_withoutOutlier"
coef <-"Last_Session_Intake"
results_LastSessionIntake_amygdala_withoutOutlier<-DEA(rse_gene_amygdala_fent_withoutOutlier, 'amygdala', formula, name, coef)
save(results_LastSessionIntake_amygdala_withoutOutlier, file = 'processed-data/05_DEA/results_LastSessionIntake_amygdala_withoutOutlier.Rdata')
length(which(results_LastSessionIntake_amygdala_withoutOutlier[[1]]$adj.P.Val<0.1))
#  0




## Plot gene expression vs total/last/slope of drug intake

plot_gene_expr_vs_intake <- function(brain_region, sample_var, gene_id){

    rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'fent', sep="_")))

    if(sample_var=='First_Hour_Infusion_Slope'){
        x_label="First hr infusion slope"
    }

    else if(sample_var=='Total_Intake'){
        x_label="Total drug intake"
    }

    else if(sample_var=='Last_Session_Intake'){
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

## Top most affected gene from DEA for 1st hr infusion slope
dea_results <- results_FirstHrIntakeSlope_habenula[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p1 <- plot_gene_expr_vs_intake('habenula', 'First_Hour_Infusion_Slope', gene_id)

## Top most affected gene from DEA for total intake
dea_results <- results_TotalIntake_habenula[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p2 <- plot_gene_expr_vs_intake('habenula', 'Total_Intake', gene_id)

## Top most affected gene from DEA for last intake
dea_results <- results_LastSessionIntake_habenula[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p3 <- plot_gene_expr_vs_intake('habenula', 'Last_Session_Intake', gene_id)

plot_grid(p1, p2, p3, nrow=1)
ggsave('plots/05_DEA/01_Modeling/geneExpr_VS_drugIntake_habenula_withOutlier.pdf', width = 20, height = 7, units = "cm")



##############################
## Amygdala fentanyl samples
##############################

## Top most affected gene from DEA for 1st hr infusion slope
dea_results <- results_FirstHrIntakeSlope_amygdala[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p1 <- plot_gene_expr_vs_intake('amygdala', 'First_Hour_Infusion_Slope', gene_id)

## Top most affected gene from DEA for total intake
dea_results <- results_TotalIntake_amygdala[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p2 <- plot_gene_expr_vs_intake('amygdala', 'Total_Intake', gene_id)

## Top most affected gene from DEA for last intake
dea_results <- results_LastSessionIntake_amygdala[[1]]
gene_id <- rownames(dea_results[order(dea_results$adj.P.Val, decreasing = FALSE), ][1,])
p3 <- plot_gene_expr_vs_intake('amygdala', 'Last_Session_Intake', gene_id)

plot_grid(p1, p2, p3, nrow=1)
ggsave('plots/05_DEA/01_Modeling/geneExpr_VS_drugIntake_amygdala_withOutlier.pdf', width = 20, height = 7, units = "cm")







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
# date     2024-04-22
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# backports              1.4.1     2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3     2015-07-28 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0    2023-10-26 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# broom                  1.0.5     2023-06-09 [1] CRAN (R 4.3.0)
# car                    3.1-2     2023-03-30 [1] CRAN (R 4.3.0)
# carData                3.0-5     2022-01-06 [1] CRAN (R 4.3.0)
# cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
# checkmate              2.3.1     2023-12-04 [1] CRAN (R 4.3.1)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# curl                   5.2.1     2024-03-01 [1] CRAN (R 4.3.1)
# data.table             1.15.2    2024-02-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0    2023-11-06 [1] Bioconductor
# digest                 0.6.34    2024-01-11 [1] CRAN (R 4.3.1)
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                * 4.0.16    2024-02-20 [1] Bioconductor 3.18 (R 4.3.2)
# evaluate               0.23      2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# foreign                0.8-86    2023-11-28 [1] CRAN (R 4.3.1)
# Formula                1.2-5     2023-02-24 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.6    2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2024-02-17 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-30 [1] Bioconductor
# gghalves               0.1.4     2022-11-20 [1] CRAN (R 4.3.0)
# ggplot2              * 3.5.0     2024-02-23 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0     2023-02-10 [1] CRAN (R 4.3.0)
# ggrepel              * 0.9.5     2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4     2022-10-13 [1] CRAN (R 4.3.0)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                  5.1-1     2023-09-12 [1] CRAN (R 4.3.0)
# htmlTable              2.4.2     2023-10-29 [1] CRAN (R 4.3.1)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.1)
# htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.3.1)
# httr                 * 1.4.7     2023-08-15 [1] CRAN (R 4.3.0)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# jsonlite             * 1.8.8     2023-12-04 [1] CRAN (R 4.3.1)
# knitr                  1.45      2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1    2023-11-02 [1] Bioconductor
# locfit                 1.5-9.9   2024-03-01 [1] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-26 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# mgcv                   1.9-1     2023-12-21 [1] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nlme                   3.1-164   2023-11-27 [1] CRAN (R 4.3.1)
# nnet                   7.3-19    2023-05-03 [1] CRAN (R 4.3.2)
# pheatmap             * 1.0.12    2019-01-04 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# pwr                    1.3-0     2020-03-17 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26      2024-03-05 [1] CRAN (R 4.3.1)
# rpart                  4.1.23    2023-12-05 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
# rstatix                0.7.2     2023-02-01 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-26 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sdamr                  0.2.0     2022-11-16 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# smplot2              * 0.1.0     2024-03-13 [1] Github (smin95/smplot2@052f4f9)
# SparseArray            1.2.4     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr                1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0    2023-11-06 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.42      2024-02-08 [1] CRAN (R 4.3.1)
# xml2                 * 1.3.6     2023-12-04 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-26 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-26 [1] Bioconductor
# zoo                    1.8-12    2023-04-13 [1] CRAN (R 4.3.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
