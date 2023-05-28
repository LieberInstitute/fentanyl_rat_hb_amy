library(here)
library(SummarizedExperiment)
library(ggplot2)
library(scater)
library(rlang)
library(smplot2)
library(cowplot)
library(Hmisc)
library(lme4)
library(variancePartition)
library(reshape2)
library(sessioninfo)



################################################################################
##                        3. Explore gene-level effects
################################################################################

load(here('processed-data/04_EDA/02_PCA/rse_gene_habenula_filt.Rdata'), verbose = TRUE)
load(here('processed-data/04_EDA/02_PCA/rse_gene_amygdala_filt.Rdata'), verbose = TRUE)


#######################   3.1 Explanatory variables   #######################

## Compute the % of gene expression variance explained by each sample variable

## Plot density function for % of variance explained
expl_var<- function(brain_region){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## % of variance in gene expression explained by each variable

    if (brain_region=='habenula'){
        variables=c("Substance", "Batch_RNA_extraction", "Total_Num_Fentanyl_Sessions", "mitoRate", "concordMapRate",
                    "overallMapRate", "totalAssignedGene", "RIN", "detected_num_genes", "library_size", "Total_RNA_amount",
                    "RNA_concentration")
    }
    else{
        variables=c("Substance", "Batch_RNA_extraction", "Batch_lib_prep", "Total_Num_Fentanyl_Sessions",
                    "mitoRate", "concordMapRate","overallMapRate", "totalAssignedGene", "RIN", "detected_num_genes",
                    "library_size", "Total_RNA_amount", "RNA_concentration")
    }

    exp_vars<-getVarianceExplained(RSE, variables=variables, exprs_values = "logcounts")

    ## Plot of explanatory variables
    ## Variables' colors
    colors=c("Substance"= 'turquoise4', "Batch_RNA_extraction"='bisque2', "Batch_lib_prep"='blueviolet',
             "Total_Num_Fentanyl_Sessions"='indianred1', 'mitoRate'='khaki3', 'totalAssignedGene'='plum2',
             'overallMapRate'='turquoise', 'concordMapRate'='lightsalmon','detected_num_genes'='skyblue2',
             'library_size'='palegreen3', 'RIN'='rosybrown3', 'Total_RNA_amount'='brown4', 'RNA_concentration'='blue3')

    ## Plot density graph for each variable
    p<-plotExplanatoryVariables(exp_vars, theme_size = 12, nvars_to_plot = Inf)
    p + scale_colour_manual(values = colors[c(variables)]) + labs(color="Variables")
    ## Save plot
    ggsave(filename = paste0('plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/ExplanatoryVars_', brain_region, '.pdf'), width = 35, height = 25, units = "cm")
    return(exp_vars)

}

## Plots and data
expl_vars_habenula <- as.data.frame(expl_var("habenula"))
expl_vars_amygdala <- as.data.frame(expl_var("amygdala"))



## Examine expression of most affected genes by each sample variable

## Plots of gene expression lognorm counts
plot_gene_expr <- function(brain_region, sample_var, gene_id){

   rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
   expl_vars <- eval(parse_expr(paste0('expl_vars_', brain_region)))

   if(sample_var=='Substance'){
        colors=c('Fentanyl'='turquoise3', 'Saline'='yellow3')
        x_label="Substance"
   }

    else if(sample_var=='Total_Num_Fentanyl_Sessions'){
        colors=c('24'='salmon', '22'='pink2')
        x_label="Total num of Fentanyl Sessions"
    }

    else if(sample_var=='Batch_RNA_extraction'){
        colors=c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2')
        x_label="RNA extraction Batch"
    }

    else if(sample_var=='Batch_lib_prep'){
        colors=c('1'='darkgoldenrod3', '2'='mediumpurple2', '3'= 'darkmagenta')
        x_label="Library preparation Batch"
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

    ## Percentage of variance explained by the variable
    percentage <- signif(expl_vars[gene_id, sample_var], digits=3)

    ## Boxplots for discrete variables
    if (sample_var=='Substance' | sample_var=='Total_Num_Fentanyl_Sessions' | sample_var=='Batch_RNA_extraction' | sample_var=='Batch_lib_prep') {
        plot <- ggplot(data = as.data.frame(data), mapping = aes(x = !! rlang::sym(sample_var),
                                                                 y = gene_expr, color = !! rlang::sym(sample_var))) +
            geom_boxplot(size = 0.25, width=0.32, color='black', outlier.color = "#FFFFFFFF") +
            geom_jitter( aes(shape=Batch_RNA_extraction), width = 0.15, alpha = 1, size = 1) +
            stat_smooth (geom="line", alpha=0.6, size=0.4, span=0.3, method = lm, aes(group=1), color='orangered3') +
            scale_color_manual(values = colors) +
            scale_shape_manual(name='Batch RNA extraction', values=c('1'=8, '2'=10, '3'=15)) +
            theme_bw() +
            guides(color="none") +
            labs(title = gene_ids,
                 subtitle = paste0("Variance explained: ", percentage, '%'),
                 y= 'lognorm counts', x = x_label) +
            theme(axis.title = element_text(size = (7)),
                  axis.text = element_text(size = (6)),
                  plot.title = element_text(hjust=0.5, size=7.5, face="bold"),
                  plot.subtitle = element_text(size = 7, color='gray40'),
                  legend.text = element_text(size=6),
                  legend.title = element_text(size=7))
    }

    ## Scatterplots for continuous variables
    else {
        colors <- c('mitoRate'='khaki3', 'totalAssignedGene'='plum2', 'overallMapRate'='turquoise', 'concordMapRate'='lightsalmon',
                    'detected_num_genes'='skyblue2', 'library_size'='palegreen3', 'RIN'='rosybrown3', 'Total_RNA_amount'='brown4',
                    'RNA_concentration'='blue3')

        plot <- ggplot(as.data.frame(data), aes(x=eval(parse_expr(sample_var)), y=gene_expr)) +
            geom_point( aes(shape=Batch_RNA_extraction), color=colors[sample_var], size=2) +
            stat_smooth (geom="line", alpha=0.4, size=0.4, span=0.25, method = lm, color='orangered3') +
            scale_shape_manual(name='Batch RNA extraction', values=c('1'=8, '2'=10, '3'=15)) +
            theme_bw() +
            guides(color="none") +
            labs(title = gene_ids,
                 subtitle = paste0("Variance explained: ", percentage, '%'),
                 y= 'lognorm counts', x = gsub('_', ' ', sample_var)) +
            theme(plot.margin=unit (c (0.4,0.1,0.4,0.1), 'cm'),
                  axis.title = element_text(size = (7)),
                  axis.text = element_text(size = (6)),
                  plot.title = element_text(hjust=0.5, size=7.5, face="bold"),
                  plot.subtitle = element_text(size = 7, color='gray40'),
                  legend.text = element_text(size=6),
                  legend.title = element_text(size=7))
    }

    return(plot)
}



gene_expr_expl_var <- function(brain_region, sample_var){

    expl_vars <- eval(parse_expr(paste0('expl_vars_', brain_region)))

    ## Order genes by % of variance explained by the sample var and extract top 6 affected genes
    top_genes_expl_vars <- expl_vars[order(expl_vars[,sample_var], decreasing =  TRUE),][1:6,]

    ## Plots of lognorm counts of those genes
    i=1
    plots=list()
    for (gene_id in rownames(top_genes_expl_vars)){
        plots[[i]] = plot_gene_expr(brain_region, sample_var, gene_id)
        i=i+1
    }

    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=2)
    ggsave(here(paste("plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/Top_affected_genes_by_", capitalize(sample_var), "_", brain_region, ".pdf", sep="")), width = 26, height = 13, units = "cm")
}


## Plots
gene_expr_expl_var('habenula', 'Substance')
gene_expr_expl_var('habenula', 'library_size')
gene_expr_expl_var('habenula', 'totalAssignedGene')
gene_expr_expl_var('habenula', 'Batch_RNA_extraction')
gene_expr_expl_var('habenula', 'RNA_concentration')
gene_expr_expl_var('habenula', 'Total_RNA_amount')
gene_expr_expl_var('habenula', 'mitoRate')

gene_expr_expl_var('amygdala', 'Substance')
gene_expr_expl_var('amygdala', 'library_size')
gene_expr_expl_var('amygdala', 'totalAssignedGene')
gene_expr_expl_var('amygdala', 'Batch_RNA_extraction')
gene_expr_expl_var('amygdala', 'RNA_concentration')
gene_expr_expl_var('amygdala', 'Total_RNA_amount')
gene_expr_expl_var('amygdala', 'mitoRate')







#######################   3.2 Variance Partition   #######################

## Fraction of variation attributable to each variable after correcting for all other variables


## 3.2.1 Canonical Correlation Analysis (CCA)

## Asses the correlation between each pair of sample variables

## Plot Heatmap of CC
plot_CCA<- function(brain_region){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## Define variables
    if (brain_region == 'habenula'){
        formula = ~ Substance + Batch_RNA_extraction + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                    library_size + detected_num_genes + RNA_concentration + Total_RNA_amount
    }
    else {
        formula = ~ Substance + Batch_RNA_extraction + Batch_lib_prep + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                    library_size + detected_num_genes + RNA_concentration + Total_RNA_amount
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
        filename=paste("plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/CCA_heatmap_", brain_region, ".pdf", sep="")
    )
}

plot_CCA('habenula')
plot_CCA('amygdala')





## 3.2.2 Model fit

## Fit a linear mixed model (LMM) that takes continuous variables as fixed effects and categorical variables as random effects

varPartAnalysis <- function(brain_region){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## Ignore genes with variance 0
    genes_var_zero<-which(apply(assays(RSE)$logcounts, 1, var)==0)
    RSE <- RSE[-genes_var_zero, ]

    ## Variables
    if (brain_region=='habenula'){
        formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + (1|Total_Num_Fentanyl_Sessions) +
                      mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                      library_size + detected_num_genes + RNA_concentration + Total_RNA_amount
    }
    else {
        formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + (1| Batch_lib_prep) + (1|Total_Num_Fentanyl_Sessions) +
                      mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                      library_size + detected_num_genes + RNA_concentration + Total_RNA_amount
    }

    ## Loop over each gene to fit model and extract variance explained by each variable
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts, formula, colData(RSE))

    # Sort variables by median fraction of variance explained
    vp <- sortCols(varPart)
    p <- plotVarPart(vp)
    ggsave(fileName,  p, width = 40, height = 20, units = "cm")
}

varPartAnalysis('habenula')
varPartAnalysis('amygdala')


## Violin plots






## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# setting  value
# version  R version 4.2.2 (2022-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-05-23
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   2.19.2 @ /private/var/folders/r6/94s0dsks4m3298b0d1mrp5tw0000gn/T/AppTranslocation/A38DEB34-4EE2-4D7E-B40B-40B091ABE956/d/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
