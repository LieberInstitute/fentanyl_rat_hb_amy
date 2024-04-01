
library(here)
library(readxl)
library(SummarizedExperiment)
library(ggplot2)
library(rlang)
library(scater)
library(cowplot)
library(Hmisc)
library(lme4)
library(variancePartition)
library(pheatmap)
library(smplot2)

library(reshape2)

library(sessioninfo)



################################################################################
##                        3. Explore gene-level effects
################################################################################

load(here('processed-data/04_EDA/02_PCA/rse_gene_habenula_filt.Rdata'), verbose = TRUE)
load(here('processed-data/04_EDA/02_PCA/rse_gene_amygdala_filt.Rdata'), verbose = TRUE)

## Load additional rat behavioral data (1st hr infusion slopes computed by simple linear regression)
covariate_data <- as.data.frame(read_excel("raw-data/covariate_sample_info.xlsx"))
colnames(covariate_data) <- gsub(' ', '_', colnames(covariate_data))
colnames(covariate_data) <- gsub('1st', 'First',  gsub("_\\((mg|#)\\)", '', colnames(covariate_data)))

## Append behavioral covariates to sample data
colData(rse_gene_habenula_filt) <- merge(colData(rse_gene_habenula_filt), covariate_data[,-c(2,3)], by='Rat_ID', sort=FALSE)
colData(rse_gene_amygdala_filt) <- merge(colData(rse_gene_amygdala_filt), covariate_data[,-c(2,3)], by='Rat_ID', sort=FALSE)

## Create table with rat behavioral data used
rat_behavioral_data <- subset(covariate_data, !Rat_ID %in% c('LgA 03', 'LgA 07', 'LgA 17'))[,-9]
write.csv(rat_behavioral_data, file="raw-data/rat_behavioral_data.csv", row.names = FALSE)



## 3.1 Analysis of explanatory variables

## 3.1.1 Computation of gene-wise expression variance percentages:
## Compute the % of gene expression variance explained by each sample variable

## Variables' colors
colors=c("Substance"= 'turquoise4', "Batch_RNA_extraction"='bisque2', "Batch_lib_prep"='blueviolet',
         "Total_Num_Fentanyl_Sessions"='indianred1', 'mitoRate'='khaki3', 'totalAssignedGene'='plum2',
         'overallMapRate'='turquoise', 'concordMapRate'='lightsalmon','detected_num_genes'='skyblue2',
         'library_size'='palegreen3', 'RIN'='rosybrown3', 'Total_RNA_amount'='brown4', 'RNA_concentration'='blue3',
         'Total_Intake'='yellow2', 'Last_Session_Intake'='pink','First_Hour_Infusion_Slope'='lightblue2')

## Plot density function for % of variance explained

expl_var<- function(brain_region, variables, all_vars){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## % of variance in gene expression explained by each variable

    if (brain_region=='habenula'){
        variables=variables[-3]
    }

    exp_vars<-getVarianceExplained(RSE, variables=variables, exprs_values = "logcounts")

    ## Plot density graph for each variable
    p<-plotExplanatoryVariables(exp_vars, theme_size = 16, nvars_to_plot = Inf)
    p + scale_colour_manual(values = colors[c(variables)]) + labs(color="Variables")
    ## Save plot
    ggsave(filename = paste0('plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/ExplanatoryVars_', brain_region, '_', all_vars,'.pdf'), width = 25, height = 20, units = "cm")
    return(exp_vars)

}

## Plots and data

## All sample variables
variables <- c("Substance", "Batch_RNA_extraction", "Batch_lib_prep", "Total_Num_Fentanyl_Sessions",
               "mitoRate", "concordMapRate","overallMapRate", "totalAssignedGene", "RIN", "detected_num_genes",
               "library_size", "Total_RNA_amount", "RNA_concentration", "Total_Intake", "Last_Session_Intake",
               "First_Hour_Infusion_Slope")

expl_vars_habenula_all <- as.data.frame(expl_var("habenula", variables, 'all'))
expl_vars_amygdala_all <- as.data.frame(expl_var("amygdala", variables, 'all'))

## Without Last_Session_Intake
variables <- variables[!variables=="Last_Session_Intake"]

expl_vars_habenula_without_last_sess_int <- as.data.frame(expl_var("habenula", variables, 'without_last_sess_int'))
expl_vars_amygdala_without_last_sess_int <- as.data.frame(expl_var("amygdala", variables, 'without_last_sess_int'))




## 3.1.2 Expression exploration of most affected genes
## Examine expression of most affected genes by each sample variable

## Plots of gene expression lognorm counts
plot_gene_expr <- function(brain_region, sample_var, gene_id){

   rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
   expl_vars <- eval(parse_expr(paste0('expl_vars_', brain_region, '_all')))

   if(sample_var=='Substance'){
        sample_colors=c('Fentanyl'='turquoise3', 'Saline'='yellow3')
        x_label="Substance"
   }

    else if(sample_var=='Total_Num_Fentanyl_Sessions'){
        sample_colors=c('24'='salmon', '22'='pink2')
        x_label="Total num of Fentanyl Sessions"
    }

    else if(sample_var=='Batch_RNA_extraction'){
        sample_colors=c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2')
        x_label="RNA extraction Batch"
    }

    else if(sample_var=='Batch_lib_prep'){
        sample_colors=c('1'='darkgoldenrod3', '2'='mediumpurple2', '3'= 'darkmagenta')
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
            scale_color_manual(values = sample_colors) +
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
        colors <- colors

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
gene_expr_expl_var('habenula', 'Total_Intake')
gene_expr_expl_var('habenula', 'Last_Session_Intake')
gene_expr_expl_var('habenula', 'First_Hour_Infusion_Slope')

gene_expr_expl_var('amygdala', 'Substance')
gene_expr_expl_var('amygdala', 'totalAssignedGene')
gene_expr_expl_var('amygdala', 'library_size')
gene_expr_expl_var('amygdala', 'Batch_RNA_extraction')
gene_expr_expl_var('amygdala', 'First_Hour_Infusion_Slope')
gene_expr_expl_var('amygdala', 'Total_RNA_amount')
gene_expr_expl_var('amygdala', 'Last_Session_Intake')
gene_expr_expl_var('amygdala', 'Total_Intake')
gene_expr_expl_var('amygdala', 'mitoRate')
gene_expr_expl_var('amygdala', 'RNA_concentration')





## 3.2 Variance Partition Analysis

## Fraction of variation attributable to each variable after correcting for all other variables


## 3.2.1 Canonical Correlation Analysis (CCA)

## Assess the correlation between each pair of sample variables

## Plot Heatmap of CC
plot_CCA<- function(brain_region){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## Define variables
    if (brain_region == 'habenula'){
        formula = ~ Substance + Batch_RNA_extraction + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                    library_size + detected_num_genes + RNA_concentration + Total_RNA_amount +
                    Total_Intake + Last_Session_Intake + First_Hour_Infusion_Slope
    }
    else {
        formula = ~ Substance + Batch_RNA_extraction + Batch_lib_prep + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                    library_size + detected_num_genes + RNA_concentration + Total_RNA_amount +
                    Total_Intake + Last_Session_Intake + First_Hour_Infusion_Slope
    }

    ## Correlations
    C=canCorPairs(formula, colData(RSE))
    ## Heatmap
    pheatmap(
        C,
        color = hcl.colors(50, "YlOrRd", rev = TRUE),
        fontsize=11,
        border_color = "black",
        height = 5.5,
        width = 6,
        filename=paste("plots/04_EDA/03_Explore_gene_level_effects/02_Var_Partition/CCA_heatmap_", brain_region, ".pdf", sep="")
    )

    return(C)
}

CCA_habenula <- plot_CCA('habenula')
CCA_amygdala <- plot_CCA('amygdala')



## Scatterplots/boxplots for each pair of correlated variables

corr_plots <- function(brain_region, sample_var1, sample_var2){

    rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
    CCA <- eval(parse_expr(paste0('CCA_', brain_region)))

    if(sample_var1=='Substance'){
        colors=c('Fentanyl'='turquoise3', 'Saline'='yellow3')
        x_label="Substance"
    }

    else if(sample_var1=='Total_Num_Fentanyl_Sessions'){
        colors=c('24'='salmon', '22'='pink2')
        x_label="Total num of Fentanyl Sessions"
    }

    else if(sample_var1=='Batch_RNA_extraction'){
        colors=c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2')
        x_label="RNA extraction Batch"
    }

    else if(sample_var1=='Batch_lib_prep'){
        colors=c('1'='darkgoldenrod3', '2'='mediumpurple2', '3'= 'darkmagenta')
        x_label="Library preparation Batch"
    }

    data <- colData(rse_gene)

    ## Boxplots for categorical variable vs continuous variable
    if (sample_var1=='Substance' | sample_var1=='Total_Num_Fentanyl_Sessions' | sample_var1=='Batch_RNA_extraction' | sample_var1=='Batch_lib_prep') {
        plot <- ggplot(data = as.data.frame(data), mapping = aes(x = !! rlang::sym(sample_var1),
                                                                 y = !! rlang::sym(sample_var2),
                                                                 color = !! rlang::sym(sample_var1))) +
            geom_boxplot(size = 0.25, width=0.32, color='black', outlier.color = "#FFFFFFFF") +
            geom_jitter( aes(shape=Batch_RNA_extraction), width = 0.15, alpha = 1, size = 1) +
            stat_smooth (geom="line", alpha=0.6, size=0.4, span=0.3, method = lm, aes(group=1), color='orangered3') +
            scale_color_manual(values = colors) +
            scale_shape_manual(name='Batch RNA extraction', values=c('1'=8, '2'=10, '3'=15)) +
            theme_bw() +
            guides(color="none") +
            labs(subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits=3)), y = gsub('_', ' ', sample_var2), x = x_label) +
            theme(axis.title = element_text(size = (7)),
                  axis.text = element_text(size = (6)),
                  plot.subtitle = element_text(size = 7, color='gray40'),
                  legend.text = element_text(size=6),
                  legend.title = element_text(size=7))
    }

    ## Scatterplots for continuous variable vs continuous variable
    else {
        ## Color and shape samples by RNA extraction batch
        colors=c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2')

        plot <- ggplot(as.data.frame(data), aes(x=eval(parse_expr(sample_var1)),
                                                y=eval(parse_expr(sample_var2)),
                                                color=Batch_RNA_extraction)) +
            geom_point( aes(shape=Batch_RNA_extraction), size=2) +
            stat_smooth (geom="line", alpha=0.4, size=0.4, span=0.25, method = lm, color='orangered3') +
            scale_shape_manual(name='Batch RNA extraction', values=c('1'=8, '2'=10, '3'=15)) +
            scale_color_manual(name='Batch RNA extraction', values = colors) +
            theme_bw() +
            labs(subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits=3)), y= gsub('_', ' ', sample_var2), x = gsub('_', ' ', sample_var1)) +
            theme(plot.margin=unit (c (0.4,0.1,0.4,0.1), 'cm'),
                  axis.title = element_text(size = (7)),
                  axis.text = element_text(size = (6)),
                  plot.subtitle = element_text(size = 7, color='gray40'),
                  legend.text = element_text(size=6),
                  legend.title = element_text(size=7))
    }

    return(plot)
}

## Multiple plots
multiple_corr_plots <- function(brain_region, sample_vars, name){

    ## Pairs of samples
    pairs <- list()
    for (i in 1:length(sample_vars)){
        pairs[[i]] <- merge(sample_vars[i], sample_vars[-c(1:i)])
    }

    k=1
    plots=list()
    for (i in 1:(length(pairs)-1)){
        for (j in 1:dim(pairs[[i]])[1]){
            plots[[k]] = corr_plots(brain_region, pairs[[i]][j, 1], pairs[[i]][j, 2])
            k=k+1
        }
    }

    if (length(sample_vars)>2){
        width = 45
        height = 20
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
                  plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
                  plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]], nrow=3)
    }
    else{
        width = 8.5
        height = 6
        plot_grid(plots[[1]])
    }

    ggsave(here(paste("plots/04_EDA/03_Explore_gene_level_effects/02_Var_Partition/Corr_", brain_region, "_", name, ".pdf", sep="")),
           width = width, height = height, units = "cm")

}


## Plots of correlated variables identified in the heatmaps

sample_vars <- c('Batch_RNA_extraction', 'library_size', 'totalAssignedGene', 'RNA_concentration', 'Total_RNA_amount', 'mitoRate')
multiple_corr_plots('habenula', sample_vars, 'Batch_RNA_extraction')

multiple_corr_plots('habenula', c('concordMapRate', 'overallMapRate'), 'concord_overallMapRates')
multiple_corr_plots('habenula', c('Substance', 'detected_num_genes'), 'Substance_detectedNumGenes')
multiple_corr_plots('habenula', c('mitoRate', 'detected_num_genes'), 'mitoRate_detectedNumGenes')
multiple_corr_plots('habenula', c('Substance', 'mitoRate'), 'Substance_mitoRate')

multiple_corr_plots('amygdala', c('library_size', 'totalAssignedGene'), 'libSize_totalAssig')
multiple_corr_plots('amygdala', c('concordMapRate', 'overallMapRate'), 'concord_overallMapRates')
multiple_corr_plots('amygdala', c('Substance', 'library_size'), 'Substance_librarySize')
multiple_corr_plots('amygdala', c('Substance', 'totalAssignedGene'), 'Substance_totalAssignedGene')
multiple_corr_plots('amygdala', c('Batch_RNA_extraction', 'Total_RNA_amount'), 'TotalRNAamount_BatchRNAextraction')
multiple_corr_plots('amygdala', c('Batch_lib_prep', 'Total_RNA_amount'), 'TotalRNAamount_BatchlibPrep')
multiple_corr_plots('amygdala', c('Batch_RNA_extraction', 'mitoRate'), 'mitoRate_BatchRNAextraction')
multiple_corr_plots('amygdala', c('Total_Num_Fentanyl_Sessions', 'RNA_concentration'), 'TotalNumFenSessions_RNAconcentration')




## 3.2.2 Model fit

## Fit a linear mixed model (LMM) that takes continuous variables as fixed effects and categorical variables as random effects

varPartAnalysis <- function(brain_region, formula, name){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## Ignore genes with variance 0
    genes_var_zero<-which(apply(assays(RSE)$logcounts, 1, var)==0)
    if (length(genes_var_zero)>0){
        RSE <- RSE[-genes_var_zero, ]
    }

    ## Loop over each gene to fit model and extract variance explained by each variable
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts, formula, colData(RSE))

    # Sort variables by median fraction of variance explained
    vp <- sortCols(varPart)
    p <- plotVarPart(vp)
    ggsave(filename=paste("plots/04_EDA/03_Explore_gene_level_effects/02_Var_Partition/VarPart_", brain_region, name, ".pdf", sep=""),
           p, width = 40, height = 20, units = "cm")
}

## Violin plots

## Habenula plots
## Define variables; random effects indicated with (1| )
formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + (1|Total_Num_Fentanyl_Sessions) +
                mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                library_size + detected_num_genes + RNA_concentration + Total_RNA_amount
varPartAnalysis('habenula', formula, '')

## Amygdala plots
formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + (1| Batch_lib_prep) + (1|Total_Num_Fentanyl_Sessions) +
                mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                library_size + detected_num_genes + RNA_concentration + Total_RNA_amount
varPartAnalysis('amygdala', formula, '')


## Plots without correlated variables

## Habenula plots without detected_num_genes, concordMapRate, library_size, totalAssignedGene, RNA_concentration and Total_RNA_amount; remove Total_Num_Fentanyl_Sessions since it doesn't contribute
formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + overallMapRate + RIN + mitoRate
varPartAnalysis('habenula', formula, '_withoutCorrVars')

## Amygdala plots without library_size, Total_RNA_amount, RNA_concentration, mitoRate, concordMapRate, detected_num_genes
## and Total_Num_Fentanyl_Sessions
formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + (1| Batch_lib_prep) +
             + overallMapRate + totalAssignedGene + RIN + mitoRate
varPartAnalysis('amygdala', formula, '_withoutCorrVars')
## Amygdala plot without totalAssignedGene as well
formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + (1| Batch_lib_prep) + overallMapRate + RIN + mitoRate
varPartAnalysis('amygdala', formula, '_withoutTotalAssignedGene')





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
