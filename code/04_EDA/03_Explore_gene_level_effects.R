
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
save(rse_gene_habenula_filt, file = 'processed-data/04_EDA/03_Explore_gene_level_effects/rse_gene_habenula_filt.Rdata')
save(rse_gene_amygdala_filt, file = 'processed-data/04_EDA/03_Explore_gene_level_effects/rse_gene_amygdala_filt.Rdata')

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

expl_var<- function(brain_region, variables, all_vars, substance){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
    file = paste0('plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/ExplanatoryVars_',
           brain_region, '_', all_vars,'.pdf')

    if(!is.null(substance)){
        RSE <- RSE[,RSE$Substance==substance]
        file = paste0('plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/ExplanatoryVars_',
                    brain_region, '_', substance, '_', all_vars,'.pdf')
    }

    ## % of variance in gene expression explained by each variable

    if (brain_region=='habenula'){
        variables=variables[-3]
    }

    exp_vars<-getVarianceExplained(RSE, variables=variables, exprs_values = "logcounts")

    ## Plot density graph for each variable
    p<-plotExplanatoryVariables(exp_vars, theme_size = 16, nvars_to_plot = Inf)
    p + scale_colour_manual(values = colors[c(variables)]) + labs(color="Variables")

    ## Save plot
    ggsave(filename = file, width = 25, height = 20, units = "cm")
    return(exp_vars)

}

## Plots and data

## All samples:
## All sample variables
variables <- names(colors)
expl_vars_habenula_all <- as.data.frame(expl_var("habenula", variables, 'all', NULL))
expl_vars_amygdala_all <- as.data.frame(expl_var("amygdala", variables, 'all', NULL))

## Without library_size and detected_num_genes
variables <- variables[!variables %in% c("library_size", "detected_num_genes")]
expl_vars_habenula_without_Lib_detectNum <- as.data.frame(expl_var("habenula", variables, 'without_Lib_detectNum', NULL))
expl_vars_amygdala_without_Lib_detectNum <- as.data.frame(expl_var("amygdala", variables, 'without_Lib_detectNum', NULL))


## Subset to fentanyl samples:
variables <- c("Total_Num_Fentanyl_Sessions", "mitoRate", "concordMapRate","overallMapRate", "totalAssignedGene",
               "RIN", "detected_num_genes", "library_size", "Total_RNA_amount", "RNA_concentration", "Total_Intake",
               "Last_Session_Intake", "First_Hour_Infusion_Slope")

expl_vars_habenula_Fent_all <- as.data.frame(expl_var("habenula", variables, 'all', 'Fentanyl'))
expl_vars_amygdala_Fent_all <- as.data.frame(expl_var("amygdala", variables, 'all', 'Fentanyl'))

variables <- variables[!variables %in% c("library_size", "detected_num_genes")]
expl_vars_habenula_Fent_without_Lib_detectNum <- as.data.frame(expl_var("habenula", variables, 'without_Lib_detectNum', 'Fentanyl'))
expl_vars_amygdala_Fent_without_Lib_detectNum <- as.data.frame(expl_var("amygdala", variables, 'without_Lib_detectNum', 'Fentanyl'))




## 3.1.2 Expression exploration of most affected genes
## Examine expression of most affected genes by each sample variable

## Plots of gene expression lognorm counts
plot_gene_expr <- function(brain_region, sample_var, gene_id, substance){

    rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    if(substance!='allSamples'){
        rse_gene <- rse_gene[,colData(rse_gene)$Substance==substance]
        expl_vars <- eval(parse_expr(paste0('expl_vars_', brain_region, '_Fent_all')))
    }
    else{
        expl_vars <- eval(parse_expr(paste0('expl_vars_', brain_region, '_all')))
    }

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


gene_expr_expl_var <- function(brain_region, sample_var, substance){

    if(substance!='allSamples'){
        expl_vars <- eval(parse_expr(paste0('expl_vars_', brain_region, '_Fent_all')))
    }
    else{
        expl_vars <- eval(parse_expr(paste0('expl_vars_', brain_region, '_all')))
    }

    ## Order genes by % of variance explained by the sample var and extract top 6 affected genes
    top_genes_expl_vars <- expl_vars[order(expl_vars[,sample_var], decreasing =  TRUE),][1:6,]

    ## Plots of lognorm counts of those genes
    i=1
    plots=list()
    for (gene_id in rownames(top_genes_expl_vars)){
        plots[[i]] = plot_gene_expr(brain_region, sample_var, gene_id, substance)
        i=i+1
    }

    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=2)
    ggsave(here(paste("plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/Top_affected_genes_by_", capitalize(sample_var), "_",
                      brain_region, '_', substance, ".pdf", sep="")), width = 26, height = 13, units = "cm")
}


## Plots for all samples:
gene_expr_expl_var('habenula', 'Substance', 'allSamples')
gene_expr_expl_var('habenula', 'library_size', 'allSamples')
gene_expr_expl_var('habenula', 'totalAssignedGene', 'allSamples')
gene_expr_expl_var('habenula', 'Batch_RNA_extraction', 'allSamples')
gene_expr_expl_var('habenula', 'RNA_concentration', 'allSamples')
gene_expr_expl_var('habenula', 'Total_RNA_amount', 'allSamples')
gene_expr_expl_var('habenula', 'mitoRate', 'allSamples')
gene_expr_expl_var('habenula', 'Total_Intake', 'allSamples')
gene_expr_expl_var('habenula', 'Last_Session_Intake', 'allSamples')
gene_expr_expl_var('habenula', 'First_Hour_Infusion_Slope', 'allSamples')

gene_expr_expl_var('amygdala', 'Substance', 'allSamples')
gene_expr_expl_var('amygdala', 'totalAssignedGene', 'allSamples')
gene_expr_expl_var('amygdala', 'library_size', 'allSamples')
gene_expr_expl_var('amygdala', 'Batch_RNA_extraction', 'allSamples')
gene_expr_expl_var('amygdala', 'First_Hour_Infusion_Slope', 'allSamples')
gene_expr_expl_var('amygdala', 'Total_RNA_amount', 'allSamples')
gene_expr_expl_var('amygdala', 'Last_Session_Intake', 'allSamples')
gene_expr_expl_var('amygdala', 'Total_Intake', 'allSamples')
gene_expr_expl_var('amygdala', 'mitoRate', 'allSamples')
gene_expr_expl_var('amygdala', 'RNA_concentration', 'allSamples')


## Plots for fentanyl samples
gene_expr_expl_var('habenula', 'detected_num_genes', 'Fentanyl')
gene_expr_expl_var('habenula', 'overallMapRate', 'Fentanyl')
gene_expr_expl_var('habenula', 'Last_Session_Intake', 'Fentanyl')
gene_expr_expl_var('habenula', 'Total_Intake', 'Fentanyl')
gene_expr_expl_var('habenula', 'First_Hour_Infusion_Slope', 'Fentanyl')

gene_expr_expl_var('amygdala', 'mitoRate', 'Fentanyl')
gene_expr_expl_var('amygdala', 'library_size', 'Fentanyl')
gene_expr_expl_var('amygdala', 'Last_Session_Intake', 'Fentanyl')
gene_expr_expl_var('amygdala', 'Total_Intake', 'Fentanyl')
gene_expr_expl_var('amygdala', 'First_Hour_Infusion_Slope', 'Fentanyl')





## 3.2 Variance Partition Analysis

## Fraction of variation attributable to each variable after correcting for all other variables


## 3.2.1 Canonical Correlation Analysis (CCA)
## Assess the correlation between each pair of sample variables

## Plot Heatmap of CC
plot_CCA<- function(brain_region, all_vars, substance){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    ## For fentanyl samples only
    if(substance!='allSamples'){
        RSE<-RSE[,colData(RSE)$Substance==substance]

        if(all_vars=='all_vars'){
            ## Define variables
            formula = ~ Total_Num_Fentanyl_Sessions + mitoRate + overallMapRate + concordMapRate + totalAssignedGene +
                RIN + library_size + detected_num_genes + RNA_concentration + Total_RNA_amount +
                Total_Intake + Last_Session_Intake + First_Hour_Infusion_Slope
        }
        else{
            formula = ~ Total_Num_Fentanyl_Sessions + mitoRate + overallMapRate + concordMapRate + totalAssignedGene +
                RIN  + RNA_concentration + Total_RNA_amount + Total_Intake + Last_Session_Intake + First_Hour_Infusion_Slope
        }

    }

    ## For both fentanyl and saline samples
    else{
        if (brain_region == 'habenula'){
            if(all_vars=='all_vars'){
                formula = ~ Substance + Batch_RNA_extraction + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                    library_size + detected_num_genes + RNA_concentration + Total_RNA_amount +
                    Total_Intake + Last_Session_Intake + First_Hour_Infusion_Slope
            }
            else{
                formula = ~ Substance + Batch_RNA_extraction + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                    RNA_concentration + Total_RNA_amount + Total_Intake + Last_Session_Intake +
                    First_Hour_Infusion_Slope
            }
        }
        else {
            if(all_vars=='all_vars'){
                formula = ~ Substance + Batch_RNA_extraction + Batch_lib_prep + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN +
                    library_size + detected_num_genes + RNA_concentration + Total_RNA_amount +
                    Total_Intake + Last_Session_Intake + First_Hour_Infusion_Slope
            }
            else{
                formula = ~ Substance + Batch_RNA_extraction + Batch_lib_prep + Total_Num_Fentanyl_Sessions +
                    mitoRate + overallMapRate + concordMapRate + totalAssignedGene + RIN + RNA_concentration +
                    Total_RNA_amount + Total_Intake + Last_Session_Intake + First_Hour_Infusion_Slope
            }
        }
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
        filename=paste("plots/04_EDA/03_Explore_gene_level_effects/02_Var_Partition/CCA_heatmap_",
                       brain_region, '_', substance, '_', all_vars, ".pdf", sep="")
    )

    return(C)
}

## All samples
CCA_habenula_all_vars_allSamples <- plot_CCA('habenula', 'all_vars', 'allSamples')
CCA_amygdala_all_vars_allSamples <- plot_CCA('amygdala', 'all_vars', 'allSamples')
CCA_habenula_without_Lib_detectNum_allSamples <- plot_CCA('habenula', 'without_Lib_detectNum', 'allSamples')
CCA_amygdala_without_Lib_detectNum_allSamples <- plot_CCA('amygdala', 'without_Lib_detectNum', 'allSamples')

## Fentanyl samples only
CCA_habenula_all_vars_Fentanyl <- plot_CCA('habenula', 'all_vars', 'Fentanyl')
CCA_amygdala_all_vars_Fentanyl <- plot_CCA('amygdala', 'all_vars', 'Fentanyl')
CCA_habenula_without_Lib_detectNum_Fentanyl <- plot_CCA('habenula', 'without_Lib_detectNum', 'Fentanyl')
CCA_amygdala_without_Lib_detectNum_Fentanyl <- plot_CCA('amygdala', 'without_Lib_detectNum', 'Fentanyl')



## Scatterplots/boxplots for each pair of correlated variables

corr_plots <- function(brain_region, sample_var1, sample_var2, substance){

    rse_gene <- eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
    CCA <- eval(parse_expr(paste0('CCA_', brain_region, '_all_vars_allSamples')))

    ## For fentanyl samples
    if(substance!='allSamples'){
        rse_gene <- rse_gene[,colData(rse_gene)$Substance==substance]
        CCA <- eval(parse_expr(paste0('CCA_', brain_region, '_all_vars_Fentanyl')))
    }

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
            geom_boxplot(size = 0.25, width=0.32, color='black', outlier.color = NA) +
            geom_jitter( aes(shape=Batch_RNA_extraction), width = 0.15, alpha = 1, size = 1) +
            stat_smooth (geom="line", alpha=0.6, size=0.4, span=0.3, method = lm, aes(group=1), color='orangered3') +
            scale_color_manual(values = colors) +
            scale_shape_manual(name='Batch RNA extraction', values=c('1'=8, '2'=10, '3'=15)) +
            theme_bw() +
            guides(color="none") +
            labs(subtitle = paste0("Corr: ", signif(CCA[sample_var1, sample_var2], digits=3)), y = gsub('_', ' ', sample_var2), x = x_label) +
            theme(axis.title = element_text(size = (8)),
                  axis.text = element_text(size = (6)),
                  plot.subtitle = element_text(size = 7, color='gray40'),
                  legend.text = element_text(size=7),
                  legend.title = element_text(size=8))
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
                  axis.title = element_text(size = (8)),
                  axis.text = element_text(size = (6)),
                  plot.subtitle = element_text(size = 7, color='gray40'),
                  legend.text = element_text(size=7),
                  legend.title = element_text(size=8))
    }

    return(plot)
}

## Multiple plots
multiple_corr_plots <- function(brain_region, sample_vars, name, substance){

    ## Pairs of samples
    pairs <- list()
    for (i in 1:length(sample_vars)){
        pairs[[i]] <- merge(sample_vars[i], sample_vars[-c(1:i)])
    }

    k=1
    plots=list()
    for (i in 1:(length(pairs)-1)){
        for (j in 1:dim(pairs[[i]])[1]){
            plots[[k]] = corr_plots(brain_region, pairs[[i]][j, 1], pairs[[i]][j, 2], substance)
            k=k+1
        }
    }

    if (length(sample_vars)==6){
        width = 45
        height = 20
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
                  plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]],
                  plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]], nrow=3)
    }
    else if (length(sample_vars)==5){
        width = 45
        height = 13
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
                  plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], nrow=3)
    }
    else if(length(sample_vars)==4){
        width = 27
        height = 12
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                  plots[[5]], plots[[6]], ncol=3)
    }
    else if(length(sample_vars)==3){
        width = 27
        height = 6
        plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=3)
    }
    else{
        width = 8.5
        height = 6
        plot_grid(plots[[1]])
    }

    ggsave(here(paste("plots/04_EDA/03_Explore_gene_level_effects/02_Var_Partition/Corr_", brain_region, "_", name,
                      '_', substance,".pdf", sep="")), width = width, height = height, units = "cm")

}


## Plots of correlated variables identified in the heatmaps for all samples

sample_vars <- c('Batch_RNA_extraction', 'library_size', 'totalAssignedGene', 'RNA_concentration',
                 'Total_RNA_amount', 'mitoRate')
multiple_corr_plots('habenula', sample_vars, 'Batch_RNA_extraction' ,'allSamples')
sample_vars <- c('Substance', 'detected_num_genes', 'First_Hour_Infusion_Slope',
                 'Total_Intake', 'Last_Session_Intake')
multiple_corr_plots('habenula', sample_vars, 'Substance' ,'allSamples')
sample_vars <- c('Total_Num_Fentanyl_Sessions', 'Total_Intake', 'Last_Session_Intake', 'First_Hour_Infusion_Slope')
multiple_corr_plots('habenula', sample_vars, 'Total_Num_Fentanyl_Sessions' ,'allSamples')
multiple_corr_plots('habenula', c('concordMapRate', 'overallMapRate'), 'concord_overallMapRates', 'allSamples')
multiple_corr_plots('habenula', c('mitoRate', 'detected_num_genes'), 'mitoRate_detectedNumGenes', 'allSamples')
multiple_corr_plots('habenula', c('Substance', 'mitoRate'), 'Substance_mitoRate', 'allSamples')

sample_vars <- c('Substance', 'totalAssignedGene', 'library_size', 'First_Hour_Infusion_Slope',
                 'Last_Session_Intake', 'Total_Intake')
multiple_corr_plots('amygdala', sample_vars, 'Substance', 'allSamples')
sample_vars <- c('Total_Num_Fentanyl_Sessions', 'Total_Intake', 'Last_Session_Intake', 'First_Hour_Infusion_Slope')
multiple_corr_plots('amygdala', sample_vars, 'Total_Num_Fentanyl_Sessions' ,'allSamples')
multiple_corr_plots('amygdala', c('concordMapRate', 'overallMapRate'), 'concord_overallMapRates', 'allSamples')
multiple_corr_plots('amygdala', c('Batch_RNA_extraction', 'mitoRate'), 'mitoRate_BatchRNAextraction', 'allSamples')
multiple_corr_plots('amygdala', c('Batch_RNA_extraction', 'Total_RNA_amount'), 'TotalRNAamount_BatchRNAextraction', 'allSamples')
multiple_corr_plots('amygdala', c('Batch_lib_prep', 'Total_RNA_amount'), 'TotalRNAamount_BatchlibPrep', 'allSamples')
multiple_corr_plots('amygdala', c('Total_Num_Fentanyl_Sessions', 'RNA_concentration'), 'TotalNumFenSessions_RNAconcentration', 'allSamples')


## For fentanyl samples
sample_vars <- c('Total_Intake', 'Total_RNA_amount', 'totalAssignedGene')
multiple_corr_plots('habenula', sample_vars, 'Total_Intake' ,'Fentanyl')
sample_vars <- c('Total_Num_Fentanyl_Sessions', 'Total_Intake', 'Last_Session_Intake', 'First_Hour_Infusion_Slope')
multiple_corr_plots('habenula', sample_vars, 'Total_Num_Fentanyl_Sessions' ,'Fentanyl')
multiple_corr_plots('habenula', c('library_size', 'totalAssignedGene'), 'lib_size_totalAssGene' ,'Fentanyl')
multiple_corr_plots('habenula', c('Total_RNA_amount', 'totalAssignedGene'), 'Total_RNA_amount_totalAssGene' ,'Fentanyl')
multiple_corr_plots('habenula', c('mitoRate', 'totalAssignedGene'), 'mitoRate_totalAssGene' ,'Fentanyl')
multiple_corr_plots('habenula', c('Total_RNA_amount', 'RNA_concentration'), 'Total_RNA_amount_RNAconc' ,'Fentanyl')
multiple_corr_plots('habenula', c('Total_RNA_amount', 'detected_num_genes'), 'Total_RNA_amount_detected_num' ,'Fentanyl')
multiple_corr_plots('habenula', c('RIN', 'detected_num_genes'), 'RIN_detected_num' ,'Fentanyl')
multiple_corr_plots('habenula', c('concordMapRate', 'overallMapRate'), 'concord_overallMapRates' ,'Fentanyl')
multiple_corr_plots('habenula', c('mitoRate', 'overallMapRate'), 'mitoRate_overallMapRate' ,'Fentanyl')
multiple_corr_plots('habenula', c('detected_num_genes', 'overallMapRate'), 'detected_num_overallMapRate' ,'Fentanyl')
multiple_corr_plots('habenula', c('concordMapRate', 'detected_num_genes'), 'concord_detected_num' ,'Fentanyl')
multiple_corr_plots('habenula', c('mitoRate', 'detected_num_genes'), 'mitoRate_detected_num' ,'Fentanyl')

sample_vars <- c('Total_Intake', 'overallMapRate', 'concordMapRate', 'detected_num_genes')
multiple_corr_plots('amygdala', sample_vars, 'Total_Intake' ,'Fentanyl')
sample_vars <- c('Total_Num_Fentanyl_Sessions', 'Total_Intake', 'Last_Session_Intake', 'First_Hour_Infusion_Slope')
multiple_corr_plots('amygdala', sample_vars, 'Total_Num_Fentanyl_Sessions' ,'Fentanyl')
multiple_corr_plots('amygdala', c('library_size', 'totalAssignedGene', 'mitoRate'), 'lib_size_totalAssGene_mito' ,'Fentanyl')
multiple_corr_plots('amygdala', c('Total_Num_Fentanyl_Sessions', 'Total_RNA_amount', 'RNA_concentration'), 'Total_Num_FentSess_TotalRNAam_RNAconc' ,'Fentanyl')




## 3.2.2 Model fit

## Fit a linear mixed model (LMM) that takes continuous variables as fixed effects and categorical variables as random effects

varPartAnalysis <- function(brain_region, formula, name, substance){

    RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))

    if(substance!='allSamples'){
        RSE <- RSE[,colData(RSE)$Substance==substance]
    }

    ## Ignore genes with variance 0
    genes_var_zero<-which(apply(assays(RSE)$logcounts, 1, var)==0)
    if (length(genes_var_zero)>0){
        RSE <- RSE[-genes_var_zero, ]
    }

    ## Loop over each gene to fit model and extract variance explained by each variable
    varPart <- fitExtractVarPartModel(assays(RSE)$logcounts, formula, colData(RSE))

    # Sort variables by median fraction of variance explained
    vp <- sortCols(varPart)
    p <- plotVarPart(vp, col = colors)
    ggsave(filename=paste("plots/04_EDA/03_Explore_gene_level_effects/02_Var_Partition/VarPart_",
                          brain_region, '_', substance, '_', name, ".pdf", sep=""),
           p, width = 35, height = 15, units = "cm")
}

## Violin plots without correlated variables

## All samples:

## Habenula covariates discarding:
##  * 1. Variables correlated with Substance: First_Hour_Infusion_Slope, Total_Intake, and Last_Session_Intake
##  * 2. Variables correlated with RNA extraction batch: mitoRate, totalAssignedGene,
##                                                       RNA_concentration, and Total_RNA_amount
##  * 3. Total num session / variables correlated with it: Total_Num_Fentanyl_Sessions
##  * 4. Variables highly correlated with any other that explains higher %s of variance: overallMapRate

## Define variables; random effects indicated with (1| )
formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + concordMapRate + RIN
varPartAnalysis('habenula', formula, 'finalVariableSet', 'allSamples')


## Amygdala covariates discarding:
##  * 1. Variables correlated with Substance: totalAssignedGene, First_Hour_Infusion_Slope,
##                                            Last_Session_Intake, and Total_Intake
##  * 2. Variables correlated with RNA extraction batch: mitoRate and Total_RNA_amount
##  * 3. Total num session / variables correlated with it: Total_Num_Fentanyl_Sessions and RNA_concentration
##  * 4. Variables highly correlated with any other that explains higher %s of variance: concordMapRate

formula <-  ~ (1|Substance) + (1|Batch_RNA_extraction) + (1| Batch_lib_prep) + overallMapRate  + RIN
varPartAnalysis('amygdala', formula, 'finalVariableSet', 'allSamples')



## Fentanyl samples only:

## Habenula covariates for Total Intake DEA:
##  * 1. Variables correlated with Total Intake: totalAssignedGene and Total_RNA_amount
##  * 2. Total num session / variables correlated with it: Total_Num_Fentanyl_Sessions and concordMapRate
##  * 3. Variables highly correlated with any other that explains higher %s of variance: mitoRate

formula <-  ~  Total_Intake + RIN + RNA_concentration + overallMapRate
varPartAnalysis('habenula', formula, 'Total_Intake_finalVariableSet', 'Fentanyl')

## Habenula covariates for Last Session Intake DEA:
##  * 1. Variables correlated with Last Session Intake: Total_Num_Fentanyl_Sessions
##  * 2. Total num session / variables correlated with it: concordMapRate
##  * 3. Variables highly correlated with any other that explains higher %s of variance: totalAssignedGene, Total_RNA_amount, and overallMapRate

formula <-  ~  Last_Session_Intake + RIN + RNA_concentration + mitoRate
varPartAnalysis('habenula', formula, 'Last_Session_Intake_finalVariableSet', 'Fentanyl')

## Habenula covariates for First hr Infusion Slope DEA:
##  * 1. Variables correlated with First hr Infusion Slope: Total_Num_Fentanyl_Sessions
##  * 2. Total num session / variables correlated with it:  concordMapRate
##  * 3. Variables highly correlated with any other that explains higher %s of variance: totalAssignedGene, Total_RNA_amount, and overallMapRate

formula <-  ~  First_Hour_Infusion_Slope + RIN + RNA_concentration + mitoRate
varPartAnalysis('habenula', formula, 'First_Hour_Infusion_Slope_finalVariableSet', 'Fentanyl')


## Amygdala covariates for Total Intake DEA:
##  * 1. Variables correlated with Total Intake: concordMapRate and overallMapRate
##  * 2. Total num session / variables correlated with it: Total_Num_Fentanyl_Sessions, RNA_concentration, and Total_RNA_amount
##  * 3. Variables highly correlated with any other that explains higher %s of variance: totalAssignedGene

formula <-  ~  Total_Intake + RIN + mitoRate
varPartAnalysis('amygdala', formula, 'Total_Intake_finalVariableSet', 'Fentanyl')

## Amygdala covariates for Last Session Intake DEA:
##  * 1. Variables correlated with Last Session Intake: RNA_concentration, Total_Num_Fentanyl_Sessions, and mitoRate
##  * 2. Total num session / variables correlated with it: Total_RNA_amount
##  * 3. Variables highly correlated with any other that explains higher %s of variance: overallMapRate

formula <-  ~  Last_Session_Intake + RIN + totalAssignedGene + concordMapRate
varPartAnalysis('amygdala', formula, 'Last_Session_Intake_finalVariableSet', 'Fentanyl')

## Amygdala covariates for First hr Infusion Slope DEA:
##  * 1. Variables correlated with First hr Infusion Slope:  RNA_concentration, Total_Num_Fentanyl_Sessions, overallMapRate,
##                                                           and concordMapRate
##  * 2. Total num session / variables correlated with it: Total_RNA_amount
##  * 3. Variables highly correlated with any other that explains higher %s of variance: totalAssignedGene

formula <-  ~  First_Hour_Infusion_Slope + RIN + mitoRate
varPartAnalysis('amygdala', formula, 'First_Hour_Infusion_Slope_finalVariableSet', 'Fentanyl')







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
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
# date     2024-04-15
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version       date (UTC) lib source
# abind                  1.4-5         2016-07-21 [1] CRAN (R 4.3.0)
# aod                    1.3.3         2023-12-13 [1] CRAN (R 4.3.1)
# backports              1.4.1         2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3         2015-07-28 [1] CRAN (R 4.3.0)
# beachmat               2.18.1        2024-02-17 [1] Bioconductor 3.18 (R 4.3.2)
# beeswarm               0.4.0         2021-06-01 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0        2023-10-26 [1] Bioconductor
# BiocGenerics         * 0.48.1        2023-11-02 [1] Bioconductor
# BiocNeighbors          1.20.2        2024-01-13 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel         * 1.36.0        2023-10-26 [1] Bioconductor
# BiocSingular           1.18.0        2023-11-06 [1] Bioconductor
# bitops                 1.0-7         2021-04-24 [1] CRAN (R 4.3.0)
# boot                   1.3-30        2024-02-26 [1] CRAN (R 4.3.1)
# broom                  1.0.5         2023-06-09 [1] CRAN (R 4.3.0)
# car                    3.1-2         2023-03-30 [1] CRAN (R 4.3.0)
# carData                3.0-5         2022-01-06 [1] CRAN (R 4.3.0)
# caTools                1.18.2        2021-03-28 [1] CRAN (R 4.3.0)
# cellranger             1.1.0         2016-07-27 [1] CRAN (R 4.3.0)
# checkmate              2.3.1         2023-12-04 [1] CRAN (R 4.3.1)
# cli                    3.6.2         2023-12-11 [1] CRAN (R 4.3.1)
# cluster                2.1.6         2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-19        2023-02-01 [1] CRAN (R 4.3.2)
# colorspace             2.1-0         2023-01-23 [1] CRAN (R 4.3.0)
# corpcor                1.6.10        2021-09-16 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.3         2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2         2022-09-29 [1] CRAN (R 4.3.0)
# data.table             1.15.2        2024-02-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0        2023-11-06 [1] Bioconductor
# DelayedMatrixStats     1.24.0        2023-11-06 [1] Bioconductor
# digest                 0.6.34        2024-01-11 [1] CRAN (R 4.3.1)
# dplyr                  1.1.4         2023-11-17 [1] CRAN (R 4.3.1)
# EnvStats               2.8.1         2023-08-22 [1] CRAN (R 4.3.0)
# evaluate               0.23          2023-11-01 [1] CRAN (R 4.3.1)
# fANCOVA                0.6-1         2020-11-13 [1] CRAN (R 4.3.0)
# fansi                  1.0.6         2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1         2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1         2023-02-24 [1] CRAN (R 4.3.0)
# foreign                0.8-86        2023-11-28 [1] CRAN (R 4.3.1)
# Formula                1.2-5         2023-02-24 [1] CRAN (R 4.3.0)
# generics               0.1.3         2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.6        2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11        2024-02-17 [1] Bioconductor
# GenomicRanges        * 1.54.1        2023-10-30 [1] Bioconductor
# ggbeeswarm             0.7.2         2023-04-29 [1] CRAN (R 4.3.0)
# gghalves               0.1.4         2022-11-20 [1] CRAN (R 4.3.0)
# ggplot2              * 3.5.0         2024-02-23 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0         2023-02-10 [1] CRAN (R 4.3.0)
# ggrepel                0.9.5         2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4         2022-10-13 [1] CRAN (R 4.3.0)
# glue                   1.7.0         2024-01-09 [1] CRAN (R 4.3.1)
# gplots                 3.1.3.1       2024-02-02 [1] CRAN (R 4.3.1)
# gridExtra              2.3           2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4         2023-08-21 [1] CRAN (R 4.3.0)
# gtools                 3.9.5         2023-11-20 [1] CRAN (R 4.3.1)
# here                 * 1.0.1         2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                * 5.1-1         2023-09-12 [1] CRAN (R 4.3.0)
# htmlTable              2.4.2         2023-10-29 [1] CRAN (R 4.3.1)
# htmltools              0.5.7         2023-11-03 [1] CRAN (R 4.3.1)
# htmlwidgets            1.6.4         2023-12-06 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0        2023-10-26 [1] Bioconductor
# irlba                  2.3.5.1       2022-10-03 [1] CRAN (R 4.3.0)
# iterators              1.0.14        2022-02-05 [1] CRAN (R 4.3.0)
# KernSmooth             2.23-22       2023-07-10 [1] CRAN (R 4.3.2)
# knitr                  1.45          2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3         2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.22-5        2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4         2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1        2023-11-02 [1] Bioconductor
# lme4                 * 1.1-35.1.9000 2024-03-12 [1] Github (lme4/lme4@d1a36a8)
# lmerTest               3.1-3         2020-10-23 [1] CRAN (R 4.3.0)
# magrittr               2.0.3         2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60.0.1    2024-01-13 [1] CRAN (R 4.3.1)
# Matrix               * 1.6-5         2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0        2023-10-26 [1] Bioconductor
# matrixStats          * 1.2.0         2023-12-11 [1] CRAN (R 4.3.1)
# minqa                  1.2.6         2023-09-11 [1] CRAN (R 4.3.0)
# munsell                0.5.0         2018-06-12 [1] CRAN (R 4.3.0)
# mvtnorm                1.2-4         2023-11-27 [1] CRAN (R 4.3.1)
# nlme                   3.1-164       2023-11-27 [1] CRAN (R 4.3.1)
# nloptr                 2.0.3         2022-05-26 [1] CRAN (R 4.3.0)
# nnet                   7.3-19        2023-05-03 [1] CRAN (R 4.3.2)
# numDeriv               2016.8-1.1    2019-06-06 [1] CRAN (R 4.3.0)
# pbkrtest               0.5.2         2023-01-19 [1] CRAN (R 4.3.0)
# pheatmap             * 1.0.12        2019-01-04 [1] CRAN (R 4.3.0)
# pillar                 1.9.0         2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3         2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9         2023-10-02 [1] CRAN (R 4.3.1)
# purrr                  1.0.2         2023-08-10 [1] CRAN (R 4.3.0)
# pwr                    1.3-0         2020-03-17 [1] CRAN (R 4.3.0)
# R6                     2.5.1         2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7         2023-12-11 [1] CRAN (R 4.3.1)
# rbibutils              2.2.16        2023-10-25 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3         2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12        2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14     2024-01-09 [1] CRAN (R 4.3.1)
# Rdpack                 2.6           2023-11-08 [1] CRAN (R 4.3.1)
# readxl               * 1.4.3         2023-07-06 [1] CRAN (R 4.3.0)
# remaCor                0.0.18        2024-02-08 [1] CRAN (R 4.3.0)
# reshape2             * 1.4.4         2020-04-09 [1] CRAN (R 4.3.0)
# RhpcBLASctl            0.23-42       2023-02-11 [1] CRAN (R 4.3.0)
# rlang                * 1.1.3         2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26          2024-03-05 [1] CRAN (R 4.3.1)
# rpart                  4.1.23        2023-12-05 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4         2023-11-05 [1] CRAN (R 4.3.1)
# rstatix                0.7.2         2023-02-01 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0        2023-07-07 [1] CRAN (R 4.3.0)
# rsvd                   1.0.5         2021-04-16 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0         2023-10-26 [1] Bioconductor
# S4Vectors            * 0.40.2        2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# ScaledMatrix           1.10.0        2023-11-06 [1] Bioconductor
# scales                 1.3.0         2023-11-28 [1] CRAN (R 4.3.1)
# scater               * 1.30.1        2023-11-16 [1] Bioconductor
# scuttle              * 1.12.0        2023-11-06 [1] Bioconductor
# sdamr                  0.2.0         2022-11-16 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2         2021-12-06 [1] CRAN (R 4.3.0)
# SingleCellExperiment * 1.24.0        2023-11-06 [1] Bioconductor
# smplot2              * 0.1.0         2024-03-13 [1] Github (smin95/smplot2@052f4f9)
# SparseArray            1.2.4         2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats      1.14.0        2023-10-26 [1] Bioconductor
# statmod                1.5.0         2023-01-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3         2023-12-11 [1] CRAN (R 4.3.1)
# stringr                1.5.1         2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0        2023-11-06 [1] Bioconductor
# systemfonts            1.0.5         2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7         2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1         2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.1         2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.0         2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4         2023-10-22 [1] CRAN (R 4.3.1)
# variancePartition    * 1.32.5        2024-02-17 [1] Bioconductor 3.18 (R 4.3.2)
# vctrs                  0.6.5         2023-12-01 [1] CRAN (R 4.3.1)
# vipor                  0.4.7         2023-12-18 [1] CRAN (R 4.3.1)
# viridis                0.6.5         2024-01-29 [1] CRAN (R 4.3.1)
# viridisLite            0.4.2         2023-05-02 [1] CRAN (R 4.3.0)
# withr                  3.0.0         2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.42          2024-02-08 [1] CRAN (R 4.3.1)
# XVector                0.42.0        2023-10-26 [1] Bioconductor
# zlibbioc               1.48.0        2023-10-26 [1] Bioconductor
# zoo                    1.8-12        2023-04-13 [1] CRAN (R 4.3.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
