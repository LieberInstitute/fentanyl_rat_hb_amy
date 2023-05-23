library(here)
library(SummarizedExperiment)
library(ggplot2)
library(scater)
library(rlang)
library(smplot2)
library(cowplot)
library(Hmisc)
library(sessioninfo)



################################################################################
##                        3. Explore gene-level effects
################################################################################

load(here('processed-data/04_EDA/02_PCA/rse_gene_habenula_filt.Rdata'), verbose = TRUE)
load(here('processed-data/04_EDA/02_PCA/rse_gene_amygdala_filt.Rdata'), verbose = TRUE)


#######################   3.1 Explanatory variables   #######################

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
    percentage <- signif(exp_vars[gene_id, sample_var], digits=3)

    ## Boxplots for discrete variables
    if (sample_var=='Substance' | sample_var=='Total_Num_Fentanyl_Sessions' | sample_var=='Batch_RNA_extraction' | sample_var=='Batch_lib_prep') {
        plot <- ggplot(data = as.data.frame(data), mapping = aes(x = !! rlang::sym(sample_var),
                                                                 y = gene_expr, color = !! rlang::sym(sample_var),
                                                                 color = !! rlang::sym(sample_var))) +
            geom_boxplot(size = 0.25, width=0.32, color='black', outlier.color = "#FFFFFFFF") +
            geom_jitter(width = 0.15, alpha = 1, size = 1) +
            scale_color_manual(values = colors) +
            sm_hgrid() +
            labs(title = gene_ids,
                 subtitle = paste0("Variance explained: ", percentage, '%'),
                 y= 'lognorm counts', x = x_label) +
            theme(axis.title = element_text(size = (7)),
                  axis.text = element_text(size = (6)),
                  plot.title = element_text(hjust=0.5, size=7.5, face="bold"),
                  plot.subtitle = element_text(size = 7, color='gray40'))
    }

    ## Scatterplots for continuous variables
    else {
        colors <- c('mitoRate'='khaki3', 'totalAssignedGene'='plum2', 'overallMapRate'='turquoise', 'concordMapRate'='lightsalmon',
                    'detected_num_genes'='skyblue2', 'library_size'='palegreen3', 'RIN'='rosybrown3', 'Total_RNA_amount'='brown4',
                    'RNA_concentration'='blue3')

        plot <- ggplot(as.data.frame(data), aes(x=eval(parse_expr(sample_var)), y=gene_expr)) +
            geom_point(color=colors[sample_var], show.legend = FALSE, size=2) +
            theme_classic(base_size = 5) +
            labs(title = gene_ids,
                 subtitle = paste0("Variance explained: ", percentage, '%'),
                 y= 'lognorm counts', x = gsub('_', ' ', sample_var)) +
            theme(plot.margin=unit (c (0.4,0.1,0.4,0.1), 'cm'),
                  axis.title = element_text(size = (7)),
                  axis.text = element_text(size = (6)),
                  plot.title = element_text(hjust=0.5, size=7.5, face="bold"),
                  plot.subtitle = element_text(size = 7, color='gray40'))
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
    ggsave(here(paste("plots/04_EDA/03_Explore_gene_level_effects/01_Expl_vars/Top_affected_genes_by_", capitalize(sample_var), "_", brain_region, ".pdf", sep="")), width = 18, height = 13, units = "cm")
}


## Plots
gene_expr_expl_var('habenula', 'Substance')
gene_expr_expl_var('habenula', 'library_size')
gene_expr_expl_var('habenula', 'totalAssignedGene')
gene_expr_expl_var('habenula', 'Batch_RNA_extraction')
gene_expr_expl_var('habenula', 'RNA_concentration')
gene_expr_expl_var('habenula', 'Total_RNA_amount')



