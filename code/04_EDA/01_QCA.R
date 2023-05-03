
library(here)
library(SummarizedExperiment)
library(readxl)
library(stringr)
library(ggplot2)
library(rlang)
library(smplot2)
library(Hmisc)
library(cowplot)
library(scater)
library(sessioninfo)


########################   Exploratory Data Analysis   #########################


################################################################################
##                        1. Quality Control Analysis
################################################################################
## (Note: genes were not filtered and counts were not normalized in these analyses)

## Load RSE object at gene level
load(here('processed-data/03_Data_preparation/rse_gene_sample_info.Rdata'), verbose=TRUE)



## 1.1 QC variable construction and correction

## Add library size for each sample
colData(rse_gene)$library_size <- apply(assay(rse_gene), 2, sum)

## Add detected number of genes (not zero-expressed genes) for each sample
colData(rse_gene)$detected_num_genes <- apply(assay(rse_gene), 2, function(x){length(x[which(x>0)])})

## Create a single and complete RIN column (add LIBD RIN if Psomagen is not available)
colData(rse_gene)$RIN <- as.numeric(apply(colData(rse_gene), 1, function(x){if (x['Psomagen_RIN']!="error"){x['Psomagen_RIN']} else {x['RIN_LIBD']}}))

## Correct data for RNA concentration
colData(rse_gene)$RNA_concentration <- as.vector(sapply(rse_gene$RNA_concentration, function(x){if (strsplit(x, '')[[1]][1]=='*'){strsplit(x, '\\*')[[1]][2]} else {x}}))
colData(rse_gene)$RNA_concentration <- as.numeric(colData(rse_gene)$RNA_concentration)

## Change numeric data to char
colData(rse_gene)$Total_Num_Fentanyl_Sessions <- as.character(colData(rse_gene)$Total_Num_Fentanyl_Sessions)
colData(rse_gene)$Num_Fentanyl_Sessions_six_hrs <- as.character(colData(rse_gene)$Num_Fentanyl_Sessions_six_hrs)





## 1.2 Evaluate QC metrics for groups of samples

## QC metrics of interest
qc_metrics <- c('mitoRate', 'overallMapRate', 'totalAssignedGene', 'concordMapRate', 'library_size', 'detected_num_genes', 'RIN', 'RNA_concentration', 'Total_RNA_amount')

## Create variable for Brain Region + Substance
colData(rse_gene)$Brain_Region_and_Substance <- apply(colData(rse_gene), 1, function(x){ capitalize(paste(x['Brain_Region'], x['Substance'])) })
save(rse_gene, file = 'processed-data/04_EDA/01_QCA/rse_gene_with_QCvars.Rdata')
## Sample variables of interest
sample_variables <- c("Brain_Region", "Substance", "Brain_Region_and_Substance", "Num_Fentanyl_Sessions_six_hrs", 'Total_Num_Fentanyl_Sessions')


## Function to create boxplots of QC metrics for groups of samples

QC_boxplots <- function(qc_metric, sample_var){

    if (sample_var=="Brain_Region"){
        colors=c('amygdala'='palegreen3', 'habenula'='orchid1')
        violin_width=1
        jitter_width=0.1
        x_label="Brain Region"
    }
    else if (sample_var=="Substance"){
        colors=c('Fentanyl'='turquoise3', 'Saline'='yellow3')
        violin_width=0.7
        jitter_width=0.1
        x_label="Substance"
    }
    else if (sample_var=="Brain_Region_and_Substance"){
        colors=c('Amygdala Fentanyl'='springgreen3' , 'Amygdala Saline'='yellowgreen', 'Habenula Fentanyl'='hotpink1', 'Habenula Saline'='violet')
        violin_width=0.7
        jitter_width=0.09
        x_label="Brain Region & Substance"
    }
    else if (sample_var=='Total_Num_Fentanyl_Sessions'){
        colors=c('24'='salmon', '22'='pink2')
        violin_width=0.7
        jitter_width=0.1
        x_label="Total Number of Fentanyl Sessions"
    }
    else if (sample_var=='Num_Fentanyl_Sessions_six_hrs'){
        colors=c('18'='dodgerblue', '16'='lightskyblue')
        violin_width=0.7
        jitter_width=0.1
        x_label="Number of 6hrs Fentanyl Sessions"
    }

    y_label=str_replace_all(qc_metric, c("_"=" "))

    data <- data.frame(colData(rse_gene))
    plot <- ggplot(data = data, mapping = aes(x = !! rlang::sym(sample_var), y = !! rlang::sym(qc_metric), color = !! rlang::sym(sample_var))) +
                geom_violin(alpha = 0, size = 0.4, color='black', width=violin_width)+
                geom_jitter(width = jitter_width, alpha = 0.7, size = 2) +
                geom_boxplot(alpha = 0, size = 0.4, width=0.1, color='black') +
                scale_color_manual(values = colors) +
                sm_hgrid() +
                labs(y= y_label, x = x_label) +
                theme(axis.title = element_text(size = (9)),
                      axis.text = element_text(size = (8)))

    return(plot)
}


## Multiple plots
for (sample_var in sample_variables){

    if (sample_var=="Brain_Region_and_Substance"){
        width=45
        height=35
    }
    else {
        width=35
        height=30
    }

    i=1
    plots = list()
    for (qc_metric in qc_metrics) {
        plots[[i]]<- QC_boxplots(qc_metric, sample_var)
        i=i+1
    }
    combine_plots(
        list(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]]),
        plotgrid.args = list(nrow = 3)
    )
    ggsave(paste("plots/04_EDA/01_QCA/QC_boxplots_", sample_var,".pdf", sep=""), width=width, height=height, units = "cm")
}





## 1.3 Compare QC metrics of samples

## Correlation betweern RNA concentration/RNA amount and the rest of QC variables

corrs <- as.data.frame(t(round(cor(data[c('mitoRate', 'totalAssignedGene', 'overallMapRate', 'concordMapRate',
                          'detected_num_genes', 'RIN', 'library_size')], data[c('RNA_concentration' ,'Total_RNA_amount')],
                   method = c("pearson")), 3)))

## Create variable with QC metrics
QC_data <- data.frame('QC_values'=c(unlist(data$mitoRate), unlist(data$overallMapRate), unlist(data$totalAssignedGene),
                                    unlist(data$concordMapRate), NA, NA, NA),
                      'QC_var_name'=c(rep('mitoRate', 33), rep('overallMapRate', 33),  rep('totalAssignedGene', 33),
                                      rep('concordMapRate', 33), 'detected_num_genes', 'RIN', 'library_size'),
                      'Brain_Region'=c(rep(data$Brain_Region, 4), NA, NA, NA),
                      'RNA_concentration'=c(rep(data$RNA_concentration, 4), NA, NA, NA),
                      'Total_RNA_amount'=c(rep(data$Total_RNA_amount, 4), NA, NA, NA))


## Plot correlations between RNA quantities and QC metrics of the samples
corr_plots <- function(var1){
    if (var1=='RNA_concentration'){
        dist=150
        dist_corr_coeff=0.07
    }

    else {
        dist=0.1
        dist_corr_coeff=0.12
    }

    values = c('mitoRate'='khaki3', 'totalAssignedGene'='plum2', 'overallMapRate'='turquoise', 'concordMapRate'='lightsalmon',
               'detected_num_genes'='skyblue2', 'library_size'='palegreen3', 'RIN'='rosybrown3')

    p1 <- ggplot(data, aes(x=eval(parse_expr(var1)), y=detected_num_genes)) +
        geom_point(aes(shape=Brain_Region), color='skyblue2', show.legend = FALSE) +
        stat_smooth (geom="line", alpha=0.4, size=1.1, span=0.1, method = lm, color='skyblue2') +
        theme_classic() +
        theme(plot.margin = unit(c(1,2,1,2), "cm")) +
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('amygdala'=3, 'habenula'=2)) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='QC metrics') +
        geom_text(x = max(eval(parse_expr(paste0('data$', var1))))-dist, y = max(data$detected_num_genes),
                  label = paste0('r=', corrs[var1, 'detected_num_genes']),
                  color = 'skyblue2', size=3)

    p2 <- ggplot(data, aes(x=eval(parse_expr(var1)), y=RIN)) +
        geom_point(aes(shape=Brain_Region), color='rosybrown3', show.legend = FALSE) +
        stat_smooth (geom="line", alpha=0.4, size=1.1, span=0.1, method = lm, color='rosybrown3') +
        theme_classic() +
        theme(plot.margin = unit(c(1,2,1,2), "cm")) +
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('amygdala'=3, 'habenula'=2)) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='') +
        geom_text(x = max(eval(parse_expr(paste0('data$', var1))))-dist, y = max(data$RIN),
                  label = paste0('r=', corrs[var1, 'RIN']),
                  color = 'rosybrown3', size=3)

    p3 <- ggplot(data, aes(x=eval(parse_expr(var1)), y=library_size)) +
        geom_point(aes(shape=Brain_Region), color='palegreen3', show.legend = FALSE) +
        stat_smooth (geom="line", alpha=0.4, size=1.1, span=0.1, method = lm, color='palegreen3') +
        theme_classic() +
        theme(plot.margin = unit(c(1,2,1,2), "cm")) +
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('amygdala'=3, 'habenula'=2)) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='') +
        geom_text(x = max(eval(parse_expr(paste0('data$', var1))))-dist, y = max(data$library_size),
                  label = paste0('r=', corrs[var1, 'library_size']),
                  color = 'palegreen3', size=3)


    p4 <- ggplot(QC_data, aes(x=eval(parse_expr(var1)), y=QC_values, color=QC_var_name)) +
        geom_point(aes(shape=Brain_Region)) +
        stat_smooth (geom="line", alpha=0.4, size=1.1, span=0.1, method = lm) +
        theme_classic() +
        theme(plot.margin = unit(c(1,0,1,0), "cm")) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='') +
        scale_color_manual(values = values) +
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('amygdala'=3, 'habenula'=2)) +
        labs(shape="Brain Region", colour="QC variables") +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$mitoRate)+.03,
                  label = paste0('r=', corrs[var1, 'mitoRate']),
                  color = 'khaki3', size=3) +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$totalAssignedGene)+.03,
                  label = paste0('r=', corrs[var1, 'totalAssignedGene']),
                  color = 'plum2', size=3) +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$overallMapRate)-0.05,
                  label = paste0('r=', corrs[var1, 'overallMapRate']),
                  color = 'turquoise', size=3) +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$concordMapRate)+dist_corr_coeff,
                  label = paste0('r=', corrs[var1, 'concordMapRate']),
                  color = 'lightsalmon', size=3)

    plot_grid(p1, p2, p3, p4, nrow=1)
    ggsave(here(paste("plots/04_EDA/01_QCA/Corr_QCmetrics_vs_", var1,".pdf", sep="")), width = 50, height = 10, units = "cm")


}

corr_plots('RNA_concentration')
corr_plots('Total_RNA_amount')





## 1.4 QC sample filetring

## Find sample outliers based on their QC metrics, separately for habenula and amygdala samples
rse_gene_habenula <- rse_gene[,which(rse_gene$Brain_Region=="habenula")]
rse_gene_amygdala <- rse_gene[,which(rse_gene$Brain_Region=="amygdala")]

## Drop samples with lower library sizes, detected number of genes, RNA concentration, total RNA amount, concordMapRate,
## overallMapRate and totalAssignedGene
## Drop samples with high mitoRates and RIN numbers

## Filter all samples together

outliers_library_size <- isOutlier(rse_gene$library_size, nmads = 3, type="lower")
outliers_detected_num <- isOutlier(rse_gene$detected_num_genes, nmads = 3, type="lower")
outliers_RNA_conc <- isOutlier(rse_gene$RNA_concentration, nmads = 3, type="lower")
outliers_RNA_amount <- isOutlier(rse_gene$Total_RNA_amount, nmads = 3, type="lower")
outliers_totalAssignedGene <- isOutlier(rse_gene$totalAssignedGene, nmads = 3, type="lower")
outliers_overallMapRate <- isOutlier(rse_gene$overallMapRate, nmads = 3, type="lower")
outliers_concordMapRate <- isOutlier(rse_gene$concordMapRate, nmads = 3, type="lower")
outliers_mito<-isOutlier(rse_gene$mitoRate, nmads = 3, type="higher")
outliers_RIN<-isOutlier(rse_gene$RIN, nmads = 3, type="higher")

not_outliers<-which(! (outliers_library_size | outliers_detected_num | outliers_RNA_conc | outliers_RNA_amount |
                       outliers_totalAssignedGene | outliers_overallMapRate | outliers_concordMapRate | outliers_mito | outliers_RIN))
rse_gene_qc<-rse_gene[,not_outliers]

## Number of samples retained
dim(rse_gene_qc)[2]
# 18

## Save data
save(rse_gene_qc, file = 'processed-data/04_EDA/01_QCA/rse_gene_qc.Rdata')



# Filter habenula samples

outliers_library_size <- isOutlier(rse_gene_habenula$library_size, nmads = 3, type="lower")
outliers_detected_num <- isOutlier(rse_gene_habenula$detected_num_genes, nmads = 3, type="lower")
outliers_RNA_conc <- isOutlier(rse_gene_habenula$RNA_concentration, nmads = 3, type="lower")
outliers_RNA_amount <- isOutlier(rse_gene_habenula$Total_RNA_amount, nmads = 3, type="lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_habenula$totalAssignedGene, nmads = 3, type="lower")
outliers_overallMapRate <- isOutlier(rse_gene_habenula$overallMapRate, nmads = 3, type="lower")
outliers_concordMapRate <- isOutlier(rse_gene_habenula$concordMapRate, nmads = 3, type="lower")
outliers_mito<-isOutlier(rse_gene_habenula$mitoRate, nmads = 3, type="higher")
outliers_RIN<-isOutlier(rse_gene_habenula$RIN, nmads = 3, type="higher")

not_outliers<-which(! (outliers_library_size | outliers_detected_num | outliers_RNA_conc | outliers_RNA_amount |
                           outliers_totalAssignedGene | outliers_overallMapRate | outliers_concordMapRate | outliers_mito | outliers_RIN))
rse_gene_habenula_qc<-rse_gene_habenula[,not_outliers]

## Number of samples retained
dim(rse_gene_habenula_qc)[2]
# 15

## Save data
save(rse_gene_habenula_qc, file = 'processed-data/04_EDA/01_QCA/rse_gene_habenula_qc.Rdata')



# Filter amygdala samples

outliers_library_size <- isOutlier(rse_gene_amygdala$library_size, nmads = 3, type="lower")
outliers_detected_num <- isOutlier(rse_gene_amygdala$detected_num_genes, nmads = 3, type="lower")
outliers_RNA_conc <- isOutlier(rse_gene_amygdala$RNA_concentration, nmads = 3, type="lower")
outliers_RNA_amount <- isOutlier(rse_gene_amygdala$Total_RNA_amount, nmads = 3, type="lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_amygdala$totalAssignedGene, nmads = 3, type="lower")
outliers_overallMapRate <- isOutlier(rse_gene_amygdala$overallMapRate, nmads = 3, type="lower")
outliers_concordMapRate <- isOutlier(rse_gene_amygdala$concordMapRate, nmads = 3, type="lower")
outliers_mito<-isOutlier(rse_gene_amygdala$mitoRate, nmads = 3, type="higher")
outliers_RIN<-isOutlier(rse_gene_amygdala$RIN, nmads = 3, type="higher")

not_outliers<-which(! (outliers_library_size | outliers_detected_num | outliers_RNA_conc | outliers_RNA_amount |
                           outliers_totalAssignedGene | outliers_overallMapRate | outliers_concordMapRate | outliers_mito | outliers_RIN))
rse_gene_amygdala_qc<-rse_gene_amygdala[,not_outliers]

## Number of samples retained
dim(rse_gene_amygdala_qc)[2]
# 14

## Save data
save(rse_gene_amygdala_qc, file = 'processed-data/04_EDA/01_QCA/rse_gene_amygdala_qc.Rdata')




## 1.4.1 Boxplots of QC metrics after sample filtering

## Boxplots
boxplots_after_QC_filtering <- function(RSE, qc_metric, sample_var){

    rse_gene<-eval(parse_expr(RSE))
    colors=c('Retained'='deepskyblue', 'Dropped'='brown2')

    if (sample_var=="Brain_Region"){
        shapes=c('amygdala'=3, 'habenula'=2)
        sample_var_label="Brain Region"
    }
    else if (sample_var=="Substance"){
        shapes=c('Fentanyl'=9, 'Saline'=5)
        sample_var_label="Substance"
    }
    else if (sample_var=="Brain_Region_and_Substance"){
        shapes=c('Amygdala Fentanyl'=15 , 'Amygdala Saline'=0, 'Habenula Fentanyl'=16, 'Habenula Saline'=1)
        sample_var_label="Brain Region & Substance"
    }
    else if (sample_var=='Total_Num_Fentanyl_Sessions'){
        shapes=c('24'=8, '22'=1)
        sample_var_label="Total Number of Fentanyl Sessions"
    }
    else if (sample_var=='Num_Fentanyl_Sessions_six_hrs'){
        shapes=c('18'=8, '16'=1)
        sample_var_label="Number of 6hrs Fentanyl Sessions"
    }

    y_label=str_replace_all(qc_metric, c("_"=" "))
    data <- data.frame(colData(rse_gene))

    ## Median of the QC var values
    median<-median(eval(parse_expr(paste("rse_gene$", qc_metric, sep=""))))
    ## Mean-absolute-deviation of the QC var values
    mad<-mad(eval(parse_expr(paste("rse_gene$", qc_metric, sep=""))))

    plot <- ggplot(data = data, mapping = aes(x = '', y = !! rlang::sym(qc_metric), color = !! rlang::sym('Retention_after_QC_filtering'))) +
        geom_jitter(width = 0.2, alpha = 1, size = 2, aes(shape=eval(parse_expr((sample_var))))) +
        geom_boxplot(alpha = 0, size = 0.3, color='black') +
        scale_color_manual(values = colors) +
        scale_shape_manual(values=shapes) +
        labs(x="", y = y_label, color='Retention after QC filtering', shape=sample_var_label) +
        sm_hgrid() +
        theme_classic() +
        ## Median line
        geom_hline(yintercept = median, size=0.5) +
        ## Line of median + 3 MADs
        geom_hline(yintercept = median+(3*mad), size=0.5, linetype=2) +
        ## Line of median - 3 MADs
        geom_hline(yintercept = median-(3*mad), size=0.5, linetype=2) +
        theme(axis.title = element_text(size = (9)),
              axis.text = element_text(size = (8)),
              legend.position="right",
              legend.text = element_text(size=8),
              legend.title = element_text(size=9))

    return(plot)
}


## Multiple plots
multiple_boxplots <- function(RSE, sample_group){
    for (sample_var in sample_variables){

        i=1
        plots = list()
        for (qc_metric in qc_metrics) {
            plots[[i]]<- boxplots_after_QC_filtering(RSE, qc_metric, sample_var)
            i=i+1
        }
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]],
            nrow = 3)
        ggsave(paste("plots/04_EDA/01_QCA/Boxplots_afterQC_filtering_", sample_var, "_", sample_group, ".pdf", sep=""), width=35, height=30, units = "cm")
    }
}


## Plots

## All samples together
## Add new variable to rse_gene with info of samples retained/dropped
rse_gene$Retention_after_QC_filtering <- as.vector(sapply(rse_gene$SAMPLE_ID, function(x){if (x %in% rse_gene_qc$SAMPLE_ID){'Retained'} else{'Dropped'}}))
multiple_boxplots('rse_gene', 'all')

## Habenula samples
rse_gene_habenula$Retention_after_QC_filtering <- as.vector(sapply(rse_gene_habenula$SAMPLE_ID, function(x){if (x %in% rse_gene_habenula_qc$SAMPLE_ID){'Retained'} else{'Dropped'}}))
multiple_boxplots('rse_gene_habenula', 'habenula')

## Amygdala samples
rse_gene_amygdala$Retention_after_QC_filtering <- as.vector(sapply(rse_gene_amygdala$SAMPLE_ID, function(x){if (x %in% rse_gene_amygdala_qc$SAMPLE_ID){'Retained'} else{'Dropped'}}))
multiple_boxplots('rse_gene_amygdala', 'amygdala')







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
# date     2023-04-20
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA



