
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
library(ggrepel)
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

## Capitalize brain region info
colData(rse_gene)$Brain_Region <- capitalize(colData(rse_gene)$Brain_Region)

## Create variable for Brain Region + Substance
colData(rse_gene)$Brain_Region_and_Substance <- apply(colData(rse_gene), 1, function(x){ capitalize(paste(x['Brain_Region'], x['Substance'])) })

rse_gene_QC_vars <- rse_gene

## Save original rse with qc variables
save(rse_gene_QC_vars, file = 'processed-data/04_EDA/01_QCA/rse_gene_with_QCvars.Rdata')
## Write csv with complete sample metadata and QC metrics (Supp Table)
sample_metadata_and_QCmetrics <- colData(rse_gene)[,c('Sample_Num', 'Rat_ID', 'Brain_Region', 'Substance', 'Sex', 'SAMPLE_ID',
                                        'Num_Fentanyl_Sessions_six_hrs', 'Total_Num_Fentanyl_Sessions', 'Batch_RNA_extraction',
                                        'Batch_lib_prep', 'Batch_seq', 'RNA_concentration',
                                        'Total_RNA_amount', 'RIN_LIBD', 'Psomagen_RIN', 'RIN',
                                        'library_size', 'detected_num_genes', 'mitoRate',
                                        'totalAssignedGene', 'overallMapRate', 'concordMapRate')]
## "error" in Psomagen_RIN as 'NA'
sample_metadata_and_QCmetrics$Psomagen_RIN <- replace(sample_metadata_and_QCmetrics$Psomagen_RIN, which(sample_metadata_and_QCmetrics$Psomagen_RIN=='error'), NA)
## Order by sample number
sample_metadata_and_QCmetrics <- sample_metadata_and_QCmetrics[order(sample_metadata_and_QCmetrics$Sample_Num),]
write.csv(sample_metadata_and_QCmetrics, file="raw-data/sample_metadata_and_QCmetrics.csv", row.names = FALSE, col.names = TRUE, sep = '\t')

## Return to rse_gene object
rse_gene <- rse_gene_QC_vars





## 1.2 Evaluate QC metrics for groups of samples

## QC metrics of interest
qc_metrics <- c('mitoRate', 'overallMapRate', 'totalAssignedGene', 'concordMapRate', 'library_size', 'detected_num_genes', 'RIN', 'RNA_concentration', 'Total_RNA_amount')

## Sample variables of interest
sample_variables <- c("Brain_Region", "Substance", "Brain_Region_and_Substance", "Num_Fentanyl_Sessions_six_hrs", 'Total_Num_Fentanyl_Sessions', 'Batch_RNA_extraction', 'Batch_lib_prep')


## Function to create boxplots of QC metrics for groups of samples

QC_boxplots <- function(qc_metric, sample_var){

    if (sample_var=="Brain_Region"){
        colors=c('Amygdala'='palegreen3', 'Habenula'='orchid1')
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
        colors=c('Amygdala\nFentanyl'='springgreen3' , 'Amygdala\nSaline'='yellowgreen', 'Habenula\nFentanyl'='hotpink1', 'Habenula\nSaline'='violet')
        violin_width=0.7
        jitter_width=0.09
        x_label="Brain Region & Substance"
        sample_var='Brain_Region_and_Substance_lab'
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
    else if (sample_var=='Batch_RNA_extraction'){
        colors=c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2')
        violin_width=0.7
        jitter_width=0.1
        x_label="RNA extraction Batch"
    }
    else if (sample_var=='Batch_lib_prep'){
        colors=c('1'='darkgoldenrod3', '2'='mediumpurple2', '3'= 'darkmagenta')
        violin_width=0.7
        jitter_width=0.1
        x_label="Library preparation Batch"
    }

    y_label=str_replace_all(qc_metric, c("_"=" "))

    data <- data.frame(colData(rse_gene))
    data$Brain_Region_and_Substance_lab <- paste(as.data.frame(strsplit(data$Brain_Region_and_Substance, ' '))[1,],
                                                 as.data.frame(strsplit(data$Brain_Region_and_Substance, ' '))[2,], sep='\n')

    plot <- ggplot(data = data, mapping = aes(x = !! rlang::sym(sample_var), y = !! rlang::sym(qc_metric), color = !! rlang::sym(sample_var))) +
                geom_violin(alpha = 0, size = 0.4, color='black', width=violin_width)+
                geom_jitter(width = jitter_width, alpha = 0.7, size = 2) +
                geom_boxplot(alpha = 0, size = 0.4, width=0.1, color='black') +
                scale_color_manual(values = colors) +
                sm_hgrid() +
                labs(y= y_label, x = x_label) +
                theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                      axis.title = element_text(size = (12)),
                      axis.text = element_text(size = (11)))

    return(plot)
}


## Multiple plots
for (sample_var in sample_variables){

    if (sample_var=="Brain_Region_and_Substance"){
        width=35.5
        height=32
    }
    else if(sample_var %in% c("Num_Fentanyl_Sessions_six_hrs", "Total_Num_Fentanyl_Sessions")){
        width=32
        height=25.5
    }
    else {
        width=27
        height=25.5
    }

    i=1
    plots = list()
    for (qc_metric in qc_metrics) {
        plots[[i]]<- QC_boxplots(qc_metric, sample_var)
        i=i+1
    }
    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], nrow = 3, align = 'vh')
    ggsave(paste("plots/04_EDA/01_QCA/QC_boxplots_", sample_var,".pdf", sep=""), width=width, height=height, units = "cm")
}





## 1.3 Compare QC metrics of samples

## Correlation between RNA concentration/RNA amount and the rest of QC variables
data <- data.frame(colData(rse_gene))
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
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('Amygdala'=3, 'Habenula'=2)) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='QC metrics') +
        geom_text(x = max(eval(parse_expr(paste0('data$', var1))))-dist, y = max(data$detected_num_genes),
                  label = paste0('r = ', corrs[var1, 'detected_num_genes']),
                  color = 'skyblue2', size=3)

    p2 <- ggplot(data, aes(x=eval(parse_expr(var1)), y=RIN)) +
        geom_point(aes(shape=Brain_Region), color='rosybrown3', show.legend = FALSE) +
        stat_smooth (geom="line", alpha=0.4, size=1.1, span=0.1, method = lm, color='rosybrown3') +
        theme_classic() +
        theme(plot.margin = unit(c(1,2,1,2), "cm")) +
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('Amygdala'=3, 'Habenula'=2)) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='') +
        geom_text(x = max(eval(parse_expr(paste0('data$', var1))))-dist, y = max(data$RIN),
                  label = paste0('r = ', corrs[var1, 'RIN']),
                  color = 'rosybrown3', size=3)

    p3 <- ggplot(data, aes(x=eval(parse_expr(var1)), y=library_size)) +
        geom_point(aes(shape=Brain_Region), color='palegreen3', show.legend = FALSE) +
        stat_smooth (geom="line", alpha=0.4, size=1.1, span=0.1, method = lm, color='palegreen3') +
        theme_classic() +
        theme(plot.margin = unit(c(1,2,1,2), "cm")) +
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('Amygdala'=3, 'Habenula'=2)) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='') +
        geom_text(x = max(eval(parse_expr(paste0('data$', var1))))-dist, y = max(data$library_size),
                  label = paste0('r = ', corrs[var1, 'library_size']),
                  color = 'palegreen3', size=3)


    p4 <- ggplot(QC_data, aes(x=eval(parse_expr(var1)), y=QC_values, color=QC_var_name)) +
        geom_point(aes(shape=Brain_Region)) +
        stat_smooth (geom="line", alpha=0.4, size=1.1, span=0.1, method = lm) +
        theme_classic() +
        theme(plot.margin = unit(c(1,0,1,0), "cm")) +
        labs(x=str_replace_all(var1, c("_"=" ")), y='') +
        scale_color_manual(values = values) +
        scale_shape_manual(labels=c("Amygdala","Habenula"), values=c('Amygdala'=3, 'Habenula'=2)) +
        labs(shape="Brain Region", colour="QC variables") +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$mitoRate)+.03,
                  label = paste0('r = ', corrs[var1, 'mitoRate']),
                  color = 'khaki3', size=3) +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$totalAssignedGene)+.03,
                  label = paste0('r = ', corrs[var1, 'totalAssignedGene']),
                  color = 'plum2', size=3) +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$overallMapRate)-0.05,
                  label = paste0('r = ', corrs[var1, 'overallMapRate']),
                  color = 'turquoise', size=3) +
        geom_text(x = max(eval(parse_expr(paste0('QC_data$', var1, '[1:132]'))))-dist, y = max(data$concordMapRate)+dist_corr_coeff,
                  label = paste0('r = ', corrs[var1, 'concordMapRate']),
                  color = 'lightsalmon', size=3)

    plot_grid(p1, p2, p3, p4, nrow=1)
    ggsave(here(paste("plots/04_EDA/01_QCA/Corr_QCmetrics_vs_", var1,".pdf", sep="")), width = 50, height = 10, units = "cm")


}

corr_plots('RNA_concentration')
corr_plots('Total_RNA_amount')




## Create correlation plots for samples within each brain region

data <- data.frame(colData(rse_gene))

## Correlation between RNA concentration/amount and the QC variables in habenula and amygdala samples
brain_regions <- c('Habenula', 'Amygdala')
RNA_vars <- c('RNA_concentration', 'Total_RNA_amount')
qc_metrics <- c('mitoRate', 'totalAssignedGene', 'overallMapRate', 'concordMapRate',
                'detected_num_genes', 'RIN', 'library_size')
corrs_brain_regions  <- lapply(brain_regions,
                               function(z){sapply(RNA_vars,
                                    function(y){sapply(qc_metrics,
                                         function(x){round(cor(data[which(data$Brain_Region==z),  x], data[which(data$Brain_Region==z), y], method = c("pearson")), 3)}
                                    )}
                               )})
corrs_habeula <- corrs_brain_regions[[1]]
corrs_amygdala <- corrs_brain_regions[[2]]


## Plots

## Function to plot RNA vars vs QC metrics for habenula and amygdala samples separately
corr_plots_brain_region <- function(qc_metric, var1){

    if (qc_metric=='detected_num_genes'){
        color1 = 'skyblue3'
        color2 = 'skyblue'
        if (var1=="RNA_concentration"){
            x_dist1=0
            x_dist2=130
            y_dist1=90
            y_dist2=-260
        }
        else {
            x_dist1=0.24
            x_dist2=0.1
            y_dist1=-370
            y_dist2=-500
        }

    }

    else if (qc_metric=='mitoRate'){
        color1 = 'khaki4'
        color2 = 'khaki2'
        if (var1=="RNA_concentration"){
            x_dist1=0
            x_dist2=100
            y_dist1=0.003
            y_dist2=-0.013
        }
        else{
            x_dist1=0
            x_dist2=0.1
            y_dist1=0.0035
            y_dist2=-0.01
        }
    }

    else if (qc_metric=='totalAssignedGene'){
        color1 = 'plum4'
        color2 = 'plum1'
        if (var1=="RNA_concentration"){
            x_dist1=0
            x_dist2=130
            y_dist1=0.03
            y_dist2=-0.11
        }
        else{
            x_dist1=0.1
            x_dist2=0.1
            y_dist1=-0.11
            y_dist2=-0.03
        }
    }

    else if (qc_metric=='overallMapRate'){
        color1 = 'turquoise4'
        color2 = 'turquoise3'
        if (var1=="RNA_concentration"){
            x_dist1=0
            x_dist2=100
            y_dist1=-0.029
            y_dist2=-0.023
        }
        else{
            x_dist1=0
            x_dist2=0.1
            y_dist1=-0.029
            y_dist2=-0.018
        }
    }

    else if (qc_metric=='concordMapRate'){
        color1 = 'lightsalmon3'
        color2 = 'lightsalmon1'
        if (var1=="RNA_concentration"){
            x_dist1=0
            x_dist2=100
            y_dist1=-0.027
            y_dist2=-0.023
        }
        else{
            x_dist1=0
            x_dist2=0.1
            y_dist1=-0.029
            y_dist2=-0.018
        }
    }

    else if (qc_metric=='RIN'){
        color1 = 'rosybrown4'
        color2 = 'rosybrown2'
        if (var1=="RNA_concentration"){
            x_dist1=0
            x_dist2=100
            y_dist1=-0.1
            y_dist2=-0.25
        }
        else{
            x_dist1=0
            x_dist2=0.1
            y_dist1=-0.1
            y_dist2=-0.34
        }
    }

    else if (qc_metric=='library_size'){
        color1 = 'palegreen4'
        color2 = 'palegreen2'
        if (var1=="RNA_concentration"){
            x_dist1=0
            x_dist2=100
            y_dist1=-5000000
            y_dist2=-4000000
        }
        else{
            x_dist1=0.1
            x_dist2=0.1
            y_dist1=-6000000
            y_dist2=-6000000
        }
    }

    plot <- ggplot(data, aes(x=eval(parse_expr(var1)), y=eval(parse_expr(qc_metric)), color=Brain_Region)) +
                    geom_point(aes(shape=Brain_Region)) +
                    stat_smooth (geom="line", alpha=0.7, size=1.1, span=0.1, method = lm, show.legend = FALSE) +
                    theme_classic() +
                    scale_color_manual(name='Brain Region', values=c('Amygdala'=color1, 'Habenula'=color2)) +
                    scale_shape_manual(name='Brain Region', labels=c("Amygdala","Habenula"), values=c('Amygdala'=3, 'Habenula'=2)) +
                    labs(x=str_replace_all(var1, c("_"=" ")), y=str_replace_all(qc_metric, c("_"=" "))) +
                    geom_text(x = max(data[which(data$Brain_Region=='Habenula'),  var1]) - x_dist1,
                              y = max(data[which(data$Brain_Region=='Habenula'),  qc_metric]) + y_dist1,
                              label = paste0('r = ', corrs_habeula[qc_metric, var1]),
                              color = color2, size=3) +
                    geom_text(x = max(data[which(data$Brain_Region=='Amygdala'),  var1]) - x_dist2,
                              y = max(data[which(data$Brain_Region=='Amygdala'),  qc_metric]) + y_dist2,
                              label = paste0('r = ', corrs_amygdala[qc_metric, var1]),
                              color = color1, size=3)
    return(plot)
}

## Multiple plots
for (var1 in RNA_vars){
    i=1
    plots <- list()
    for (qc_metric in qc_metrics){

        plots[[i]] <- corr_plots_brain_region(qc_metric, var1)
        i=i+1
    }

    plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]])
    ggsave(here(paste("plots/04_EDA/01_QCA/Corr_BrainRegion_QCmetrics_vs_", var1,".pdf", sep="")), width = 40, height = 30, units = "cm")
}





## 1.4 QC-based sample filtering

## Find sample outliers based on their QC metrics, separately for habenula and amygdala samples
rse_gene_habenula <- rse_gene[,which(rse_gene$Brain_Region=="Habenula")]
rse_gene_amygdala <- rse_gene[,which(rse_gene$Brain_Region=="Amygdala")]

## Drop samples with lower library sizes, detected number of genes, RIN numbers, RNA concentration, total RNA amount, concordMapRate,
## overallMapRate and totalAssignedGene
## Drop samples with high mitoRates

## Filter all samples together

outliers_library_size <- isOutlier(rse_gene$library_size, nmads = 3, type="lower")
outliers_detected_num_genes <- isOutlier(rse_gene$detected_num_genes, nmads = 3, type="lower")
outliers_RNA_concentration <- isOutlier(rse_gene$RNA_concentration, nmads = 3, type="lower")
outliers_Total_RNA_amount <- isOutlier(rse_gene$Total_RNA_amount, nmads = 3, type="lower")
outliers_totalAssignedGene <- isOutlier(rse_gene$totalAssignedGene, nmads = 3, type="lower")
outliers_overallMapRate <- isOutlier(rse_gene$overallMapRate, nmads = 3, type="lower")
outliers_concordMapRate <- isOutlier(rse_gene$concordMapRate, nmads = 3, type="lower")
outliers_RIN<-isOutlier(rse_gene$RIN, nmads = 3, type="lower")
outliers_mitoRate<-isOutlier(rse_gene$mitoRate, nmads = 3, type="higher")

## Add QC metrics for which samples had outlier values
colData(rse_gene) <- cbind(colData(rse_gene), sapply(paste0('outliers_', qc_metrics), function(x){eval(parse_expr(x))}))

## Subset to samples without outliers
not_outliers<-which(! (outliers_library_size | outliers_detected_num_genes | outliers_RNA_concentration | outliers_Total_RNA_amount |
                       outliers_totalAssignedGene | outliers_overallMapRate | outliers_concordMapRate | outliers_mitoRate | outliers_RIN))
rse_gene_qc<-rse_gene[,not_outliers]
save(rse_gene_qc, file='processed-data/04_EDA/01_QCA/rse_gene_qc.Rdata')
## Number of samples retained
dim(rse_gene_qc)[2]
# 18



# Filter habenula samples

outliers_library_size <- isOutlier(rse_gene_habenula$library_size, nmads = 3, type="lower")
outliers_detected_num_genes <- isOutlier(rse_gene_habenula$detected_num_genes, nmads = 3, type="lower")
outliers_RNA_concentration <- isOutlier(rse_gene_habenula$RNA_concentration, nmads = 3, type="lower")
outliers_Total_RNA_amount <- isOutlier(rse_gene_habenula$Total_RNA_amount, nmads = 3, type="lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_habenula$totalAssignedGene, nmads = 3, type="lower")
outliers_overallMapRate <- isOutlier(rse_gene_habenula$overallMapRate, nmads = 3, type="lower")
outliers_concordMapRate <- isOutlier(rse_gene_habenula$concordMapRate, nmads = 3, type="lower")
outliers_RIN<-isOutlier(rse_gene_habenula$RIN, nmads = 3, type="lower")
outliers_mitoRate<-isOutlier(rse_gene_habenula$mitoRate, nmads = 3, type="higher")

colData(rse_gene_habenula) <- cbind(colData(rse_gene_habenula), sapply(paste0('outliers_', qc_metrics), function(x){eval(parse_expr(x))}))

not_outliers<-which(! (outliers_library_size | outliers_detected_num_genes | outliers_RNA_concentration | outliers_Total_RNA_amount |
                           outliers_totalAssignedGene | outliers_overallMapRate | outliers_concordMapRate | outliers_mitoRate | outliers_RIN))
rse_gene_habenula_qc<-rse_gene_habenula[,not_outliers]
save(rse_gene_habenula_qc, file='processed-data/04_EDA/01_QCA/rse_gene_habenula_qc.Rdata')
## Number of samples retained
dim(rse_gene_habenula_qc)[2]
# 15



# Filter amygdala samples

outliers_library_size <- isOutlier(rse_gene_amygdala$library_size, nmads = 3, type="lower")
outliers_detected_num_genes <- isOutlier(rse_gene_amygdala$detected_num_genes, nmads = 3, type="lower")
outliers_RNA_concentration <- isOutlier(rse_gene_amygdala$RNA_concentration, nmads = 3, type="lower")
outliers_Total_RNA_amount <- isOutlier(rse_gene_amygdala$Total_RNA_amount, nmads = 3, type="lower")
outliers_totalAssignedGene <- isOutlier(rse_gene_amygdala$totalAssignedGene, nmads = 3, type="lower")
outliers_overallMapRate <- isOutlier(rse_gene_amygdala$overallMapRate, nmads = 3, type="lower")
outliers_concordMapRate <- isOutlier(rse_gene_amygdala$concordMapRate, nmads = 3, type="lower")
outliers_RIN<-isOutlier(rse_gene_amygdala$RIN, nmads = 3, type="lower")
outliers_mitoRate<-isOutlier(rse_gene_amygdala$mitoRate, nmads = 3, type="higher")

colData(rse_gene_amygdala) <- cbind(colData(rse_gene_amygdala), sapply(paste0('outliers_', qc_metrics), function(x){eval(parse_expr(x))}))

not_outliers<-which(! (outliers_library_size | outliers_detected_num_genes | outliers_RNA_concentration | outliers_Total_RNA_amount |
                           outliers_totalAssignedGene | outliers_overallMapRate | outliers_concordMapRate | outliers_mitoRate | outliers_RIN))
rse_gene_amygdala_qc<-rse_gene_amygdala[,not_outliers]
save(rse_gene_amygdala_qc, file='processed-data/04_EDA/01_QCA/rse_gene_amygdala_qc.Rdata')
## Number of samples retained
dim(rse_gene_amygdala_qc)[2]
# 14



## 1.4.1 Boxplots of QC metrics after sample filtering

## Boxplots
boxplots_after_QC_filtering <- function(RSE_obj, qc_metric, sample_var){

    RSE<-eval(parse_expr(RSE_obj))
    colors=c('FALSE'='deepskyblue', 'TRUE'='brown2')

    if (sample_var=="Brain_Region"){
        shapes=c('Amygdala'=3, 'Habenula'=2)
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
    else if (sample_var=='Batch_RNA_extraction'){
        shapes=c('1'=2, '2'=3, '3'=12)
        sample_var_label="RNA extraction Batch"
    }
    else if (sample_var=="Batch_lib_prep"){
        shapes=c('1'=2, '2'=3, '3'=13)
        sample_var_label="Library preparation Batch"
    }

    y_label=str_replace_all(qc_metric, c("_"=" "))
    data <- data.frame(colData(RSE))

    ## Median of the QC var values
    median<-median(eval(parse_expr(paste("RSE$", qc_metric, sep=""))))
    ## Mean-absolute-deviation of the QC var values
    mad<-mad(eval(parse_expr(paste("RSE$", qc_metric, sep=""))))

    ## Label outlier sample (if any)
    if (RSE_obj!='rse_gene'){
        data$outlier_label <- apply(data, 1, function(x){if(x[paste0('outliers_', qc_metric)]==TRUE){x['SAMPLE_ID']}
            else {NA}})
    }
    else{
        data$outlier_label <- rep(NA, dim(rse_gene)[2])
    }


    ## Jitter position
    pos <- position_jitter(width = 0.2, seed = 2)

    plot <- ggplot(data = data, mapping = aes(x = '', y = !! rlang::sym(qc_metric),
                                              color = !! rlang::sym(paste0('outliers_', qc_metric)), label=outlier_label)) +
        geom_jitter(alpha = 1, size = 2.5, position = pos, aes(shape=eval(parse_expr((sample_var))))) +
        geom_boxplot(alpha = 0, size = 0.3, color='black') +
        scale_color_manual(values = colors) +
        scale_shape_manual(values = shapes) +
        labs(x="", y = y_label, color='Outlier QC metric', shape=sample_var_label) +
        sm_hgrid() +
        theme_classic() +
        ## Median line
        geom_hline(yintercept = median, size=0.5) +
        ## Line of median + 3 MADs
        geom_hline(yintercept = median+(3*mad), size=0.5, linetype=2) +
        ## Line of median - 3 MADs
        geom_hline(yintercept = median-(3*mad), size=0.5, linetype=2) +
        ## Labels of outlier samples
        geom_label_repel(color='gray30', size=3.2, max.overlaps = Inf,
                         box.padding = 0.7,
                         show.legend=FALSE,
                         position = pos,
                         min.segment.length = 0) +
        theme(plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm"),
              axis.title = element_text(size = (12)),
              axis.text = element_text(size = (11)),
              legend.position="right",
              legend.text = element_text(size=11),
              legend.title = element_text(size=12))

    return(plot)
}


## Multiple plots
multiple_boxplots <- function(RSE, sample_group){
    for (sample_var in sample_variables){

        if(sample_var %in% c("Total_Num_Fentanyl_Sessions", "Num_Fentanyl_Sessions_six_hrs")){
            height=23
            width=40
        }
        else if(sample_var=="Brain_Region_and_Substance"){
            height=23
            width=38
        }
        else{
            height=23
            width=37
        }

        i=1
        plots = list()
        for (qc_metric in qc_metrics) {
            plots[[i]]<- boxplots_after_QC_filtering(RSE, qc_metric, sample_var)
            i=i+1
        }
        plot_grid(plots[[1]], plots[[2]], plots[[3]],
                  plots[[4]], plots[[5]], plots[[6]],
                  plots[[7]], plots[[8]], plots[[9]],
            ncol = 3, align = 'vh')
        ggsave(paste("plots/04_EDA/01_QCA/Boxplots_afterQC_filtering_", sample_var, "_", sample_group, ".pdf", sep=""), width=width, height=height, units = "cm")
    }
}


## Plots

## All samples together
## Add new variables to rse_gene with info of samples retained/dropped
rse_gene$Retention_after_QC_filtering <- as.vector(sapply(rse_gene$SAMPLE_ID, function(x){if (x %in% rse_gene_qc$SAMPLE_ID){'Retained'} else{'Dropped'}}))
rse_gene$Retention_sample_label <- c(rep(NA, ncol(rse_gene)))
rse_gene_complete <- rse_gene
save(rse_gene_complete, file="processed-data/04_EDA/01_QCA/rse_gene_complete.Rdata")
multiple_boxplots('rse_gene', 'all')

## Habenula samples
rse_gene_habenula$Retention_after_QC_filtering <- as.vector(sapply(rse_gene_habenula$SAMPLE_ID, function(x){if (x %in% rse_gene_habenula_qc$SAMPLE_ID){'Retained'} else{'Dropped'}}))
rse_gene_habenula$Retention_sample_label <- as.vector(sapply(rse_gene_habenula$SAMPLE_ID, function(x){if (x %in% rse_gene_habenula_qc$SAMPLE_ID){NA} else{x}}))
rse_gene_habenula_complete <- rse_gene_habenula
save(rse_gene_habenula_complete, file="processed-data/04_EDA/01_QCA/rse_gene_habenula_complete.Rdata")
multiple_boxplots('rse_gene_habenula', 'habenula')

## Amygdala samples
rse_gene_amygdala$Retention_after_QC_filtering <- as.vector(sapply(rse_gene_amygdala$SAMPLE_ID, function(x){if (x %in% rse_gene_amygdala_qc$SAMPLE_ID){'Retained'} else{'Dropped'}}))
rse_gene_amygdala$Retention_sample_label <- as.vector(sapply(rse_gene_amygdala$SAMPLE_ID, function(x){if (x %in% rse_gene_amygdala_qc$SAMPLE_ID){NA} else{x}}))
rse_gene_amygdala_complete <- rse_gene_amygdala
save(rse_gene_amygdala_complete, file="processed-data/04_EDA/01_QCA/rse_gene_amygdala_complete.Rdata")
multiple_boxplots('rse_gene_amygdala', 'amygdala')







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
# date     2024-03-06
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# ! package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# backports              1.4.1     2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3     2015-07-28 [1] CRAN (R 4.3.0)
# V beachmat               2.18.0    2024-02-17 [1] Bioconductor 3.18 (R 4.3.2) (on disk 2.18.1)
# beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0    2023-10-26 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# BiocNeighbors          1.20.2    2024-01-13 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           1.36.0    2023-10-26 [1] Bioconductor
# BiocSingular           1.18.0    2023-11-06 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# broom                  1.0.5     2023-06-09 [1] CRAN (R 4.3.0)
# car                    3.1-2     2023-03-30 [1] CRAN (R 4.3.0)
# carData                3.0-5     2022-01-06 [1] CRAN (R 4.3.0)
# cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
# checkmate              2.3.1     2023-12-04 [1] CRAN (R 4.3.1)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.2)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# V data.table             1.15.0    2024-02-29 [1] CRAN (R 4.3.1) (on disk 1.15.2)
# DelayedArray           0.28.0    2023-11-06 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-11-06 [1] Bioconductor
# digest                 0.6.34    2024-01-11 [1] CRAN (R 4.3.1)
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# V edgeR                  4.0.15    2024-02-20 [1] Bioconductor 3.18 (R 4.3.2) (on disk 4.0.16)
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
# ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.3.0)
# gghalves               0.1.4     2022-11-20 [1] CRAN (R 4.3.0)
# V ggplot2              * 3.4.4     2024-02-23 [1] CRAN (R 4.3.1) (on disk 3.5.0)
# ggpubr                 0.6.0     2023-02-10 [1] CRAN (R 4.3.0)
# ggrepel              * 0.9.5     2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4     2022-10-13 [1] CRAN (R 4.3.0)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                * 5.1-1     2023-09-12 [1] CRAN (R 4.3.0)
# htmlTable              2.4.2     2023-10-29 [1] CRAN (R 4.3.1)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.1)
# htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.0)
# knitr                  1.45      2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                  3.58.1    2023-11-02 [1] Bioconductor
# V locfit                 1.5-9.8   2024-03-01 [1] CRAN (R 4.3.1) (on disk 1.5.9.9)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-26 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# nnet                   7.3-19    2023-05-03 [1] CRAN (R 4.3.2)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# pwr                    1.3-0     2020-03-17 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.1)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
# rlang                * 1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# V rmarkdown              2.25      2024-03-05 [1] CRAN (R 4.3.1) (on disk 2.26)
# rpart                  4.1.23    2023-12-05 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
# rstatix                0.7.2     2023-02-01 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-26 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# ScaledMatrix           1.10.0    2023-11-06 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# scater               * 1.30.1    2023-11-16 [1] Bioconductor
# scuttle              * 1.12.0    2023-11-06 [1] Bioconductor
# sdamr                  0.2.0     2022-11-16 [1] CRAN (R 4.3.0)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# SingleCellExperiment * 1.24.0    2023-11-06 [1] Bioconductor
# smplot2              * 0.1.0     2023-06-07 [1] Github (smin95/smplot2@836f909)
# SparseArray            1.2.4     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats      1.14.0    2023-10-26 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1     2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0    2023-11-06 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.1     2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.1)
# vipor                  0.4.7     2023-12-18 [1] CRAN (R 4.3.1)
# viridis                0.6.5     2024-01-29 [1] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.42      2024-02-08 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-26 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-26 [1] Bioconductor
# zoo                    1.8-12    2023-04-13 [1] CRAN (R 4.3.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# V ── Loaded and on-disk version mismatch.
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────



