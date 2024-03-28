
library(here)
library(SummarizedExperiment)
library(jaffelab)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(cowplot)
library(rlang)
library(stringr)
library(smplot2)
library(sessioninfo)



################################################################################
##                        2. Explore sample-level effects
################################################################################
## (Note: genes were filtered and counts normalized in this analysis)

load(here('processed-data/03_Data_preparation/rse_gene_filt.Rdata'), verbose=TRUE)
load(here('processed-data/04_EDA/01_QCA/rse_gene_habenula_complete.Rdata'), verbose = TRUE)
load(here('processed-data/04_EDA/01_QCA/rse_gene_amygdala_complete.Rdata'), verbose = TRUE)

## Habenula filtered data
rse_gene_habenula_filt <- rse_gene_filt[,colData(rse_gene_filt)$Brain_Region=='Habenula']
## Add info for outlier samples after QC sample filtering of habenula samples
rse_gene_habenula_filt$Retention_after_QC_filtering <- rse_gene_habenula_complete$Retention_after_QC_filtering
rse_gene_habenula_filt$Retention_sample_label <- rse_gene_habenula_complete$Retention_sample_label

## Amygdala filtered data
rse_gene_amygdala_filt <- rse_gene_filt[,colData(rse_gene_filt)$Brain_Region=='Amygdala']
## Add info for outlier samples after QC sample filtering of amygdala samples
rse_gene_amygdala_filt$Retention_after_QC_filtering <- rse_gene_amygdala_complete$Retention_after_QC_filtering
rse_gene_amygdala_filt$Retention_sample_label <- rse_gene_amygdala_complete$Retention_sample_label



## Principal Component Analysis: explore samples' gene expression variation

## 2.1 PCA data obtention

## Function to obtain PCs

PCA <- function(brain_region){
    if (is.null(brain_region)){
        RSE <- rse_gene_filt
    }
    else{
        RSE<-eval(parse_expr(paste("rse_gene", brain_region, 'filt', sep="_")))
    }
    pca<-prcomp(t(assays(RSE)$logcounts))
    ## % of the variance explained by each PC
    pca_vars<- getPcaVars(pca)
    pca_vars_labs<- paste0(
        "PC", seq(along = pca_vars), ": ",
        pca_vars, "% Var Expl")

    ## Join PCs and sample data
    pca_data<-cbind(pca$x,colData(RSE))
    pca_data<-as.data.frame(pca_data)
    return(list(pca_data, pca_vars_labs))
}





## 2.2 PCA visualization

## Define colors per sample variable
colors <- list('Brain_Region'=c('Amygdala'='palegreen3', 'Habenula'='orchid1'),
               'Substance'=c('Fentanyl'='turquoise3', 'Saline'='yellow3'),
               'Brain_Region_and_Substance'=c('Amygdala Fentanyl'='springgreen3' , 'Amygdala Saline'='yellowgreen',
                                              'Habenula Fentanyl'='hotpink1', 'Habenula Saline'='violet'),
               'Total_Num_Fentanyl_Sessions'=c('24'='salmon', '22'='pink2'),
               'Num_Fentanyl_Sessions_six_hrs'=c('18'='dodgerblue', '16'='lightskyblue'),
               'Batch_RNA_extraction'= c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2'),
               'Batch_lib_prep'=c('1'='darkgoldenrod3', '2'='mediumpurple2', '3'= 'darkmagenta'))


## Colors to highlight outlier and segregated (standout) samples
rare_and_poorQC_samples_colors_hab <- c("5_F_LHb_13"="magenta", "3_F_LHb_09"='blue',
                                        "1_F_LHb_01"='darkorange',  "18_S_LHb_20"='darkgreen')
rare_and_poorQC_samples_colors_amyg <- c("33_S_Amyg_20"="orange1", "10_S_Amyg_06"="royalblue3",
                                         "34_S_Amyg_22"="springgreen4", "14_S_Amyg_14"='orchid1')

## PCx vs PCy plots

PCx_vs_PCy <- function (PCx, PCy, pca_data, pca_vars_labs, sample_var, brain_region, sample_labels, sample_shape) {

    ## Sample shapes
    if(sample_shape=='RNAbatch'){
        shape_name = 'Batch RNA extraction'
        shapes <- c('1'=8, '2'=10, '3'=15)
        shape='Batch_RNA_extraction'
    }
    else{
        shape_name = NA
        shapes <- 1
        shape='NULL'
    }

    ## Colors of squares around samples to highlight
    if (is.null(brain_region)){
        sample_square_colors <- c(NA, NA, NA, NA)
    }
    else if(brain_region=='habenula'){
        sample_square_colors <- rare_and_poorQC_samples_colors_hab
    }
    else if(brain_region=='amygdala'){
        sample_square_colors <- rare_and_poorQC_samples_colors_amyg
    }

    ## Sample labels and colors (outliers only or + rare ones)
    if (sample_labels=='outliers_only'){
        pca_data$sample_labels <- pca_data$Retention_sample_label
        ## Outliers in red labels
        pca_data$sample_label_colors <- sapply(pca_data$Retention_sample_label, function(x){if(!is.na(x)){'red'}
                                                                                            else {NA}})
    }
    else{
        pca_data$sample_labels <- pca_data$outlier_or_rare_samples_labels
        ## Red labels for outliers, gray for segregated ones
        pca_data$sample_label_colors <- sapply(pca_data$outlier_or_rare_samples_labels,
                                               function(x){if(!is.na(x) & x %in% pca_data$Retention_sample_label){'red'}
                                                           else if (!is.na(x)) {'gray35'}
                                                           else {NA}})
    }

    ## Plot margin
    if(sample_var=='Brain_Region'){margin=0.7}
    else if(sample_var=='Substance'){
        if (is.null(brain_region)){margin=0.7}
        else {margin=0.9}
    }
    else if(sample_var=='Brain_Region_and_Substance'){margin=0.1}
    else if(sample_var=='Total_Num_Fentanyl_Sessions'){margin=0.1}
    else if(sample_var=='Num_Fentanyl_Sessions_six_hrs'){margin=0.1}
    else if(sample_var=='Batch_RNA_extraction'){
        if (is.null(brain_region)){margin=0.7}
        else {margin=0.9}
    }
    else if(sample_var=='Batch_lib_prep'){
        if (is.null(brain_region)){margin=0.7}
        else {margin=0.9}
    }

    plot <- ggplot(data=pca_data,
                  aes(x=eval(parse_expr(PCx)),y=eval(parse_expr(PCy)),
                  color=eval(parse_expr(sample_var)),
                  label=sample_labels)) +
            theme_classic() +
            # Add squares around outlier and rare samples
            geom_point(data=subset(pca_data, !is.na(sample_labels)), aes(color=sample_labels),
                       pch = 0, size=5.6, stroke = 1) +
            scale_color_manual(values = sample_square_colors) +
            guides(color='none') +
            new_scale_color() +
            geom_point(aes(color=eval(parse_expr(sample_var)), shape=eval(parse_expr(shape))), size=3.5) +
            scale_color_manual(values = colors[[sample_var]]) +
            scale_shape_manual(name = shape_name, values=shapes) +
            ## Labels of outlier and rare samples
            geom_label_repel(label=pca_data$sample_labels, color=pca_data$sample_label_colors,
                             size=3.2,
                             max.overlaps = Inf,
                             box.padding = 0.7,
                             show.legend=FALSE,
                             min.segment.length = 0) +
            labs(x= pca_vars_labs[strtoi(gsub("PC","", PCx))], y = pca_vars_labs[strtoi(gsub("PC","", PCy))],
                 color=str_replace_all(sample_var, c("_"=" "))) +
            theme(legend.position="right",
                  plot.margin=unit (c (1,margin,1,margin), 'cm'),
                  axis.title = element_text(size = (12)),
                  axis.text = element_text(size = (11)),
                  legend.text = element_text(size=11),
                  legend.title = element_text(size=12))

    return(plot)
}


## All PCA plots

plot_PCAs<-function(brain_region, sample_labels, sample_shape){

    ## PC data
    pca<-PCA(brain_region)
    pca_data <- pca[[1]]
    pca_vars_labs<-pca[[2]]

    ## Plots for ALL samples
    if (is.null(brain_region)){
        for (PCs in list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"), c("PC3", "PC4"), c("PC5", "PC6"))){
            plots<-list()
            i=1
            for (sample_var in c("Brain_Region", "Substance", "Brain_Region_and_Substance",
                                 "Batch_RNA_extraction", "Batch_lib_prep",
                                 "Num_Fentanyl_Sessions_six_hrs", "Total_Num_Fentanyl_Sessions")){
                p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars_labs, sample_var, brain_region, sample_labels, sample_shape)
                plots[[i]]=p
                i=i+1
            }
            plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow = 2, align = 'vh')
            ## Save plots
            ggsave(paste("plots/04_EDA/02_PCA/",PCs[1],"_vs_",PCs[2],"_all_samples_", sample_shape, ".pdf", sep=""),
                   width = 60, height = 19, units = "cm")
        }
    }
    ## For habenula and amygdala separately
    else {
        for (PCs in list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"), c("PC3", "PC4"), c("PC5", "PC6"))){
            plots<-list()
            i=1
            for (sample_var in c("Substance", "Batch_RNA_extraction", "Batch_lib_prep",
                                 "Total_Num_Fentanyl_Sessions", "Num_Fentanyl_Sessions_six_hrs")){
                p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars_labs, sample_var, brain_region, sample_labels, sample_shape)
                plots[[i]]=p
                i=i+1
            }
            plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], nrow = 2, align = 'vh')

            ## Save plots
            ggsave(paste("plots/04_EDA/02_PCA/",PCs[1],"_vs_",PCs[2],"_", brain_region, "_", sample_labels,
                          "_", sample_shape,".pdf", sep=""),
                   width = 47, height = 19, units = "cm")

        }
    }
}


## Plots
plot_PCAs(NULL, 'none', 'RNAbatch')
plot_PCAs(NULL, 'none', 'no_shape')

plot_PCAs('habenula', 'outliers_only', 'RNAbatch')
plot_PCAs('habenula', 'outliers_only', 'no_shape')
plot_PCAs('habenula', 'outliers_plus_rare', 'RNAbatch')
plot_PCAs('habenula', 'outliers_plus_rare', 'no_shape')

plot_PCAs('amygdala', 'outliers_only', 'RNAbatch')
plot_PCAs('amygdala', 'outliers_only', 'no_shape')
plot_PCAs('amygdala', 'outliers_plus_rare', 'RNAbatch')
plot_PCAs('amygdala', 'outliers_plus_rare', 'no_shape')





## 2.3 Manual sample filtering

## 2.3.1 Identify rare samples in PCA plots

#####################
##  Amygdala plots
#####################
amyg_pca_data <- PCA('amygdala')[[1]]

#############  A) Saline samples within fentanyl samples' group:  ##############

amyg_sal_pca_data <- amyg_pca_data[which(amyg_pca_data$Substance=='Saline'),]

## -> One is the "34_S_Amyg_22" outlier sample: saline sample with the second lowest PC2 value
amyg_sal_pca_data[order(amyg_sal_pca_data$PC2), 'SAMPLE_ID'][2]
#  "34_S_Amyg_22"

## -> The other is the "14_S_Amyg_14" sample: saline sample with the lowest value in PC2
amyg_sal_pca_data[order(amyg_sal_pca_data$PC2), 'SAMPLE_ID'][1]
#  "14_S_Amyg_14"


##################  B) Sample with the highest value in PC6:  ##################

## -> It is the "34_S_Amyg_22" outlier sample
amyg_pca_data[which.max(amyg_pca_data$PC6), 'SAMPLE_ID']
#  "34_S_Amyg_22"


## All rare amygdala samples
rare_amyg_samples <- c("34_S_Amyg_22", "14_S_Amyg_14")

## Add sample ID label for rare/outlier samples
rse_gene_amygdala_filt$outlier_or_rare_samples_labels <- sapply(rse_gene_amygdala_filt$SAMPLE_ID,
                                                                function(x){if(x %in% rse_gene_amygdala_filt$Retention_sample_label | x %in% rare_amyg_samples){x}
                                                                    else {NA}})
save(rse_gene_amygdala_filt, file='processed-data/04_EDA/02_PCA/rse_gene_amygdala_filt.Rdata')


#####################
## Habenula samples
#####################
hab_pca_data <- PCA('habenula')[[1]]

################  A) Samples that have the lowest PC2 values:  #################

## -> One is the "5_F_LHb_13" outlier sample
hab_pca_data[order(hab_pca_data$PC2), 'SAMPLE_ID'][1]
#  "5_F_LHb_13"

## -> The second is the "3_F_LHb_09" sample
hab_pca_data[order(hab_pca_data$PC2), 'SAMPLE_ID'][2]
#  "3_F_LHb_09"


########  B) Fentanyl samples that appear within saline samples' group:  ########

## -> One is the "5_F_LHb_13" outlier sample: sample with the lowest PC4 value
hab_pca_data[which.min(hab_pca_data$PC4), 'SAMPLE_ID']
#  "5_F_LHb_13"

hab_fenta_pca_data <- hab_pca_data[which(hab_pca_data$Substance=='Fentanyl'),]
## -> Another is the "1_F_LHb_01" sample: fentanyl sample with the third highest PC3 value
hab_fenta_pca_data[order(hab_fenta_pca_data$PC3, decreasing = TRUE), 'SAMPLE_ID'][3]
#  "1_F_LHb_01"


##################  C) Saline samples within fentanyl group:  ##################

hab_sal_pca_data <- hab_pca_data[which(hab_pca_data$Substance=="Saline"),]
## -> It is the "18_S_LHb_20" sample: the saline sample with the highest PC4 value
hab_sal_pca_data[which.max(hab_sal_pca_data$PC4), 'SAMPLE_ID']
#  "18_S_LHb_20"


## All rare habenula samples
rare_hab_samples <- c("5_F_LHb_13", "3_F_LHb_09", "1_F_LHb_01", "18_S_LHb_20")

## Add sample ID label for rare/outlier samples
rse_gene_habenula_filt$outlier_or_rare_samples_labels <- sapply(rse_gene_habenula_filt$SAMPLE_ID,
                                                                function(x){if(x %in% rse_gene_habenula_filt$Retention_sample_label | x %in% rare_hab_samples){x}
                                                                            else {NA}})
save(rse_gene_habenula_filt, file='processed-data/04_EDA/02_PCA/rse_gene_habenula_filt.Rdata')

## Create same variables for rse_gene_filt
rse_gene_filt$outlier_or_rare_samples_labels <- rep(NA, dim(rse_gene_filt)[2])




## 2.3.2 Explore QC metrics of rare samples

## Plot all QC metrics for a rare/outlier sample

qc_metrics <- c('mitoRate', 'overallMapRate', 'totalAssignedGene', 'concordMapRate', 'library_size', 'detected_num_genes', 'RIN', 'RNA_concentration', 'Total_RNA_amount')

rare_and_poorQC_samples_colors <- c(rare_and_poorQC_samples_colors_hab, rare_and_poorQC_samples_colors_amyg)

QC_boxplot <- function(qc_metric, sample_var, brain_region, rare_sample_ID, sample_shape){

    ## Brain region data
    data <- as.data.frame(colData(eval(parse_expr(paste('rse_gene', brain_region, 'filt', sep='_')))))

    ## Sample shapes
    if(sample_shape=='RNAbatch'){
        shape_name = 'Batch RNA extraction'
        shapes <- c('1'=8, '2'=10, '3'=15)
        shape='Batch_RNA_extraction'
    }
    else{
        shape_name = NA
        shapes <- 1
        shape='NULL'
    }

    ## Add colors for rare/outlier samples from each brain region
    data$outlier_or_rare_samples_colors <- sapply(as.vector(data$outlier_or_rare_samples_labels),
                                                  function(x){if(!is.na(x) && x==rare_sample_ID){rare_and_poorQC_samples_colors[rare_sample_ID]}
                                                      else if(!is.na(x)){'gray35'}
                                                      else {NA}})

    if (sample_var=="Brain_Region"){
        violin_width=1
        jitter_width=0.1
        x_label="Brain Region"
    }
    else if (sample_var=="Substance"){
        violin_width=0.7
        jitter_width=0.1
        x_label="Substance"
    }
    else if (sample_var=="Brain_Region_and_Substance"){
        violin_width=0.7
        jitter_width=0.09
        x_label="Brain Region & Substance"
        data$Brain_Region_and_Substance <- gsub(' ', '\n', data$Brain_Region_and_Substance)
        names(colors$Brain_Region_and_Substance) <- gsub(' ', '\n', names(colors$Brain_Region_and_Substance))
    }
    else if (sample_var=='Total_Num_Fentanyl_Sessions'){
        violin_width=0.7
        jitter_width=0.1
        x_label="Total Number of Fentanyl Sessions"
    }
    else if (sample_var=='Num_Fentanyl_Sessions_six_hrs'){
        violin_width=0.7
        jitter_width=0.1
        x_label="Number of 6hrs Fentanyl Sessions"
    }
    else if (sample_var=='Batch_RNA_extraction'){
        violin_width=0.7
        jitter_width=0.1
        x_label="RNA extraction Batch"
    }
    else if (sample_var=='Batch_lib_prep'){
        violin_width=0.7
        jitter_width=0.1
        x_label="Library preparation Batch"
    }

    y_label=str_replace_all(qc_metric, c("_"=" "))

    pos <- position_jitter(width = 0.2, seed = 2)
    plot <- ggplot(data = data, mapping = aes(x = !! rlang::sym(sample_var),
                                              y = !! rlang::sym(qc_metric),
                                              color = !! rlang::sym(sample_var),
                                              label = outlier_or_rare_samples_labels)) +
        geom_violin(alpha = 0, size = 0.4, color='gray50', width=violin_width)+
        geom_jitter(alpha = 2, size = 2.4,  position = pos, aes(shape=eval(parse_expr(shape)))) +
        geom_boxplot(alpha = 0, size = 0.4, width=0.08, color='gray50') +
        sm_hgrid(legends = TRUE) +
        scale_color_manual(values = colors[[sample_var]]) +
        guides(color="none") +
        scale_shape_manual(name = shape_name, values=shapes) +
        geom_label_repel(label=data$outlier_or_rare_samples_labels,
                         color=data$outlier_or_rare_samples_colors,
                         size=3.2, max.overlaps = Inf,
                         box.padding = 0.7,
                         position = pos,
                         min.segment.length = 0) +
        labs(y= y_label, x = x_label) +
        theme(plot.margin=unit (c(0.6, 0.6, 0.6, 0.6), 'cm'),
              legend.position="right",
              axis.title = element_text(size = (12)),
              axis.text = element_text(size = (11)),
              legend.text = element_text(size=11),
              legend.title = element_text(size=12))

    return(plot)
}

## Plot multiple QC metrics per sample
multiple_QC_boxplots <- function(brain_region, rare_sample_ID){

    plots <- list()
    for(i in 1:length(qc_metrics)){

        ## For substance and without shape
        p <- QC_boxplot(qc_metrics[i], 'Substance', brain_region, rare_sample_ID, 'no_shape')
        plots[[i]] <- p
    }

    plot_grid(plotlist = plots, ncol=3, align = 'vh')
    ggsave(paste0("plots/04_EDA/02_PCA/", rare_sample_ID, "_QC_boxplots_", brain_region, "_no_shape.pdf"), width = 23, height = 23, units = "cm")
}


#######################
##  Habenula samples
#######################

########### Sample "5_F_LHb_13" ###########
## Sample "5_F_LHb_13" has the least number of expressed genes
multiple_QC_boxplots('habenula', '5_F_LHb_13')
colData(rse_gene_habenula_filt)[which.min(rse_gene_habenula_filt$detected_num_genes), 'SAMPLE_ID']
#  "5_F_LHb_13"

########### Sample "1_F_LHb_01" ###########
## Sample "1_F_LHb_01" has the lowest rate of concordant reads
multiple_QC_boxplots('habenula', '1_F_LHb_01')
colData(rse_gene_habenula_filt)[which.min(rse_gene_habenula_filt$concordMapRate), 'SAMPLE_ID']
#  "1_F_LHb_01"

########### Sample "3_F_LHb_09" ###########
## Sample "3_F_LHb_09" has the least amount of RNA
multiple_QC_boxplots('habenula', '3_F_LHb_09')
colData(rse_gene_habenula_filt)[which.min(rse_gene_habenula_filt$Total_RNA_amount), 'SAMPLE_ID']
#  "3_F_LHb_09"

########### Sample "18_S_LHb_20" ###########
## Sample "18_S_LHb_20" has the highest RNA concentration and RIN, and the largest number of expressed genes
multiple_QC_boxplots('habenula', '18_S_LHb_20')
colData(rse_gene_habenula_filt)[which.max(rse_gene_habenula_filt$RNA_concentration), 'SAMPLE_ID']
#  "18_S_LHb_20"
colData(rse_gene_habenula_filt)[which.max(rse_gene_habenula_filt$detected_num_genes), 'SAMPLE_ID']
#   "18_S_LHb_20"
colData(rse_gene_habenula_filt)[which.max(rse_gene_habenula_filt$RIN), 'SAMPLE_ID']
#   "18_S_LHb_20"


#######################
##  Amygdala samples
#######################

########### Sample "34_S_Amyg_22" ###########

## Sample "34_S_Amyg_22" has the biggest library size, the highest proportion of reads assigned to genes and
## the highest fraction of reads that mapped to the mitochondrial chr
multiple_QC_boxplots('amygdala', '34_S_Amyg_22')
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$library_size), 'SAMPLE_ID']
#  "34_S_Amyg_22"
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$totalAssignedGene), 'SAMPLE_ID']
#  "34_S_Amyg_22"
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$mitoRate), 'SAMPLE_ID']
#  "34_S_Amyg_22"

########### Sample "33_S_Amyg_20" ###########

## Sample "33_S_Amyg_20" has the least number of expressed genes but the highest amount of RNA and the highest
## fractions of both, concordant reads and reads that mapped to the reference genome
multiple_QC_boxplots('amygdala', '33_S_Amyg_20')
colData(rse_gene_amygdala_filt)[which.min(rse_gene_amygdala_filt$detected_num_genes), 'SAMPLE_ID']
#  "33_S_Amyg_20"
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$Total_RNA_amount), 'SAMPLE_ID']
#  "33_S_Amyg_20"
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$overallMapRate), 'SAMPLE_ID']
#  "33_S_Amyg_20"
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$concordMapRate), 'SAMPLE_ID']
#  "33_S_Amyg_20"

########### Sample "14_S_Amyg_14" ###########

## Sample "14_S_Amyg_14" has the lowest fraction of mitochondrial reads
multiple_QC_boxplots('amygdala', '14_S_Amyg_14')
colData(rse_gene_amygdala_filt)[which.min(rse_gene_amygdala_filt$mitoRate), 'SAMPLE_ID']
#  "14_S_Amyg_14

########### Sample "10_S_Amyg_06" ###########

## Sample "10_S_Amyg_06" has the lowest rates for concordant and overall mapping reads
multiple_QC_boxplots('amygdala', '10_S_Amyg_06')
colData(rse_gene_amygdala_filt)[which.min(rse_gene_amygdala_filt$overallMapRate), 'SAMPLE_ID']
#  "10_S_Amyg_06"
colData(rse_gene_amygdala_filt)[which.min(rse_gene_amygdala_filt$concordMapRate), 'SAMPLE_ID']
#  "10_S_Amyg_06"




## 2.3.3 Remove rare/outlier samples

# amyg_samples_to_remove <- c("34_S_Amyg_22", "14_S_Amyg_14")
# hab_samples_to_remove <- c("5_F_LHb_13")
# rse_gene_amygdala_filt <- rse_gene_amygdala_filt[,-which(rse_gene_amygdala_filt$SAMPLE_ID %in% amyg_samples_to_remove)]
# rse_gene_habenula_filt <- rse_gene_habenula_filt[,-which(rse_gene_habenula_filt$SAMPLE_ID %in% hab_samples_to_remove)]





## 2.4 Explore differences within fentanyl and saline sample groups

## PC boxplots

PC_boxplot <- function(PC, sample_var, brain_region){

    ## PC data
    pca<-PCA(brain_region)
    pca_data <- pca[[1]]
    pca_vars_labs<-pca[[2]]

    ## Add rare samples' labels and colors
    pca_data$outlier_or_rare_samples_labels <- sapply(pca_data$SAMPLE_ID,
                                                      function(x){if(x %in% names(rare_and_poorQC_samples_colors)){x}
                                                                  else{NA}})
    pca_data$outlier_or_rare_samples_colors <- sapply(pca_data$outlier_or_rare_samples_labels,
                                                      function(x){if(!is.na(x) && x %in% names(table(pca_data$Retention_sample_label))){'red'}
                                                                  else if(!is.na(x)){'gray35'}
                                                                  else{NA}})

    if(sample_var=='Substance'){
        x_label='Substance'
    }

    else if(sample_var=='Total_Num_Fentanyl_Sessions'){
        x_label="Total num of Fentanyl Sessions"
    }

    else if(sample_var=='Batch_RNA_extraction'){
        x_label="RNA extraction Batch"
    }

    else if(sample_var=='Batch_lib_prep'){
        x_label="Library preparation Batch"
    }

    pos <- position_jitter(width = 0.2, seed = 2)

    plot <- ggplot(data = pca_data, mapping = aes(x = !! rlang::sym(sample_var),y = !! rlang::sym(PC),
                                                  color = Substance, label=outlier_or_rare_samples_labels)) +
        geom_boxplot(size = 0.35, width=0.32, color='black', outlier.color = "#FFFFFFFF") +
        geom_jitter(alpha = 1, size = 1.4, position = pos) +
        scale_color_manual(values = c('Fentanyl'='turquoise3', 'Saline'='yellow3')) +
        labs(y= pca_vars_labs[strtoi(gsub("PC","", PC))], x = x_label, color='Substance') +
        sm_hgrid(legends = TRUE) +
        geom_label_repel(color=pca_data$outlier_or_rare_samples_colors, size=2.5, max.overlaps = Inf,
                         box.padding = 0.7,
                         position = pos,
                         min.segment.length = 0) +
        theme(plot.margin=unit (c (0.5, 0.5, 0.5, 0.5), 'cm'),
              axis.title = element_text(size = (12)),
              axis.text = element_text(size = (11)),
              legend.text = element_text(size=11),
              legend.title = element_text(size=12))

    return(plot)
}

## Multiple plots
multiple_PC_boxplots <- function(brain_region){
    ## Habenula
    if (brain_region=='habenula'){

        sample_vars <- c('Substance', 'Total_Num_Fentanyl_Sessions' ,'Batch_RNA_extraction')
        PCs <- c('PC1', 'PC3', 'PC4')

        i=1
        plots=list()

        for (PC in PCs){
            for (sample_var in sample_vars){
                plots[[i]] <- PC_boxplot(PC, sample_var, 'habenula')
                i=i+1
            }
        }
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], nrow=3)
        ggsave(here(paste("plots/04_EDA/02_PCA/Boxplots_PCs_", brain_region, ".pdf", sep="")), width = 32, height = 24, units = "cm")
    }

    ## Amygdala
    else {
        sample_vars <- c('Substance', 'Total_Num_Fentanyl_Sessions' ,'Batch_RNA_extraction', 'Batch_lib_prep')
        PCs <- c('PC1', 'PC2')

        i=1
        plots=list()

        for (PC in PCs){
            for (sample_var in sample_vars){
                plots[[i]] <- PC_boxplot(PC, sample_var, brain_region)
                i=i+1
            }
        }
        plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], nrow=2)
        ggsave(here(paste("plots/04_EDA/02_PCA/Boxplots_PCs_", brain_region, ".pdf", sep="")), width = 44, height = 16, units = "cm")
    }

}

multiple_PC_boxplots('habenula')
multiple_PC_boxplots('amygdala')







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
# user   system  elapsed
# 186.09    58.89 87926.95
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
# date     2024-03-08
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# abind                  1.4-5      2016-07-21 [1] CRAN (R 4.3.0)
# backports              1.4.1      2021-12-13 [1] CRAN (R 4.3.0)
# base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0     2023-10-26 [1] Bioconductor
# BiocGenerics         * 0.48.1     2023-11-02 [1] Bioconductor
# bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.3.0)
# broom                  1.0.5      2023-06-09 [1] CRAN (R 4.3.0)
# car                    3.1-2      2023-03-30 [1] CRAN (R 4.3.0)
# carData                3.0-5      2022-01-06 [1] CRAN (R 4.3.0)
# checkmate              2.3.1      2023-12-04 [1] CRAN (R 4.3.1)
# cli                    3.6.2      2023-12-11 [1] CRAN (R 4.3.1)
# cluster                2.1.6      2023-12-01 [1] CRAN (R 4.3.1)
# colorspace             2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.3      2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# data.table             1.15.2     2024-02-29 [1] CRAN (R 4.3.1)
# DelayedArray           0.28.0     2023-11-06 [1] Bioconductor
# digest                 0.6.34     2024-01-11 [1] CRAN (R 4.3.1)
# dplyr                  1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# edgeR                  4.0.16     2024-02-20 [1] Bioconductor 3.18 (R 4.3.2)
# evaluate               0.23       2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# foreign                0.8-86     2023-11-28 [1] CRAN (R 4.3.1)
# Formula                1.2-5      2023-02-24 [1] CRAN (R 4.3.0)
# fs                     1.6.3      2023-07-20 [1] CRAN (R 4.3.0)
# gargle                 1.5.2      2023-07-20 [1] CRAN (R 4.3.0)
# generics               0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.6     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11     2024-02-17 [1] Bioconductor
# GenomicRanges        * 1.54.1     2023-10-30 [1] Bioconductor
# gghalves               0.1.4      2022-11-20 [1] CRAN (R 4.3.0)
# ggnewscale           * 0.4.10     2024-02-08 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.0      2024-02-23 [1] CRAN (R 4.3.1)
# ggpubr                 0.6.0      2023-02-10 [1] CRAN (R 4.3.0)
# ggrepel              * 0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# ggsignif               0.6.4      2022-10-13 [1] CRAN (R 4.3.0)
# glue                   1.7.0      2024-01-09 [1] CRAN (R 4.3.1)
# googledrive            2.1.1      2023-06-11 [1] CRAN (R 4.3.0)
# gridExtra              2.3        2017-09-09 [1] CRAN (R 4.3.0)
# gtable                 0.3.4      2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.3.0)
# Hmisc                  5.1-1      2023-09-12 [1] CRAN (R 4.3.0)
# htmlTable              2.4.2      2023-10-29 [1] CRAN (R 4.3.1)
# htmltools              0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
# htmlwidgets            1.6.4      2023-12-06 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0     2023-10-26 [1] Bioconductor
# jaffelab             * 0.99.32    2023-05-28 [1] Github (LieberInstitute/jaffelab@21e6574)
# knitr                  1.45       2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3      2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# limma                  3.58.1     2023-11-02 [1] Bioconductor
# locfit                 1.5-9.9    2024-03-01 [1] CRAN (R 4.3.1)
# magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# MASS                   7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.1)
# Matrix                 1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0     2023-10-26 [1] Bioconductor
# matrixStats          * 1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# nlme                   3.1-164    2023-11-27 [1] CRAN (R 4.3.1)
# nnet                   7.3-19     2023-05-03 [1] CRAN (R 4.3.2)
# pillar                 1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# purrr                  1.0.2      2023-08-10 [1] CRAN (R 4.3.0)
# pwr                    1.3-0      2020-03-17 [1] CRAN (R 4.3.0)
# R6                     2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.3.0)
# ragg                   1.2.7      2023-12-11 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12     2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.3.1)
# rlang                * 1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26       2024-03-05 [1] CRAN (R 4.3.1)
# rpart                  4.1.23     2023-12-05 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
# rstatix                0.7.2      2023-02-01 [1] CRAN (R 4.3.0)
# rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0      2023-10-26 [1] Bioconductor
# S4Vectors            * 0.40.2     2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# sdamr                  0.2.0      2022-11-16 [1] CRAN (R 4.3.0)
# segmented              2.0-3      2024-02-16 [1] CRAN (R 4.3.2)
# sessioninfo          * 1.2.2      2021-12-06 [1] CRAN (R 4.3.0)
# smplot2              * 0.1.0      2023-06-07 [1] Github (smin95/smplot2@836f909)
# SparseArray            1.2.4      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# statmod                1.5.0      2023-01-06 [1] CRAN (R 4.3.0)
# stringi                1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# SummarizedExperiment * 1.32.0     2023-11-06 [1] Bioconductor
# systemfonts            1.0.5      2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  1.3.1      2024-01-24 [1] CRAN (R 4.3.1)
# tidyselect             1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# withr                  3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# xfun                   0.42       2024-02-08 [1] CRAN (R 4.3.1)
# XVector                0.42.0     2023-10-26 [1] Bioconductor
# zlibbioc               1.48.0     2023-10-26 [1] Bioconductor
# zoo                    1.8-12     2023-04-13 [1] CRAN (R 4.3.0)
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

