
library(here)
library(SummarizedExperiment)
library(jaffelab)
library(ggplot2)
library(ggrepel)
library(stringr)
library(cowplot)
library(rlang)
library(sessioninfo)



########################   Dimensionality Reduction   #########################


################################################################################
##                        1. Principal Component Analysis
################################################################################
## (Note: genes are filtered and counts normalized in these analyses)

load(here('processed-data/03_Data_preparation/rse_gene_filt.Rdata'), verbose=TRUE)
load(here('processed-data/04_EDA/01_QCA/rse_gene_habenula_complete.Rdata'), verbose = TRUE)
load(here('processed-data/04_EDA/01_QCA/rse_gene_amygdala_complete.Rdata'), verbose = TRUE)

## Habenula filtered data
rse_gene_habenula_filt <- rse_gene_filt[,colData(rse_gene_filt)$Brain_Region=='Habenula']
## Add info for outlier samples after QC sample filtering of habenula samples
rse_gene_habenula_filt$Retention_after_QC_filtering <- rse_gene_habenula_complete$Retention_after_QC_filtering
rse_gene_habenula_filt$Retention_sample_label <- rse_gene_habenula_complete$Retention_sample_label
save(rse_gene_habenula_filt, file='processed-data/04_EDA/02_PCA/rse_gene_habenula_filt.Rdata')

## Amygdala filtered data
rse_gene_amygdala_filt <- rse_gene_filt[,colData(rse_gene_filt)$Brain_Region=='Amygdala']
## Add info for outlier samples after QC sample filtering of amygdala samples
rse_gene_amygdala_filt$Retention_after_QC_filtering <- rse_gene_amygdala_complete$Retention_after_QC_filtering
rse_gene_amygdala_filt$Retention_sample_label <- rse_gene_amygdala_complete$Retention_sample_label
save(rse_gene_amygdala_filt, file='processed-data/04_EDA/02_PCA/rse_gene_amygdala_filt.Rdata')



### Explore samples' gene expression variation


## 1.1 Generate PCA data

## Function to obtain PC
PCA<-function(brain_region){
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



## 1.2 Generate PCA plots

## PCx vs PCy Plots

PCx_vs_PCy <- function (PCx, PCy, pca_data, pca_vars_labs, sample_var, brain_region) {

    if(sample_var=='Brain_Region'){
        colors=c('Amygdala'='palegreen3', 'Habenula'='orchid1')
        legend_height = 1
        legend_width = 1
        legend_text_size = 8
        legend_title_size = 10
        margin=1.5
    }

    else if(sample_var=='Substance'){
        colors=c('Fentanyl'='turquoise3', 'Saline'='yellow3')
        legend_height = 1
        legend_width = 1
        legend_text_size = 8
        legend_title_size = 10
        margin=1.5
    }

    else if(sample_var=='Brain_Region_and_Substance'){
        colors=c('Amygdala Fentanyl'='springgreen3' , 'Amygdala Saline'='yellowgreen', 'Habenula Fentanyl'='hotpink1', 'Habenula Saline'='violet')
        legend_height = 1
        legend_width = 1
        legend_text_size = 8
        legend_title_size = 10
        margin=0.3
    }

    else if(sample_var=='Total_Num_Fentanyl_Sessions'){
        colors=c('24'='salmon', '22'='pink2')
        legend_height = 1
        legend_width = 1
        legend_text_size = 8
        legend_title_size = 10
        margin=0.2
    }

    else if(sample_var=='Num_Fentanyl_Sessions_six_hrs'){
        colors=c('18'='dodgerblue', '16'='lightskyblue')
        legend_height = 1
        legend_width = 1
        legend_text_size = 8
        legend_title_size = 10
        margin=0.1
    }

    else if(sample_var=='Batch'){
        colors=c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2' , '4'= 'navajowhite2')
        legend_height = 1
        legend_width = 1
        legend_text_size = 8
        legend_title_size = 10
        margin=1.5
    }

    plot <- ggplot(data=pca_data,
                aes(x=eval(parse_expr(PCx)),y=eval(parse_expr(PCy)),
                color=eval(parse_expr(sample_var)),
                label=outlier_or_rare_samples_labels) ) +
            theme_classic() +
            theme(legend.position="right", plot.margin=unit (c (1,margin,1,margin), 'cm'),
                  legend.text = element_text(size=legend_text_size),
                  legend.title = element_text(size=legend_title_size)) +
            geom_point(size=2) +
            scale_color_manual(values = colors) +
            ## Labels of removed samples
            geom_label_repel(color=pca_data$outlier_or_rare_samples_colors, size=2, max.overlaps = Inf,
                             box.padding = 0.7,
                             show.legend=FALSE,
                             min.segment.length = 0) +
            labs(x= pca_vars_labs[strtoi(gsub("PC","", PCx))], y = pca_vars_labs[strtoi(gsub("PC","", PCy))],
                 color=str_replace_all(sample_var, c("_"=" ")))
    return(plot)
}


## All PCA plots
plot_PCAs<-function(brain_region){

    ## PC data
    pca<-PCA(brain_region)
    pca_data <- pca[[1]]
    pca_vars_labs<-pca[[2]]

    ## Plots for ALL SAMPLES
    if (is.null(brain_region)){
        for (PCs in list(c("PC1", "PC2"), c("PC3", "PC4"), c("PC5", "PC6"))){
            plots<-list()
            i=1
            for (sample_var in c("Brain_Region", "Substance", "Brain_Region_and_Substance",
                                 "Total_Num_Fentanyl_Sessions", "Num_Fentanyl_Sessions_six_hrs", "Batch")){
                p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars_labs, sample_var, brain_region)
                plots[[i]]=p
                i=i+1
            }
            plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow = 2)
            ## Save plots
            ggsave(paste("plots/04_EDA/02_PCA/",PCs[1],"_vs_",PCs[2],"_all_samples.pdf", sep=""),
                   width = 40, height = 19, units = "cm")
        }
    }
    ## For habenula and amygdala separately
    else {
        for (PCs in list(c("PC1", "PC2"), c("PC3", "PC4"), c("PC5", "PC6"))){
            plots<-list()
            i=1
            for (sample_var in c("Substance", "Total_Num_Fentanyl_Sessions", "Num_Fentanyl_Sessions_six_hrs", "Batch")){
                p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars_labs, sample_var, brain_region)
                plots[[i]]=p
                i=i+1
            }
            plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], nrow = 2)
            ## Save plots
            ggsave(paste("plots/04_EDA/02_PCA/",PCs[1],"_vs_",PCs[2],"_", brain_region,".pdf", sep=""),
                   width = 28, height = 20, units = "cm")
        }
    }
}


## Plots
plot_PCAs(NULL)
plot_PCAs('habenula')
plot_PCAs('amygdala')





## 1.3 Manual sample filtering

## 1.3.1 Identify rare samples in PCA plots

#####################
##  Amygdala plots
#####################
amyg_pca_data <- PCA('amygdala')[[1]]

######## 1. Saline samples within fentanyl samples' group:  ########

amyg_sal_pca_data <- amyg_pca_data[which(amyg_pca_data$Substance=='Saline'),]

## -> One is the "34_S_Amyg_22" outlier sample: saline sample with the second lowest PC2 value
amyg_sal_pca_data[order(amyg_sal_pca_data$PC2), 'SAMPLE_ID'][2]
#  "34_S_Amyg_22"

## -> The other is the "14_S_Amyg_14" sample: saline sample with the lowest value in PC2
amyg_sal_pca_data[order(amyg_sal_pca_data$PC2), 'SAMPLE_ID'][1]
#  "14_S_Amyg_14"


######## 2. Sample with the highest value in PC4:  ########

## -> It is the "16_S_Amyg_18" saline sample
amyg_pca_data[which.max(amyg_pca_data$PC4), 'SAMPLE_ID']
#  "16_S_Amyg_18"


######## 3. Sample with the lowest value in PC6:  ########

## -> It is the "34_S_Amyg_22" outlier sample
amyg_pca_data[which.min(amyg_pca_data$PC6), 'SAMPLE_ID']
#  "34_S_Amyg_22"


## All rare amygdala samples
rare_amyg_samples <- c("34_S_Amyg_22", "14_S_Amyg_14", "16_S_Amyg_18")

## Add sample ID label for rare/outlier samples
rse_gene_amygdala_filt$outlier_or_rare_samples_labels <- sapply(rse_gene_amygdala_filt$SAMPLE_ID,
                                                                function(x){if(x %in% rse_gene_amygdala_filt$Retention_sample_label | x %in% rare_amyg_samples){x}
                                                                    else {NA}})
## Labels' colors
rse_gene_amygdala_filt$outlier_or_rare_samples_colors <- sapply(rse_gene_amygdala_filt$SAMPLE_ID,
                                                                function(x){if(x %in% rse_gene_amygdala_filt$Retention_sample_label){'gray30'}
                                                                    else if (x %in% rare_amyg_samples){'mistyrose3'}
                                                                    else{NA}})


#####################
## Habenula samples
#####################
hab_pca_data <- PCA('habenula')[[1]]

######## 1. Samples that have the lowest PC2 values:  ########

## -> One is the "5_F_LHb_13" outlier sample
hab_pca_data[order(hab_pca_data$PC2), 'SAMPLE_ID'][1]
#  "5_F_LHb_13"

## -> The second is the "3_F_LHb_09" sample
hab_pca_data[order(hab_pca_data$PC2), 'SAMPLE_ID'][2]
#  "3_F_LHb_09"


######## 2. Fentanyl samples that appear within saline samples' group:  ########

## -> One is the "5_F_LHb_13" outlier sample: sample with the highest PC4 value
hab_pca_data[which.max(hab_pca_data$PC4), 'SAMPLE_ID']
#  "5_F_LHb_13"

hab_fenta_pca_data <- hab_pca_data[which(hab_pca_data$Substance=='Fentanyl'),]
## -> Another is the "1_F_LHb_01" sample: fentanyl sample with the third lowest PC3 value
hab_fenta_pca_data[order(hab_fenta_pca_data$PC3), 'SAMPLE_ID'][3]
#  "1_F_LHb_01"


######## 3. Saline samples within fentanyl group:  ########

hab_sal_pca_data <- hab_pca_data[which(hab_pca_data$Substance=="Saline"),]
## -> It is the "18_S_LHb_20" sample: the saline sample with the lowest PC4 value
hab_sal_pca_data[which.min(hab_sal_pca_data$PC4), 'SAMPLE_ID']
#  "18_S_LHb_20"


## All rare habenula samples
rare_hab_samples <- c("5_F_LHb_13", "3_F_LHb_09", "1_F_LHb_01", "18_S_LHb_20")

## Add sample ID label for rare/outlier samples
rse_gene_habenula_filt$outlier_or_rare_samples_labels <- sapply(rse_gene_habenula_filt$SAMPLE_ID,
                                                                function(x){if(x %in% rse_gene_habenula_filt$Retention_sample_label | x %in% rare_hab_samples){x}
                                                                            else {NA}})
## Labels' colors
rse_gene_habenula_filt$outlier_or_rare_samples_colors <- sapply(rse_gene_habenula_filt$SAMPLE_ID,
                                                                function(x){if(x %in% rse_gene_habenula_filt$Retention_sample_label){'gray30'}
                                                                            else if (x %in% rare_hab_samples){'mistyrose3'}
                                                                            else{NA}})


## Create same variables for rse_gene
rse_gene_filt$outlier_or_rare_samples_labels <- rep(NA, dim(rse_gene_filt)[2])
rse_gene_filt$outlier_or_rare_samples_colors <- rep(NA, dim(rse_gene_filt)[2])







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
# date     2023-05-10
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   2.19.2 @ /private/var/folders/r6/94s0dsks4m3298b0d1mrp5tw0000gn/T/AppTranslocation/8639B884-1A6F-4EA4-9186-4C76049940ED/d/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)

