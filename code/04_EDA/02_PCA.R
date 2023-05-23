
library(here)
library(SummarizedExperiment)
library(jaffelab)
library(ggplot2)
library(ggrepel)
library(stringr)
library(smplot2)
library(cowplot)
library(rlang)
library(sessioninfo)



################################################################################
##                        2. Explore sample-level effects
################################################################################

#######################   Principal Component Analysis   #######################

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


## 2.1 PCA data obtention

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



## 2.2 PCA visualization

## PCx vs PCy Plots

PCx_vs_PCy <- function (PCx, PCy, pca_data, pca_vars_labs, sample_var, brain_region) {

    legend_height = 1
    legend_width = 1
    legend_text_size = 8
    legend_title_size = 10

    if(sample_var=='Brain_Region'){
        colors=c('Amygdala'='palegreen3', 'Habenula'='orchid1')
        margin=0.7
    }

    else if(sample_var=='Substance'){
        colors=c('Fentanyl'='turquoise3', 'Saline'='yellow3')
        if (is.null(brain_region)){
            margin=0.7
        }
        else {
            margin=0.9
        }
    }

    else if(sample_var=='Brain_Region_and_Substance'){
        colors=c('Amygdala Fentanyl'='springgreen3' , 'Amygdala Saline'='yellowgreen', 'Habenula Fentanyl'='hotpink1', 'Habenula Saline'='violet')
        margin=0.1
    }

    else if(sample_var=='Total_Num_Fentanyl_Sessions'){
        colors=c('24'='salmon', '22'='pink2')
        margin=0.1
    }

    else if(sample_var=='Num_Fentanyl_Sessions_six_hrs'){
        colors=c('18'='dodgerblue', '16'='lightskyblue')
        margin=0.1
    }

    else if(sample_var=='Batch_RNA_extraction'){
        colors=c('1'='darksalmon', '2'='darkseagreen3', '3'= 'lightsteelblue2')
        if (is.null(brain_region)){
            margin=0.7
        }
        else {
            margin=0.9
        }
    }
    else if(sample_var=='Batch_lib_prep'){
        colors=c('1'='darkgoldenrod3', '2'='mediumpurple2', '3'= 'darkmagenta')
        if (is.null(brain_region)){
            margin=0.7
        }
        else {
            margin=0.9
        }
    }


    plot <- ggplot(data=pca_data,
                aes(x=eval(parse_expr(PCx)),y=eval(parse_expr(PCy)),
                color=eval(parse_expr(sample_var)),
                label=outlier_or_rare_samples_labels) ) +
            theme_classic() +
            theme(legend.position="right", plot.margin=unit (c (1,margin,1,margin), 'cm'),
                  legend.text = element_text(size=legend_text_size),
                  legend.title = element_text(size=legend_title_size)) +
            geom_point(aes(shape=Batch_RNA_extraction), size=3) +
            scale_color_manual(values = colors) +
            scale_shape_manual(name='Batch RNA extraction', values=c('1'=8, '2'=10, '3'=15)) +
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
plot_PCAs<-function(brain_region, filename){

    ## PC data
    pca<-PCA(brain_region)
    pca_data <- pca[[1]]
    pca_vars_labs<-pca[[2]]

    ## Plots for ALL SAMPLES
    if (is.null(brain_region)){
        for (PCs in list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"), c("PC3", "PC4"), c("PC5", "PC6"))){
            plots<-list()
            i=1
            for (sample_var in c("Brain_Region", "Substance", "Brain_Region_and_Substance",
                                 "Batch_RNA_extraction", "Batch_lib_prep",
                                 "Num_Fentanyl_Sessions_six_hrs", "Total_Num_Fentanyl_Sessions")){
                p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars_labs, sample_var, brain_region)
                plots[[i]]=p
                i=i+1
            }
            plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], nrow = 2)
            ## Save plots
            ggsave(paste("plots/04_EDA/02_PCA/",PCs[1],"_vs_",PCs[2],"_all_samples.pdf", sep=""),
                   width = 52, height = 19, units = "cm")
        }
    }
    ## For habenula and amygdala separately
    else {
        for (PCs in list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"), c("PC3", "PC4"), c("PC5", "PC6"))){
            plots<-list()
            i=1
            for (sample_var in c("Substance", "Batch_RNA_extraction", "Batch_lib_prep",
                                 "Total_Num_Fentanyl_Sessions", "Num_Fentanyl_Sessions_six_hrs")){
                p<-PCx_vs_PCy(PCs[1], PCs[2], pca_data, pca_vars_labs, sample_var, brain_region)
                plots[[i]]=p
                i=i+1
            }
            plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], nrow = 2)

            ## Save plots
            if (is.null(filename)){
                ggsave(paste("plots/04_EDA/02_PCA/",PCs[1],"_vs_",PCs[2],"_", brain_region,".pdf", sep=""),
                       width = 42, height = 20, units = "cm")
            }
            else {
                ggsave(paste("plots/04_EDA/02_PCA/New_",PCs[1],"_vs_",PCs[2],"_", brain_region,".pdf", sep=""),
                       width = 42, height = 20, units = "cm")
            }

        }
    }
}


## Plots
plot_PCAs(NULL, NULL)
plot_PCAs('habenula', NULL)
plot_PCAs('amygdala', NULL)





## 2.3 Manual sample filtering

## 2.3.1 Identify rare samples in PCA plots

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
                                                                    else if (x %in% rare_amyg_samples){'yellow4'}
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
                                                                            else if (x %in% rare_hab_samples){'yellow4'}
                                                                            else{NA}})

## Create same variables for rse_gene_filt
rse_gene_filt$outlier_or_rare_samples_labels <- rep(NA, dim(rse_gene_filt)[2])
rse_gene_filt$outlier_or_rare_samples_colors <- rep(NA, dim(rse_gene_filt)[2])




## 2.3.2 Explore QC metrics of rare samples

#######################
##  Habenula samples
#######################

########### Sample "5_F_LHb_13" ###########

## Sample "5_F_LHb_13" has the least number of expressed genes
colData(rse_gene_habenula_filt)[which.min(rse_gene_habenula_filt$detected_num_genes), 'SAMPLE_ID']
#  "5_F_LHb_13"

########### Sample "1_F_LHb_01" ###########

## Sample "1_F_LHb_01" has the lowest rate of concordant reads
colData(rse_gene_habenula_filt)[which.min(rse_gene_habenula_filt$concordMapRate), 'SAMPLE_ID']
#  "1_F_LHb_01"

########### Sample "3_F_LHb_09" ###########

## Sample "3_F_LHb_09" has the least amount of RNA
colData(rse_gene_habenula_filt)[which.min(rse_gene_habenula_filt$Total_RNA_amount), 'SAMPLE_ID']
#  "3_F_LHb_09"

########### Sample "18_S_LHb_20" ###########

## Sample "18_S_LHb_20" has the highest RNA concentration and RIN, and the largest number of expressed genes
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
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$library_size), 'SAMPLE_ID']
#  "34_S_Amyg_22"
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$totalAssignedGene), 'SAMPLE_ID']
#  "34_S_Amyg_22"
colData(rse_gene_amygdala_filt)[which.max(rse_gene_amygdala_filt$mitoRate), 'SAMPLE_ID']
#  "34_S_Amyg_22"

########### Sample "33_S_Amyg_20" ###########

## Sample "33_S_Amyg_20" has the least number of expressed genes but the highest amount of RNA and the highest
## fractions of both, concordant reads and reads that mapped to the reference genome
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
colData(rse_gene_amygdala_filt)[which.min(rse_gene_amygdala_filt$mitoRate), 'SAMPLE_ID']
#  "14_S_Amyg_14

########### Sample "10_S_Amyg_06" ###########

## Sample "10_S_Amyg_06" has the lowest rates for concordant and overall mapping reads
colData(rse_gene_amygdala_filt)[which.min(rse_gene_amygdala_filt$overallMapRate), 'SAMPLE_ID']
#  "10_S_Amyg_06"
colData(rse_gene_amygdala_filt)[which.min(rse_gene_amygdala_filt$concordMapRate), 'SAMPLE_ID']
#  "10_S_Amyg_06"

########### Sample "16_S_Amyg_18" ###########

## No special QC information identified



## Plot QC metrics of rare samples

QC_boxplot <- function(qc_metric, sample_var, brain_region, rare_sample_ID){

    data<-as.data.frame(colData(rse_gene_filt))
    data_brain_region <- as.data.frame(colData(eval(parse_expr(paste('rse_gene', brain_region, 'filt', sep='_')))))
    ## Add labels and colors for rare/outlier samples from each brain region
    data$outlier_or_rare_samples_labels <- apply(data, 1, function(x){if(x['SAMPLE_ID'] %in% data_brain_region$SAMPLE_ID)
                                                                        {data_brain_region[which(data_brain_region$SAMPLE_ID==x['SAMPLE_ID']), 'outlier_or_rare_samples_labels']}
                                                                     else {NA} })
    data$outlier_or_rare_samples_colors <- apply(data, 1, function(x){if(x['SAMPLE_ID']==rare_sample_ID){'red'}
                                                                     else if(x['SAMPLE_ID'] %in% data_brain_region$SAMPLE_ID){'gray55'}
                                                                     else {NA}})

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

    pos <- position_jitter(width = 0.2, seed = 2)
    plot <- ggplot(data = data, mapping = aes(x = !! rlang::sym(sample_var),
                                              y = !! rlang::sym(qc_metric),
                                              color = !! rlang::sym(sample_var),
                                              label = outlier_or_rare_samples_labels)) +
        geom_violin(alpha = 0, size = 0.4, color='gray80', width=violin_width)+
        geom_jitter(alpha = 1, size = 2.4,  position = pos, aes(shape=Batch_RNA_extraction)) +
        geom_boxplot(alpha = 0, size = 0.4, width=0.08, color='gray80') +
        sm_hgrid(legends = TRUE) +
        scale_color_manual(values = colors) +
        guides(color="none") +
        ## Shape for RNA extraction batch
        scale_shape_manual(name='Batch RNA extraction', labels=c("1","2","3"), values=c('1'=8, '2'=10, '3'=15)) +
        geom_label_repel(color=data$outlier_or_rare_samples_colors, size=2.4, max.overlaps = Inf,
                         box.padding = 0.7,
                         position = pos,
                         min.segment.length = 0) +
        labs(y= y_label, x = x_label) +
        theme(axis.title = element_text(size = (9)),
              axis.text = element_text(size = (8)))

    return(plot)
}


#######################
##  Habenula samples
#######################

########### Sample "5_F_LHb_13" ###########
## Sample "5_F_LHb_13" has the least number of expressed genes
QC_boxplot('detected_num_genes', 'Brain_Region_and_Substance', 'habenula', "5_F_LHb_13")
ggsave("plots/04_EDA/02_PCA/5_F_LHb_13_QC_boxplot_habenula.pdf", width = 20, height = 15, units = "cm")

########### Sample "1_F_LHb_01" ###########
## Sample "1_F_LHb_01" has the lowest rate of concordant reads
QC_boxplot('concordMapRate', 'Brain_Region_and_Substance', 'habenula', "1_F_LHb_01")
ggsave("plots/04_EDA/02_PCA/1_F_LHb_01_QC_boxplot_habenula.pdf", width = 20, height = 15, units = "cm")

########### Sample "3_F_LHb_09" ###########

## Sample "3_F_LHb_09" has the least amount of RNA
QC_boxplot('Total_RNA_amount', 'Brain_Region_and_Substance', 'habenula', "3_F_LHb_09")
ggsave("plots/04_EDA/02_PCA/3_F_LHb_09_QC_boxplot_habenula.pdf", width = 20, height = 15, units = "cm")

########### Sample "18_S_LHb_20" ###########

## Sample "18_S_LHb_20" has the highest RNA concentration and RIN, and the largest number of expressed genes
p1 <- QC_boxplot('RNA_concentration', 'Brain_Region_and_Substance', 'habenula', "18_S_LHb_20")
p2 <- QC_boxplot('RIN', 'Brain_Region_and_Substance', 'habenula', "18_S_LHb_20")
p3 <- QC_boxplot('detected_num_genes', 'Brain_Region_and_Substance', 'habenula', "18_S_LHb_20")
plot_grid(p1, p2, p3, ncol=3)
ggsave("plots/04_EDA/02_PCA/18_S_LHb_20_QC_boxplot_habenula.pdf", width = 60, height = 15, units = "cm")



#######################
##  Amygdala samples
#######################

########### Sample "34_S_Amyg_22" ###########

## Sample "34_S_Amyg_22" has the biggest library size, the highest proportion of reads assigned to genes and
## the highest fraction of reads that mapped to the mitochondrial chr
p1 <- QC_boxplot('library_size', 'Brain_Region_and_Substance', 'amygdala', "34_S_Amyg_22")
p2 <- QC_boxplot('totalAssignedGene', 'Brain_Region_and_Substance', 'amygdala', "34_S_Amyg_22")
p3 <- QC_boxplot('mitoRate', 'Brain_Region_and_Substance', 'amygdala', "34_S_Amyg_22")
plot_grid(p1, p2, p3, ncol=3)
ggsave("plots/04_EDA/02_PCA/34_S_Amyg_22_QC_boxplot_amygdala.pdf", width = 60, height = 15, units = "cm")

########### Sample "33_S_Amyg_20" ###########

## Sample "33_S_Amyg_20" has the least number of expressed genes but the highest amount of RNA and the highest
## fractions of both, concordant reads and reads that mapped to the reference genome
p1 <- QC_boxplot('detected_num_genes', 'Brain_Region_and_Substance', 'amygdala', "33_S_Amyg_20")
p2 <- QC_boxplot('Total_RNA_amount', 'Brain_Region_and_Substance', 'amygdala', "33_S_Amyg_20")
p3 <- QC_boxplot('concordMapRate', 'Brain_Region_and_Substance', 'amygdala', "33_S_Amyg_20")
p4 <- QC_boxplot('overallMapRate', 'Brain_Region_and_Substance', 'amygdala', "33_S_Amyg_20")
plot_grid(p1, p2, p3, p4, ncol=2)
ggsave("plots/04_EDA/02_PCA/33_S_Amyg_20_QC_boxplot_amygdala.pdf", width = 40, height = 30, units = "cm")

########### Sample "14_S_Amyg_14" ###########

## Sample "14_S_Amyg_14" has the lowest fraction of mitochondrial reads
QC_boxplot('mitoRate', 'Brain_Region_and_Substance', 'amygdala', "14_S_Amyg_14")
ggsave("plots/04_EDA/02_PCA/14_S_Amyg_14_QC_boxplot_amygdala.pdf", width = 20, height = 15, units = "cm")

########### Sample "10_S_Amyg_06" ###########

## Sample "10_S_Amyg_06" has the lowest rates for concordant and overall mapping reads
p1 <- QC_boxplot('concordMapRate', 'Brain_Region_and_Substance', 'amygdala', "10_S_Amyg_06")
p2 <- QC_boxplot('overallMapRate', 'Brain_Region_and_Substance', 'amygdala', "10_S_Amyg_06")
plot_grid(p1, p2, ncol=2)
ggsave("plots/04_EDA/02_PCA/10_S_Amyg_06_QC_boxplot_amygdala.pdf", width = 40, height = 15, units = "cm")





## 2.3.3 Remove rare/outlier samples

# amyg_samples_to_remove <- c("34_S_Amyg_22", "14_S_Amyg_14")
# hab_samples_to_remove <- c("5_F_LHb_13")
# rse_gene_amygdala_filt <- rse_gene_amygdala_filt[,-which(rse_gene_amygdala_filt$SAMPLE_ID %in% amyg_samples_to_remove)]
# rse_gene_habenula_filt <- rse_gene_habenula_filt[,-which(rse_gene_habenula_filt$SAMPLE_ID %in% hab_samples_to_remove)]
#
# ## New PCA plots
# plot_PCAs('habenula', 'new_plots')
# plot_PCAs('amygdala', 'new_plots')







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

