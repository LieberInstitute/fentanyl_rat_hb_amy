
library(here)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(GGally)
library(RRHO2)
library(ComplexHeatmap)
library(circlize)
library(sessioninfo)


################################################################################
##                                2. Comparisons
################################################################################

## Habenula results
load(here("processed-data/05_DEA/results_Substance_uncorr_vars_habenula.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_FirstHrIntakeSlope_habenula.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_TotalIntake_habenula.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_LastSessionIntake_habenula.Rdata"), verbose = TRUE)
load(here('processed-data/05_DEA/de_genes_Substance_habenula.Rdata'), verbose = TRUE)

## Amygdala results
load(here("processed-data/05_DEA/results_Substance_uncorr_vars_amygdala.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_FirstHrIntakeSlope_amygdala.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_TotalIntake_amygdala.Rdata"), verbose = TRUE)
load(here("processed-data/05_DEA/results_LastSessionIntake_amygdala.Rdata"), verbose = TRUE)
load(here('processed-data/05_DEA/de_genes_Substance_amygdala.Rdata'), verbose = TRUE)


## Plotting functions:
## Function to add gene DE info with respect to two groups
add_DE_info <-function(t_stats_1, t_stats_2, name_1, name_2) {

    DE<-vector()

    for (i in 1:dim(t_stats_1)[1]) {
        ## DE genes in both DGE analyses
        if (t_stats_1$adj.P.Val[i]<0.05 && t_stats_2$adj.P.Val[i]<0.05) {
            DE<-append(DE, "sig Both")
        }
        ## DE genes in only one DEA
        else if (t_stats_1$adj.P.Val[i]<0.05 && !t_stats_2$adj.P.Val[i]<0.05) {
            DE<-append(DE, paste("sig", name_1))
        }

        else if (t_stats_2$adj.P.Val[i]<0.05 && !t_stats_1$adj.P.Val[i]<0.05) {
            DE<-append(DE, paste("sig", name_2))
        }
        ## No DE genes in neither group
        else {
            DE<-append(DE, "None")
        }
    }
    return(DE)
}


## Compare t-stats of genes from different DGE analyses
t_stat_plot <- function(t_stats_1, t_stats_2, name_1, name_2, brain_region){

    ## Spearman correlation coeff
    rho <- cor(t_stats_1$t, t_stats_2$t, method = "spearman")
    rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))

    ## Colors and transparency
    cols <- c("darkorange3", "palegreen3","orchid2", "darkgrey")
    names(cols) <- c("sig Both", paste0("sig ", name_1), paste0("sig ", name_2), "None")
    alphas <- c( 1, 1, 1, 0.5)
    names(alphas) <- names(cols)

    ## Shared genes to show
    up_shared <- c("Dio2", "Kiaa0408L", "Hspa5", "Ap3s1", "Kcnj3")
    down_shared <- c("Ccdc187", "Itpkb", "Flnc", "Itgal", "Sox8")

    ## Merge t-stats
    t_stats<-data.frame(symbol = t_stats_1$Symbol, t1=t_stats_1$t, t2=t_stats_2$t)
    ## Add DE info for the genes in both DEAs
    t_stats$DEG<-add_DE_info(t_stats_1, t_stats_2, name_1, name_2)
    t_stats$DEG <- factor(t_stats$DEG, levels=names(cols))

    plot <- ggplot(t_stats, aes(x = t1, y = t2, color=DEG, alpha=DEG)) +
        geom_point(size = 1) +
        scale_color_manual(values = cols) +
        scale_alpha_manual(values = alphas) +
        labs(x = paste("t-stats", name_1),
             y = paste("t-stats", name_2),
             subtitle = rho_anno,
             color = "Differential expression",
             parse = T) +
        guides(alpha = 'none', color = guide_legend(override.aes = list(size=1.3))) +
        theme_bw() +
        theme(plot.margin = unit(c(1,1,1,1), "cm"),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.text = element_text(size=11),
              legend.title = element_text(size=12))


    if(brain_region == "hab_amy"){
        plot <- plot + geom_text_repel(data = subset(t_stats, symbol %in% c(up_shared, down_shared)),
                                       aes(label = symbol),
                                       size=3,
                                       color='black',
                                       alpha = 1,
                                       max.overlaps = Inf,
                                       box.padding = 0.15,
                                       segment.size = unit(0.35, 'mm'),
                                       segment.alpha = 0.4,
                                       show.legend=FALSE)
    }
    plot

    name1 <- gsub(' ', '_', name_1)
    name2 <- gsub(' ', '_', name_2)
    ggsave(paste0('plots/05_DEA/02_Comparisons/t_stats_', name1, '_vs_', name2, '_', brain_region, '.pdf'), width = 6.6, height = 4.5)
}



## 2.1 Comparison of gene DE signal for substance and behavior in habenula

## DE t-stats of the genes for substance and behavioral covariates
t_stats_Substance_habenula <- results_Substance_uncorr_vars_habenula[[1]]
t_stats_FirstHrIntakeSlope_habenula <- results_FirstHrIntakeSlope_habenula[[1]]
t_stats_TotalIntake_habenula <- results_TotalIntake_habenula[[1]]
t_stats_LastSessionIntake_habenula <- results_LastSessionIntake_habenula[[1]]

## Df with gene metadata and DE statistics
df_habenula <- cbind(t_stats_Substance_habenula[,c("seqnames", "start", "end", "width", "strand", "Length", "ensemblID",
                                                   "EntrezID", "Symbol", "meanExprs", "logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_FirstHrIntakeSlope_habenula[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_TotalIntake_habenula[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_LastSessionIntake_habenula[,c("logFC", "t", "P.Value", "adj.P.Val")])

colnames(df_habenula)[c(1,11:26)] <- c("chr", paste(c("logFC", "t", "P.Value", "adj.P.Val"),
                                           rep(c('Substance', 'FirstHrIntakeSlope', 'TotalIntake', 'LastSessionIntake'), c(4,4,4,4)),
                                           sep='_'))
df_habenula$EntrezID <- as.character(df_habenula$EntrezID)

## Pairwise correlations between t-stats
ggpairs(df_habenula, columns = c("t_Substance", "t_FirstHrIntakeSlope", "t_TotalIntake", "t_LastSessionIntake")) + theme_bw()
ggsave(here('plots/05_DEA/02_Comparisons/t_stats_pairs_habenula.pdf'))


## Compare DE signal for substance vs behavior in habenula
t_stat_plot(t_stats_Substance_habenula, t_stats_FirstHrIntakeSlope_habenula,
            'Substance', 'First hr infusion slope', 'habenula')
t_stat_plot(t_stats_Substance_habenula, t_stats_TotalIntake_habenula,
            'Substance', 'Total intake', 'habenula')
t_stat_plot(t_stats_Substance_habenula, t_stats_LastSessionIntake_habenula,
            'Substance', 'Last session Intake', 'habenula')

## Compare DE signal for behavioral covariates in habenula
t_stat_plot(t_stats_FirstHrIntakeSlope_habenula, t_stats_TotalIntake_habenula,
            'First hr infusion slope', 'Total intake', 'habenula')
t_stat_plot(t_stats_FirstHrIntakeSlope_habenula, t_stats_LastSessionIntake_habenula,
            'First hr infusion slope', 'Last session Intake', 'habenula')
t_stat_plot(t_stats_TotalIntake_habenula, t_stats_LastSessionIntake_habenula,
            'Total intake', 'Last session Intake', 'habenula')


#===============================================================================
#                  Rank-rank hypergeometric test (RRHO2)
#===============================================================================
## Create lists of genes with signed pval (-log10(p)*sign(logFC))
gene_list_Substance_Hb <- data.frame(Genes = t_stats_Substance_habenula$ensemblID,
                                     DDE = -log10(t_stats_Substance_habenula$P.Value)*
                                          sign(t_stats_Substance_habenula$logFC),
                                     stringsAsFactors = FALSE)

gene_list_FirstHrIntakeSlope_Hb <- data.frame(Genes = t_stats_FirstHrIntakeSlope_habenula$ensemblID,
                                              DDE = -log10(t_stats_FirstHrIntakeSlope_habenula$P.Value)*
                                                     sign(t_stats_FirstHrIntakeSlope_habenula$logFC),
                                              stringsAsFactors = FALSE)

gene_list_TotalIntake_Hb <- data.frame(Genes = t_stats_TotalIntake_habenula$ensemblID,
                                       DDE = -log10(t_stats_TotalIntake_habenula$P.Value)*
                                              sign(t_stats_TotalIntake_habenula$logFC),
                                       stringsAsFactors = FALSE)

gene_list_LastSessionIntake_Hb <- data.frame(Genes = t_stats_LastSessionIntake_habenula$ensemblID,
                                             DDE = -log10(t_stats_LastSessionIntake_habenula$P.Value)*
                                                    sign(t_stats_LastSessionIntake_habenula$logFC),
                                             stringsAsFactors = FALSE)

#-------------------------------------------------------------------------------
## Create RRHO2 object and run their function
RRHO_obj <-  RRHO2_initialize(gene_list_Substance_Hb, gene_list_FirstHrIntakeSlope_Hb,
                              labels = c("Substance, Hb", "First Hour Intake Slope, Hb"),
                              log10.ind = TRUE, stepsize = stepsize)
## Visualize the heatmap
RRHO2_heatmap(RRHO_obj)
mat <- RRHO_obj$hypermat
#-------------------------------------------------------------------------------

## Calculate overlap signif manually based on their source code in:
## https://rdrr.io/github/RRHO2/RRHO2/src/R/RRHO2_initialize.R.
## I don't trust their results!

## Compute overlap signif between each pair top gene sets in condition 1 and 2
overlap_hyper <- function(a,b) {

    count<-as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b])))
    log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)

    return(c(counts=count,
             log.pval=as.numeric(log.pval)
    ))
}

## Compute overlap signif between each pair top gene sets in condition 1 and 2,
## inverting ranking of genes in condition 1 (to assess enrichment of top + genes in 2 among top - in 1)
overlap_hyper_inverted_sample1 <- function(a,b) {

    count<-as.integer(sum(as.numeric(rev(sample1)[1:a] %in% sample2[1:b])))
    log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)

    return(c(counts=count,
             log.pval=as.numeric(log.pval)
    ))
}

## Manual implementation of RRHO
RRHO_manual <- function(list1, list2, title1, title2, filename){

    ## Rank genes by signed pvalue
    list1 <- list1[order(list1[, 2], decreasing = TRUE), ]
    list2 <- list2[order(list2[, 2], decreasing = TRUE), ]

    nlist1 <- length(list1[, 1])
    nlist2 <- length(list2[, 1])

    sample1 = list1[, 1]
    sample2 = list2[, 1]
    stepsize = ceiling(sqrt(dim(list1)[1]))

    ## Gene universe size
    n <- length(sample1)

    ## Assess enrichment among each pair of top gene sets
    indexes<- expand.grid(i=seq(1,n,by=stepsize), j=seq(1,n,by=stepsize))
    overlaps<- apply(indexes, 1, function(x) c(x['i'], x['j'], overlap_hyper(x['i'], x['j'])))
    ov <- as.data.frame(t(overlaps))

    ## Signif overlaps
    head(ov[which(ov$log.pval>2.995732), ])
    #         i     j counts log.pval signs
    # 14838 261 14951    243 3.538316     1
    # 14967 261 15081    246 4.286130     1
    tail(ov[which(ov$log.pval>2.995732), ])
    #         i     j counts log.pval signs
    # 16257 261 16381    260 3.373324     1
    # 16258 391 16381    389 4.131022     1

    ## Martrix with -log(pvals)
    m <- matrix(ov$log.pval, ncol = sqrt(dim(indexes)[1]), byrow = F)
    rownames(m) <- colnames(m) <- 1:sqrt(dim(indexes)[1])

    ## Assess enrichment among each pair of top gene sets with inverted gene ranking in condition 1
    overlaps_rev <- apply(indexes, 1, function(x) c(x['i'], x['j'], overlap_hyper_inverted_sample1(x['i'], x['j'])))
    ov_rev <- as.data.frame(t(overlaps_rev))
    m_rev <- matrix(ov_rev$log.pval, ncol = sqrt(dim(indexes)[1]), byrow = F)
    rownames(m_rev) <- colnames(m_rev) <- 1:sqrt(dim(indexes)[1])

    ## Signif enrichments
    dim(ov_rev[which(ov_rev$log.pval>2.995732), ])
    # [1] 15610     5

    ## Final matrix with enrich pvals based on quadrant
    stepList1 <- seq(1, nlist1, stepsize)
    stepList2 <- seq(1, nlist2, stepsize)

    len1 <- length(stepList1)
    len2 <- length(stepList2)

    boundary1 <- sum(list1[stepList1,2] > 0)
    boundary2 <- sum(list2[stepList2,2] > 0)

    hypermat <- matrix(NA, nrow=nrow(m), ncol=ncol(m))

    ## Top + genes in 1 among Top + genes in 2
    hypermat[1:boundary1, 1:boundary2] <- m[1:boundary1,1:boundary2]
    ## Top - genes in 1 among Top - genes in 2
    hypermat[(boundary1+1):len1, (boundary2+1):len2] <- m[(boundary1+1):len1, (boundary2+1):len2]
    ## Top + genes in 1 among Top - genes in 2
    hypermat[1:boundary1, (boundary2+1):len2] <- m_rev[len1:(len1 - boundary1 + 1), (boundary2+1):len2]
    ## Top - genes in 1 among Top + genes in 2
    hypermat[(boundary1+1):len1, 1:boundary2] <- m_rev[(len1 - boundary1):1, 1:boundary2]

    ## Plot heatmap
    col_fun <- colorRamp2(
        breaks = c(min(hypermat), quantile(hypermat, 1:9/10), max(hypermat)),
        colors = c("#FFFFF0", "lightyellow", "#FFEC8B", "#FFD700", "#FFB90F",
                   "#FFA500", "#FF7B00FF", "#FF7B00FF", "#FF0000FF", "#CD0000",
                   "#8B1A1A"))

    ## Labels for top up and down genes contrasted
    labs1 <- c(paste0("top up ", seq(1, nlist1, stepsize)[1:boundary1]),
               paste0("top down ", N - (seq(1, nlist1, stepsize)[(boundary1+1):nrow(hypermat)])))

    labs2 <- c(paste0("top up ", seq(1, nlist2, stepsize)[1:boundary2]),
               paste0("top down ", N - (seq(1, nlist1, stepsize)[(boundary2+1):ncol(hypermat)])))

    tolabel1 <- sort(unique(c(1, which(seq_len(nrow(hypermat)) %% 5 == 0), nrow(hypermat))))

    tolabel2 <- sort(unique(c(1, which(seq_len(ncol(hypermat)) %% 5 == 0), ncol(hypermat))))

    labs1_some <- rep("", nrow(hypermat))
    labs2_some <- rep("", ncol(hypermat))

    labs1_some[tolabel1] <- labs1[tolabel1]
    labs2_some[tolabel2] <- labs2[tolabel2]


    # Draw heatmap
    h <- Heatmap(
        hypermat,
        name = "-log(p)",
        col = col_fun,

        # preserve order
        cluster_rows = FALSE,
        cluster_columns = FALSE,

        # cleaner labels
        show_row_names = T,
        show_column_names = T,

        # labels every 10
        row_labels = labs1_some,
        column_labels = labs2_some,

        row_names_gp = grid::gpar(fontsize = 6),
        column_names_gp = grid::gpar(fontsize = 6),

        # aesthetics
        border = TRUE,
        na_col = "grey95", row_names_side = "left",

        heatmap_legend_param = list(
            at = c(0, signif(max(hypermat)/2, 2),
                   # signif(median(hypermat), 2),
                   signif(max(hypermat), 2)),
            title = "-log(p)",
            legend_height = unit(4, "cm")
        )
    )

    png(paste0(here("plots/05_DEA/02_Comparisons/RRHO_"), filename,  ".png"),
        width = 1500, height = 1300, res = 300)

    pdf(paste0(here("plots/05_DEA/02_Comparisons/RRHO_"), filename,  ".pdf"),
        width = 4, height = 3.5)

    draw(h)

    # Add grid elements
    grid.text(title2,
              x = unit(0.5, "npc"),
              y = unit(0, "npc") + unit(4, "mm"),
              gp = gpar(fontsize = 10, fontface = "bold"))

    grid.text(title1,
              x = unit(0, "mm"),
              y = unit(0.5, "npc"),
              rot = 90,
              gp = gpar(fontsize = 10, fontface = "bold"))

    # Close device
    dev.off()

}

## Run analysis
RRHO_manual(list1 = gene_list_Substance_Hb, list2 = gene_list_LastSessionIntake_Hb,
            title1 = "Substance Hb", title2 = "Last session intake Hb",
            filename = "Substance_vs_Last_Session_Intake_habenula")

RRHO_manual(list1 = gene_list_Substance_Hb, list2 = gene_list_TotalIntake_Hb,
            title1 = "Substance Hb", title2 = "Total intake Hb",
            filename = "Substance_vs_Total_Intake_habenula")

RRHO_manual(list1 = gene_list_Substance_Hb, list2 = gene_list_FirstHrIntakeSlope_Hb,
            title1 = "Substance Hb", title2 = "First hr intake slope Hb",
            filename = "Substance_vs_First_hr_infusion_slope_habenula")

RRHO_manual(list1 = gene_list_TotalIntake_Hb, list2 = gene_list_LastSessionIntake_Hb,
            title1 = "Total intake Hb", title2 = "Last session intake Hb",
            filename = "Total_Intake_vs_Last_Session_Intake_habenula")

RRHO_manual(list1 = gene_list_TotalIntake_Hb, list2 = gene_list_FirstHrIntakeSlope_Hb,
            title1 = "Total intake Hb", title2 = "First hr intake slope Hb",
            filename = "Total_Intake_vs_First_hr_infusion_slope_habenula")

RRHO_manual(list1 = gene_list_FirstHrIntakeSlope_Hb, list2 = gene_list_LastSessionIntake_Hb,
            title1 = "First hr intake slope Hb",  title2 = "Last session intake Hb",
            filename = "First_hr_infusion_slope_vs_Last_Session_Intake_habenula")
#===============================================================================



## 2.2 Comparison of gene DE signal for substance and behavior in amygdala

t_stats_Substance_amygdala <- results_Substance_uncorr_vars_amygdala[[1]]
t_stats_FirstHrIntakeSlope_amygdala <- results_FirstHrIntakeSlope_amygdala[[1]]
t_stats_TotalIntake_amygdala <- results_TotalIntake_amygdala[[1]]
t_stats_LastSessionIntake_amygdala <- results_LastSessionIntake_amygdala[[1]]

df_amygdala <- cbind(t_stats_Substance_amygdala[,c("seqnames", "start", "end", "width", "strand", "Length", "ensemblID",
                                                   "EntrezID", "Symbol", "meanExprs", "logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_FirstHrIntakeSlope_amygdala[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_TotalIntake_amygdala[,c("logFC", "t", "P.Value", "adj.P.Val")],
                     t_stats_LastSessionIntake_amygdala[,c("logFC", "t", "P.Value", "adj.P.Val")])

colnames(df_amygdala)[c(1,11:26)]<- c("chr", paste(c("logFC", "t", "P.Value", "adj.P.Val"),
                                                  rep(c('Substance', 'FirstHrIntakeSlope', 'TotalIntake', 'LastSessionIntake'),
                                                      c(4,4,4,4)),sep='_'))
df_amygdala$EntrezID <- as.character(df_amygdala$EntrezID)

## Pairwise correlations between t-stats
ggpairs(df_amygdala, columns = c("t_Substance", "t_FirstHrIntakeSlope", "t_TotalIntake", "t_LastSessionIntake")) + theme_bw()
ggsave(here('plots/05_DEA/02_Comparisons/t_stats_pairs_amygdala.pdf'))


## t-stats comparisons in amygdala
t_stat_plot(t_stats_Substance_amygdala, t_stats_FirstHrIntakeSlope_amygdala,
            'Substance', 'First hr infusion slope', 'amygdala')
t_stat_plot(t_stats_Substance_amygdala, t_stats_TotalIntake_amygdala,
            'Substance', 'Total intake', 'amygdala')
t_stat_plot(t_stats_Substance_amygdala, t_stats_LastSessionIntake_amygdala,
            'Substance', 'Last session Intake', 'amygdala')
t_stat_plot(t_stats_FirstHrIntakeSlope_amygdala, t_stats_TotalIntake_amygdala,
            'First hr infusion slope', 'Total intake', 'amygdala')
t_stat_plot(t_stats_FirstHrIntakeSlope_amygdala, t_stats_LastSessionIntake_amygdala,
            'First hr infusion slope', 'Last session Intake', 'amygdala')
t_stat_plot(t_stats_TotalIntake_amygdala, t_stats_LastSessionIntake_amygdala,
            'Total intake', 'Last session Intake', 'amygdala')


## RRHO analysis:
gene_list_Substance_Amy <- data.frame(Genes = t_stats_Substance_amygdala$ensemblID,
                                      DDE = -log10(t_stats_Substance_amygdala$P.Value)*
                                            sign(t_stats_Substance_amygdala$logFC),
                                            stringsAsFactors = FALSE)

gene_list_FirstHrIntakeSlope_Amy <- data.frame(Genes = t_stats_FirstHrIntakeSlope_amygdala$ensemblID,
                                               DDE = -log10(t_stats_FirstHrIntakeSlope_amygdala$P.Value)*
                                                   sign(t_stats_FirstHrIntakeSlope_amygdala$logFC),
                                                   stringsAsFactors = FALSE)

gene_list_TotalIntake_Amy <- data.frame(Genes = t_stats_TotalIntake_amygdala$ensemblID,
                                        DDE = -log10(t_stats_TotalIntake_amygdala$P.Value)*
                                           sign(t_stats_TotalIntake_amygdala$logFC),
                                            stringsAsFactors = FALSE)

gene_list_LastSessionIntake_Amy <- data.frame(Genes = t_stats_LastSessionIntake_amygdala$ensemblID,
                                              DDE = -log10(t_stats_LastSessionIntake_amygdala$P.Value)*
                                                 sign(t_stats_LastSessionIntake_amygdala$logFC),
                                                 stringsAsFactors = FALSE)

RRHO_manual(list1 = gene_list_Substance_Amy, list2 = gene_list_LastSessionIntake_Amy,
            title1 = "Substance Amyg", title2 = "Last session intake Amyg",
            filename = "Substance_vs_Last_Session_Intake_amygdala")

RRHO_manual(list1 = gene_list_Substance_Amy, list2 = gene_list_TotalIntake_Amy,
            title1 = "Substance Amyg", title2 = "Total intake Amyg",
            filename = "Substance_vs_Total_Intake_amygdala")

RRHO_manual(list1 = gene_list_Substance_Amy, list2 = gene_list_FirstHrIntakeSlope_Amy,
            title1 = "Substance Amyg", title2 = "First hr intake slope Amyg",
            filename = "Substance_vs_First_hr_infusion_slope_amygdala")

RRHO_manual(list1 = gene_list_TotalIntake_Amy, list2 = gene_list_LastSessionIntake_Amy,
            title1 = "Total intake Amyg", title2 = "Last session intake Amyg",
            filename = "Total_Intake_vs_Last_Session_Intake_amygdala")

RRHO_manual(list1 = gene_list_TotalIntake_Amy, list2 = gene_list_FirstHrIntakeSlope_Amy,
            title1 = "Total intake Amyg", title2 = "First hr intake slope Amyg",
            filename = "Total_Intake_vs_First_hr_infusion_slope_amygdala")

RRHO_manual(list1 = gene_list_FirstHrIntakeSlope_Amy, list2 = gene_list_LastSessionIntake_Amy,
            title1 = "First hr intake slope Amyg",  title2 = "Last session intake Amyg",
            filename = "First_hr_infusion_slope_vs_Last_Session_Intake_amygdala")



## 2.3 Comparison of gene DE signal for substance and behavior in habenula and amygdala

colnames(df_habenula)[11:26] <- paste0(colnames(df_habenula)[11:26], '_habenula')
colnames(df_amygdala)[11:26] <- paste0(colnames(df_amygdala)[11:26], '_amygdala')
df_DEstats_complete <- cbind(df_habenula, df_amygdala[,11:26])

## Create supp table
save(df_DEstats_complete, file = "processed-data/05_DEA/DEAs_results_all_genes_hab_amyg.Rdata")
write.table(df_DEstats_complete, "processed-data/Supplementary_Tables/TableS10_DEAs_results_all_genes_hab_amyg.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')

## Scatter plots
ggpairs(df_DEstats_complete,
        columns = paste0(c("t_Substance", "t_FirstHrIntakeSlope", "t_TotalIntake", "t_LastSessionIntake"), rep(c('_habenula', '_amygdala'), c(4,4)))) + theme_bw() +
    theme(strip.text.x = element_text(size = 10),
          strip.text.y = element_text(size = 10))
ggsave(here('plots/05_DEA/02_Comparisons/t_stats_pairs_hab_amyg.pdf'), height = 18, width = 20)


## Number of common habenula and amygdala DEGs for substance
common_genes <- intersect(de_genes_habenula$ensemblID, de_genes_amygdala$ensemblID)
common_DEGs_hab_amyg <- de_genes_habenula[which(de_genes_habenula$ensemblID %in% common_genes), ]
dim(common_DEGs_hab_amyg)
# [1] 106  19

## Data frame with common genes
write.table(common_DEGs_hab_amyg, "processed-data/05_DEA/common_DEGs_hab_amyg.tsv", row.names = FALSE, col.names = TRUE, sep = '\t')

## Compare DE signal for substance in habenula vs amygdala
t_stat_plot(t_stats_Substance_habenula, t_stats_Substance_amygdala,
            'Substance habenula', 'Substance amygdala', 'hab_amy')

## RRHO analysis:
RRHO_manual(list1 = gene_list_Substance_Hb, list2 = gene_list_Substance_Amy,
            title1 = "Substance Hb", title2 = "Substance Amyg",
            filename = "Substance_habenula_vs_Substance_amygdala")






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
# date     2024-04-23
# rstudio  2023.12.1+402 Ocean Storm (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
# Biobase              * 2.62.0    2023-10-26 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-02 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# cowplot              * 1.1.3     2024-01-22 [1] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-11-06 [1] Bioconductor
# digest                 0.6.34    2024-01-11 [1] CRAN (R 4.3.1)
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.1)
# evaluate               0.23      2023-11-01 [1] CRAN (R 4.3.1)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
# fastmap                1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.6    2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11    2024-02-17 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-30 [1] Bioconductor
# GGally               * 2.2.1     2024-02-14 [1] CRAN (R 4.3.1)
# ggplot2              * 3.5.0     2024-02-23 [1] CRAN (R 4.3.1)
# ggstats                0.6.0     2024-04-05 [1] CRAN (R 4.3.1)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.0)
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                    1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmltools              0.5.7     2023-11-03 [1] CRAN (R 4.3.1)
# IRanges              * 2.36.0    2023-10-26 [1] Bioconductor
# knitr                  1.45      2023-10-30 [1] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.0)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.1)
# limma                * 3.58.1    2023-11-02 [1] Bioconductor
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-5     2024-01-11 [1] CRAN (R 4.3.1)
# MatrixGenerics       * 1.14.0    2023-10-26 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.3.1)
# progress               1.2.3     2023-12-06 [1] CRAN (R 4.3.1)
# purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.1)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.1)
# RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.1)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.1)
# rmarkdown              2.26      2024-03-05 [1] CRAN (R 4.3.1)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.1)
# rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# S4Arrays               1.2.0     2023-10-26 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.1)
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# SparseArray            1.2.4     2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.0)
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
# XVector                0.42.0    2023-10-26 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-26 [1] Bioconductor
#
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
