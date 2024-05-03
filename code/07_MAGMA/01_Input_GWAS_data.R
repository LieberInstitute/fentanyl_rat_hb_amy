
library(here)
library(data.table)
library(biomaRt)
library(sessioninfo)
library(org.Hs.eg.db)


####################   Generalized Gene-Set Analysis of GWAS Data (MAGMA)   ######################

# Note: this code is based on https://github.com/LieberInstitute/Habenula_Pilot/tree/master/code/13_MAGMA
# Input files were originally generated there by Louise Huuki-Myers (@lahuuki)


################################################################################
##                       1. GWAS input data preparation
################################################################################


## Function to extract positions of human SNPs if missing in GWAS data
## (Not needed)
SNP_loc_extraction <- function(snp_list){

    ## Divide SNP list into 1000 groups
    snp_lists <- split(snp_list, cut(seq_along(snp_list), 1000, labels=FALSE))

    output <- vector()
    for (snps in snp_lists){
        ## Mart
        mart <- useEnsembl(biomart="ENSEMBL_MART_SNP",
                           host="https://grch37.ensembl.org",  ## Use Grch37 assembly
                           dataset="hsapiens_snp") ## human SNPs

        message(Sys.time(), " - Searching SNP position")
        ## Obtain ID, chr, and start position of SNPs
        snploc <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start'),
                        values = snps,
                        filters ="snp_filter",
                        mart = mart)
        snploc <- as.data.frame(snploc)
        colnames(snploc) <- c('SNP', 'CHR', 'BP')

        output <- rbind(output, snploc)
    }

    return(output)
}


# ------------------------------------------------------------------------------
#                                 GWAS SZC data
# ------------------------------------------------------------------------------

#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

## Internal paths (data is in JHPCE)
snploc_szc = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.snploc")))

dim(snploc_szc)
# [1] 7659767      3

## Explore data
head(snploc_szc)
#          SNP CHR        BP
# 1 rs62513865   8 101592213
# 2 rs79643588   8 106973048
# 3 rs17396518   8 108690829
# 4   rs983166   8 108681675
# 5 rs28842593   8 103044620
# 6  rs7014597   8 104152280

## Confirm no repearted variant IDs
which(duplicated(snploc_szc$ID))
# integer(0)

## Check no variants in sexual chrs
unique(snploc_szc$CHR)
# [1]  8  2 16  7  6 15 14 12 17  1  3 11 22 20 21  9  4 13  5 10 18 19

## Remove NA
snploc_szc <- subset(snploc_szc, !is.na(SNP) & !is.na(CHR) & !is.na(BP))
dim(snploc_szc)
# [1] 7659767      3

## Save
save(snploc_szc, file=here('processed-data/07_MAGMA/Input_GWAS_data/SCZ/SCZ_PGC3_wave3.european.autosome.public.v3.snploc'))


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

snp_pval_szc = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.pval")))

dim(snp_pval_szc)
# [1] 7659767      3

## Explore data
head(snp_pval_szc)
#          SNP      P      N
# 1 rs62513865 0.4847 130644
# 2 rs79643588 0.5605 130644
# 3 rs17396518 0.8145 130644
# 4   rs983166 0.5704 130644
# 5 rs28842593 0.7488 130644
# 6  rs7014597 0.5034 130644

## Check we have the same SNPs in snploc and snp_pval files
setdiff(snploc_szc$SNP, snp_pval_szc$SNP)
# character(0)

## Remove NA
snp_pval_szc <- subset(snp_pval_szc, !is.na(SNP) & !is.na(P) & !is.na(N))
dim(snp_pval_szc)
# [1] 7659767      3

## Save
save(snp_pval_szc, file='processed-data/07_MAGMA/Input_GWAS_data/SCZ/SCZ_PGC3_wave3.european.autosome.public.v3.pval')



# ------------------------------------------------------------------------------
#                                 GWAS MDD data
# ------------------------------------------------------------------------------

#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

snploc_mdd2019 = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/mdd2019edinburgh/PGC_UKB_depression_genome-wide.snploc")))

dim(snploc_mdd2019)
# [1] 8481694       3

head(snploc_mdd2019)
#          SNP CHR        BP
# 1  rs2326918   6 130840091
# 2  rs7929618  11 134897743
# 3 rs66941928   3 176666749
# 4  rs7190157  16   8600861
# 5 rs12364336  11 100009976
# 6  rs6977693   7 145771806

## Remove repeated SNPs
rep_SNPs <- snploc_mdd2019[which(duplicated(snploc_mdd2019$SNP)), 'SNP']
snploc_mdd2019[which(snploc_mdd2019$SNP %in% rep_SNPs),]
#                 SNP CHR        BP
# 1929499 rs188230556   7  45554166
# 2254826 rs188230556   7  45554166
# 3634584 rs116109352   2 213600534
# 4332918  rs75488681  12 109823850
# 6044110  rs75488681  12 109823850
# 6199830 rs116109352   2 213600534

snploc_mdd2019 <- snploc_mdd2019[-which(duplicated(snploc_mdd2019$SNP)), ]
dim(snploc_mdd2019)
# [1] 8481691       3

## Check no variants in sexual chrs
unique(snploc_mdd2019$CHR)
# [1]  6 11  3 16  7  1 12 14  2 17 10 18 13  9  5  4 21  8 19 15 20 22

## Remove NA
snploc_mdd2019 <- subset(snploc_mdd2019, !is.na(SNP) & !is.na(CHR) & !is.na(BP))
dim(snploc_mdd2019)
# [1] 8481691       3


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

snp_pval_mdd2019 = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/mdd2019edinburgh/PGC_UKB_depression_genome-wide.pval")))

dim(snp_pval_mdd2019)
# [1] 8483301       2

## Explore data
head(snp_pval_mdd2019)
#          SNP         P
# 1  rs2326918 0.0756100
# 2  rs7929618 0.0004804
# 3 rs66941928 0.9502000
# 4  rs7190157 0.5992000
# 5 rs12364336 0.2450000
# 6  rs6977693 0.1442000

## No repeated rsIDs
which(duplicated(snp_pval_mdd2019$SNP))
# integer(0)

## Subset to common SNPs in snploc and snp_pval files
snploc_mdd2019 <- snploc_mdd2019[which(snploc_mdd2019$SNP %in% snp_pval_mdd2019$SNP), ]
snp_pval_mdd2019 <- snp_pval_mdd2019[which(snp_pval_mdd2019$SNP %in% snploc_mdd2019$SNP),]
dim(snploc_mdd2019)
# [1] 8481690       3
dim(snp_pval_mdd2019)
# [1] 8481690       3

setdiff(snploc_mdd2019$SNP, snp_pval_mdd2019$SNP)
# character(0)

## Remove NA
snp_pval_mdd2019 <- subset(snp_pval_mdd2019, !is.na(SNP) & !is.na(P))
dim(snp_pval_mdd2019)
# [1] 8481690       2

## Save
save(snploc_mdd2019, file=here('processed-data/07_MAGMA/Input_GWAS_data/MDD_2019/MDD_PGC_UKB_depression_genome-wide.snploc'))
save(snp_pval_mdd2019, file="processed-data/07_MAGMA/Input_GWAS_data/MDD_2019/MDD_PGC_UKB_depression_genome-wide.pval")



# ------------------------------------------------------------------------------
#                                 GWAS Panic data
# ------------------------------------------------------------------------------

#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

snploc_panic = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/panic2019/pgc-panic2019.snploc")))

dim(snploc_panic)
# [1] 10151624        3

head(snploc_panic)
#               SNP CHR      BP
# 1      rs11250701  10 1689546
# 2 chr10_2622752_D  10 2622752
# 3       rs7085086  10  151476
# 4     rs113494187  10 1593759
# 5     rs117915320  10 1708106
# 6     rs182753344  10  790310

## Delete not valid CHR entries
chrs <- as.character(1:22)
snploc_panic <- snploc_panic[which(snploc_panic$CHR %in% chrs), ]

unique(snploc_panic$CHR)
# [1] "10" "1"  "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "2"  "21" "22"
# [16] "3"  "4"  "5"  "6"  "7"  "8"  "9"
dim(snploc_panic)
# [1] 10151300        3


## Only unique SNPs
rep_SNPs <- snploc_panic[which(duplicated(snploc_panic$SNP)), 'SNP']
snploc_panic[which(snploc_panic$SNP %in% rep_SNPs),]
#                      SNP CHR        BP
# 58433         rs74733400  10  12000000
# 68879         rs74733400  10  12000000
# 823717       rs188690587   1  93000000
# 833311       rs188690587   1  93000000
# 968325  chr1_105000000_I   1 105000000
# 1081973 chr1_105000000_I   1 105000000
# 1801583         rs855274  12  48000000
# 1811723         rs855274  12  48000000
# 2879007 chr14_81000000_D  14  81000000
# 2889922 chr14_81000000_D  14  81000000
# 4905653       rs61573637   2  72000000
# 4914297       rs61573637   2  72000000

snploc_panic <- snploc_panic[-which(duplicated(snploc_panic$SNP)), ]
dim(snploc_panic)
# [1] 10151294        3

which(duplicated(snploc_panic$SNP))
# integer(0)


## Remove NA
snploc_panic <- subset(snploc_panic, !is.na(SNP) & !is.na(CHR) & !is.na(BP))
dim(snploc_panic)
# [1] 10151294        3

## Save
save(snploc_panic, file=here('processed-data/07_MAGMA/Input_GWAS_data/Panic/Panic_PGC_panic2019.snploc'))


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

snp_pval_panic = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/panic2019/pgc-panic2019.pval")))

dim(snp_pval_panic)
# [1] 10151300        3

head(snp_pval_panic)
#               SNP       P    N
# 1      rs11250701 0.53280 9907
# 2 chr10_2622752_D 0.26150 9907
# 3       rs7085086 0.28970 9907
# 4     rs113494187 0.05117 9907
# 5     rs117915320 0.10030 9907
# 6     rs182753344 0.74840 9907

## Keep duplicates
# ----------------------------------------------------------------|
#    Caution: p-values are not the same for repeated SNP IDs!!!   |
#    They may represent different variants in the same position   |
# ----------------------------------------------------------------|
rep_SNPs <- snp_pval_panic[which(duplicated(snp_pval_panic$SNP)), 'SNP']
snp_pval_panic[which(snp_pval_panic$SNP %in% rep_SNPs),]
#                      SNP       P    N
# 58433         rs74733400 0.85430 9907
# 68879         rs74733400 0.89080 9907
# 823717       rs188690587 0.39910 9907
# 833311       rs188690587 0.32750 9907
# 968325  chr1_105000000_I 0.64240 9907
# 1081973 chr1_105000000_I 0.63300 9907
# 1801583         rs855274 0.88540 9907
# 1811723         rs855274 0.92810 9907
# 2879007 chr14_81000000_D 0.67070 9907
# 2889922 chr14_81000000_D 0.88160 9907
# 4905653       rs61573637 0.04726 9907
# 4914297       rs61573637 0.06132 9907

## Same SNPs
setdiff(snploc_panic$SNP, snp_pval_panic$SNP)
# character(0)

## Remove NA
snp_pval_panic <- subset(snp_pval_panic, !is.na(SNP) & !is.na(P))
dim(snp_pval_panic)
# [1] 10151300        3

## Save
save(snp_pval_panic, file="processed-data/07_MAGMA/Input_GWAS_data/Panic/Panic_PGC_panic2019.pval")



# ------------------------------------------------------------------------------
#                                 GWAS SUD data
# ------------------------------------------------------------------------------

#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

snploc_SUD = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/sud2020op/opi.DEPvEXP_EUR.noAF.snploc")))

dim(snploc_SUD)
# [1] 4187212       3

head(snploc_SUD)
#          SNP CHR        BP
# 1 rs10868284   9  87977652
# 2 rs62291883   3 185313842
# 3 rs61956327  12  55828708
# 4 rs60994383  17  56247101
# 5 rs12531896   7 115132342
# 6 rs35515951   3  97907440

## Remove repeated SNPs
rep_SNPs <- snploc_SUD[which(duplicated(snploc_SUD$SNP)), 'SNP']
head(snploc_SUD[which(snploc_SUD$SNP %in% rep_SNPs),])
#                SNP CHR        BP
# 64312   rs58887011   4  25086197
# 95726   rs58996925   1  56267033
# 210428 rs112812624   6  91998732
# 416482  rs58644696   1  41350292
# 497789 rs114379693   1 165629897
# 601798 rs112957741  13  57866839

snploc_SUD <- snploc_SUD[-which(duplicated(snploc_SUD$SNP)), ]
dim(snploc_SUD)
# [1] 4187189       3

## Check no variants in sexual chrs
unique(snploc_SUD$CHR)
# [1]  9  3 12 17  7  4  6  2  8 10  5 13  1 11 18 22 15 19 20 14 21 16

## Remove NAs
snploc_SUD <- subset(snploc_SUD, !is.na(SNP) & !is.na(CHR) & !is.na(BP))
dim(snploc_SUD)
# [1] 4187189       3


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

snp_pval_SUD = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/sud2020op/opi.DEPvEXP_EUR.noAF.pval")))

dim(snp_pval_SUD)
# [1] 4211587       3

head(snp_pval_SUD)
#          SNP      P    N
# 1 rs10868284 0.1075 5709
# 2 rs62291883 0.1134 6072
# 3 rs61956327 0.1458 5260
# 4 rs60994383 0.2846 5140
# 5 rs12531896 0.1857 6115
# 6 rs35515951 0.2337 6062

## Unique SNPs
which(duplicated(snp_pval_SUD$SNP))
# integer(0)

## Subset to common SNPs in snploc and snp_pval
snploc_SUD <- snploc_SUD[which(snploc_SUD$SNP %in% snp_pval_SUD$SNP), ]
snp_pval_SUD <- snp_pval_SUD[which(snp_pval_SUD$SNP %in% snploc_SUD$SNP),]
dim(snploc_SUD)
# [1] 4187118       3
dim(snp_pval_SUD)
# [1] 4187118       3

setdiff(snploc_SUD$SNP, snp_pval_SUD$SNP)
# character(0)

## Remove NA
snp_pval_SUD <- subset(snp_pval_SUD, !is.na(SNP) & !is.na(P))
dim(snp_pval_SUD)
# [1] 4187118       3

## Save
save(snploc_SUD, file=here('processed-data/07_MAGMA/Input_GWAS_data/SUD/SUD_DEPvEXP_EUR.noAF.snploc'))
save(snp_pval_SUD, file=here('processed-data/07_MAGMA/Input_GWAS_data/SUD/SUD_DEPvEXP_EUR.noAF.pval'))



# ------------------------------------------------------------------------------
#                                 GWAS MDD data (2)
# ------------------------------------------------------------------------------

#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

snploc_MDD = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122.snploc")))

dim(snploc_MDD)
# [1] 11700     3

head(snploc_MDD)
#          SNP CHR       BP
# 1 rs12137092   1  4666965
# 2  rs3003457   1 17181435
# 3   rs754171   1 17196325
# 4  rs2501818   1 17192406
# 5  rs2977232   1 17183742
# 6  rs3104441   1 17191335

## No repeated SNPs
which(duplicated(snploc_MDD$SNP))
# integer(0)

## Check no variants in sexual chrs
unique(snploc_MDD$CHR)
# [1]  1 10 11 12 13 14 15 16 17  2  3  4  5  6 18 19 20 21 22  7  8  9

## Remove NAs
snploc_MDD <- subset(snploc_MDD, !is.na(SNP) & !is.na(CHR) & !is.na(BP))
dim(snploc_MDD)
# [1] 11700     3

## Save
save(snploc_MDD, file=here('processed-data/07_MAGMA/Input_GWAS_data/MDD/MDD.phs001672.pha005122.snploc'))


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

snp_pval_MDD = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/MDD/MDD.phs001672.pha005122.pval")))

dim(snp_pval_MDD)
# [1] 1700     2

head(snp_pval_MDD)
#          SNP         P
# 1 rs12137092 4.145e-05
# 2  rs3003457 6.309e-05
# 3   rs754171 8.449e-05
# 4  rs2501818 4.033e-05
# 5  rs2977232 5.965e-05
# 6  rs3104441 5.213e-05

## Unique SNPs
which(duplicated(snp_pval_MDD$SNP))
# integer(0)

## Common SNPs in snploc and snp_pval
setdiff(snploc_MDD$SNP, snp_pval_MDD$SNP)
# character(0)

## Remove NA
snp_pval_MDD <- subset(snp_pval_MDD, !is.na(SNP) & !is.na(P))
dim(snp_pval_MDD)
# [1] 11700     2

## Save
save(snp_pval_MDD, file=here('processed-data/07_MAGMA/Input_GWAS_data/MDD/MDD.phs001672.pha005122.pval'))



# ------------------------------------------------------------------------------
#                                 GWAS OUD data
# ------------------------------------------------------------------------------

## Input GWAS data
gwas_OUD <- fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/OUD/OUD.phs001672.pha004954.txt"), skip=17)

dim(gwas_OUD)
# [1] 1322   14

head(gwas_OUD)
#        SNP ID Marker Type   P-value Chr ID Chr Position Rank Allele1 Allele2
# 1: rs12585674          NA 2.674e-05     13    100763239   NA       C       G
# 2: rs16957671          NA 5.659e-05     13    100765829   NA       T       C
# 3: rs17074629          NA 6.999e-05     13     30477600   NA       A       G
# 4: rs17088924          NA 2.402e-05     13     64605453   NA       T       C
# 5: rs17579424          NA 8.402e-05     13    100756534   NA       A       G
# 6: rs56397918          NA 8.784e-05     13     50835248   NA       C       G
# Sample size Effect Allele Odds ratio CI low CI high Bin ID
# 1:       76287             C         NA     NA      NA   1096
# 2:       79474             T         NA     NA      NA   1096
# 3:       79634             A         NA     NA      NA   1061
# 4:       77428             T         NA     NA      NA   1078
# 5:       79102             A         NA     NA      NA   1096
# 6:       77965             C         NA     NA      NA   1071


#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

## Columns of interest
snploc_OUD <- gwas_OUD[, c("SNP ID", "Chr ID", "Chr Position")]
colnames(snploc_OUD) <- c('SNP', 'CHR', 'BP')

head(snploc_OUD)
#           SNP    CHR           BP
# 1: rs12585674     13    100763239
# 2: rs16957671     13    100765829
# 3: rs17074629     13     30477600
# 4: rs17088924     13     64605453
# 5: rs17579424     13    100756534
# 6: rs56397918     13     50835248

## No repeated SNPs
which(duplicated(snploc_OUD$SNP))
# integer(0)

## Check no variants in sexual chrs
unique(snploc_OUD$CHR)
#  [1] 13 14 15 16 17 18  1  2  3  4  5  6  7  8  9 10 11 12 19 20 22

## Remove NAs
snploc_OUD <- subset(snploc_OUD, !is.na(SNP) & !is.na(CHR) & !is.na(BP))
dim(snploc_OUD)
# [1] 1322    3

snploc_OUD <- as.data.frame(snploc_OUD)

## Save
save(snploc_OUD, file=here('processed-data/07_MAGMA/Input_GWAS_data/OUD/OUD.phs001672.pha004954.snploc'))


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

snp_pval_OUD <-  gwas_OUD[, c("SNP ID", "P-value")]
colnames(snp_pval_OUD) <- c('SNP', 'P')

head(snp_pval_OUD)
#           SNP         P
# 1: rs12585674 2.674e-05
# 2: rs16957671 5.659e-05
# 3: rs17074629 6.999e-05
# 4: rs17088924 2.402e-05
# 5: rs17579424 8.402e-05
# 6: rs56397918 8.784e-05

## Common SNPs in snploc and snp_pval
setdiff(snploc_OUD$SNP, snp_pval_OUD$SNP)
# character(0)

## Remove NAs
snp_pval_OUD <- subset(snp_pval_OUD, !is.na(SNP) & !is.na(P))
dim(snp_pval_OUD)
# [1] 1322    3

snp_pval_OUD <- as.data.frame(snp_pval_OUD)

## Save
save(snp_pval_OUD, file=here('processed-data/07_MAGMA/Input_GWAS_data/OUD/OUD.phs001672.pha004954.pval'))

