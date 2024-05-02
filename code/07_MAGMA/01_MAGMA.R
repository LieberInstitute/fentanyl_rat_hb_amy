
library(here)
library(data.table)
library(biomaRt)
library(sessioninfo)
library(org.Hs.eg.db)


####################   Generalized Gene-Set Analysis of GWAS Data (MAGMA)   ######################

# Note: this code is based on https://github.com/LieberInstitute/Habenula_Pilot/tree/master/code/13_MAGMA


################################################################################
##                       1. GWAS input data preparation
################################################################################


## Function to extract positions of human SNPs if missing in GWAS data
SNP_loc_extraction <- function(snp_list){

    output <- vector()
    i=1
    for (snps in snp_list[i:i+999]){
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
        i=i+1000
    }


    return(output)
}


# ------------------------------------------------------------------------------
#                                 GWAS SZC data
# ------------------------------------------------------------------------------

## Internal paths (data is in JHPCE)
gwas_szc = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz")))

dim(gwas_szc)
# [1] 7659767      14

## Explore data
head(gwas_szc)
# CHROM         ID       POS     A1     A2  FCAS  FCON IMPINFO         BETA
# <int>     <char>     <int> <char> <char> <num> <num>   <num>        <num>
#  8    rs62513865 101592213      C      T 0.930 0.927   0.963  0.011997738
#  8    rs79643588 106973048      G      A 0.907 0.906   0.997 -0.008596847
#  8    rs17396518 108690829      T      G 0.565 0.566   0.985 -0.002102208
#  8      rs983166 108681675      A      C 0.564 0.563   0.988  0.004897985
#  8    rs28842593 103044620      T      C 0.840 0.840   0.948 -0.003897586
#  8     rs7014597 104152280      G      C 0.841 0.838   0.994  0.007898723
#   SE   PVAL  NCAS  NCON     NEFF
# <num>  <num> <int> <int>    <num>
# 0.0171 0.4847 53386 77258 58749.13
# 0.0148 0.5605 53386 77258 58749.13
# 0.0087 0.8145 53386 77258 58749.13
# 0.0087 0.5704 53386 77258 58749.13
# 0.0121 0.7488 53386 77258 58749.13
# 0.0117 0.5034 53386 77258 58749.13

## Confirm no repearted variants
which(duplicated(gwas_szc$ID))
# integer(0)

## Check no variants in sexual chrs
unique(gwas_szc$CHROM)
# [1]  8  2 16  7  6 15 14 12 17  1  3 11 22 20 21  9  4 13  5 10 18 19


#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

snploc_szc <- gwas_szc[,c('ID', 'CHROM', 'POS')]
## (col names for MAGMA)
colnames(snploc_szc) <- c('SNP', 'CHR', 'BP')

head(snploc_szc)
#        SNP   CHR        BP
# rs62513865     8 101592213
# rs79643588     8 106973048
# rs17396518     8 108690829
#   rs983166     8 108681675
# rs28842593     8 103044620
#  rs7014597     8 104152280

## Remove NA
snploc_szc <- subset(snploc_szc, !is.na(SNP) & !is.na(CHR) & !is.na(BP))

write.table(snploc_szc,
            file = here("processed-data/07_MAGMA/Input_GWAS_data/SCZ", "SCZ_PGC3_wave3.european.autosome.public.v3.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F)


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

gwas_szc$N <- gwas_szc$NCAS + gwas_szc$NCON
snp_pval_szc <- gwas_szc[, c('ID', 'PVAL','N')]
colnames(snp_pval_szc) <- c('SNP', 'P', 'N')

head(snp_pval_szc)
#        SNP      P      N
# rs62513865 0.4847 130644
# rs79643588 0.5605 130644
# rs17396518 0.8145 130644
#   rs983166 0.5704 130644
# rs28842593 0.7488 130644
#  rs7014597 0.5034 130644

## Remove NA
snp_pval_szc <- subset(snp_pval_szc, !is.na(SNP) & !is.na(P) & !is.na(N))

write.table(snp_pval_szc,
            file = here("processed-data/07_MAGMA/Input_GWAS_data/SCZ", "SCZ_PGC3_wave3.european.autosome.public.v3.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F
)



# ------------------------------------------------------------------------------
#                                 GWAS MDD data
# ------------------------------------------------------------------------------

gwas_mdd = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/mdd2019edinburgh/PGC_UKB_depression_genome-wide.txt")))

dim(gwas_mdd)
# [1] 8483301       7

head(gwas_mdd)
#    MarkerName A1 A2   Freq   LogOR StdErrLogOR         P
#     rs2326918  a  g 0.8452  0.0106      0.0060 0.0756100
#     rs7929618  c  g 0.1314 -0.0224      0.0064 0.0004804
#    rs66941928  t  c 0.8031  0.0003      0.0055 0.9502000
#     rs7190157  a  c 0.3517  0.0024      0.0045 0.5992000
#    rs12364336  a  g 0.8685  0.0075      0.0064 0.2450000
#     rs6977693  t  c 0.8544  0.0089      0.0061 0.1442000

colnames(gwas_mdd)[1] <- 'SNP'

which(duplicated(gwas_mdd$SNP))
# integer(0)


#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

## No SNP chr/location provided:
## Add position and chr for MDD SNPs based on SZC data
common_MDD_SNPs <- snploc_szc[which(snploc_szc$SNP %in% gwas_mdd$SNP), ]
dim(common_MDD_SNPs)
# [1] 7280414       3

## Obtain position of MDD SNPs not present in SZC data
MDD_snps_no_loc <- gwas_mdd[which(! gwas_mdd$SNP %in% common_MDD_SNPs$SNP), 'SNP']
length(MDD_snps_no_loc)
# [1] 1202887

## Extract location for these SNPs (HEREEEE)
MDD_snps_pos <- SNP_loc_extraction(MDD_snps_no_loc)
slurmjobs::job_single(name = "00_snp_lookup", memory = "25G", create_shell = TRUE, command = "Rscript 00_snp_lookup.R")

## Bind all SNPs
snploc_mdd <- rbind(common_MDD_SNPs, MDD_snps_pos)
dim(snploc_mdd)

## Check no variants in sexual chrs
unique(snploc_mdd$CHR)


write.table(snploc_mdd,
            file = here("processed-data/07_MAGMA/Input_GWAS_data/MDD", "MDD_PGC_UKB_depression_genome-wide.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F)


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

snp_pval_mdd <- gwas_mdd[, c('SNP', 'P')]
head(snp_pval_mdd)
#        SNP      P      N
# rs62513865 0.4847 130644
# rs79643588 0.5605 130644
# rs17396518 0.8145 130644
#   rs983166 0.5704 130644
# rs28842593 0.7488 130644
#  rs7014597 0.5034 130644

write.table(snp_pval_mdd,
            file = here("processed-data/07_MAGMA/Input_GWAS_data/MDD", "MDD_PGC_UKB_depression_genome-wide.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F)



# ------------------------------------------------------------------------------
#                                 GWAS Panic data
# ------------------------------------------------------------------------------

gwas_panic = as.data.frame(fread(here("/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/panic2019/pgc-panic2019.vcf.tsv.gz")))

dim(gwas_panic)
# [1]

head(gwas_panic)

which(duplicated(gwas_panic$ID))
# integer(0)

## Check no variants in sexual chrs
unique(gwas_panic$CHROM)

#  - - - - - - - - - - - - - - A) SNP location file  - - - - - - - - - - - - - -

snploc_panic <- gwas_panic[,c('ID', '#CHROM', 'POS')]
colnames(snploc_panic) <- c('SNP', 'CHR', 'BP')

head(snploc_panic)

write.table(snploc_panic,
            file = here("processed-data/07_MAGMA/Input_GWAS_data/Panic", "Panic_pgc-panic2019.snploc"),
            sep = "\t", col.names = T, row.names = F, quote = F)


#  - - - - - - - - - - - - - - - B) SNP p-values - - - - - - - - - - - - - - - -

gwas_panic$N <- gwas_panic$NCAS + gwas_panic$NCON
snp_pval_panic <- gwas_panic[, c('ID', 'PVAL','N')]
colnames(snp_pval_panic) <- c('SNP', 'P', 'N')

head(snp_pval_panic)

write.table(snp_pval_panic,
            file = here("processed-data/07_MAGMA/Input_GWAS_data/Panic", "Panic_pgc-panic2019.pval"),
            sep = "\t", col.names = T, row.names = F, quote = F)






