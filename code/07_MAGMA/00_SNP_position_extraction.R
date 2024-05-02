
library(here)
library(data.table)
library(biomaRt)
library(sessioninfo)

## Load input MDD & OUD SNPs (with missing location info)
load(here('processed-data/07_MAGMA/Input_GWAS_data/MDD/MDD_snps_no_loc.txt'), verbose = TRUE)
# load(here('processed-data/07_MAGMA/Input_GWAS_data/OUD/OUD_snps_no_loc.txt'), verbose = TRUE)


## Function to extract positions of human SNPs if missing in GWAS data
SNP_loc_extraction <- function(snp_list){

    ## Mart
    mart <- useEnsembl(biomart="ENSEMBL_MART_SNP",
                       host="https://grch37.ensembl.org",  ## Use Grch37 assembly
                       dataset="hsapiens_snp") ## human SNPs

    message(Sys.time(), " - Searching SNP position")
    ## Obtain ID, chr, and start position of SNPs
    snploc <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start'),
                    values = snp_list,
                    filters ="snp_filter",
                    mart = mart)
    snploc <- as.data.frame(snploc)
    colnames(snploc) <- c('SNP', 'CHR', 'BP')

    return(snploc)
}


## MDD SNPs
MDD_SNPs_extracted_pos <- SNP_loc_extraction(MDD_snps_no_loc)
save(MDD_SNPs_extracted_pos, file='processed-data/07_MAGMA/MDD/MDD_SNPs_extracted_pos.snploc')



## OUD SNPs

