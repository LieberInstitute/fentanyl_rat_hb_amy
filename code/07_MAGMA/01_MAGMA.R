library(here)




#### GWAS SZC Data ####
gwas = fread(here("processed-data", "13_MAGMA","GWAS", "SCZ", "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"))
dim(gwas)
# [1] 7659767      14
