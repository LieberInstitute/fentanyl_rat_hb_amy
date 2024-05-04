#!/bin/bash -l
#SBATCH --output=logs/03_MAGMA_scz2022.txt
#SBATCH --error=logs/03_MAGMA_scz2022.txt
#SBATCH --partition=shared
#SBATCH --job-name=03_MAGMA_scz2022
#SBATCH --mem=25GB

################################################################################
##                                4. Run MAGMA
################################################################################

echo "**** MAGMA running ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## load modules
module load magma/1.10
## list modules
module list


# ------------------------------------------------------------------------------
#                                 GWAS SCZ data
# ------------------------------------------------------------------------------

ANNOT_PREFIX="../../processed-data/07_MAGMA/Input_GWAS_data/SCZ/SCZ_PGC3_wave3.european.autosome.public.v3"
OUTPUT_PREFIX="../../processed-data/07_MAGMA/Output/SCZ/SCZ_PGC3_wave3.european.autosome.public.v3"

## Step 1: Annotate SNPs onto genes 
## (with GRCh37 (h19) human genome reference build)
magma --annotate --snp-loc $ANNOT_PREFIX.snploc\
	--gene-loc ../../processed-data/07_MAGMA/geneloc/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc\
	--out $OUTPUT_PREFIX



#magma --annotate --snp-loc /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.snploc\
	#--gene-loc /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/geneloc/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc\
	#--out $OUTPUT_PREFIX



## Step 2 Gene Analysis - SNP p-values (model SNPwise-mean)
## Use same reference dataset as in HumanPilot project
magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval $ANNOT_PREFIX.pval ncol=N\
	--gene-annot $OUTPUT_PREFIX.genes.annot\
	--out $OUTPUT_PREFIX


## Step 3 Gene Set Analysis 
magma --gene-results $OUTPUT_PREFIX.genes.raw\
	--set-annot ../../processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/07_MAGMA/Output/SCZ/SCZ_MAGMA

echo "**** Job ends ****"
date
