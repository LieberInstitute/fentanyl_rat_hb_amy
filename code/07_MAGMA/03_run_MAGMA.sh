#!/bin/bash -l
#SBATCH --output=logs/03_MAGMA_scz2022.txt
#SBATCH --error=logs/03_MAGMA_scz2022.txt
#SBATCH --partition=shared
#SBATCH --job-name=03_MAGMA_scz2022
#SBATCH --mem=25GB

################################################################################
##                                3. Run MAGMA
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

echo "**** SCZ analysis ends ****"
date


# ------------------------------------------------------------------------------
#                                 GWAS MDD data
# ------------------------------------------------------------------------------

ANNOT_PREFIX="../../processed-data/07_MAGMA/Input_GWAS_data/MDD_2019/MDD_PGC_UKB_depression_genome-wide"
OUTPUT_PREFIX="../../processed-data/07_MAGMA/Output/MDD_2019/MDD_PGC_UKB_depression_genome-wide"

## Step 1: Annotate SNPs onto genes
## (GRCh37 (h19) human build)
magma --annotate --snp-loc $ANNOT_PREFIX.snploc\
	--gene-loc ../../processed-data/07_MAGMA/geneloc/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc\
	--out $OUTPUT_PREFIX

## Step 2 Gene Analysis - SNP p-values (model SNPwise-mean)
magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval $ANNOT_PREFIX.pval N=807553\
	--gene-annot $OUTPUT_PREFIX.genes.annot\
	--out $OUTPUT_PREFIX

## Step 3 Gene Set Analysis
magma --gene-results $OUTPUT_PREFIX.genes.raw\
	--set-annot ../../processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/07_MAGMA/Output/MDD_2019/MDD_2019_MAGMA

echo "**** MDD analysis ends ****"
date


# ------------------------------------------------------------------------------
#                                 GWAS Panic data
# ------------------------------------------------------------------------------

ANNOT_PREFIX="../../processed-data/07_MAGMA/Input_GWAS_data/Panic/Panic_PGC_panic2019"
OUTPUT_PREFIX="../../processed-data/07_MAGMA/Output/Panic/Panic_PGC_panic2019"

## Step 1: Annotate SNPs onto genes
## (GRCh37 (h19) human build)
magma --annotate --snp-loc $ANNOT_PREFIX.snploc\
	--gene-loc ../../processed-data/07_MAGMA/geneloc/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc\
	--out $OUTPUT_PREFIX

## Step 2 Gene Analysis - SNP p-values (model SNPwise-mean)
magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval $ANNOT_PREFIX.pval ncol=N\
	--gene-annot $OUTPUT_PREFIX.genes.annot\
	--out $OUTPUT_PREFIX

## Step 3 Gene Set Analysis
magma --gene-results $OUTPUT_PREFIX.genes.raw\
	--set-annot ../../processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/07_MAGMA/Output/Panic/Panic_MAGMA

echo "**** Panic analysis ends ****"
date


# ------------------------------------------------------------------------------
#                                 GWAS SUD data
# ------------------------------------------------------------------------------

ANNOT_PREFIX="../../processed-data/07_MAGMA/Input_GWAS_data/SUD/SUD_DEPvEXP_EUR.noAF"
OUTPUT_PREFIX="../../processed-data/07_MAGMA/Output/SUD/SUD_DEPvEXP_EUR.noAF"

## Step 1: Annotate SNPs onto genes
## (GRCh37 (h19) human build)
magma --annotate --snp-loc $ANNOT_PREFIX.snploc\
	--gene-loc ../../processed-data/07_MAGMA/geneloc/GRCh38-ensembl93_to_hg19-lifted_30k-expressing-GENES.gene.loc\
	--out $OUTPUT_PREFIX

## Step 2 Gene Analysis - SNP p-values (model SNPwise-mean)
magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval $ANNOT_PREFIX.pval ncol=N\
	--gene-annot $OUTPUT_PREFIX.genes.annot\
	--out $OUTPUT_PREFIX

## Step 3 Gene Set Analysis
magma --gene-results $OUTPUT_PREFIX.genes.raw\
	--set-annot ../../processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/07_MAGMA/Output/SUD/SUD_MAGMA

echo "**** SUD analysis ends ****"
date


# ------------------------------------------------------------------------------
#                                 GWAS MDD data
# ------------------------------------------------------------------------------

ANNOT_PREFIX="../../processed-data/07_MAGMA/Input_GWAS_data/MDD/MDD.phs001672.pha005122"
OUTPUT_PREFIX="../../processed-data/07_MAGMA/Output/MDD/MDD.phs001672.pha005122"

## Step 1: Annotate SNPs onto genes
## (with GRCh38 human genome reference build)
magma --annotate --snp-loc $ANNOT_PREFIX.snploc\
	--gene-loc ../../processed-data/07_MAGMA/geneloc/GRCh38_Ensembl-93_GENES_chr-x-y-mt.gene.loc\
	--out $OUTPUT_PREFIX

## Step 2 Gene Analysis - SNP p-values (model SNPwise-mean)
magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval $ANNOT_PREFIX.pval N=1154267\
	--gene-annot $OUTPUT_PREFIX.genes.annot\
	--out $OUTPUT_PREFIX

## Step 3 Gene Set Analysis
magma --gene-results $OUTPUT_PREFIX.genes.raw\
	--set-annot ../../processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/07_MAGMA/Output/MDD/MDD_MAGMA

echo "**** MDD analysis ends ****"
date


# ------------------------------------------------------------------------------
#                                 GWAS OUD data
# ------------------------------------------------------------------------------

ANNOT_PREFIX="../../processed-data/07_MAGMA/Input_GWAS_data/OUD/OUD.phs001672.pha004954"
OUTPUT_PREFIX="../../processed-data/07_MAGMA/Output/OUD/OUD.phs001672.pha004954"

## Step 1: Annotate SNPs onto genes
## (with GRCh38 human genome reference build)
magma --annotate --snp-loc $ANNOT_PREFIX.snploc\
	--gene-loc ../../processed-data/07_MAGMA/geneloc/GRCh38_Ensembl-93_GENES_chr-x-y-mt.gene.loc\
	--out $OUTPUT_PREFIX

## Step 2 Gene Analysis - SNP p-values (model SNPwise-mean)
magma --bfile /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur\
	--pval $ANNOT_PREFIX.pval N=1154267\
	--gene-annot $OUTPUT_PREFIX.genes.annot\
	--out $OUTPUT_PREFIX

## Step 3 Gene Set Analysis
magma --gene-results $OUTPUT_PREFIX.genes.raw\
	--set-annot ../../processed-data/07_MAGMA/Input_Gene_Sets/gene_sets.txt gene-col=Gene set-col=Set\
	--out ../../processed-data/07_MAGMA/Output/OUD/OUD_MAGMA

echo "**** OUD analysis ends ****"
date


