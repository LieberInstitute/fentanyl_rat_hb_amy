#!/bin/bash
#$ -cwd
#$ -N "create_rse_tx"
#$ -o ../../processed-data/02_SPEAQeasy/03-create_rse_tx.log
#$ -e ../../processed-data/02_SPEAQeasy/03-create_rse_tx.log
#$ -l mf=20G,h_vmem=20G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.1.x # for compatibility with SPEAQeasy
Rscript 03-create_rse_tx.R

echo "**** Job ends ****"
date
