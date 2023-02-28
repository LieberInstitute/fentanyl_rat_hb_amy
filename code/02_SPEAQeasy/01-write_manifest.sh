#!/bin/bash
#$ -cwd
#$ -N "write_manifest"
#$ -o ../../processed-data/02_SPEAQeasy/01-write_manifest.log
#$ -e ../../processed-data/02_SPEAQeasy/01-write_manifest.log
#$ -l mf=5G,h_vmem=5G

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 01-write_manifest.R

echo "**** Job ends ****"
date
