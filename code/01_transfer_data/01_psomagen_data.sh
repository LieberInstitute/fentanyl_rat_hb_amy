#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=300G
#$ -N psomagen_data
#$ -o logs/psomagen_data.$TASK_ID.txt
#$ -e logs/psomagen_data.$TASK_ID.txt
#$ -m e
#$ -t 1-2
#$ -tc 2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module list

## Create output directory
mkdir -p /dcs04/lieber/marmaypag/fentanylRat_LIBD4205/fentanyl_rat_hb_amy/raw-data/FASTQ/psomagen

## Locate directory
ORIGINALDIR=$(awk "NR==${SGE_TASK_ID}" 01_psomagen_directories.txt)
echo "Processing sample ${ORIGINALDIR}"
date

BASEDIR=$(basename ${ORIGINALDIR})
ORIGINALHOME=$(dirname ${ORIGINALDIR})

## Determine amount of data to transfer
du -sk --apparent-size ${ORIGINALDIR}/ | awk '{$1=$1/(1024^3); print $1, "TB";}'

## List owners of all the files in the original path
find ${ORIGINALDIR}/ -exec ls -l {} \; | grep -v total | tr -s ' ' | cut -d ' ' -f3,9-

## Copy from /dcs04/lieber/data/transfers/psomagen to dcs04/marmaypag
rsync -rltgvh --chown=:lieber_lcolladotor ${ORIGINALDIR}/ /dcs04/lieber/marmaypag/fentanylRat_LIBD4205/fentanyl_rat_hb_amy/raw-data/FASTQ/psomagen/${BASEDIR}/

## Update permissions
sh /dcs04/lieber/lcolladotor/_jhpce_org_LIBD001/update_permissions_spatialteam.sh /dcs04/lieber/marmaypag/fentanylRat_LIBD4205/fentanyl_rat_hb_amy/raw-data/FASTQ/psomagen

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
