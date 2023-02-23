#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=300G
#$ -N psomagen_data
#$ -o logs/psomagen_data.txt
#$ -e logs/psomagen_data.txt
#$ -m e

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

## Locate directory
ORIGINALDIR=$(awk "NR==${SGE_TASK_ID}" qSVA_dcl01_dirs.txt)
echo "Processing sample ${ORIGINALDIR}"
date

BASEDIR=$(basename ${ORIGINALDIR})
ORIGINALHOME=$(dirname ${ORIGINALDIR})

## Determine amount of data to transfer
du -sk --apparent-size ${ORIGINALDIR}/ | awk '{$1=$1/(1024^3); print $1, "TB";}'

## List owners of all the files in the original path
find ${ORIGINALDIR}/ -exec ls -l {} \; | grep -v total | tr -s ' ' | cut -d ' ' -f3,9-

## Copy from dcl01 to dcs04
rsync -rltgvh --chown=:lieber_lcolladotor ${ORIGINALDIR}/ /dcs04/lieber/lcolladotor/qSVA_LIBD3080/${BASEDIR}/


echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
