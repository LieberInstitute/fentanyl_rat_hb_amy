#!/bin/bash
#$ -l mem_free=40G,h_vmem=40G,h_fsize=800G
#$ -o ../../processed-data/02_SPEAQeasy/SPEAQeasy_output.log
#$ -e ../../processed-data/02_SPEAQeasy/SPEAQeasy_output.log
#$ -N run_pipeline
#$ -cwd

REPO_DIR=$(git rev-parse --show-toplevel)
SPEAQEASY_DIR=$REPO_DIR/code/02_SPEAQeasy/SPEAQeasy
PROCESSED_DIR=$REPO_DIR/processed-data/02_SPEAQeasy

module load nextflow/20.01.0
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $SPEAQEASY_DIR/main.nf \
    --sample "paired" \
    --reference "rat" \
    --strand "forward" \
    --input "$PROCESSED_DIR" \
    --output "$PROCESSED_DIR/out" \
    -w "$PROCESSED_DIR/work" \
    --annotation "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation" \
    -with-report "$PROCESSED_DIR/02_run_pipeline.html" \
    -profile jhpce

#   Log successful runs on non-test data in a central location. Please adjust
#   the log path here if it is changed at the top!
bash $SPEAQEASY_DIR/scripts/track_runs.sh $PROCESSED_DIR/SPEAQeasy_output.log

#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash $SPEAQEASY_DIR/scripts/generate_logs.sh $PROCESSED_DIR/SPEAQeasy_output.log
