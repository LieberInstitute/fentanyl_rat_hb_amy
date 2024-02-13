This directory contains code used to run SPEAQeasy on the raw FASTQ sequences to obtain RSE objects for downstream analyses.

SPEAQeasy was run with default settings as of commit `87ba0b4`, with the exception of using ENSEMBL release 109 as annotation (i.e. `ensembl_version_rat = "109"`).

- `01-write_manifest.*`: Produce the `samples.manifest` file, a required input to SPEAQeasy.
- `02-run_pipeline.sh`: The shell script invoking SPEAQeasy itself.
- `03-create_rse_tx.*`: A script to manually run the `CountObjects` process of SPEAQeasy, required because a bug at the time used "primary" instead of "main" annotation, and didn't initially produce the important `rse_tx` object.
