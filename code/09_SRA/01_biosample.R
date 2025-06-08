library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('processed-data', '02_SPEAQeasy', 'samples.manifest')
pheno_path = here('raw-data', 'sample_metadata_and_QCmetrics.csv')
out_path = here('processed-data', '09_SRA', 'biosample.tsv')

#   Read in the manifest just to make sure the phenotype table includes the
#   expected set of IDs (though we'll only really use the phenotype table)
man_df = read.table(
        man_path, col.names = c('r1', 'md1', 'r2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    select(sample_name)

pheno_df = read_csv(pheno_path, show_col_types = FALSE) |>
    rename(sample_name = SAMPLE_ID) |>
    mutate(
        organism = 'rattus norvegicus',
        isolate = sprintf(
            '%s substance group (sample %s)', Substance, sample_name
        ),
        geo_loc_name = 'United States: Baltimore, MD',
        sex = tolower(Sex),
        biomaterial_provider = 'Lieber Institute for Brain Development: 855 North Wolfe Street, Suite 300, 3rd Floor, Baltimore, MD 21205',
        treatment = Substance
    )

#   Phenotype data and manifest should include same set of sample IDs
stopifnot(setequal(pheno_df$sample_name, man_df$sample_name))

write_tsv(pheno_df, out_path)

session_info()
