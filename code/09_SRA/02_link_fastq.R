library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('processed-data', '02_SPEAQeasy', 'samples.manifest')
dest_dir = here('raw-data', 'FASTQ', 'flat_dir')

dir.create(dest_dir, showWarnings = FALSE)

#   Read in manifest from SPEAQeasy
man_df = read.table(
        man_path,
        col.names = c('filename', 'md1', 'filename2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    select(filename, filename2, sample_name)

#   Assume one pair of files per sample, each ending in '.fastq.gz'
stopifnot(nrow(man_df) == length(unique(man_df$sample_name)))
stopifnot(all(grepl('\\.fastq\\.gz$', man_df$filename)))
stopifnot(all(grepl('\\.fastq\\.gz$', man_df$filename2)))

#   Link original files to the new directory
file.symlink(
    man_df$filename,
    file.path(dest_dir, sprintf('%s_1.fastq.gz', man_df$sample_name))
) |>
    all() |>
    stopifnot()
file.symlink(
    man_df$filename2,
    file.path(dest_dir, sprintf('%s_2.fastq.gz', man_df$sample_name))
) |>
    all() |>
    stopifnot()

session_info()
