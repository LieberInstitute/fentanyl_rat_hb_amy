library('here')
library('jaffelab')
library('sessioninfo')
library('readxl')
library('tidyverse')

man_path = here('processed-data', '02_SPEAQeasy', 'samples.manifest')
sample_info_path = here(
    'raw-data', 'FentanylvsSaline_SelfAdministration_RNAextraction.xlsx'
)
fastq_dir = here('raw-data', 'FASTQ', 'psomagen')

r1 = list.files(
    fastq_dir, pattern = paste0('.*_1\\.fastq\\.gz$'), full.names = TRUE,
    recursive = TRUE
)

r2 = list.files(
    fastq_dir, pattern = paste0('.*_2\\.fastq\\.gz$'), full.names = TRUE,
    recursive = TRUE
)

sample_info = read_excel(
    sample_info_path, sheet = "SubjectInformation for Analysis"
) |>
    filter(!is.na(!!as.symbol('Sample No.')))

#   Use underscores in sample IDs instead of hyphens or spaces
sample_info['Tissue Punch Label'] = gsub(
    ' |-', '_', sample_info[['Tissue Punch Label']]
)

sample_ids = sample_info[
    match(ss(basename(r1), '_'), sample_info[['Sample No.']]),
] |>
    pull('Tissue Punch Label')

man = paste(r1, 0, r2, 0, sample_ids, sep = '\t')

writeLines(man, con = man_path)

session_info()
