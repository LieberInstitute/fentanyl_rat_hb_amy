library('here')
library('jaffelab')
library('sessioninfo')
library('readxl')
library('tidyverse')

man_path = here('processed-data', '02_SPEAQeasy', 'samples.manifest')
sample_info_path = here(
    'raw-data', 'FentanylvsSaline_SelfAdministration_RNAextraction.xlsx'
)

get_fastq_paths = function(read_num) {
    return(
        c(
            list.files(
                here('raw-data', 'FASTQ', 'psomagen', 'AN00012355'),
                pattern = paste0('.*_', read_num, '\\.fastq\\.gz$'),
                full.names = TRUE
            ),
            list.files(
                here('raw-data', 'FASTQ', 'psomagen', 'AN00012356'),
                pattern = paste0('.*_', read_num, '\\.fastq\\.gz$'),
                full.names = TRUE
            )
        )
    )
}

sample_info = read_excel(
    sample_info_path, sheet = "SubjectInformation for Analysis"
)

#   Use underscores in sample IDs instead of hyphens or spaces
sample_info['Tissue Punch Label'] = gsub(
    ' |-', '_', sample_info[['Tissue Punch Label']]
)

r1 = get_fastq_paths(1)
sample_ids = sample_info[
        match(ss(basename(r1), '_'), sample_info[['Sample No.']]),
    ] |>
    pull('Tissue Punch Label')
    
man = paste(r1, 0, get_fastq_paths(2), 0, sample_ids, sep = '\t')

writeLines(man, con = man_path)

session_info()
