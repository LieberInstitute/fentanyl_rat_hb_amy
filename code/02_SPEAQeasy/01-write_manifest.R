library('here')
library('jaffelab')
library('sessioninfo')

man_path = here('processed-data', '02_SPEAQeasy', 'samples.manifest')

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

dir.create(dirname(man_path), showWarnings = FALSE)

r1 = get_fastq_paths(1)
sample_ids = ss(basename(r1), '_')
man = paste(r1, 0, get_fastq_paths(2), 0, sample_ids, sep = '\t')

writeLines(man, con = man_path)

session_info()
