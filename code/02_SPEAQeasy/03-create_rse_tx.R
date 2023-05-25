#   SPEAQeasy was run before a bug fix, and as a result the 'rse_tx' output
#   object was initially not created, and the reference transcripts FASTA
#   included some 'primary' transcripts, despite 'main' being specified via the
#   configuration variable 'anno_build'.
#   
#   Here we manually run a version of 'create_count_objects.R' that uses just
#   the 'main' transcripts of the FASTA and simply builds and writes the
#   'rse_tx' object.

## Required libraries
library("derfinder")
library("BiocParallel")
library("Biostrings")
library("GenomicRanges")
library("GenomicFeatures")
library("jaffelab")
library("getopt")
library("rafalib")
library("devtools")
library("SummarizedExperiment")
library("DelayedArray")
library("matrixStats")
library("plyr")
library("rtracklayer")
library("here")

#   Where the 'CountObjects' process was originally run
work_dir = here(
    'processed-data', '02_SPEAQeasy', 'work', 'a8',
    '134546d24ba3ee0c36879442289d20'
)

out_dir = here('processed-data', '02_SPEAQeasy', 'out', 'count_objects')

#   Command-line options normally passed via the CountObjects process
opt = list(
    'anno_suffix' = 'rat_ensembl_109_main',
    'organism' = 'rat',
    'experiment' = 'Jlab_experiment',
    'prefix' = '',
    'paired' = TRUE,
    'ercc' = FALSE,
    'cores' = 1,
    'stranded' = 'reverse',
    'salmon' = FALSE,
    'star' = FALSE,
    'output' = '/dcs04/lieber/marmaypag/fentanylRat_LIBD4205/fentanyl_rat_hb_amy/processed-data/02_SPEAQeasy/out'
)

#  Reference-specific libraries and sequence names
library("org.Rn.eg.db")
chr_names <- c(1:20, "X", "Y")
mito_chr <- "MT"

#  The default "prefix" is empty (i.e. "")
if (nchar(opt$prefix) > 0) {
    EXPNAME <- paste0(opt$experiment, "_", opt$prefix)
} else {
    EXPNAME <- opt$experiment
}

## read in pheno
manifest <- read.table(
    file.path(work_dir, "samples_complete.manifest"),
    sep = " ", header = FALSE, stringsAsFactors = FALSE
)
metrics <- data.frame(
    "SAMPLE_ID" = manifest[, ncol(manifest) - 1],
    "strandness" = manifest[, ncol(manifest)],
    stringsAsFactors = FALSE
)
metrics$SAMPLE_ID <- as.character(metrics$SAMPLE_ID)

N <- length(metrics$SAMPLE_ID)


###############################################################################
#  Parse FastQC logs and summaries
###############################################################################

fastqcdata <- c(
    "SeqLength", "percentGC", "phred1", "phred2", "phred3", "phred4",
    "phredGT30", "phredGT35", "Adapter1", "Adapter2", "Adapter3"
)
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))

#  Possible columns we expect in FastQC summaries
fastqc_stat_names <- c(
    "Basic Statistics",
    "Per base sequence quality",
    "Per tile sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
    "Kmer Content"
)

#  For choosing the post-trimming FastQC files for any samples that
#  were trimmed. If 'only_check_trim', then return a logical instead,
#  indicating if a sample was trimmed.
get_file <- function(id, suffix, read, only_check_trim = FALSE) {
    if (file.exists(paste0(id, read, "_trimmed", suffix))) {
        if (only_check_trim) {
            return(TRUE)
        } else {
            return(paste0(id, read, "_trimmed", suffix))
        }
    } else {
        if (only_check_trim) {
            return(FALSE)
        } else {
            return(paste0(id, read, "_untrimmed", suffix))
        }
    }
}

#  Parse a vector of FastQC summaries, containing each sample in the order
#  present in 'metrics'. Return a matrix of statistics
parse_summary <- function(summary_paths, metrics) {
    temp_rows <- lapply(summary_paths, function(x) {
        temp <- readLines(file.path(work_dir, x))
        vals <- ss(temp, "\t")
        names(vals) <- ss(temp, "\t", 2)
        
        stopifnot(all(names(vals) %in% fastqc_stat_names))
        vals <- vals[fastqc_stat_names]
    })
    
    actual_colnames <- gsub(" ", "_", tolower(fastqc_stat_names))
    stat_mat <- matrix(unlist(temp_rows),
                       ncol = length(fastqc_stat_names),
                       byrow = TRUE,
                       dimnames = list(metrics$SAMPLE_ID, actual_colnames)
    )
    
    return(stat_mat)
}

#  Parse a vector of FastQC "data reports", containing each sample in the order
#  present in 'metrics'. Return a matrix of statistics
parse_data <- function(data_path, metrics) {
    R <- lapply(data_path, function(x) {
        scan(
            file.path(work_dir, x), what = "character", sep = "\n",
            quiet = TRUE, strip = TRUE
        )
    })
    names(R) <- metrics$SAMPLE_ID
    
    ## Split list into sublists of metric categories
    zz <- lapply(R, function(x) splitAt(x, which(x == ">>END_MODULE") + 1))
    
    # sequence length
    seqlen <- sapply(zz, function(x) {
        index <- grep(">>Basic Statistics", x)
        stopifnot(ss(x[[index]][9], "\t") == "Sequence length")
        return(ss(x[[index]][9], "\t", 2))
    })
    
    # percent GC
    gcp <- sapply(zz, function(x) {
        index <- grep(">>Basic Statistics", x)
        stopifnot(ss(x[[index]][10], "\t") == "%GC")
        return(ss(x[[index]][10], "\t", 2))
    })
    
    # median phred scores (at roughly 1/4, 1/2, 3/4, and end of seq length)
    # get positions
    index <- grep(">>Per base sequence quality", zz[[1]])
    len <- round((length(zz[[1]][[index]]) - 3) / 4)
    pos <- c(len + 3, 2 * len + 3, 3 * len + 3, length(zz[[1]][[index]]) - 1)
    nameSuf <- ss(zz[[1]][[index]][pos], "\t", 1)
    fastqcdata[3:6] <- paste0("phred", nameSuf)
    
    phred <- lapply(zz, function(x) {
        index <- grep(">>Per base sequence quality", x)
        stopifnot(ss(x[[index]][2], "\t", 3) == "Median")
        return(ss(x[[index]][pos], "\t", 3))
    })
    phred <- matrix(unlist(phred), ncol = 4, byrow = TRUE)
    
    # proportion of reads above phred 30 and 35
    sc <- lapply(zz, function(x) {
        return(x[[grep(">>Per sequence quality scores", x)]])
    })
    
    phred2 <- lapply(sc, function(x) x[3:(length(x) - 1)])
    phred2 <- lapply(phred2, function(x) {
        data.frame(score = ss(x, "\t", 1), count = ss(x, "\t", 2))
    })
    phred2 <- lapply(phred2, function(x) {
        data.frame(x, cumulRev = rev(cumsum(rev(as.numeric(levels(x$count))[x$count]))))
    })
    phred2 <- lapply(phred2, function(x) {
        data.frame(x, prop = x$cumulRev / x$cumulRev[1])
    })
    phred2 <- lapply(phred2, function(x) c(x[which(x$score %in% c(30, 35)), "prop"], 0, 0)[1:2])
    phred2 <- matrix(unlist(phred2), ncol = 2, byrow = TRUE)
    
    # Illumina adapter content (at roughly 1/2, 3/4, and end of seq length)
    # get positions
    ac <- lapply(zz, function(x) {
        return(x[[grep(">>Adapter Content", x)]])
    })
    
    len <- round((length(ac[[1]]) - 3) / 5)
    pos <- c(3 * len + 2, 4 * len + 2, length(ac[[1]]) - 1)
    nameSuf <- ss(ac[[1]][pos], "\t", 1)
    fastqcdata[9:11] <- paste0("Adapter", nameSuf)
    
    adap <- lapply(ac, function(x) x[pos])
    adap <- lapply(adap, function(x) ss(x, "\t", 2))
    adap <- matrix(as.numeric(unlist(adap)), ncol = 3, byrow = TRUE)
    
    combined <- data.frame(SeqLen = unlist(seqlen), GCprec = unlist(gcp), phred, phred2, adap)
    rownames(combined) <- NULL
    
    return(list(combined, fastqcdata))
}

#  Parse FastQC summaries and "data reports"; append to 'metrics'
if (opt$paired) {
    summary_paths_1 <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "_1")
    summary_paths_2 <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "_2")
    
    #  Get stat values for each mate
    stat_mat <- parse_summary(summary_paths_1, metrics)
    stat_mat_2 <- parse_summary(summary_paths_2, metrics)
    
    #  Combine into a single matrix and simplify some values
    stat_mat[] <- paste(stat_mat, stat_mat_2, sep = "/")
    stat_mat[stat_mat == "PASS/PASS"] <- "PASS"
    stat_mat[stat_mat == "WARN/WARN"] <- "WARN"
    stat_mat[stat_mat == "FAIL/FAIL"] <- "FAIL"
    stat_mat[stat_mat == "NA/NA"] <- "NA"
    
    metrics <- cbind(metrics, stat_mat)
    
    for (i in 1:2) {
        data_path <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_fastqc_data.txt", read = paste0("_", i))
        temp <- parse_data(data_path, metrics)
        combined <- temp[[1]]
        fastqcdata <- temp[[2]]
        names(combined) <- paste0(fastqcdata, "_R", i)
        metrics <- cbind(metrics, combined)
    }
} else {
    summary_paths <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "")
    metrics <- cbind(metrics, parse_summary(summary_paths, metrics))
    
    data_paths <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_fastqc_data.txt", read = "")
    temp <- parse_data(data_paths, metrics)
    combined <- temp[[1]]
    fastqcdata <- temp[[2]]
    rownames(combined) <- NULL
    names(combined) <- fastqcdata
    
    metrics <- cbind(metrics, combined)
}

#  Add a metric indicating if each sample was trimmed
if (opt$paired) {
    metrics$trimmed <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "_1", only_check_trim = TRUE)
} else {
    metrics$trimmed <- sapply(metrics$SAMPLE_ID, get_file, suffix = "_summary.txt", read = "", only_check_trim = TRUE)
}

sampIDs <- as.vector(metrics$SAMPLE_ID)

###############################################################################
#  Read in transcript pseudo-alignment stats
###############################################################################

gencodeGTF <- import(
    con = list.files(work_dir, pattern = ".*\\.gtf", full.names = TRUE),
    format = "gtf"
)

txMap <- NULL
txRR <- gencodeGTF[which(gencodeGTF$type == "transcript")]
names(txRR) <- txRR$transcript_id

if (opt$salmon) {
    #----------------------------------------------------------------------
    #  Salmon quantification
    #----------------------------------------------------------------------
    
    ## observed tpm and number of reads
    txTpm <- bplapply(sampIDs, function(x) {
        read.table(file.path(work_dir, paste0(x, "_quant.sf")), header = TRUE)$TPM
    },
    BPPARAM = MulticoreParam(opt$cores)
    )
    txTpm <- do.call(cbind, txTpm)
    
    txNumReads <- bplapply(sampIDs, function(x) {
        read.table(file.path(work_dir, paste0(x, "_quant.sf")), header = TRUE)$NumReads
    },
    BPPARAM = MulticoreParam(opt$cores)
    )
    txNumReads <- do.call(cbind, txNumReads)
    
    ## get names of transcripts
    txNames <- read.table(
        file.path(work_dir, paste0(sampIDs[1], "_quant.sf")),
        header = TRUE
    )$Name
    txNames <- as.character(txNames)
} else {
    #----------------------------------------------------------------------
    #  Kallisto quantification
    #----------------------------------------------------------------------
    
    txMatrices <- bplapply(
        sampIDs,
        function(x) {
            read.table(
                file.path(work_dir, paste0(x, "_abundance.tsv")),
                header = TRUE
            )
        },
        BPPARAM = MulticoreParam(opt$cores)
    )
    
    txTpm <- do.call(cbind, lapply(txMatrices, function(x) x$tpm))
    txNumReads <- do.call(cbind, lapply(txMatrices, function(x) x$est_counts))
    txNames <- as.character(txMatrices[[1]]$target_id)
    rm(txMatrices)
}

colnames(txTpm) <- colnames(txNumReads) <- sampIDs

#   These organisms use a FASTA downloaded from Gencode
if (opt$organism %in% c("hg19", "hg38", "mm10")) {
    ## expected format of transcripts fasta header (and thus txMatrices and txNames):
    # ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|lncRNA|
    #       1                  2                3                   4                  5           6      7    8
    #      t_id             gene_id           o_g_id              o_t_id              t_name     gsym   tlen  biotype
    txMap <- t(ss(txNames, "\\|", c(1, 7, 2, 6, 8)))
    txMap <- as.data.frame(txMap)
    colnames(txMap) <- c(
        "gencodeTx", "txLength", "gencodeID", "Symbol", "gene_type"
    )
    
    if (!all(txMap$gencodeTx %in% names(txRR))) {
        stop("Some transcripts do not appear to have corresponding annotation. If using custom annotation, please ensure the GTF has all transcripts present in the FASTA")
    }
    txRR <- txRR[txMap$gencodeTx]
    
    rownames(txMap) <- rownames(txTpm) <- rownames(txNumReads) <- txMap$gencodeTx
} else { # rat; just take transcript ID from txNames, e.g. ENST00000456328.2
    #   Transcript names in the GTF exclude the transcript version suffix, which
    #   is part of the names from the transcripts FASTA. Use the GTF's
    #   convention
    txNames <- ss(txNames, "\\.")
    
    #   Here we circumvent the issue with the transcripts FASTA including some
    #   'primary' transcripts by subsetting to those seen in the GTF ('main'
    #   transcripts only)
    in_gtf = txNames %in% names(txRR)
    txNames = txNames[in_gtf]
    txTpm = txTpm[in_gtf,]
    txNumReads = txNumReads[in_gtf,]
    # if (!all(txNames %in% names(txRR))) {
    #     stop("Some transcripts do not appear to have corresponding annotation. If using custom annotation, please ensure the GTF has all transcripts present in the FASTA")
    # }
    
    txRR <- txRR[txNames]
    
    ## get transcript length from genomic ranges, summing up exon lengths
    # dtex <- as.data.table(subset(gencodeGTF, type=='exon'))
    # setkey(dtex, transcript_id, start)
    # dtx <- dtex[, .(numexons=.N, tlen=sum(end-start+1)), by=transcript_id]
    # setkey(dtx, transcript_id)
    # stopifnot(identical(names(txRR), dtx[txNames]$transcript_id))
    # txMap <- data.frame(
    #     gencodeTx=txRR$transcript_id, txLength=dtx[txNames]$tlen,
    #     gencodeID=txRR$gene_id, Symbol=txRR$gene_name
    # )
    
    rownames(txTpm) <- rownames(txNumReads) <- txRR$transcript_id
}

rm(txNames)

###############################################################################
#  Process alignment logs from HISAT2 or STAR, as applicable
###############################################################################

#  Return a desired value from a HISAT2 alignment summary.
#
#  Here 'log_lines' is a whitespace-stripped character vector containing each
#  line in the summary; 'phrase' is a string that should be present in the
#  lines of interest; 'get_perc' is a logical determining whether to extract
#  the percentage enclosed in parenthesis (rather than the first integer
#  present in the line); 'num_hits' is the number of times 'phrase' should be
#  present in a line.
parse_hisat_line <- function(log_lines, phrase, get_perc, num_hits = 1) {
    index <- grep(phrase, log_lines)
    if (length(index) != num_hits) {
        stop(
            paste0(
                "Phrase '", phrase, "' in HISAT2 log was expected ", num_hits,
                " time(s) but was observed ", length(index), " times."
            )
        )
    }
    
    if (get_perc) {
        #  Grab a percentage enclosed in parenthesis
        val <- ss(ss(log_lines[index], "\\(", 2), "%")
    } else {
        #  Grab a raw integer (the first number in every line)
        val <- gsub("%", "", ss(log_lines[index], " "))
    }
    
    val <- sum(as.numeric(val))
    
    return(val)
}

hisatStats <- function(logFile) {
    y <- scan(logFile,
              what = "character", sep = "\n",
              quiet = TRUE, strip = TRUE
    )
    if (opt$paired) {
        perc_paired <- parse_hisat_line(y, "%) were paired; of these:", TRUE)
        if (perc_paired == 100) {
            ## 100% of reads paired
            reads <- 2 * parse_hisat_line(y, "%) were paired; of these:", FALSE)
            unaligned <- parse_hisat_line(y, "%) aligned 0 times", FALSE)
        } else {
            ## Combo of paired and unpaired (from trimming w/ '--keep_unpaired')
            reads <- 2 * parse_hisat_line(y, "%) were paired; of these:", FALSE) +
                parse_hisat_line(y, "%) were unpaired; of these:", FALSE)
            unaligned <- parse_hisat_line(y, "%) aligned 0 times", FALSE, 2)
        }
    } else {
        reads <- parse_hisat_line(y, "reads; of these:", FALSE)
        unaligned <- parse_hisat_line(y, "%) aligned 0 times", FALSE)
    }
    
    o <- data.frame(
        numReads = reads,
        numMapped = reads - unaligned,
        numUnmapped = unaligned,
        overallMapRate = parse_hisat_line(y, "% overall alignment rate", FALSE) / 100
    )
    
    if (opt$paired) {
        o$concordMapRate <- (parse_hisat_line(y, "%) aligned concordantly exactly 1 time", TRUE) +
                                 parse_hisat_line(y, "%) aligned concordantly >1 times", TRUE)) / 100
    }
    
    return(o)
}

#  Return a 1-row data frame of STAR alignment metrics, given a single
#  sample ID and a character vector of metric names to extract. Note this
#  vector is "fixed" (the function is hardcoded to work for a particular
#  value of 'metric_names').
starStats <- function(id, metric_names) {
    #  Infer log path from sample ID, then read
    star_log <- readLines(paste0(id, "_STAR_alignment.log"))
    
    #  Infer whether trimming was performed from whether a post-trimming FastQC
    #  report exists for this sample ID
    if (opt$paired) {
        is_trimmed <- file.exists(paste0(id, "_1_trimmed_summary.txt"))
    } else {
        is_trimmed <- file.exists(paste0(id, "_trimmed_summary.txt"))
    }
    
    #  Grab the lines including each metric of interest
    key_lines <- star_log[sapply(metric_names, function(n) grep(n, star_log))]
    
    #  Extract the numeric value for each metric in those lines
    metric_values <- as.numeric(ss(key_lines, "\t", 2))
    
    o <- data.frame(
        "trimmed" = is_trimmed,
        "numReads" = metric_values[1],
        "numMapped" = sum(metric_values[5:6]),
        "numUnmapped" = sum(metric_values[2:4]),
        "overallMapRate" = sum(metric_values[5:6]) / metric_values[1]
    )
    
    return(o)
}

#  Extract alignment stats, and add to 'metrics'
if (opt$star) {
    metric_names <- c(
        "Number of input reads",
        "Number of reads unmapped: too many mismatches",
        "Number of reads unmapped: too short",
        "Number of reads unmapped: other",
        "Uniquely mapped reads number",
        "Number of reads mapped to multiple loci"
    )
    
    alignment_stats <- lapply(metrics$SAMPLE_ID, function(id) starStats(id, metric_names))
    alignment_stats <- do.call(rbind, alignment_stats)
} else { # HISAT2 is used for alignment
    logFiles <- file.path(work_dir, paste0(metrics$SAMPLE_ID, "_align_summary.txt"))
    names(logFiles) <- metrics$SAMPLE_ID
    
    alignment_stats <- do.call(rbind, lapply(logFiles, hisatStats))
}

metrics <- cbind(metrics, alignment_stats)

###############################################################################
#  Process BAM files
###############################################################################

### confirm total mapping
bamFile <- file.path(work_dir, paste0(metrics$SAMPLE_ID, "_sorted.bam"))
metrics$bamFile <- file.path(opt$output, "alignment", "bam_sort", bamFile)

metrics$totalMapped <- unlist(bplapply(bamFile, getTotalMapped,
                                       chrs = chr_names,
                                       BPPARAM = MulticoreParam(opt$cores)
))
metrics$mitoMapped <- unlist(bplapply(bamFile, getTotalMapped,
                                      chrs = mito_chr,
                                      BPPARAM = MulticoreParam(opt$cores)
))
metrics$mitoRate <- metrics$mitoMapped / (metrics$mitoMapped + metrics$totalMapped)

###############################################################################
#  Get reference data regarding genes and exons from GTF
###############################################################################

gencodeGENES <- gencodeGTF[gencodeGTF$type == "gene", ]
names(gencodeGENES) <- gencodeGENES$gene_id

if (opt$organism == "mm10") {
    gencodeEXONS <- as.data.frame(gencodeGTF)[which(gencodeGTF$type == "exon"), c("seqnames", "start", "end", "exon_id")]
    names(gencodeEXONS) <- c("Chr", "Start", "End", "exon_gencodeID")
} else if (opt$organism == "rat") {
    gencodeEXONS <- as.data.frame(gencodeGTF)[which(gencodeGTF$type == "exon"), c("seqnames", "start", "end", "exon_id")]
    names(gencodeEXONS) <- c("Chr", "Start", "End", "exon_ensemblID")
} else { # human
    gencodeEXONS <- as.data.frame(gencodeGTF)[which(gencodeGTF$type == "exon"), c("seqnames", "start", "end", "gene_id", "exon_id")]
    names(gencodeEXONS) <- c("Chr", "Start", "End", "gene_id", "exon_gencodeID")
}

if (opt$organism %in% c("hg19", "hg38")) {
    ### exons in PAR regions
    par_y <- grep("PAR_Y", gencodeEXONS$gene_id)
    gencodeEXONS$exon_gencodeID[par_y] <- paste0(gencodeEXONS$exon_gencodeID[par_y], "_PAR_Y")
    gencodeEXONS[, -match("gene_id", colnames(gencodeEXONS))]
}


###############################################################################
#  Read in gene counts
###############################################################################

### read in annotation ##
geneFn <- file.path(
    work_dir, paste0(metrics$SAMPLE_ID, "_", opt$anno_suffix, "_Genes.counts")
)
names(geneFn) <- metrics$SAMPLE_ID

geneMap <- read.delim(geneFn[1], skip = 1, as.is = TRUE)[
    , c("Geneid", "Length")
]

## organize gene map
indices <- match(geneMap$Geneid, gencodeGENES$gene_id)
if (any(is.na(indices))) {
    stop("Not all genes observed in FeatureCounts output are in GTF.")
}

#   Read in gene coordinates and strand from GTF
geneMap$Chr <- as.character(seqnames(gencodeGENES)[indices])
geneMap$Start <- start(gencodeGENES)[indices]
geneMap$End <- end(gencodeGENES)[indices]
geneMap$Strand <- as.character(strand(gencodeGENES)[indices])
rownames(geneMap) <- geneMap$Geneid

if (opt$organism == "rat") {
    geneMap$ensemblID <- geneMap$Geneid
} else {
    geneMap$gencodeID <- geneMap$Geneid
    geneMap$ensemblID <- ss(geneMap$Geneid, "\\.")
    geneMap$gene_type <- gencodeGENES$gene_type[indices]
}
geneMap$Geneid <- NULL

#  Get the 'Symbol' and 'EntrezID' columns
if (opt$organism %in% c("hg19", "hg38")) {
    geneMap$Symbol <- gencodeGENES$gene_name[indices]
    geneMap$EntrezID <- mapIds(org.Hs.eg.db, geneMap$Symbol, "ENTREZID", "SYMBOL")
} else if (opt$organism == "mm10") {
    temp <- gencodeGENES$gene_name[indices]
    geneMap$EntrezID <- mapIds(org.Mm.eg.db, temp, "ENTREZID", "SYMBOL")
    geneMap$Symbol <- mapIds(org.Mm.eg.db, temp, "MGI", "SYMBOL")
} else { # 'rat'
    geneMap$Symbol <- gencodeGENES$gene_name[indices]
    geneMap$EntrezID <- mapIds(org.Rn.eg.db, geneMap$Symbol, "ENTREZID", "SYMBOL")
}

## counts
geneCountList <- bplapply(geneFn,
                          function(x) {
                              cat(".")
                              read.delim(pipe(paste("cut -f7", x)), as.is = TRUE, skip = 1)[, 1]
                          },
                          BPPARAM = MulticoreParam(opt$cores)
)
geneCounts <- do.call("cbind", geneCountList)
rownames(geneCounts) <- rownames(geneMap)
geneCounts <- geneCounts[, metrics$SAMPLE_ID] # put in order

# number of reads assigned
geneStatList <- lapply(paste0(geneFn, ".summary"),
                       read.delim,
                       row.names = 1
)
geneStats <- do.call("cbind", geneStatList)
colnames(geneStats) <- metrics$SAMPLE_ID
metrics$totalAssignedGene <- as.numeric(geneStats[1, ] / colSums(geneStats))

#  Add all the other stats from featureCounts at the gene level
geneStats_t <- t(geneStats)
colnames(geneStats_t) <- paste0("gene_", colnames(geneStats_t))
metrics <- cbind(metrics, geneStats_t)

###############################################################################
#  Create transcript RangedSummarizedExperiment object
###############################################################################

#   SPEAQeasy settings were written to a CSV, which will be read in here and
#   used to populate the metadata of each RSE
meta_csv <- read.csv(file.path(work_dir, "params.csv"), header = FALSE)
rse_meta <- meta_csv[, 2]
names(rse_meta) <- meta_csv[, 1]
rse_meta <- list("SPEAQeasy_settings" = as.list(rse_meta))

rse_tx <- SummarizedExperiment(
    assays = list("counts" = txNumReads, "tpm" = txTpm),
    colData = metrics, rowRanges = txRR, metadata = rse_meta
)

save(
    rse_tx,
    file = file.path(
        out_dir, paste0("rse_tx_", EXPNAME, "_n", N, ".Rdata")
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
