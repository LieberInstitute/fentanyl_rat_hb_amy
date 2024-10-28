library(tidyverse)
library(here)
library(sessioninfo)

man_path = here('processed-data', '02_SPEAQeasy', 'samples.manifest')
pheno_path = here('raw-data', 'sample_metadata_and_QCmetrics.csv')
out_path = here('processed-data', '09_SRA', 'metadata.tsv')

design_description_text = paste(
    "Drugs",
    "Fentanyl Citrate (Cayman Chemical) was diluted in 0.9% sterile saline at 10 ug/mL. Brevital sodium (Henry Schein) was diluted to 10 mg/mL with 0.9% sterile saline.",
    "Surgeries",
    "Rats were anesthetized with isoflurane (3-5% induction, 1-2% maintenance) and administered the analgesic rimadyl (5 mg/kg, SC, Zoetis) and the antibiotic cefazolin (70 mg/kg, SC, West-Ward) prior to catheter implantation. A silastic jugular catheter, constructed as described previously (Thomsen & Caine, 2005), was threaded under the skin from a subcutaneous base in the mid-scapular region and inserted into the right jugular vein (Thomsen & Caine, 2005). Catheters were flushed daily with ~0.1 mL of sterile saline solution containing 10 mg/mL gentamicin sulfate (VetOne) and 100 IU/mL Heparin (Sagent). Once behavioral training began, catheters were flushed before and after self-administration (and daily on rest days), and were assessed daily for blood return. In the event that blood return was absent, rats were administered 0.1 mL of brevital sodium (10 mg/mL in sterile saline, Henry Schein). If the rat did not become ataxic within 10 seconds, the catheter was considered not patent and the rat was removed from the study.",
    "Self-Administration training",
    "After recovery, rats were trained in standard Med Associates operant chambers housed inside sound-attenuating chambers (Med Associates, St Albans, VT, USA). Operant chambers were fitted with two retractable levers (active and inactive), stimulus lights above each lever, speakers, and a house light. Boxes were controlled using Med-PC IV software (Med Associates, St Albans, VT, USA).",
    "Rats were trained to lever-press on a fixed-ratio 1 schedule, initially in short-access 2-hour sessions, and then they progressed to long-access 6-hour sessions (~6 days per week). Both levers were extended for the entire duration of the sessions. The stimulus light above the active lever indicated drug (or saline) availability. For the experimental group of rats, active lever presses resulted in infusion of 0.5 µg fentanyl dissolved in 50 µl sterile saline delivered over 2.8 seconds (Fragale et al., 2020). Infusion was accompanied by a 2.8 second tone, and followed by a 20 second timeout during which the stimulus light and house light were extinguished. After 20 seconds elapsed, the house light and the stimulus light were re-illuminated to indicate availability of drug (or saline). Training occurred for 22-24 sessions, after which the rats’ brains were harvested and rats were assessed for self-administration performance.",
    "Brain tissue extraction and bulk RNA-sequencing",
    "Rats were rapidly decapitated without anesthesia 60-90 minutes after the last long-access SA session, and extracted brains were fresh frozen in isopentane and stored at −80 °C. For tissue punching, brains were sectioned into ~2 mm coronal slabs using a rat brain matrix (Stainless Steel Alto Coronal, 1.0 mm matrix, Small Rat 175-300gm), and bilateral tissue punches were collected from the habenula (#39443001RM, Leica Biosystems, nominal diameter 1.25 mm) and amygdala (nominal diameter 1.25 mm) on petri dishes placed in a dry ice bucket. Punched tissue was ejected into tubes sitting in dry ice, and sample tubes were stored at −80 °C until ready for RNA processing. Total RNA was isolated by using TRIzol Reagent homogenization (#15596018, Invitrogen) and chloroform layer separation. RNA was then purified using an RNeasy Micro (#74004, Qiagen) kit with an RNase-Free DNase step (Mat. No. 1023460, Qiagen) according to manufacturer’s instructions. RNA concentration and purity were measured using a NanoDrop Eight (ThermoScientific). RNA quality control assays were performed on the Agilent 2100 Bioanalyzer, and the RNA integrity number for all samples (±SEM) ranged from 7.2 to 8.6 (mean = 7.9 ± 0.06).",
    "Ribosomal RNA depletion and library preparation (Illumina Ribo-Zero) was performed. ERCC spike-ins were included. Specifically, Illumina Stranded Total RNA Prep with Ribo-Zero Plus Kit was used for library preparation from the habenula samples (10 ng of RNA input), and TruSeq Stranded Total RNA with Ribo-Zero Human/Mouse/Rat Gold Kit was used for library preparation from the amygdala samples (100 ng of RNA input). RNA-sequencing was carried out by Psomagen with the following configuration: ~80M paired reads per sample on an Illumina NovaSeq 6000 S4 (150 bp PE)."
)

meta_df = read.table(
        man_path,
        col.names = c('filename', 'md1', 'filename2', 'md2', 'sample_name')
    ) |>
    as_tibble() |>
    select(filename, filename2, sample_name) |>
    left_join(
        read_csv(pheno_path, show_col_types = FALSE) |>
            rename(sample_name = SAMPLE_ID),
        by = 'sample_name'
    ) |>
    mutate(
        library_ID = sample_name,
        title = sprintf(
            'Bulk RNA-seq of rat %s: %s-treated group',
            tolower(Brain_Region), tolower(Substance)
        ),
        library_strategy = 'RNA-Seq',
        library_source = 'TRANSCRIPTOMIC',
        library_selection = 'Inverse rRNA',
        library_layout = 'paired',
        platform = 'ILLUMINA',
        instrument_model = 'Illumina NovaSeq 6000',
        design_description = factor(design_description_text),
        filetype = 'fastq'
    ) |>
    #   Re-order columns to match SRA's expectations
    select(
        sample_name, library_ID, title, library_strategy,
        library_source, library_selection, library_layout, platform,
        instrument_model, design_description, filetype, filename,
        filename2
    )

#   Also check disk usage (to make sure we comply with SRA requirements)
disk_usage = sapply(
    c(meta_df$filename, meta_df$filename2),
    function(x) {
        as.integer(
            system(sprintf('du -k $(readlink %s) | cut -f 1', x), intern = TRUE)
        )
    }
)
message(
    sprintf(
        "FASTQ files occupy a total of %sGB.", round(sum(disk_usage) / 1e6, 1)
    )
)

#   Individual files must be less than 100GB
stopifnot(all(disk_usage < 100e6))

write_tsv(meta_df, out_path)

session_info()
