# =============================================================================
# 0_Download_data.R
# =============================================================================
#
# START HERE — This is the first script to run for reproducing the analysis.
#
# Where things live:
#   - CODE  -> this GitHub repository (the analysis scripts, notebooks, and
#              Stan models you already have after cloning).
#   - DATA  -> downloaded by this script from two external sources:
#                * Zenodo (record 20754663) : processed + raw data, intermediate
#                                             model outputs, and figures.
#                * NCBI SRA (PRJNA1426049)  : raw sequencing reads.
#
# This script fetches the data (the code is already present from GitHub) and
# wires everything into the project. It performs three main steps:
#
#   1. Download the data archive from Zenodo (record 20754663)
#      - Code_and_raw_data.zip is unpacked into the project, populating the
#        data/, and Intermediate_data/ directories.
#
#   2. Download raw sequencing data from NCBI SRA (BioProject PRJNA1426049)
#      - Fetches paired-end FASTQ files via ENA and saves to SRA/fastq/
#      - Run metadata is saved to SRA/metadata/
#      - All FASTQ files are concatenated into SRA/combined_R2.fastq.gz
#        and decompressed to SRA/combined_R2.fastq
#      - Helper scripts live in code/sra_python/
#
#   3. Open analysis notebooks (Quarto .qmd files) for the next steps:
#      - 1_Run_QM_qPCR.qmd        : Joint QM-qPCR model of hake survey samples
#      - 2_sdmTMB_smooths_13sp.qmd : Spatial smooths of joint model output
#      - 3_All_Figures.qmd         : Generate all manuscript figures
#
# Prerequisites:
#   - R packages: here, jsonlite, R.utils (auto-installed if missing)
#   - Python 3 with internet access (for SRA download scripts)
#   - ~1 hour timeout is set for large file downloads
#
# Project directory structure:
#   project_root/
#   ├── code/              Analysis scripts and Stan models   (from GitHub)
#   ├── data/              Raw qPCR, metabarcoding, metadata   (from Zenodo)
#   ├── Intermediate_data/ Model outputs and processed data (from Zenodo)
#   ├── raw_plots/         Figures created with other software (from Zenodo)
#   └── SRA/               Raw sequencing data                 (from NCBI SRA)
#       ├── fastq/         Per-run FASTQ files
#       ├── metadata/      SRA run metadata
#       ├── combined_R2.fastq.gz  Concatenated reverse reads (compressed)
#       └── combined_R2.fastq     Concatenated reverse reads (decompressed)
#
# =============================================================================
library(here)


# --- Download files from Zenodo --------------------------------------------------------------

options(timeout = 3600) 
record_id <- "20754663"

meta <- jsonlite::fromJSON(paste0("https://zenodo.org/api/records/", record_id))

files <- meta$files

df <- data.frame(
  key     = files$key,
  size_Mb = round(files$size / 1024^2, 2),
  size_Gb = round(files$size / 1024^3, 2),
  link    = files$links$self)

print(df, right = FALSE);rm(df)

# This will download the files into the project directory
download.file(
  url = files$links$self[1],   # first file; adjust index as needed
  destfile = files$key[1],
  mode = "wb")


# --- Unzip files from Zenodo -----------------------------------------------------------------
tmp <- tempfile()
unzip("Data.zip", exdir = tmp, unzip = "unzip")

unlink(file.path(tmp, "__MACOSX"), recursive = TRUE)
inner <- file.path(tmp)
files <- list.files(inner, all.files = TRUE, no.. = TRUE, full.names = TRUE)
file.copy(files, ".", recursive = TRUE)

unlink(tmp, recursive = TRUE)   # clean up the temp copy


# --- Open one of the analysis notebooks -------------------------------------
# file.edit(here("Code", "1_Run_QM_qPCR.qmd"))
# file.edit(here("Code", "2_sdmTMB_smooths_13sp.qmd"))
file.edit(here("Code", "3_All_Figures.qmd"))


# --- Download raw sequencing data from the SRA (optional) ------------------------------------
# Download all SRA data
system2("python3", here("code", "sra_python", "download_sra.py"))

# Concatenate all SRA fastq data
system2("python3", here("code", "sra_python", "concatenate_fastq.py"))

# Unzipt the SRA fastq concatenated data
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
R.utils::gunzip(here('SRA','combined_R2.fastq.gz'), remove = FALSE)

