# =============================================================================
# 00_Download_data.R
# =============================================================================
#
# START HERE — This is the first script to run for reproducing the analysis.
#
# Where things live:
#   - CODE  -> this GitHub repository (Code/: analysis scripts and Stan model
#              you already have after cloning).
#   - DATA  -> downloaded by this script from Zenodo (record 21141367):
#              the empirical + simulated input data and all simulation/model
#              output, packaged as Smart-eDNA-sampling-data.zip.
#
# Because the archive is ~19 GB, Zenodo's web/API uploaders can't reliably
# push it as a single file, so it is stored on Zenodo as 10 sequential parts
# (Smart-eDNA-sampling-data.zip.part-aa ... part-aj). This script performs
# four steps:
#
#   1. Download every part from the Zenodo record and verify each one's
#      checksum against what Zenodo reports (catches a corrupted/partial
#      download before wasting time reassembling it).
#   2. Reassemble the parts into a single Smart-eDNA-sampling-data.zip.
#   3. Unzip the archive into the project, populating Data/ and Output/.
#   4. Open the next script in the pipeline (03_Plots.R by default — see
#      below for why).
#
# Prerequisites:
#   - R packages: here, jsonlite (auto-installed if missing)
#   - ~2 hour timeout is set per file for the large download
#
# Project directory structure:
#   project_root/
#   ├── Code/     Analysis scripts and Stan model      (from GitHub)
#   ├── Data/     Input eDNA data                       (from Zenodo)
#   ├── Output/   Simulation and model output            (from Zenodo)
#   └── Plots/    Manuscript figures, created by 03_Plots.R
#
# =============================================================================
library(here)


# --- Download files from Zenodo ---------------------------------------------

options(timeout = 7200)
record_id <- "21141367"

# Public, unauthenticated endpoint — works once the record is published.
# If you need to pull from the record while it is still an unpublished draft,
# set a ZENODO_TOKEN environment variable (never hardcode a token here) and
# use the deposit API instead:
#   token <- Sys.getenv("ZENODO_TOKEN")
#   meta  <- jsonlite::fromJSON(paste0(
#     "https://zenodo.org/api/deposit/depositions/", record_id,
#     "?access_token=", token))
#   files <- data.frame(key = meta$files$filename, size = meta$files$filesize,
#                        checksum = meta$files$checksum,
#                        link = meta$files$links$download)
meta  <- jsonlite::fromJSON(paste0("https://zenodo.org/api/records/", record_id))
files <- meta$files

df <- data.frame(
  key     = files$key,
  size_Mb = round(files$size / 1024^2, 2),
  size_Gb = round(files$size / 1024^3, 2),
  link    = files$links$self)

print(df, right = FALSE); rm(df)

# Download every part (the archive is split into 10 parts; adjust if the
# record ever changes to a single file)
for (i in seq_along(files$key)) {
  cat("Downloading", files$key[i], "\n")
  download.file(
    url      = files$links$self[i],
    destfile = files$key[i],
    mode     = "wb")
}

# Verify each downloaded part against the checksum Zenodo reports, so a
# truncated/corrupted download is caught now rather than after reassembly
local_md5    <- tools::md5sum(files$key)
expected_md5 <- sub("^md5:", "", files$checksum)
if (!all(unname(local_md5) == expected_md5)) {
  stop("Checksum mismatch on: ",
       paste(files$key[unname(local_md5) != expected_md5], collapse = ", "),
       " — delete the file(s) and re-run this script.")
}
cat("All", length(files$key), "parts downloaded and verified.\n")


# --- Reassemble the parts into a single zip ----------------------------------

zip_name <- "Smart-eDNA-sampling-data.zip"
parts    <- sort(files$key[grepl("\\.part-", files$key)])

con_out <- file(zip_name, "wb")
for (p in parts) {
  con_in <- file(p, "rb")
  while (length(chunk <- readBin(con_in, "raw", 1024^2 * 64)) > 0) {
    writeBin(chunk, con_out)
  }
  close(con_in)
}
close(con_out)

unlink(parts)  # the individual parts are no longer needed once reassembled


# --- Unzip files from Zenodo --------------------------------------------------

tmp <- tempfile()
unzip(zip_name, exdir = tmp, unzip = "unzip")

unlink(file.path(tmp, "__MACOSX"), recursive = TRUE)
inner <- file.path(tmp)
extracted <- list.files(inner, all.files = TRUE, no.. = TRUE, full.names = TRUE)
file.copy(extracted, ".", recursive = TRUE)

unlink(tmp, recursive = TRUE)   # clean up the temp copy
# unlink(zip_name)               # uncomment to remove the ~19 GB zip once unpacked


# --- Open the next script in the pipeline -------------------------------------
# Output/Complete_db/ (just unpacked) already has everything 03_Plots.R needs,
# so most users can skip straight to reproducing the figures. Only re-run
# 01_GP_sim.R / 02_Data_manip.R if you want to regenerate the simulations and
# intermediate output from scratch (compute-heavy; originally run on an HPC
# cluster).
file.edit(here("Code", "01_GP_sim.R"))
file.edit(here("Code", "02_Data_manip.R"))
file.edit(here("Code", "03_Plots.R"))
