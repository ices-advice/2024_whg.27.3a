library(DATRAS)
library(maps); library(mapdata)
library(mgcv)
library(surveyIndex)

source("00_functions.R")

## No correction - WGNSSK 2024 (1983-2023)
##    Data
##    DATRAS exchange files for all years
##    TNG (2004-2023) and TOR (2008-2023) Danish survey data from Rie
outdir <- "surveyindex_no_correction_2024"
years <- 1983:2023
if (! dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
odir <- setwd(outdir)

source("../01_preprocess_data.R")
calc_Index(dat = "WhitingData_subset_noCorrection.RData", output_folder = ".")
make_all_plots(dat = "WhitingData_subset_noCorrection.RData", wd = ".", years = years)
run_retro_lo(dat = "WhitingData_subset_noCorrection.RData", "dat_mqs_grid.RData", wd = ".")