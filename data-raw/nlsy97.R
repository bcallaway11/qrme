# =============================================================================
# Title: Build nlsy97 package dataset
# Description: Loads the NLSY97 father-son analysis sample used in Callaway,
#   Li, Murtazashvili, and Tsyawo (2026) and saves it as package data. The
#   raw .dta file is not distributed with the package; this script lives in
#   data-raw/ for provenance only.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

library(haven)

raw_file <- file.path(
  "~", "Dropbox", "QR Measurement Error", "application", "nlsy97-data",
  "Father_SonPairs1997_2015NLS_20191206.dta"
)

nlsy97 <- as.data.frame(read_dta(raw_file))

# Apply the standard sample restriction used in the paper
nlsy97 <- subset(nlsy97, !is.na(educF))

# Add log-income variables used throughout the paper
nlsy97$lci <- log(nlsy97$incS)   # log child income
nlsy97$lpi <- log(nlsy97$incF)   # log parent (father) income

# Drop the raw id columns — not needed for the analysis
nlsy97$id_hh    <- NULL
nlsy97$id_youth <- NULL

usethis::use_data(nlsy97, overwrite = TRUE, compress = "xz")
