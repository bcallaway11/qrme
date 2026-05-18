# =============================================================================
# Title: Build nlsy97_tsme_fit package dataset
# Description: Trims the full tsme() result from the NLSY97 intergenerational
#   mobility application (Frank copula, Laplace ME, n_mix=1, n_boot=200) to
#   keep only the summary outputs needed for print/summary/autoplot. The large
#   rqs objects (nome_qyx, nome_qtx, qrytx) are dropped to keep the dataset
#   under 1 MB.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

raw_file <- file.path(
  "~", "Dropbox", "QR Measurement Error", "application", "results",
  "tsme_main_2026-05-16_17-14-03.rds"
)

obj <- readRDS(raw_file)
fit <- obj$fit

# Drop the large rqs objects — not needed for print/summary/autoplot
fit$nome_qyx <- NULL
fit$nome_qtx <- NULL
fit$qrytx    <- NULL

# Ensure the S3 class is set (may be missing if saved before class() was added)
class(fit) <- "tsme"

nlsy97_tsme_fit <- fit

usethis::use_data(nlsy97_tsme_fit, overwrite = TRUE, compress = "xz")
