# 30: item and test information for the 8, from a graded response model. EXPLORATORY,
# NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 17 November 2026 (placement card S3-26, from the Fable review),
# after the registered analyses and 10-29 had been seen.
#
# A unidimensional graded response model (Samejima) of the 8, reverse-keyed statements
# recoded (7 - x), marginal ML (mirt, EM), all 720 participants (missing answers left
# missing). The agreement factor is not modelled, so this is a check beside the CFA, not
# a replacement for it. Also fitted to the 7 (no U4).
# Reported: discrimination (a) and thresholds (b1-b5) per statement; item and test
# information on theta from -3 to +3; the theta range where test information >= 5
# (reliability >= .80 on that scale) and its peak; marginal reliability; M2-based fit
# (RMSEA, CFI) where mirt can compute it; the theta of the observed-score percentiles
# (EAP scores) so the information can be read against the sample's spread.
# Writes aggregates only:
#   output/real/30_items.csv, 30_information.csv, 30_summary.csv, 30_information.png,
#   30_exploratory_irt.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/30_exploratory_irt.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "mirt", "ggplot2")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
source("study2b/analysis/config.R")
stopifnot(SF2B$mode == "real")
T0 <- Sys.time()
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; X <- A$scale_raw
SETS <- list(eight = D$sf8, seven = setdiff(D$sf8, "sf2b_07"))
stopifnot(D$source[["sf2b_07"]] == "U4")
TH <- seq(-3, 3, by = .05)

ITEMS <- list(); INFO <- list(); SUMM <- list()
for (s in names(SETS)) {
  it <- SETS[[s]]
  Y <- as.data.frame(X[, it]); for (v in intersect(it, D$reversed)) Y[[v]] <- 7 - Y[[v]]
  m <- mirt(Y, 1, itemtype = "graded", verbose = FALSE, technical = list(NCYCLES = 2000))
  stopifnot(extract.mirt(m, "converged"))
  cf <- coef(m, IRTpars = TRUE, simplify = TRUE)$items
  ITEMS[[s]] <- data.table(set = s, item = it, source = D$source[it], reversed = it %in% D$reversed, cf)
  ii <- sapply(seq_along(it), function(j) iteminfo(extract.item(m, j), Theta = matrix(TH)))
  colnames(ii) <- it
  ti <- rowSums(ii)
  INFO[[s]] <- data.table(set = s, theta = TH, test = ti, ii)
  fs <- fscores(m, method = "EAP", full.scores = TRUE)[, 1]
  m2 <- tryCatch(M2(m, na.rm = TRUE), error = function(e) NULL)
  above5 <- TH[ti >= 5]
  SUMM[[s]] <- data.table(set = s, marginal_reliability = unname(marginal_rxx(m)),
    peak_theta = TH[which.max(ti)], peak_information = max(ti),
    theta_info_ge5_from = if (length(above5)) min(above5) else NA_real_, theta_info_ge5_to = if (length(above5)) max(above5) else NA_real_,
    eap_p05 = quantile(fs, .05), eap_p50 = quantile(fs, .50), eap_p95 = quantile(fs, .95),
    share_eap_outside_info_ge5 = if (length(above5)) mean(fs < min(above5) | fs > max(above5)) else NA_real_,
    m2_rmsea = if (!is.null(m2)) m2$RMSEA else NA_real_, m2_cfi = if (!is.null(m2)) m2$CFI else NA_real_)
  sf2b_log(s, " done")
}
items <- rbindlist(ITEMS, fill = TRUE); info <- rbindlist(INFO, fill = TRUE); summ <- rbindlist(SUMM)

## Figure: test information, the 8 and the 7, with the sample's 5th-95th EAP percentiles ---------
p <- ggplot(info, aes(theta, test, linetype = set)) + geom_line() +
  geom_hline(yintercept = 5, colour = "grey60") +
  annotate("rect", xmin = summ[set == "eight", eap_p05], xmax = summ[set == "eight", eap_p95], ymin = -Inf, ymax = Inf, alpha = .08) +
  labs(x = "Latent trait (theta)", y = "Test information", linetype = NULL,
       caption = "Grey line: information 5 (reliability .80). Shaded: the sample's 5th-95th percentile of theta (the 8).") +
  theme_minimal(base_size = 11)
ggsave(sf2b_out("30_information.png"), p, width = 6.5, height = 4, dpi = 200)

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(items)), sf2b_out("30_items.csv")); fwrite(num(copy(info)), sf2b_out("30_information.csv")); fwrite(num(copy(summ)), sf2b_out("30_summary.csv"))
write_json(list(
  analysis = "30: graded response model of the 8 and the 7; item and test information. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 17 November 2026 after 10-29 were seen (card S3-26).",
  model = "Unidimensional graded response model, mirt EM, reverse-keyed statements recoded 7 - x; agreement factor not modelled.",
  items = num(copy(items)), summary = num(copy(summ)), mirt_version = as.character(packageVersion("mirt")),
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("30_exploratory_irt.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(num(copy(items), 3)); print(num(copy(summ), 3))
print(num(copy(info), 2)[theta %in% c(-3, -2, -1, 0, 1, 2, 3)])
sf2b_log("30 done")
