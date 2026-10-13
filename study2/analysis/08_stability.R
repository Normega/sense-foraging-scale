# 08: H7, two-week stability (Section 5, analysis 6), in the follow-up pairs.
#
# For each observed score: ICC(2,1), absolute agreement (irr::icc, two-way, single),
# with a 95% CI; SEM = SD_baseline x sqrt(1 - ICC); MDC95 = 1.96 x sqrt(2) x SEM, in
# raw-score units. H7 holds for Access, Valuing and savouring, and Access + Valuing,
# EACH separately, when the ICC point estimate is >= .70. Understanding and Letting be
# are reported with their ICC and baseline omega (the H1 model's), no threshold.
# PREREG: the decision uses all pairs (out-of-window completions are "kept and
# flagged"); the ICC without them is reported beside it, and without participants
# flagged at baseline (sensitivity). SD_baseline is the baseline SD
# in the same pairs; the full-sample baseline SD version is reported too.
# Also: responders vs non-responders on baseline scores (invited non-responders when the
# invitation list from draw_retest_sample.R is present, otherwise every other analysable
# participant), and the agreement index's stability (short-form pairs only: the
# follow-up has no candidates).

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; P2 <- A$retest$pairs; T1 <- A$retest$t1; T2 <- A$retest$t2
H1 <- readRDS(sf2_derived("h1"))
fi <- function(f) D$sf_items[D$sf_factor == f]
SC <- list(Access = fi("Access"), Valuing = fi("Valuing"), AccessValuing = c(fi("Access"), fi("Valuing")),
           Understanding = fi("Understanding"), LettingBe = fi("LettingBe"))
sc <- function(X) as.data.table(lapply(SC, function(it) sf2_score(X, it, D$sf_rev)))
S1 <- sc(T1); S2 <- sc(T2)
sfp <- SF2_AGREE_PAIRS[grepl("^sf6_", reversed) & grepl("^sf6_", forward)]
agree <- function(X) { both <- !is.na(X[, sfp$reversed]) & !is.na(X[, sfp$forward])
  v <- (X[, sfp$reversed] + X[, sfp$forward]); v[!both] <- NA; ifelse(rowSums(both) > 0, rowSums(v, na.rm = TRUE) / (2 * rowSums(both)), NA) }
S1$agreement_sf <- agree(T1); S2$agreement_sf <- agree(T2)

icc_row <- function(s, keep, label) {
  x <- cbind(S1[[s]], S2[[s]])[keep, , drop = FALSE]; x <- x[complete.cases(x), , drop = FALSE]
  if (nrow(x) < 10) return(data.table(score = s, pairs = label, n = nrow(x)))
  ic <- irr::icc(x, model = "twoway", type = "agreement", unit = "single")
  sdb <- sd(x[, 1]); sem <- sdb * sqrt(1 - ic$value)
  sdf <- sd(A$scores[[if (s == "agreement_sf") "agreement_sf" else s]], na.rm = TRUE)
  data.table(score = s, pairs = label, n = nrow(x), icc = ic$value, lower = ic$lbound, upper = ic$ubound,
             sd_baseline = sdb, sem = sem, mdc95 = 1.96 * sqrt(2) * sem,
             sd_baseline_full_sample = sdf, sem_full_sd = sdf * sqrt(1 - ic$value), mdc95_full_sd = 1.96 * sqrt(2) * sdf * sqrt(1 - ic$value),
             mean_change = mean(x[, 2] - x[, 1]))
}
all_p <- rep(TRUE, nrow(P2)); inw <- P2$in_window
unfl <- !A$meta$flagged[match(P2$participant, A$meta$participant)]
ICC <- rbindlist(c(lapply(c(names(SC), "agreement_sf"), icc_row, keep = all_p, label = "all pairs"),
                   lapply(c(names(SC), "agreement_sf"), icc_row, keep = inw, label = "within 14-21 days"),
                   lapply(c(names(SC), "agreement_sf"), icc_row, keep = unfl, label = "without flagged participants")), fill = TRUE)
om <- setNames(H1$rule$omega, H1$rule$factor)
ICC[, baseline_omega := om[score]]
ICC[, h7 := fifelse(score %in% c("Access", "Valuing", "AccessValuing") & pairs == "all pairs", icc >= .70, NA)]
fwrite(ICC[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("08_stability.csv"))

# Responders vs non-responders on baseline scores.
M <- A$meta; S <- A$scores
inv <- A$retest$invited
comp_set <- if (!is.null(inv)) M$participant %in% inv else M$reached_end_validity
resp <- M$participant %in% P2$participant
RR <- rbindlist(lapply(names(SC), function(s) {
  a <- S[[s]][resp & comp_set]; b <- S[[s]][!resp & comp_set]; a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 2 || length(b) < 2) return(data.table(score = s, n_responders = length(a), n_comparison = length(b)))
  t <- t.test(a, b); sp <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  data.table(score = s, comparison = if (!is.null(inv)) "invited non-responders" else "all other completers",
             n_responders = length(a), n_comparison = length(b), mean_responders = mean(a), mean_comparison = mean(b),
             d = (mean(a) - mean(b)) / sp, p = t$p.value)
}), fill = TRUE)
fwrite(RR[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("08_responders.csv"))

dec <- list(n_pairs = nrow(P2), n_outside_window = sum(!P2$in_window),
            days = as.list(round(quantile(P2$days, c(0, .25, .5, .75, 1)), 1)),
            h7 = ICC[pairs == "all pairs" & score %in% c("Access", "Valuing", "AccessValuing"), .(score, n, icc, lower, upper, supported = h7)],
            invited_list = !is.null(inv))
write_json(dec, sf2_out("08_h7_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
print(ICC[, .(score, pairs, n, icc = round(icc, 3), lower = round(lower, 3), upper = round(upper, 3), mdc95 = round(mdc95, 3), h7)])
print(RR)
