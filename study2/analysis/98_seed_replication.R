# 98 (simulated mode, optional, not part of run_all.R): does recovery hold across
# simulated samples, not just the one 99 checks? Repeats 00_simulate and 01_prepare
# with SF2_REPLICATIONS seeds (default 8), each in its own data directory under
# <data>/replications/, and for each sample fits the H1 model and computes the
# quantities whose single-sample checks are interval coverage: the H6 d's, the H7
# ICCs and the following-week coefficients. Writes 98_seed_replication.csv (one row
# per seed and quantity) and 98_seed_replication_summary.csv (coverage per quantity).
#   SF2_MODE=simulated Rscript study2/analysis/98_seed_replication.R

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
stopifnot(SF2$mode == "simulated")
K <- as.integer(Sys.getenv("SF2_REPLICATIONS", "8"))
seeds <- 20261100 + seq_len(K)
root <- sf2_path("replications")

one <- function(seed) {
  dd <- file.path(root, paste0("seed_", seed))
  Sys.setenv(SF2_MODE = "simulated", SF2_DATA_DIR = dd, SF2_SEED_SIM = seed, SF2_CORES = "1")
  for (s in c("00_simulate", "01_prepare")) {
    st <- system2(file.path(R.home("bin"), "Rscript.exe"), file.path("study2", "analysis", paste0(s, ".R")), stdout = FALSE, stderr = FALSE)
    if (st != 0) return(data.table(seed = seed, quantity = s, error = "script failed"))
  }
  tr <- readRDS(file.path(dd, "derived", "simulation_truth_simulated.rds"))
  A <- readRDS(file.path(dd, "derived", "sf2_simulated.rds"))
  D <- A$defs; M <- A$meta; S <- A$scores; F4 <- D$factors4; items <- D$sf_items
  out <- list()
  # H1: loadings against the estimand, the random-intercept variance
  member <- sf2_member(items, D$sf_factor, F4)
  f <- sf2_sem(sf2_esem_syntax(items, F4, TRUE), as.data.frame(A$sf_raw), items, target = sf2_target(member))
  K4 <- sf2_orient(sf2_loadings(f, F4, items), member, D$sf_rev)$keyed
  ri <- sf2_ri(f)
  out$h1 <- data.table(seed = seed, quantity = c("loadings: min Tucker congruence", "loadings: RMSE", "RI variance"),
                       estimate = c(min(diag(psych::factor.congruence(K4, tr$keyed_std[items, F4]))), sqrt(mean((K4 - tr$keyed_std[items, F4])^2)), ri$variance),
                       lower = c(NA, NA, ri$variance - 1.96 * ri$se), upper = c(NA, NA, ri$variance + 1.96 * ri$se),
                       truth = c(1, 0, tr$ri_variance))
  # H6 d
  g <- ifelse(M$practice_group %in% "dedicated", 1L, ifelse(M$practice_group %in% "never", 0L, NA_integer_))
  out$h6 <- rbindlist(lapply(names(tr$h6_d), function(s) { x <- sf2_kg(S[[s]], g, data.table(posting = M$posting_contemplative))
    data.table(seed = seed, quantity = paste("H6 d", s), estimate = x$d, lower = x$d_lower, upper = x$d_upper, truth = tr$h6_d[[s]]) }))
  # H7 ICC
  P2 <- A$retest$pairs; fi <- function(f) items[D$sf_factor == f]
  SC <- list(Access = fi("Access"), Valuing = fi("Valuing"), AccessValuing = c(fi("Access"), fi("Valuing")),
             Understanding = fi("Understanding"), LettingBe = fi("LettingBe"))
  out$h7 <- rbindlist(lapply(names(tr$icc), function(s) {
    x <- cbind(sf2_score(A$retest$t1, SC[[s]], D$sf_rev), sf2_score(A$retest$t2, SC[[s]], D$sf_rev)); x <- x[complete.cases(x), ]
    ic <- irr::icc(x, "twoway", "agreement", "single")
    data.table(seed = seed, quantity = paste("H7 ICC", s), estimate = ic$value, lower = ic$lbound, upper = ic$ubound, truth = tr$icc[[s]]) }))
  # the following week
  av1 <- S$AccessValuing[match(P2$participant, M$participant)]
  out$fw <- rbindlist(lapply(names(tr$following_week), function(v) {
    d <- data.table(y2 = A$retest$past_week_t2[, v], av = av1, y1 = A$retest$past_week_t1[, v]); d <- d[complete.cases(d)]
    m <- lm(scale(y2) ~ scale(av) + scale(y1), d); ci <- confint(m)["scale(av)", ]
    data.table(seed = seed, quantity = paste("following week", v), estimate = coef(m)[["scale(av)"]], lower = ci[[1]], upper = ci[[2]],
               truth = tr$following_week[[v]]) }))
  unlink(file.path(dd, "raw"), recursive = TRUE)   # the simulated extractions are large; the derived files stay
  rbindlist(out, fill = TRUE)
}

t0 <- Sys.time()
res <- rbindlist(sf2_par_lapply(seeds, one, data = list(root = root), cores = min(K, SF2$cores)), fill = TRUE)
res[, covered := lower <= truth & upper >= truth]
fwrite(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)], sf2_out("98_seed_replication.csv"))
summ <- res[, .(seeds = .N, mean_estimate = round(mean(estimate), 3), truth = round(mean(truth), 3),
                coverage = if (all(is.na(covered))) NA_real_ else round(mean(covered, na.rm = TRUE), 3),
                range = paste(round(range(estimate), 3), collapse = " to ")), by = quantity]
fwrite(summ, sf2_out("98_seed_replication_summary.csv"))
print(summ, nrows = 100)
cov <- res[!is.na(covered), mean(covered)]
sf2_log(K, " seeds in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min; interval coverage over all quantities and seeds: ",
        round(cov, 3), " (", sum(res$covered, na.rm = TRUE), " of ", sum(!is.na(res$covered)), ")")
