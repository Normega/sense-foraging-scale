# 03: H1, structure (Section 5, model 1; decisions 17, 33).
#
# The registered rule on the 8 (confirmatory) and the 10 (reported beside it): one
# content factor plus the agreement factor; supported when CFI >= .95, AND local fit (no
# residual correlation above .10 in absolute value, or freeing the largest changes the
# content omega by < .03 and every standardized content loading by < .05), AND the
# unbiased SRMR close-fit test (H0: SRMR <= .05) is not rejected at .05.
# Reported, not tested: RMSEA and its 90% CI, model-specific cutoffs, loadings, the
# agreement factor's variance and share, omega, alpha, and the model without the
# agreement factor. Also: the 10 with the three exploratory statements (their loadings
# and agreement loadings), and the 8 without flagged participants.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; X <- A$scale_raw; M <- A$meta

h1 <- list(sf8 = sf2b_h1(X, D$sf8, D$reversed), sf10 = sf2b_h1(X, D$sf10, D$reversed),
           sf8_unflagged = sf2b_h1(X[!M$flagged, ], D$sf8, D$reversed))
for (k in names(h1)) sf2b_log("H1 ", k, ": ", if (isTRUE(h1[[k]]$supported)) "supported" else "not supported",
                              " (CFI ", r3(h1[[k]]$stats$cfi.scaled), ", max residual ", r3(h1[[k]]$max_residual),
                              ", close-fit p ", r3(h1[[k]]$closefit_p), ")")

# Comparison models: without the agreement factor; the 13 statements.
fit_of <- function(items, ri) sf2b_try(sf2b_sem(sf2b_scale_syntax(items, ri), as.data.frame(X[, items]), ordered = items))
noRI <- list(sf8 = fit_of(D$sf8, FALSE), sf10 = fit_of(D$sf10, FALSE))
f13 <- fit_of(D$items, TRUE)

fit_row <- function(label, f) c(list(model = label, converged = !sf2b_failed(f)), sf2b_fit_stats(f))
fits <- rbindlist(list(
  fit_row("the 8, with agreement (H1)", h1$sf8$fit), fit_row("the 8, without agreement", noRI$sf8),
  fit_row("the 10, with agreement", h1$sf10$fit), fit_row("the 10, without agreement", noRI$sf10),
  fit_row("the 13 statements, with agreement", f13), fit_row("the 8, without flagged participants", h1$sf8_unflagged$fit)), fill = TRUE)
fwrite(fits, sf2b_out("03_fit.csv"))

loadings <- rbindlist(list(
  cbind(form = "8", h1$sf8$loadings), cbind(form = "10", h1$sf10$loadings),
  if (!sf2b_failed(f13)) { s <- lavInspect(f13, "std")
    data.table(form = "13", item = D$items, content = s$lambda[D$items, "S"], agreement = s$lambda[D$items, "RI"],
               keyed_content = s$lambda[D$items, "S"] * ifelse(D$items %in% D$reversed, -1, 1),
               residual_variance = diag(s$theta)[D$items]) }), fill = TRUE)
loadings[, `:=`(source = D$source[item], role = D$role[item])]
num_cols <- names(loadings)[vapply(loadings, is.numeric, TRUE)]
loadings[, (num_cols) := lapply(.SD, r3), .SDcols = num_cols]
fwrite(loadings, sf2b_out("03_loadings.csv"))
fwrite(as.data.table(round(h1$sf8$residuals, 3), keep.rownames = "item"), sf2b_out("03_residuals_8.csv"))
fwrite(as.data.table(round(h1$sf10$residuals, 3), keep.rownames = "item"), sf2b_out("03_residuals_10.csv"))

# PREREG: model-specific cutoffs (reported only): the fit indices' distribution at this N in
# data generated from the fitted model itself (the model true), and from it with two
# residual correlations of .30 planted (the power simulation's "small flaw"). Cutoff =
# the 5th percentile of CFI and the 95th of RMSEA and SRMR under each.
dfi <- function(h, items, reps = as.integer(Sys.getenv("SF2B_DFI_REPS", "200"))) {
  if (!isTRUE(h$converged) || reps < 1) return(NULL)
  f <- h$fit; std <- lavInspect(f, "std")
  P <- cov2cor(lavInspect(f, "implied")$cov[items, items]); Dg <- sqrt(diag(lavInspect(f, "implied")$cov[items, items]))
  th <- lavInspect(f, "th"); TH <- lapply(items, function(v) sort(th[startsWith(names(th), paste0(v, "|"))] / Dg[v]))
  res <- 1 - std$lambda[items, "S"]^2 - std$lambda[items, "RI"]^2
  Pf <- P; for (p in list(c(1, 2), c(3, 4))) { i <- p[1]; j <- p[2]; Pf[i, j] <- Pf[j, i] <- P[i, j] + .30 * sqrt(res[i] * res[j]) }
  one <- function(b) {
    set.seed(seed + b)
    Pp <- if (b > reps) Pf else P
    Z <- MASS::mvrnorm(n, rep(0, length(items)), Pp)
    Y <- vapply(seq_along(items), function(j) findInterval(Z[, j], TH[[j]]) + 1L, integer(n)); colnames(Y) <- items
    g <- sf2b_try(sf2b_sem(sf2b_scale_syntax(items), as.data.frame(Y), ordered = items))
    if (sf2b_failed(g)) return(c(pop = if (b > reps) 2 else 1, cfi = NA, rmsea = NA, srmr = NA))
    fm <- fitMeasures(g, c("cfi.scaled", "rmsea.scaled", "srmr"))
    c(pop = if (b > reps) 2 else 1, cfi = fm[[1]], rmsea = fm[[2]], srmr = fm[[3]])
  }
  r <- do.call(rbind, sf2b_par_lapply(seq_len(2 * reps), one, data = list(P = P, Pf = Pf, TH = TH, items = items, n = nrow(X),
                                                                           reps = reps, seed = SF2B$seed + 300, one = NULL)))
  r <- as.data.table(r)
  r[, .(population = c("model true", "two residual correlations of .30")[pop[1]], reps = .N, converged = sum(!is.na(cfi)),
        cfi_p05 = r3(quantile(cfi, .05, na.rm = TRUE)), rmsea_p95 = r3(quantile(rmsea, .95, na.rm = TRUE)),
        srmr_p95 = r3(quantile(srmr, .95, na.rm = TRUE))), by = pop][, pop := NULL][]
}
cut8 <- dfi(h1$sf8, D$sf8)
if (!is.null(cut8)) fwrite(cut8, sf2b_out("03_model_specific_cutoffs_8.csv"))

alpha8 <- function(items) { Z <- X[, items]; rv <- intersect(items, D$reversed); Z[, rv] <- 7 - Z[, rv]
  suppressWarnings(psych::alpha(Z[complete.cases(Z), ], check.keys = FALSE, warnings = FALSE)$total$raw_alpha) }
summ <- function(h, items) {
  if (!isTRUE(h$converged)) return(list(converged = FALSE, error = h$error, supported = FALSE))
  list(converged = TRUE, supported = h$supported, n = h$stats$ntotal,
       cfi = h$stats$cfi.scaled, cfi_ok = h$cfi_ok,
       max_residual = r3(h$max_residual), max_residual_pair = h$max_residual_pair,
       freeing_it_changes_omega_by = r3(h$d_omega), freeing_it_changes_a_loading_by_up_to = r3(h$d_lambda), local_ok = h$local_ok,
       usrmr = r3(h$usrmr), usrmr_ci = r3(h$usrmr_ci), closefit_p = r3(h$closefit_p), close_ok = h$close_ok,
       rmsea = h$stats$rmsea.scaled, rmsea_ci = c(h$stats$rmsea.ci.lower.scaled, h$stats$rmsea.ci.upper.scaled),
       srmr = h$stats$srmr, chisq = h$stats$chisq.scaled, df = h$stats$df, p = h$stats$pvalue.scaled,
       omega = r3(h$omega), alpha = r3(alpha8(items)),
       agreement_variance = r3(h$ri_variance), agreement_se = r3(h$ri_se), agreement_p = signif(h$ri_p, 3),
       agreement_share_of_item_variance = r3(h$ri_share))
}
out <- list(h1 = summ(h1$sf8, D$sf8), the_10 = summ(h1$sf10, D$sf10), h1_unflagged = summ(h1$sf8_unflagged, D$sf8),
            n_flagged_left_out = sum(M$flagged),
            exploratory_statements = loadings[form == "13" & role == "exploratory", .(item, source, content, keyed_content, agreement)])
write_json(out, sf2b_out("03_h1_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
saveRDS(list(h1 = h1, noRI = noRI, f13 = f13), sf2b_derived("structure"))
print(fits[, .(model, converged, cfi.scaled, rmsea.scaled, srmr)]); print(cut8)
