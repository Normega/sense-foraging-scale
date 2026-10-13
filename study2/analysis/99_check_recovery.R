# 99 (simulated mode only): does the pipeline recover the model the data came from,
# and did every registered exclusion and flag fire on the cases built to trigger it?
#
# Recovery is judged against the ESTIMAND (00_simulate.R): with cross-loadings in the
# generating model, the target-rotated ESEM converges to a rotation of the generating
# factors, not to the generating factors themselves; the generating values are shown
# beside it. Tolerances are sampling-based: each check states its rule.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
stopifnot(SF2$mode == "simulated")
tr <- readRDS(sf2_derived("simulation_truth")); A <- readRDS(sf2_derived("sf2"))
H1 <- readRDS(sf2_derived("h1")); D <- A$defs; F4 <- D$factors4; items <- D$sf_items
res <- list(); add <- function(area, check, value, rule, pass) res[[length(res) + 1]] <<- data.table(area = area, check = check,
  value = as.character(value), rule = rule, pass = isTRUE(pass))

# ---- recovery -------------------------------------------------------------------------------------------
L <- fread(sf2_out("03_h1_loadings.csv")); Lk <- as.matrix(L[, paste0(F4, "_keyed"), with = FALSE]); rownames(Lk) <- L$item
Tt <- tr$keyed_std[items, F4]; G <- tr$keyed_std_generating[items, F4]
cong <- diag(psych::factor.congruence(Lk, Tt)); dmax <- max(abs(Lk - Tt)); rmse <- sqrt(mean((Lk - Tt)^2))
add("loadings", "Tucker congruence with the estimand, per factor", paste(round(cong, 3), collapse = " "), ">= .98 each", all(cong >= .98))
# A loading's SE here is about .04, so the largest of 80 absolute errors is expected near .12.
add("loadings", "RMSE vs estimand (vs generating)", paste0(round(rmse, 3), " (", round(sqrt(mean((Lk - G)^2)), 3), ")"), "< .05", rmse < .05)
add("loadings", "largest |estimate - estimand| (80 loadings)", paste0(round(dmax, 3), " (", names(which.max(apply(abs(Lk - Tt), 1, max))), ")"),
    "< .15", dmax < .15)
Phi_e <- as.matrix(fread(sf2_out("03_h1_factor_correlations.csv"))[, -1]); dimnames(Phi_e) <- list(F4, F4)
dphi <- max(abs(Phi_e - tr$factor_cor[F4, F4]))
add("factor correlations", "largest |estimate - estimand|", round(dphi, 3), "< .08", dphi < .08)
add("factor correlations", "estimate | estimand | generating (A-V, A-L, U-V)",
    paste(sapply(list(c("Access", "Valuing"), c("Access", "LettingBe"), c("Understanding", "Valuing")), function(p)
      paste(round(Phi_e[p[1], p[2]], 2), round(tr$factor_cor[p[1], p[2]], 2), round(tr$factor_cor_generating[p[1], p[2]], 2), sep = "|")), collapse = "  "), "reported", TRUE)
ri <- H1$decision$ri
add("random intercept", "variance (theta metric): estimate vs true .22", paste0(round(ri$variance, 3), " (SE ", round(ri$se, 3), ")"),
    "within 2.5 SE", abs(ri$variance - tr$ri_variance) < 2.5 * ri$se)
add("random intercept", "share of item variance, %", round(ri$mean_pct_item_variance, 1), "~10", TRUE)
om <- setNames(H1$rule$omega, H1$rule$factor); dom <- max(abs(om[F4] - tr$omega[F4]))
add("omega", "model-based omega per factor vs estimand", paste(round(om[F4], 3), collapse = " "), "|diff| < .05 each", dom < .05)

LV <- fread(sf2_out("04_latent_validity.csv"))
m <- merge(LV, tr$latent_observed[factor %in% F4], by = c("score", "factor"), suffixes = c("", "_true"))
cov_lv <- mean(m$lower <= m$r_true & m$upper >= m$r_true); mad_lv <- mean(abs(m$r - m$r_true))
add("validity (latent)", paste0("95% CI coverage of the estimand (", nrow(m), " cells)"), round(cov_lv, 3), ">= .90", cov_lv >= .90)
add("validity (latent)", "mean |estimate - estimand|; vs generating", paste(round(mad_lv, 3), "/", round(mean(abs(m$r - m$r_generating)), 3)), "< .04 (estimand)", mad_lv < .04)
OB <- fread(sf2_out("04_observed_validity.csv"))
mo <- merge(OB, tr$observed_observed, by = c("score", "factor"), suffixes = c("", "_true"))
cov_ob <- mean(mo$lower <= mo$r_true & mo$upper >= mo$r_true)
add("validity (observed)", paste0("95% CI coverage of the population r (", nrow(mo), " cells)"), round(cov_ob, 3), ">= .90", cov_ob >= .90)
VD <- fromJSON(sf2_out("04_validity_decisions.json"))
for (k in seq_len(nrow(VD$h3d))) { f <- VD$h3d$factor[k]
  add("validity (H3d)", paste0(f, ": multiple R, estimate [bootstrap CI] | estimand | generating"),
      sprintf("%.3f [%.3f, %.3f] | %.3f | %.3f", VD$h3d$R[k], VD$h3d$lower[k], VD$h3d$upper[k], tr$h3d[[f]], tr$h3d_generating[[f]]),
      "CI covers the estimand", VD$h3d$lower[k] <= tr$h3d[[f]] && VD$h3d$upper[k] >= tr$h3d[[f]]) }
add("validity (H3b)", "Valuing: r(SBI) - r(Open-Mindedness), estimate [bootstrap CI] | estimand",
    sprintf("%.3f [%.3f, %.3f] | %.3f", VD$h3b$difference, VD$h3b$lower, VD$h3b$upper, tr$h3b), "CI covers the estimand",
    VD$h3b$lower <= tr$h3b && VD$h3b$upper >= tr$h3b)

ICC <- fread(sf2_out("08_stability.csv"))[pairs == "all pairs"]
for (s in names(tr$icc)) { x <- ICC[score == s]
  add("stability", paste0("ICC ", s, ": estimate [CI] vs population"), sprintf("%.3f [%.3f, %.3f] vs %.3f", x$icc, x$lower, x$upper, tr$icc[[s]]),
      "CI covers", x$lower <= tr$icc[[s]] && x$upper >= tr$icc[[s]]) }
x <- ICC[score == "agreement_sf"]
add("stability", "ICC agreement index (short form)", sprintf("%.3f [%.3f, %.3f] vs %.3f", x$icc, x$lower, x$upper, tr$icc_agreement),
    "CI covers", x$lower <= tr$icc_agreement && x$upper >= tr$icc_agreement)
KG <- fread(sf2_out("07_known_groups.csv"))[analysis == "primary"]
for (s in names(tr$h6_d)) { x <- KG[score == s]
  add("known groups", paste0("d ", s, ": estimate [CI] vs population"), sprintf("%.3f [%.3f, %.3f] vs %.3f", x$d, x$d_lower, x$d_upper, tr$h6_d[[s]]),
      "CI covers", x$d_lower <= tr$h6_d[[s]] && x$d_upper >= tr$h6_d[[s]]) }

FW <- fread(sf2_out("13_following_week.csv"))
for (v in names(tr$following_week)) { x <- FW[item == v]
  add("following week", paste0(v, ": standardized coefficient of baseline Access + Valuing, estimate [CI] vs population"),
      sprintf("%.3f [%.3f, %.3f] vs %.3f (n %d)", x$beta, x$beta_lower, x$beta_upper, tr$following_week[[v]], x$n), "CI covers",
      x$beta_lower <= tr$following_week[[v]] && x$beta_upper >= tr$following_week[[v]]) }

# Interval checks taken together: with 95% intervals about 1 in 20 should miss.
iv <- rbindlist(res)[area %in% c("stability", "known groups", "following week", "validity (H3d)", "validity (H3b)")]
add("intervals", "H3b, H3d, stability, known-groups and following-week CIs covering the population value", paste(sum(iv$pass), "of", nrow(iv)),
    ">= 80% (about 95% expected)", mean(iv$pass) >= .80)

# Step 7: the selection the rules make on the population loadings (true PNA rates) vs the half's.
S7 <- readRDS(sf2_derived("step7"))
K39 <- tr$keyed_std_39; fo <- S7$factor_of; rev39 <- S7$rev39
exp_sel <- rbindlist(lapply(F4, function(f) {
  cur <- items[D$sf_factor == f]; cands <- setdiff(rownames(K39)[fo[rownames(K39)] == f], items); out <- list()
  prim <- K39[cands, f]; oth <- apply(abs(K39[cands, setdiff(F4, f), drop = FALSE]), 1, max)
  ok <- prim >= .50 & oth < .30 & tr$pna_rate[cands] <= .05
  for (r in 1:2) { weak <- cur[which.min(K39[cur, f])]; hit <- FALSE
    for (cc in cands[order(-prim)]) { if (!ok[cc] || cc %in% unlist(out)) next
      ns <- c(setdiff(cur, weak), cc)
      if (sf2_omega(K39[ns, f]) > sf2_omega(K39[cur, f]) && any(ns %in% rev39) && any(!ns %in% rev39)) { out[[r]] <- cc; cur <- ns; hit <- TRUE; break } }
    if (!hit) break }
  data.table(factor = f, expected = paste(unlist(out), collapse = " "))
}))
got <- S7$replacements[, .(selected = paste(candidate, collapse = " ")), by = factor]
cmp7 <- merge(exp_sel, got, by = "factor", all.x = TRUE); cmp7[is.na(selected), selected := ""]
for (k in seq_len(nrow(cmp7))) add("step 7 (informative)", paste0(cmp7$factor[k], ": selected (half) vs expected (population)"),
  paste0("[", cmp7$selected[k], "] vs [", cmp7$expected[k], "]"), "same set (informative: sampling can change it)",
  setequal(strsplit(cmp7$selected[k], " ")[[1]], strsplit(cmp7$expected[k], " ")[[1]]))

# ---- every exclusion and flag fired on the planted cases ------------------------------------------------
P <- tr$planted; MA <- A$meta_all
x <- merge(P, MA, by = "participant", all.x = TRUE)
chk <- function(label, rows, cond) { n <- nrow(rows); k <- sum(cond, na.rm = TRUE)
  add("exclusions and flags", label, paste0(k, " of ", n), "all", n > 0 && k == n) }
chk("staff test enrollments excluded", x[role == "staff_test"], x[role == "staff_test", excl_test & !analysable])
chk("no short form: excluded", x[grepl("^no_sf", role)], x[grepl("^no_sf", role), excl_no_sf & !analysable])
chk("failed both checks: excluded", x[attention == "fail_both"], x[attention == "fail_both", excl_attention & !analysable])
chk("failed one of two checks: flagged, kept", x[attention == "fail_one"], x[attention == "fail_one", flag_one_check & analysable & checks_reached == 2])
chk("failed the only check reached: flagged, kept", x[attention == "fail_only_reached"],
    x[attention == "fail_only_reached", flag_one_check & analysable & checks_reached == 1 & !excl_attention])
chk("never reached a check: not flagged, kept", x[attention == "never_reached"], x[attention == "never_reached", checks_reached == 0 & !flag_one_check & analysable])
chk("speeders flagged, kept", x[flag == "speeder"], x[flag == "speeder", flag_speeder & analysable])
chk("longstrings flagged, kept", x[flag == "longstring"], x[flag == "longstring", flag_longstring & analysable])
chk("> 4 'Prefer not to answer' flagged, kept", x[flag == "pna_many"], x[flag == "pna_many", flag_pna & analysable])
chk("partial sessions kept for H1, out of the validity sample", x[grepl("^partial", role)],
    x[grepl("^partial", role), analysable & !reached_end_validity])
chk("completers not excluded: in the validity sample", x[role == "completer" & attention != "fail_both"],
    x[role == "completer" & attention != "fail_both", analysable & reached_end_validity])
add("exclusions and flags", "unplanted cases excluded (should be 0)",
    sum(x[role == "completer" & attention != "fail_both", !analysable]), "0", sum(x[role == "completer" & attention != "fail_both", !analysable]) == 0)
add("exclusions and flags", "unplanted cases flagged (natural speeders/longstrings/checks)",
    sum(x[flag == "" & !attention %in% c("fail_one", "fail_only_reached") & role != "staff_test", flagged], na.rm = TRUE), "reported", TRUE)
# Repeat submissions: the flagged resubmission dropped; for the re-entry, the FIRST submission used.
raw <- fromJSON(sf2_raw_file(), simplifyVector = FALSE)$rows
sub_rows <- function(p, slug) Filter(function(r) identical(r$participant, p) && identical(r$slug, slug) && identical(r$source, "questionnaire"), raw)
pr <- P[resubmit == TRUE, participant]; rr <- sub_rows(pr, "sf-short-20")
first_sf <- vapply(items, function(i) { v <- rr[[1]]$responses[[i]]; if (is.numeric(v)) v else NA_real_ }, 0)
add("exclusions and flags", "flagged resubmission dropped (2 rows, 1 flagged; first used)",
    paste(length(rr), "rows,", sum(vapply(rr, function(r) !is.null(r$resubmission_of), TRUE)), "flagged"), "2 rows, 1 flagged, values equal",
    length(rr) == 2 && identical(unname(first_sf), unname(A$sf_raw[pr, items])))
pe <- P[reentry == TRUE, participant]; re <- sub_rows(pe, "sf2-sfs-candidates")
re <- re[order(vapply(re, `[[`, "", "received_at"))]
first_c <- vapply(D$cand_items, function(i) { v <- re[[1]]$responses[[i]]; if (is.numeric(v)) v else NA_real_ }, 0)
add("exclusions and flags", "later re-entry with different answers: first submission used",
    paste(length(re), "unflagged rows"), "first submission's answers in the data", length(re) == 2 && identical(unname(first_c), unname(A$cand_raw[pe, D$cand_items])))
INV <- readRDS(sf2_derived("invariance"))
rel <- vapply(names(INV)[grepl("posting", names(INV))], function(k) paste(unlist(INV[[k]]$free), collapse = " "), "")
add("invariance", paste0("planted intercept non-invariance (", tr$dif$item, ", posting): released by"),
    paste(names(rel)[grepl(tr$dif$item, rel)], collapse = "; "), "at least one posting model", any(grepl(tr$dif$item, rel)))
add("invariance", "items released where no non-invariance was planted (gender, country, over time)",
    paste(unlist(lapply(INV[!grepl("posting", names(INV))], function(r) unlist(r$free))), collapse = " "), "reported", TRUE)
cr <- A$country_rule
add("country", "two groups (six listed countries vs elsewhere), both >= 200 -> analysed",
    paste(names(cr$groups), unlist(cr$groups), sep = " = ", collapse = ", "), "analysed",
    isTRUE(cr$analysed) && setequal(names(cr$groups), c("six_countries", "other")))
pw2 <- A$retest$past_week_t2
add("follow-up", "past-week items read from the follow-up", paste(colSums(!is.na(pw2)), collapse = " "), "about 300 each",
    all(colSums(!is.na(pw2)) > 250))
R2 <- A$retest$enrollments; tt <- tr$retest
add("follow-up", "retest staff test excluded", sum(R2$excl_test), "1", sum(R2$excl_test) == 1)
add("follow-up", "unlinked retest code excluded", sum(R2$excl_unlinked), "1", sum(R2$excl_unlinked) == 1)
add("follow-up", "pairs", nrow(A$retest$pairs), paste(sum(tt$linked), "planted"), nrow(A$retest$pairs) == sum(tt$linked))
out_true <- tt[linked == TRUE & (days < 14 | days >= 22), participant]
add("follow-up", "outside 14-21 days flagged and kept", paste(sum(!A$retest$pairs$in_window), "flagged;", length(out_true), "planted"),
    "same set", setequal(A$retest$pairs[in_window == FALSE, participant], out_true))

R <- rbindlist(res)
fwrite(R, sf2_out("99_recovery.csv"))
print(R[, .(area, check, value, rule, pass)], nrows = 200)
Rc <- R[!grepl("informative", area)]
sf2_log("recovery checks: ", sum(Rc$pass), " of ", nrow(Rc), " pass (step 7 comparisons informative, not counted)",
        if (any(!Rc$pass)) paste0("; failing: ", paste(Rc[pass == FALSE, check], collapse = "; ")) else "")
