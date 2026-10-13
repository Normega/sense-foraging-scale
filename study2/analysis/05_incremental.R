# 05: H5a and H5b, incremental validity (Section 5, analysis 3), in the validity sample
# (everyone who reached the end of the validity block).
#
# H5a (wellbeing). For each outcome (life satisfaction, PHQ-4 total, PSS-4), a
# hierarchical regression on observed scores, complete cases:
#   Step 1  age, gender, the five FFMQ-15 facets, the eight MAIA-2 subscales
#   Step 2  + Access and Valuing and savouring          (the test)
#   Step 3  + Understanding the modes and Letting be     (secondary)
# Supported for an outcome when the Step 2 delta-R^2 is at least .02 AND its F-test is
# significant after Holm correction across the three outcomes; H5a holds if that is
# true for at least one outcome. Delta-R^2 and standardized coefficients are reported.
# Sensitivity: spirituality and religiosity added to Step 1 (same rule).
# PREREG: gender enters as a factor (man = reference, woman, other = non-binary or
# self-described); "prefer not to say" is missing, so those cases leave the complete-
# case sample. All steps of one regression use the same cases (complete on every
# variable of its last step).
# Latent version, alongside: the H1 model with the outcome and every covariate as
# single-indicator latents (uncorrelated with the RI); the outcome regressed on the
# latent Access and Valuing factors and the covariates; delta-R^2 against the same
# model with those two paths fixed at 0, and a Wald test of the two paths.
#
# H5b (practice). A proportional-odds (ordinal logistic) regression of practice
# frequency (five levels, SF2_FREQ_CODE; no regular practice = the lowest):
#   Step 1  age, gender, posting, spirituality, religiosity, the five FFMQ-15 facets,
#           the eight MAIA-2 subscales
#   Step 2  + the Access + Valuing score
# Supported when the Access + Valuing coefficient is positive and the 1-df likelihood-
# ratio test gives p < .05. Reported: the odds ratio per SD (of the score, in the
# analysed cases) and the change in McFadden's R^2 (1 - logLik / logLik of the
# thresholds-only model). Sensitivity: Access and Valuing entered separately; and
# dedicated vs never-practitioners by Firth-penalized logistic regression (logistf)
# with the same covariates.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
D <- A$defs; items <- D$sf_items; F4 <- D$factors4
VS <- A$meta$reached_end_validity
S <- A$scores[VS]; Mv <- A$meta[VS]
FFMQ <- names(D$validity$ffmq15$scores); MAIA <- names(D$validity$maia2$scores)
OUT <- c(lifesat = "lifesat", phq4 = "phq4_total", pss4 = "pss4")
R <- data.table(S[, c("Access", "Valuing", "Understanding", "LettingBe", FFMQ, MAIA, OUT), with = FALSE],
                age = Mv$age, gender = factor(Mv$gender3, levels = c("man", "woman", "other")))

zs <- function(x) (x - mean(x)) / sd(x)
R[, `:=`(spiritual = Mv$spiritual, religious = Mv$religious)]
h5a <- function(extra = character(), label = "primary") rbindlist(lapply(names(OUT), function(o) {
  y <- OUT[[o]]
  d <- R[complete.cases(R[, c(y, "age", "gender", extra, FFMQ, MAIA, "Access", "Valuing", "Understanding", "LettingBe"), with = FALSE])]
  d <- droplevels(d)
  f1 <- reformulate(c("age", "gender", extra, FFMQ, MAIA), y)
  m1 <- lm(f1, d); m2 <- update(m1, . ~ . + Access + Valuing); m3 <- update(m2, . ~ . + Understanding + LettingBe)
  a12 <- anova(m1, m2); a23 <- anova(m2, m3)
  ds <- copy(d); for (v in c(y, "age", extra, FFMQ, MAIA, "Access", "Valuing", "Understanding", "LettingBe")) set(ds, j = v, value = zs(ds[[v]]))
  b2 <- coef(lm(formula(m2), ds)); b3 <- coef(lm(formula(m3), ds))
  data.table(analysis = label, outcome = o, n = nrow(d), r2_step1 = summary(m1)$r.squared, r2_step2 = summary(m2)$r.squared,
             delta_r2 = summary(m2)$r.squared - summary(m1)$r.squared, F = a12$F[2], df1 = a12$Df[2], df2 = a12$Res.Df[2],
             p = a12[["Pr(>F)"]][2], beta_access = b2[["Access"]], beta_valuing = b2[["Valuing"]],
             delta_r2_step3 = summary(m3)$r.squared - summary(m2)$r.squared, p_step3 = a23[["Pr(>F)"]][2],
             beta_understanding = b3[["Understanding"]], beta_lettingbe = b3[["LettingBe"]])
}))
res <- rbind(h5a(), h5a(c("spiritual", "religious"), "sensitivity: + spirituality, religiosity"))
res[, p_holm := p.adjust(p, "holm"), by = analysis]
res[, supported := delta_r2 >= .02 & p_holm < .05]

# ---- latent version -------------------------------------------------------------------------------
member <- sf2_member(items, D$sf_factor, F4); tg <- sf2_target(member)
cov_obs <- c("age", "g_woman", "g_other", FFMQ, MAIA)
Lx <- data.table(age = Mv$age, g_woman = as.numeric(Mv$gender3 == "woman"), g_other = as.numeric(Mv$gender3 == "other"),
                 S[, c(FFMQ, MAIA, OUT), with = FALSE])
num_cols <- c("age", FFMQ, MAIA, OUT)
for (v in num_cols) set(Lx, j = v, value = as.numeric(scale(Lx[[v]])))
dat <- cbind(as.data.frame(A$sf_raw[VS, ]), as.data.frame(Lx))
latent_fit <- sf2_closure(function(y, restricted) {
  syn <- c(sf2_esem_syntax(items, F4, TRUE), sf2_observed_syntax(cov_obs, F4, "zero"),
           paste0("L_", y, " =~ 1*", y), paste0(y, " ~~ 0*", y), paste0("L_", y, " ~~ 0*RI"),
           paste0("L_", y, " ~ ", if (restricted) "0*Access + 0*Valuing" else "b_acc*Access + b_val*Valuing", " + ",
                  paste0("L_", cov_obs, collapse = " + ")))
  sf2_try(sf2_sem(syn, dat[, c(items, cov_obs, y)], items, target = tg))
}, items = items, F4 = F4, cov_obs = cov_obs, dat = dat, tg = tg)
lat <- rbindlist(sf2_par_lapply(names(OUT), function(o) {
  y <- OUT[[o]]; full <- latent_fit(y, FALSE); rest <- latent_fit(y, TRUE)
  if (sf2_failed(full) || sf2_failed(rest)) return(data.table(outcome = o, error = paste(sf2_msg(full), sf2_msg(rest))))
  ori <- sf2_orient(sf2_loadings(full, F4, items), member, rev_items)$orient
  ss <- as.data.table(standardizedSolution(full))[op == "~" & lhs == paste0("L_", y) & rhs %in% c("Access", "Valuing")]
  w <- lavTestWald(full, "b_acc == 0; b_val == 0")
  r2f <- lavInspect(full, "r2")[[paste0("L_", y)]]; r2r <- lavInspect(rest, "r2")[[paste0("L_", y)]]
  data.table(outcome = o, r2_restricted = r2r, r2_full = r2f, delta_r2 = r2f - r2r,
             beta_access = ss[rhs == "Access", est.std] * ori[["Access"]], beta_valuing = ss[rhs == "Valuing", est.std] * ori[["Valuing"]],
             wald_chisq = w$stat, wald_df = w$df, wald_p = w$p.value)
}, data = list(latent_fit = latent_fit, OUT = OUT, items = items, F4 = F4, member = member, rev_items = D$sf_rev)), fill = TRUE)

# ---- H5b: practice frequency --------------------------------------------------------------------------
Q <- data.table(freq = factor(Mv$freq_code, levels = 0:4, ordered = TRUE), age = Mv$age,
                gender = factor(Mv$gender3, levels = c("man", "woman", "other")), posting = Mv$posting_contemplative,
                spiritual = Mv$spiritual, religious = Mv$religious, S[, c(FFMQ, MAIA, "AccessValuing", "Access", "Valuing"), with = FALSE],
                dedicated = ifelse(Mv$practice_group %in% "dedicated", 1L, ifelse(Mv$practice_group %in% "never", 0L, NA_integer_)))
cov5b <- c("age", "gender", "posting", "spiritual", "religious", FFMQ, MAIA)
polr_fit <- function(f, d) tryCatch(MASS::polr(f, data = d, Hess = TRUE, method = "logistic"), error = function(e) e)
mcf <- function(m, m0) 1 - as.numeric(logLik(m)) / as.numeric(logLik(m0))
h5b_step <- function(add, label) {
  d <- droplevels(Q[complete.cases(Q[, c("freq", cov5b, add), with = FALSE])])
  m0 <- polr_fit(freq ~ 1, d); m1 <- polr_fit(reformulate(cov5b, "freq"), d); m2 <- polr_fit(reformulate(c(cov5b, add), "freq"), d)
  if (inherits(m1, "error") || inherits(m2, "error") || inherits(m0, "error"))
    return(list(analysis = label, error = paste(sapply(list(m0, m1, m2), function(m) if (inherits(m, "error")) conditionMessage(m) else ""), collapse = " ")))
  lr <- 2 * (as.numeric(logLik(m2)) - as.numeric(logLik(m1))); df <- length(add)
  ct <- coef(summary(m2))
  terms <- rbindlist(lapply(add, function(v) { b <- ct[v, "Value"]; se <- ct[v, "Std. Error"]; sdv <- sd(d[[v]])
    data.table(term = v, coef = b, se = se, z = b / se, p_wald = 2 * pnorm(-abs(b / se)), sd = sdv,
               or_per_sd = exp(b * sdv), or_lower = exp((b - 1.96 * se) * sdv), or_upper = exp((b + 1.96 * se) * sdv)) }))
  list(analysis = label, n = nrow(d), n_per_level = as.list(table(d$freq)), lr_chisq = lr, lr_df = df,
       lr_p = pchisq(lr, df, lower.tail = FALSE), mcfadden_step1 = mcf(m1, m0), mcfadden_step2 = mcf(m2, m0),
       delta_mcfadden = mcf(m2, m0) - mcf(m1, m0), terms = terms)
}
h5b <- h5b_step("AccessValuing", "primary: Access + Valuing")
h5b$supported <- is.null(h5b$error) && isTRUE(h5b$terms$coef > 0 && h5b$lr_p < .05)
h5b_sep <- h5b_step(c("Access", "Valuing"), "sensitivity: Access and Valuing separately")
# Firth-penalized logistic regression: dedicated vs never, the same covariates + Access + Valuing.
QF <- Q[!is.na(dedicated)]
dF <- droplevels(QF[complete.cases(QF[, c("dedicated", cov5b, "AccessValuing"), with = FALSE])])
firth <- tryCatch({
  m <- logistf::logistf(reformulate(c(cov5b, "AccessValuing"), "dedicated"), data = dF)
  i <- which(names(coef(m)) == "AccessValuing"); sdv <- sd(dF$AccessValuing)
  list(analysis = "sensitivity: dedicated vs never, Firth logistic", n_dedicated = sum(dF$dedicated), n_never = sum(dF$dedicated == 0),
       coef = unname(coef(m)[i]), or_per_sd = exp(unname(coef(m)[i]) * sdv), or_lower = exp(m$ci.lower[[i]] * sdv),
       or_upper = exp(m$ci.upper[[i]] * sdv), p_penalized_lr = unname(m$prob[[i]]), ci_method = "profile penalized likelihood")
}, error = function(e) list(analysis = "sensitivity: dedicated vs never, Firth logistic", error = conditionMessage(e)))

fwrite(res[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2_out("05_incremental_observed.csv"))
fwrite(lat[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)], sf2_out("05_incremental_latent.csv"))
prim <- res[analysis == "primary"]
dec <- list(h5a_supported = any(prim$supported), h5a = prim[, .(outcome, n, delta_r2, p_holm, supported)],
            h5a_sensitivity_supported = any(res[analysis != "primary", supported]),
            h5b_supported = h5b$supported, h5b = h5b, h5b_sensitivity = list(separately = h5b_sep, firth = firth),
            observed = res, latent = lat)
write_json(dec, sf2_out("05_h5_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
print(res[, .(analysis, outcome, n, delta_r2 = round(delta_r2, 4), p_holm = signif(p_holm, 3), supported)]); print(lat)
str(h5b[c("n", "lr_chisq", "lr_p", "delta_mcfadden")]); print(h5b$terms); print(h5b_sep$terms); str(firth)
sf2_log("H5a: ", if (dec$h5a_supported) "supported" else "not supported", "; H5b: ", if (isTRUE(dec$h5b_supported)) "supported" else "not supported")
