# 12: Access + Valuing as the scored scale. EXPLORATORY and NOT PREREGISTERED (2026-10-12), after the
# independent review: Norm decided the PSY240 trial's primary outcome is Access + Valuing (Understanding
# and Letting be secondary, no 20-item total). What the two subscales look like on Study 1's data:
# reliability, structure, score distributions, error of measurement, and what they relate to.
# In-sample: the items were chosen on these data, so every number here is optimistic.
#
# Two forms: AV10 = the 10 Access and Valuing items of the 20-item short form (what the trial gives);
# AV14 = the proposed 14 (11_followup's max-reliability sets: 8 Access, 6 Valuing), which Study 2 can
# score from the 20 plus the five re-run Study 1 items.
# Writes 12_access_valuing_study1.json.

source("study1/analysis/config.R")
suppressMessages({ library(psych); library(data.table); library(jsonlite); library(lavaan) })
stopifnot(SF1$mode == "real")
set.seed(20261012)

A <- readRDS(sf1_path("derived", "sf1_real.rds"))
m <- as.data.table(A$meta); keep <- !m$excluded_attention %in% TRUE
k <- A$items_keyed[keep, ]; raw <- A$items_raw[keep, ]; m <- m[keep]
stopifnot(identical(as.character(m$participant), rownames(k)))
rev <- A$definition$reversed; pairs <- A$definition$pairs

# Background answers the registered pipeline does not use (as in 11).
raw_file <- tail(sort(list.files(sf1_path("raw"), "^extract_.*\\.json$", full.names = TRUE)), 1)
rows <- fromJSON(raw_file, simplifyVector = FALSE)$rows
bgr <- Filter(function(x) identical(x$slug, "sf-background") && is.null(x$resubmission_of), rows)
bgr <- bgr[order(vapply(bgr, `[[`, "", "received_at"))]; bgr <- bgr[!duplicated(vapply(bgr, `[[`, "", "participant"))]
num1 <- function(v) if (is.numeric(v)) v else NA_real_
bg <- rbindlist(lapply(bgr, function(x) { b <- x$responses; y <- b$bg_practice_years
  data.table(participant = x$participant, spiritual = num1(b$bg_spiritual), religious = num1(b$bg_religious),
             years = if (is.list(y) && identical(y$option_id, "years")) num1(y$value) else NA_real_) }))
m <- merge(m, bg, by = "participant", all.x = TRUE, sort = FALSE)
stopifnot(identical(as.character(m$participant), rownames(k)))
freq_map <- c(none = 0, lt_weekly = 0.5, `1_2_week` = 1.5, `3_5_week` = 4, daily = 6.5)
m[, freq := ifelse(practice_regular == "no", 0, freq_map[as.character(practice_freq)])]
m[, log_years := ifelse(practice_regular == "no", 0, log1p(pmin(years, 40)))]
m[, dedicated := ifelse(practice_regular == "no", 0,
                        ifelse(practice_freq %in% c("3_5_week", "daily") & !is.na(years) & years >= 1 & practice_training == "yes", 1, NA))]
acq <- rowMeans(raw[, c(pairs[, 1], pairs[, 2])], na.rm = TRUE)   # the agreement index, as in 11

ids <- function(v) paste0("sf6_", v)
FORMS <- list(AV10 = list(Access = c(35, 36, 9, 8, 7), Valuing = c(12, 18, 11, 27, 44)),
              AV14 = list(Access = c(9, 7, 6, 2, 8, 14, 36, 35), Valuing = c(18, 12, 25, 27, 28, 44)))
r3 <- function(x) round(x, 3)

# ---- scores, reliability, distribution, error of measurement ---------------------------------------
score <- function(v) rowMeans(k[, ids(v), drop = FALSE], na.rm = TRUE)
desc <- function(v) {
  X <- k[, ids(v)]; s <- score(v); s <- s[!is.na(s)]; a <- suppressWarnings(alpha(X, warnings = FALSE))$total
  pc <- suppressWarnings(polychoric(X)$rho); p <- ncol(pc)
  sem <- sd(s) * sqrt(1 - a$raw_alpha)
  list(n_items = length(v), n_reversed = sum(ids(v) %in% rev), alpha = r3(a$raw_alpha),
       ordinal_alpha = r3((p / (p - 1)) * (1 - p / sum(pc))), mean_r = r3(a$average_r),
       mean = r3(mean(s)), sd = r3(sd(s)), skew = r3(skew(s)), pct_at_or_above_5.5 = round(100 * mean(s >= 5.5), 1),
       pct_at_6 = round(100 * mean(s >= 6), 1), pct_at_or_below_2 = round(100 * mean(s <= 2), 1),
       sem = r3(sem), mdc95 = r3(1.96 * sqrt(2) * sem))
}

# ---- structure: CFA on unrecoded items, WLSMV, pairwise, with and without the agreement factor -----
cfa_fit <- function(form, one_factor = FALSE, ri = TRUE) {
  it <- ids(unlist(form)); fn <- if (one_factor) "SF" else names(form)
  ms <- if (one_factor) paste0("SF =~ ", paste(it, collapse = " + ")) else
    paste(sapply(names(form), function(f) paste0(f, " =~ ", paste(ids(form[[f]]), collapse = " + "))), collapse = "\n")
  if (ri) ms <- paste(ms, paste0("RI =~ ", paste0("1*", it, collapse = " + ")), "RI ~~ NA*RI",
                      paste0("RI ~~ 0*", fn, collapse = "\n"), sep = "\n")
  fit <- cfa(ms, data = as.data.frame(raw[, it]), ordered = it, estimator = "WLSMV", parameterization = "theta",
             missing = "pairwise", std.lv = TRUE)
  fm <- fitMeasures(fit, c("ntotal", "chisq.scaled", "df", "cfi.scaled", "tli.scaled", "rmsea.scaled",
                           "rmsea.ci.upper.scaled", "srmr", "cfi.robust", "rmsea.robust"))
  s <- standardizedSolution(fit); L <- s[s$op == "=~" & s$lhs %in% fn, ]
  # Keyed loadings (reversed items flipped, each factor oriented so its keyed loadings sum positive) and
  # model-based omega on latent responses, content only.
  L$keyed <- ifelse(L$rhs %in% rev, -L$est.std, L$est.std)
  o <- sapply(fn, function(f) sign(sum(L$keyed[L$lhs == f]))); L$keyed <- L$keyed * o[L$lhs]
  om <- sapply(fn, function(f) { l <- L$est.std[L$lhs == f]; sum(abs(l))^2 / (sum(abs(l))^2 + sum(1 - l^2)) })
  ri_l <- s[s$op == "=~" & s$lhs == "RI", "est.std"]
  list(fit = fit, out = list(fit = as.list(r3(fm)),
    loadings = data.table(factor = L$lhs, item = as.integer(sub("sf6_", "", L$rhs)), keyed_loading = r3(L$keyed)),
    omega = as.list(r3(om)),
    latent_r = if (!one_factor) {
      p <- lavInspect(fit, "std")$psi[fn, fn]; r3(p[1, 2] * prod(o)) } else NULL,
    agreement_pct_item_variance = if (ri) round(100 * mean(ri_l^2), 1) else NULL))
}

out <- list(generated = format(Sys.time(), "%Y-%m-%d"), n = nrow(k),
            note = "Exploratory, not preregistered (2026-10-12). In-sample: the items were chosen on these data.")
for (fn in names(FORMS)) {
  F <- FORMS[[fn]]
  two <- cfa_fit(F); one <- cfa_fit(F, one_factor = TRUE); two_nori <- cfa_fit(F, ri = FALSE)
  dif <- tryCatch(lavTestLRT(one$fit, two$fit), error = function(e) NULL)
  sc <- data.table(Access = score(F$Access), Valuing = score(F$Valuing)); sc[, Combined := (Access + Valuing) / 2]
  out[[fn]] <- list(
    items = lapply(F, function(v) v),
    scales = c(lapply(F, desc), list(Combined = desc(unlist(F)))),
    observed_r = r3(cor(sc$Access, sc$Valuing, use = "pair")),
    two_factor = two$out, one_factor = one$out, two_factor_no_agreement = two_nori$out,
    one_vs_two = if (!is.null(dif)) list(chisq_diff = r3(dif[2, "Chisq diff"]), df = dif[2, "Df diff"], p = signif(dif[2, "Pr(>Chisq)"], 3)) else NULL)

  # ---- what the scores relate to ------------------------------------------------------------------
  rel <- function(x) lapply(sc, function(s) r3(cor(s, x, use = "pair")))
  out[[fn]]$correlates <- list(agreement_index = rel(acq), spiritual = rel(m$spiritual), religious = rel(m$religious),
    practice_freq = rel(m$freq), practice_log_years = rel(m$log_years), age = rel(as.numeric(m$age)))
  # Dedicated practitioners vs never, raw and adjusted for spirituality and religiosity (and, separately,
  # for the agreement index). d = group coefficient / pooled SD.
  kg <- function(s, adj) {
    d <- data.table(s = s, g = m$dedicated, spiritual = m$spiritual, religious = m$religious, acq = acq)[!is.na(g)]
    d <- d[complete.cases(d)]
    sp <- sqrt(((sum(d$g == 1) - 1) * var(d$s[d$g == 1]) + (sum(d$g == 0) - 1) * var(d$s[d$g == 0])) / (nrow(d) - 2))
    b <- coef(lm(as.formula(paste("s ~ g", adj)), d))["g"]
    list(n1 = sum(d$g == 1), n0 = sum(d$g == 0), d = r3(b / sp))
  }
  out[[fn]]$dedicated_vs_never <- lapply(sc, function(s) list(raw = kg(s, ""), adj_spiritual_religious = kg(s, "+ spiritual + religious"),
                                                              adj_agreement = kg(s, "+ acq")))
}

f <- sf1_out("12_access_valuing_study1.json")
write_json(out, f, auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null", dataframe = "rows")
sf1_log("access + valuing: ", f)
