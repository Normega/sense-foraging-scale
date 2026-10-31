# 21: exploratory, not preregistered (Norm's decisions, 2026-10-31): the forms carried into
# Study 2b, described in the 2a data. Real mode only, after 12 and 18.
#   - the 8: 12's best 8 (three reverse-keyed statements, all naming a mode; Norm kept the
#     third rather than swapping in U4, 17-20)
#   - the 10: the 8 plus A9 and sfc_24 (18's best mode-free pair)
#   - the three plain-worded reverse-keyed statements carried as exploratory: U4, N6, sfc_31
# For each form: content omega and agreement (12's G model), alpha overall and by vignette
# group, SEM and smallest detectable change, practitioner d, correlations; a one-factor CFA
# with the random intercept (the 2b model); and a preview of 2b's invariance test across
# people who did and did not tell the modes apart (configural, loadings, thresholds).
# In-sample throughout: the statements were chosen on these data.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; K <- A$items_keyed; X <- A$items_raw; M <- A$meta; items <- def$items; rv <- def$reversed
F12 <- readRDS(sf2a_derived("followup12"))
out21 <- function(f) sf2a_out(paste0("21_", f))
num3 <- function(T) T[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]
id <- function(src) items[match(src, def$source[items])]

L <- lavInspect(F12$G, "std")$lambda[items, c("G", "RI", "W")]; s_all <- ifelse(items %in% rv, -1, 1)
if (sum(s_all * L[, "G"]) < 0) L[, "G"] <- -L[, "G"]
set_omega <- function(set) { s <- ifelse(set %in% rv, -1, 1); l <- L[set, , drop = FALSE]
  g <- sum(s * l[, "G"])^2; a <- sum(s * l[, "RI"])^2; w <- sum(s * l[, "W"])^2; e <- sum(1 - rowSums(l^2))
  c(content = g / (g + a + w + e), agreement = a / (g + a + w + e)) }
alpha_fast <- function(Z) { S <- cov(Z, use = "pairwise.complete.obs"); k <- ncol(S); k / (k - 1) * (1 - sum(diag(S)) / sum(S)) }
grp <- M$practice_group
dval <- function(y) { a <- y[grp %in% "current"]; b <- y[grp %in% "never"]
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2) }
trial10 <- rowMeans(cbind(K, A$end)[, unname(def$av10)], na.rm = TRUE)
ta <- M$tells_apart %in% TRUE

f8 <- fread(sf2a_out("12_strongest_scale.csv"))[n == 8 & min_reversed == 0, strsplit(ids, " ")[[1]]]
f10 <- c(f8, id(c("A9", "sfc_24")))
expl <- id(c("U4", "N6", "sfc_31"))
forms <- list(`8` = f8, `10` = f10)
sc8 <- rowMeans(K[, f8], na.rm = TRUE)

describe <- function(set, label) {
  om <- set_omega(set); sc <- rowMeans(K[, set], na.rm = TRUE); a <- alpha_fast(K[, set]); sem <- sd(sc, na.rm = TRUE) * sqrt(1 - a)
  data.table(form = label, n = length(set), n_reversed = sum(set %in% rv), n_names_mode = sum(set %in% def$names_mode),
             omega_content = om[["content"]], share_agreement = om[["agreement"]], alpha = a,
             alpha_told_apart = alpha_fast(K[ta, set]), alpha_not = alpha_fast(K[!ta, set]),
             mean = mean(sc, na.rm = TRUE), sd = sd(sc, na.rm = TRUE), sem = sem, mdc95 = 1.96 * sqrt(2) * sem,
             pct_at_or_above_5.5 = 100 * mean(sc >= 5.5, na.rm = TRUE), pct_at_or_below_2 = 100 * mean(sc <= 2, na.rm = TRUE),
             d_practice = dval(sc), r_8 = cor(sc, sc8, use = "pair"), r_trial10 = cor(sc, trial10, use = "pair"),
             r_ius = cor(sc, M$ius_total, use = "pair"), r_self = cor(sc, M$vig_self, use = "pair"),
             r_separation = cor(sc, M$separation, use = "pair"))
}
desc <- rbindlist(lapply(names(forms), function(f) describe(forms[[f]], f)))
fwrite(num3(desc), out21("forms.csv"))

stm <- data.table(item = c(f10, expl), source = def$source[c(f10, expl)], text = def$text[c(f10, expl)],
                  role = c(rep("8", 8), "10", "10", rep("exploratory", 3)),
                  hypothesized = def$factor[c(f10, expl)], reversed = c(f10, expl) %in% rv,
                  names_mode = c(f10, expl) %in% def$names_mode, G_keyed = s_all[match(c(f10, expl), items)] * L[c(f10, expl), "G"],
                  RI = L[c(f10, expl), "RI"])
d02 <- fread(sf2a_out("02_item_descriptives.csv"))
stm[, agree_pct := d02$agree_pct[match(item, d02$item)]]
stm[, mean := d02$mean[match(item, d02$item)]]

# one-factor CFA with the random intercept: the 2b model
cfa_syntax <- function(set) c(paste0("S =~ NA*", paste(set, collapse = " + ")), "S ~~ 1*S",
                              paste0("RI =~ ", paste0("1*", set, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI")
fs <- c("chisq.scaled", "df", "cfi.scaled", "tli.scaled", "rmsea.scaled", "rmsea.ci.upper.scaled", "srmr")
cfa <- rbindlist(lapply(names(forms), function(f) {
  set <- forms[[f]]; fit <- sf2a_sem(cfa_syntax(set), X[, set], ordered = set)
  ri <- sf2a_ri(fit)
  l <- lavInspect(fit, "std")$lambda[set, "S"]; l <- l * ifelse(set %in% rv, -1, 1); if (sum(l) < 0) l <- -l
  stm[match(set, item), (paste0("cfa_", f)) := round(l, 3)]
  data.table(form = f, t(round(fitMeasures(fit, fs), 3)), agreement_pct_item_variance = ri$mean_pct_item_variance,
             agreement_p = ri$p, min_loading = round(min(l), 3), max_loading = round(max(l), 3))
}))
fwrite(cfa, out21("cfa.csv"))
# the exploratory three added to the 10, to see how they load beside it
set13 <- c(f10, expl); fit13 <- sf2a_sem(cfa_syntax(set13), X[, set13], ordered = set13)
l13 <- lavInspect(fit13, "std")$lambda[set13, "S"] * ifelse(set13 %in% rv, -1, 1); if (sum(l13) < 0) l13 <- -l13
stm[, cfa_with_exploratory := round(l13[item], 3)]
fwrite(num3(stm), out21("statements.csv"))

# invariance preview across the vignette groups (2b's test), on the 8 and the 10
Xg <- as.data.frame(X); Xg$vig_group <- ifelse(ta, "told apart", "not")
inv <- rbindlist(lapply(names(forms), function(f) {
  set <- forms[[f]]; D <- Xg[, c(set, "vig_group")]
  run <- function(eq) sf2a_sem(cfa_syntax(set), D, ordered = set, group = "vig_group", group.equal = eq)
  fits <- list(configural = run(character()), loadings = run("loadings"), thresholds = run(c("loadings", "thresholds")))
  rbindlist(lapply(names(fits), function(nm) { x <- fits[[nm]]
    data.table(form = f, step = nm, converged = !sf2a_failed(x), t(round(fitMeasures(x, c("chisq.scaled", "df", "cfi.scaled", "rmsea.scaled", "srmr")), 3)),
               p_vs_previous = if (nm == "configural") NA_real_ else
                 tryCatch(signif(lavTestLRT(fits[[which(names(fits) == nm) - 1]], x)[2, "Pr(>Chisq)"], 3), error = function(e) NA_real_)) }))
}))
fwrite(inv, out21("invariance_preview.csv"))
print(num3(desc)); print(cfa); print(num3(stm)[, !"text"]); print(inv)
