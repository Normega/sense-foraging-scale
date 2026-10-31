# 14: exploratory, not preregistered (Norm's question, 2026-10-31): for each hypothesized
# factor, is there a 3-statement subscale that holds together beyond the general factor
# and is not carried by agreement? Real mode only, after 12.
#
# From 12's G model (G + the random intercept RI + the wording factor W, orthogonal):
#   - every 3-statement set within one hypothesized factor (one statement per look-alike
#     pair) gets a specific factor S from the keyed residual correlations G leaves among
#     its statements (exactly identified with three: lambda_1 = sqrt(r12 r13 / r23), ...;
#     all three residuals must be positive)
#   - the keyed sum score's variance splits into G, S, agreement (RI), wording (W) and
#     error, as in 12's set_omega with S added:
#       omega_sub      (G + S) / total   reliability of the subscale as content
#       omega_specific S / total         what it measures that the general factor does not
#       share_agreement RI / total
#   - the rule: the largest omega_specific with share_agreement <= .05; validated by
#     choosing on one random half (G refitted there) and evaluating on the other
#   - the chosen five sets confirmed in two fitted models: G + RI + W + five orthogonal
#     specific factors on all 58 statements; and the 15 statements alone as five
#     correlated factors plus RI (and W when a set has reversed statements)

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; X <- A$items_raw; items <- def$items; rv <- def$reversed
F12 <- readRDS(sf2a_derived("followup12"))
out14 <- function(f) sf2a_out(paste0("14_", f))
set.seed(SF2A$seed + 14)
agree_max <- .05

g_syntax <- function(items, wording = TRUE) {
  rw <- intersect(items, rv)
  c(paste0("G =~ NA*", items[1], " + ", paste(items[-1], collapse = " + ")), "G ~~ 1*G",
    paste0("RI =~ ", paste0("1*", items, collapse = " + ")), "RI ~~ RI", "G ~~ 0*RI",
    if (wording) c(paste0("W =~ NA*", rw[1], " + ", paste(rw[-1], collapse = " + ")), "W ~~ 1*W", "W ~~ 0*G + 0*RI"))
}
g_parts <- function(fit) {
  L <- lavInspect(fit, "std")$lambda[items, , drop = FALSE]
  s <- ifelse(items %in% rv, -1, 1)
  lg <- L[, "G"]; if (sum(s * lg) < 0) lg <- -lg
  data.table(item = items, s = s, G = lg, RI = L[, "RI"], W = L[, "W"], theta = 1 - lg^2 - L[, "RI"]^2 - L[, "W"]^2)
}
keyed_resid <- function(fit) { r <- lavResiduals(fit, type = "cor")$cov[items, items]; s <- ifelse(items %in% rv, -1, 1); r * outer(s, s) }
partner <- function(i) { p <- def$pairs[def$pairs[, "a"] == i | def$pairs[, "b"] == i, , drop = FALSE]; setdiff(c(p), i) }

triple_stats <- function(P, Rk, t) {
  r12 <- Rk[t[1], t[2]]; r13 <- Rk[t[1], t[3]]; r23 <- Rk[t[2], t[3]]
  if (min(r12, r13, r23) <= 0) return(NULL)
  ls <- sqrt(c(r12 * r13 / r23, r12 * r23 / r13, r13 * r23 / r12))
  p <- P[match(t, item)]
  if (any(p$theta - ls^2 <= 0)) return(NULL)                     # a Heywood case
  g <- sum(p$s * p$G)^2; sp <- sum(ls)^2; a <- sum(p$s * p$RI)^2; w <- sum(p$s * p$W)^2; e <- sum(p$theta - ls^2)
  tot <- g + sp + a + w + e
  list(omega_sub = (g + sp) / tot, omega_specific = sp / tot, share_g = g / tot, share_agreement = a / tot,
       share_wording = w / tot, min_resid = min(r12, r13, r23), s_load = ls)
}
all_triples <- function(P, Rk) rbindlist(lapply(def$factors, function(f) {
  pool <- items[def$factor[items] == f]
  cmb <- combn(pool, 3, simplify = FALSE)
  cmb <- Filter(function(t) !any(t %in% unlist(lapply(t, partner))), cmb)
  rbindlist(lapply(cmb, function(t) { z <- triple_stats(P, Rk, t); if (is.null(z)) return(NULL)
    data.table(factor = f, ids = paste(t, collapse = " "), omega_sub = z$omega_sub, omega_specific = z$omega_specific,
               share_g = z$share_g, share_agreement = z$share_agreement, share_wording = z$share_wording,
               min_resid = z$min_resid, s_loadings = paste(sprintf("%.2f", z$s_load), collapse = " ")) }))
}))
pick <- function(T) T[share_agreement <= agree_max][order(-omega_specific), .SD[1], by = factor]

# ---- 1. every triple, full sample ------------------------------------------------------------
PG <- g_parts(F12$G); RK <- keyed_resid(F12$G)
TT <- all_triples(PG, RK)
ann <- function(T) T[, `:=`(statements = vapply(strsplit(ids, " "), function(t) paste(def$source[t], collapse = " "), ""),
                            n_reversed = vapply(strsplit(ids, " "), function(t) sum(t %in% rv), 0L),
                            n_names_mode = vapply(strsplit(ids, " "), function(t) sum(t %in% def$names_mode), 0L))]
ann(TT)
num3 <- function(T) T[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)]
fwrite(num3(TT[order(factor, -omega_specific)]), out14("all_triples.csv"))
chosen <- pick(TT)
best_sub <- TT[share_agreement <= agree_max][order(-omega_sub), .SD[1], by = factor]

# ---- 2. half-split validation ------------------------------------------------------------------
library(parallel)
halves <- replicate(5, sample(nrow(X), nrow(X) %/% 2), simplify = FALSE)
jobs <- unlist(lapply(seq_along(halves), function(h) list(list(h = h, side = "a", rows = halves[[h]]),
                                                         list(h = h, side = "b", rows = setdiff(seq_len(nrow(X)), halves[[h]])))),
               recursive = FALSE)
cl <- makeCluster(min(length(jobs), detectCores() - 2))
clusterCall(cl, setwd, getwd())
clusterEvalQ(cl, { Sys.setenv(SF2A_MODE = "real"); source("study2a/analysis/config.R"); source("study2a/analysis/sf2a_models.R") })
clusterExport(cl, c("g_syntax", "rv", "items", "X"), envir = environment())
hf <- parLapply(cl, jobs, function(j) list(h = j$h, side = j$side,
                                           fit = tryCatch(sf2a_sem(g_syntax(items), X[j$rows, ], ordered = items), error = function(e) e)))
stopCluster(cl)
cv <- rbindlist(lapply(seq_along(halves), function(h) {
  fa <- Filter(function(z) z$h == h && z$side == "a", hf)[[1]]$fit; fb <- Filter(function(z) z$h == h && z$side == "b", hf)[[1]]$fit
  if (sf2a_failed(fa) || sf2a_failed(fb)) return(NULL)
  ca <- pick(all_triples(g_parts(fa), keyed_resid(fa))); Pb <- g_parts(fb); Rb <- keyed_resid(fb)
  rbindlist(lapply(seq_len(nrow(ca)), function(i) { t <- strsplit(ca$ids[i], " ")[[1]]; z <- triple_stats(Pb, Rb, t)
    data.table(half = h, factor = ca$factor[i], ids = ca$ids[i], same_as_full = ca$ids[i] == chosen[factor == ca$factor[i], ids],
               overlap_full = length(intersect(t, strsplit(chosen[factor == ca$factor[i], ids], " ")[[1]])),
               omega_specific_b = if (is.null(z)) 0 else z$omega_specific, omega_sub_b = if (is.null(z)) NA else z$omega_sub,
               share_agreement_b = if (is.null(z)) NA else z$share_agreement) }))
}))
fwrite(num3(cv), out14("half_split.csv"))
cvs <- cv[, .(halves = .N, mean_overlap_with_full = mean(overlap_full), omega_specific_validated = mean(omega_specific_b),
              omega_sub_validated = mean(omega_sub_b, na.rm = TRUE)), by = factor]

# ---- 3. the chosen sets in fitted models -------------------------------------------------------
sets <- setNames(lapply(strsplit(chosen$ids, " "), identity), chosen$factor)
sname <- function(f) paste0("S_", gsub("[^A-Za-z]", "", f))
bif <- c(g_syntax(items), unlist(lapply(names(sets), function(f) c(
  paste0(sname(f), " =~ NA*", paste(sets[[f]], collapse = " + ")), paste0(sname(f), " ~~ 1*", sname(f)),
  paste0(sname(f), " ~~ 0*G + 0*RI + 0*W"))) ),
  unlist(combn(vapply(names(sets), sname, ""), 2, function(p) paste0(p[1], " ~~ 0*", p[2]), simplify = FALSE)))
fb <- sf2a_sem(bif, X, ordered = items)
it15 <- unlist(sets); rw15 <- intersect(it15, rv)
cor5 <- c(unlist(lapply(names(sets), function(f) c(paste0(sname(f), " =~ NA*", paste(sets[[f]], collapse = " + ")),
                                                    paste0(sname(f), " ~~ 1*", sname(f))))),
          paste0("RI =~ ", paste0("1*", it15, collapse = " + ")), "RI ~~ RI",
          paste0("RI ~~ ", paste0("0*", vapply(names(sets), sname, ""), collapse = " + ")),
          if (length(rw15) >= 2) c(paste0("W =~ NA*", paste(rw15, collapse = " + ")), "W ~~ 1*W",
                                   paste0("W ~~ 0*RI + ", paste0("0*", vapply(names(sets), sname, ""), collapse = " + "))))
fc <- sf2a_sem(cor5, X[, it15], ordered = it15)

keyedL <- function(fit, lv) { L <- lavInspect(fit, "std")$lambda; v <- L[, lv]; v <- v[v != 0]
  k <- v * ifelse(names(v) %in% rv, -1, 1); if (sum(k) < 0) k <- -k; k }
conf <- rbindlist(lapply(names(sets), function(f) {
  kb <- if (!sf2a_failed(fb)) keyedL(fb, sname(f)) else NA; kc <- if (!sf2a_failed(fc)) keyedL(fc, sname(f)) else NA
  data.table(factor = f, item = sets[[f]], source = def$source[sets[[f]]], keying = ifelse(sets[[f]] %in% rv, "reverse", "forward"),
             names_mode = sets[[f]] %in% def$names_mode, G_in_G_model = round(PG$G[match(sets[[f]], PG$item)], 3),
             specific_in_bifactor = round(kb[sets[[f]]], 3), loading_in_5factor = round(kc[sets[[f]]], 3),
             text = def$text[sets[[f]]])
}))
fwrite(conf, out14("chosen_loadings.csv"))
fits <- rbind(data.table(model = "bifactor: G + RI + W + 5 specific (58)", converged = !sf2a_failed(fb), improper = sf2a_improper(fb),
                         as.data.table(sf2a_fit_stats(fb))),
              data.table(model = "5 correlated factors + RI (15)", converged = !sf2a_failed(fc), improper = sf2a_improper(fc),
                         as.data.table(sf2a_fit_stats(fc))), fill = TRUE)
fwrite(fits, out14("fit.csv"))
if (!sf2a_failed(fc)) {
  ps <- lavInspect(fc, "std")$psi; nm <- vapply(names(sets), sname, ""); ps <- ps[nm, nm]
  o <- vapply(names(sets), function(f) { v <- lavInspect(fc, "std")$lambda[, sname(f)]; v <- v[v != 0]
    if (sum(v * ifelse(names(v) %in% rv, -1, 1)) < 0) -1 else 1 }, 1)
  ps <- ps * outer(o, o); dimnames(ps) <- list(names(sets), names(sets))
  fwrite(data.table(factor = names(sets), round(ps, 3)), out14("factor_correlations.csv"))
}

# plain scores: alpha, correlation with the RI factor score, and with 12's best-8
sc <- sapply(sets, function(t) rowMeans(A$items_keyed[, t], na.rm = TRUE))
ri_score <- lavPredict(F12$G)[, "RI"]
b12 <- fread(sf2a_out("12_strongest_scale.csv"))[n == 8 & min_reversed == 0, strsplit(ids, " ")[[1]]]
b8 <- rowMeans(A$items_keyed[, b12], na.rm = TRUE)
grp <- A$meta$practice_group
dval <- function(y) { a <- y[grp %in% "current"]; b <- y[grp %in% "never"]
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2) }
summ <- merge(chosen[, .(factor, statements, n_reversed, n_names_mode, omega_sub, omega_specific, share_g, share_agreement)], cvs, by = "factor")
summ[, `:=`(alpha = vapply(factor, function(f) psych::alpha(A$items_keyed[, sets[[f]]], warnings = FALSE)$total$raw_alpha, 0),
            r_RI_score = vapply(factor, function(f) cor(sc[, f], ri_score, use = "pair"), 0),
            r_best8 = vapply(factor, function(f) cor(sc[, f], b8, use = "pair"), 0),
            r_ius = vapply(factor, function(f) cor(sc[, f], A$meta$ius_total, use = "pair"), 0),
            d_practice = vapply(factor, function(f) dval(sc[, f]), 0))]
fwrite(num3(summ), out14("chosen_summary.csv"))
fwrite(num3(best_sub[, .(factor, statements, n_reversed, omega_sub, omega_specific, share_agreement)]), out14("best_omega_sub.csv"))

print(num3(summ)); print(fits); if (file.exists(out14("factor_correlations.csv"))) print(fread(out14("factor_correlations.csv")))
print(conf[, !"text"]); print(num3(best_sub[, .(factor, statements, omega_sub, omega_specific, share_agreement)]))
print(TT[, .(n_triples = .N, n_low_agreement = sum(share_agreement <= agree_max), max_specific = round(max(omega_specific), 3)), by = factor])
