# 12: follow-up, EXPLORATORY and NOT PREREGISTERED (Norm's questions, 2026-10-30,
# after the registered analysis): what does a plain EFA look like, and what is the
# strongest scale these 58 statements can make once overall agreement and
# negative wording are controlled? Real mode only.
#
#   1. scree: eigenvalues of the polychoric matrix (keyed) against parallel
#      analysis' 95th percentile (from 04); plain EFAs, k = 1-7 (minres, oblimin,
#      on that matrix): what each factor is made of
#   2. method-controlled models on the UNRECODED statements (WLSMV, pairwise):
#        G       one content factor, plus the random intercept (agreement) and a
#                wording factor W on the 19 reverse-keyed statements, all
#                orthogonal (a bifactor-style model with two method factors)
#        G (RI)  the same without W, to see what W takes
#        k = 2-4 rotated content factors (geomin) with the RI and W
#   3. the strongest scale: statements chosen to maximize CONTENT omega, the share
#      of a keyed sum score's variance that is G (from G's standardized
#      loadings: var = (sum s*lG)^2 + (sum s*lRI)^2 + (sum_rev s*lW)^2 + sum theta,
#      s = +1 forward, -1 reversed; agreement cancels in a balanced set, which is
#      the point of balancing). Greedy, sizes 6-12, one statement per look-alike
#      pair, with and without a minimum of reverse-keyed statements; validated by
#      selecting on one random half and evaluating on the other.

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; X <- A$items_raw; items <- def$items; rv <- def$reversed
R <- readRDS(sf2a_derived("polychoric_keyed"))
out12 <- function(f) sf2a_out(paste0("12_", f))
set.seed(SF2A$seed)

# ---- 1. scree and plain EFAs ---------------------------------------------------------------
ev <- eigen(R, only.values = TRUE)$values
pa <- fread(sf2a_out("04_parallel_analysis.csv"))
png(out12("scree.png"), width = 1400, height = 900, res = 150)
n_show <- 20
plot(seq_len(n_show), ev[seq_len(n_show)], type = "b", pch = 19, xlab = "Factor", ylab = "Eigenvalue",
     main = "Study 2a: scree, 58 statements (polychoric, reverse-keyed recoded)", las = 1)
lines(pa$factor[seq_len(n_show)], pa$random_95[seq_len(n_show)], lty = 2, col = "grey40")
abline(h = 1, lty = 3, col = "grey60")
legend("topright", c("real data (component eigenvalues)", "parallel analysis, 95th pct (factor eigenvalues)", "eigenvalue 1"),
       lty = c(1, 2, 3), pch = c(19, NA, NA), col = c("black", "grey40", "grey60"), bty = "n")
dev.off()
fwrite(data.table(factor = seq_along(ev), eigenvalue = round(ev, 3), pct = round(100 * ev / sum(ev), 1),
                  cum_pct = round(100 * cumsum(ev) / sum(ev), 1)), out12("eigenvalues.csv"))

plain <- rbindlist(lapply(1:7, function(k) {
  f <- suppressWarnings(psych::fa(R, nfactors = k, n.obs = nrow(X), fm = "minres", rotate = "oblimin"))
  L <- unclass(f$loadings)
  # orient each factor so its defining statements' loadings sum positive
  prim <- apply(abs(L), 1, which.max)
  o <- vapply(seq_len(k), function(j) if (sum(L[prim == j, j]) < 0) -1 else 1, 1); L <- sweep(L, 2, o, `*`)
  fwrite(data.table(item = rownames(L), source = def$source[rownames(L)], hypothesized = def$factor[rownames(L)],
                    keying = ifelse(rownames(L) %in% rv, "reverse", "forward"), round(L, 3)),
         out12(sprintf("plain_efa_k%d_loadings.csv", k)))
  load <- L[cbind(seq_len(nrow(L)), prim)]
  rbindlist(lapply(seq_len(k), function(j) {
    on <- rownames(L)[prim == j & load >= .40]
    comp <- table(factor(def$factor[on], levels = def$factors))
    top <- on[order(-L[on, j])][seq_len(min(4, length(on)))]
    data.table(k = k, factor = j, n_ge_40 = length(on), n_reversed = sum(on %in% rv),
               composition = paste(sprintf("%s %d", names(comp), comp)[comp > 0], collapse = ", "),
               top = paste(def$source[top], collapse = " "),
               top_text = paste(substr(def$text[top], 1, 60), collapse = " | "))
  }))
}))
fwrite(plain, out12("plain_efa_factors.csv"))

# ---- 2. method-controlled models -----------------------------------------------------------
g_syntax <- function(items, wording = TRUE) {
  rw <- intersect(items, rv)
  c(paste0("G =~ NA*", items[1], " + ", paste(items[-1], collapse = " + ")), "G ~~ 1*G",
    paste0("RI =~ ", paste0("1*", items, collapse = " + ")), "RI ~~ RI", "G ~~ 0*RI",
    if (wording) c(paste0("W =~ NA*", rw[1], " + ", paste(rw[-1], collapse = " + ")), "W ~~ 1*W", "W ~~ 0*G + 0*RI"))
}
fit_g <- function(Xs, items = def$items, wording = TRUE) sf2a_sem(g_syntax(items, wording), Xs[, items], ordered = items)
# Standardized loadings of G, RI, W per statement, G oriented so keyed loadings sum positive.
g_parts <- function(fit, items = def$items) {
  L <- lavInspect(fit, "std")$lambda[items, , drop = FALSE]
  s <- ifelse(items %in% rv, -1, 1)
  lg <- L[, "G"]; if (sum(s * lg) < 0) lg <- -lg
  data.table(item = items, s = s, G = lg, RI = L[, "RI"], W = if ("W" %in% colnames(L)) L[, "W"] else 0,
             theta = 1 - lg^2 - L[, "RI"]^2 - (if ("W" %in% colnames(L)) L[, "W"]^2 else 0))
}
# Content omega and the agreement and wording shares of a keyed sum score.
set_omega <- function(P, set) {
  p <- P[item %in% set]
  g <- sum(p$s * p$G)^2; a <- sum(p$s * p$RI)^2; w <- sum(p$s * p$W)^2; e <- sum(p$theta)
  tot <- g + a + w + e
  c(omega_content = g / tot, share_agreement = a / tot, share_wording = w / tot, omega_total = (g + a + w) / tot)
}

library(parallel)
cl <- makeCluster(min(12, detectCores() - 2))
clusterCall(cl, setwd, getwd())
clusterEvalQ(cl, { Sys.setenv(SF2A_MODE = "real"); source("study2a/analysis/config.R"); source("study2a/analysis/sf2a_models.R") })
clusterExport(cl, c("g_syntax", "fit_g", "rv", "def", "X"), envir = environment())

halves <- replicate(5, sample(nrow(X), nrow(X) %/% 2), simplify = FALSE)
jobs <- c(list(list(what = "G", rows = NULL, wording = TRUE), list(what = "G_noW", rows = NULL, wording = FALSE)),
          lapply(2:4, function(k) list(what = paste0("esem_k", k), k = k)),
          unlist(lapply(seq_along(halves), function(h) list(list(what = paste0("G_half", h, "a"), rows = halves[[h]], wording = TRUE),
                                                           list(what = paste0("G_half", h, "b"), rows = setdiff(seq_len(nrow(X)), halves[[h]]), wording = TRUE))),
                 recursive = FALSE))
sf2a_log("fitting ", length(jobs), " models on ", length(cl), " workers")
fits <- parLapply(cl, jobs, function(j) {
  t0 <- Sys.time()
  f <- tryCatch(if (!is.null(j$k)) {
    fn <- paste0("F", seq_len(j$k))
    sf2a_sem(sf2a_esem_syntax(def$items, fn, ri = TRUE, wording = rv), X, ordered = def$items)
  } else fit_g(if (is.null(j$rows)) X else X[j$rows, ], wording = j$wording), error = function(e) e)
  list(what = j$what, fit = f, mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
})
stopCluster(cl)
names(fits) <- vapply(fits, `[[`, "", "what")
for (f in fits) sf2a_log(f$what, ": ", if (sf2a_failed(f$fit)) paste("FAILED", sf2a_msg(f$fit)) else "ok", " (", round(f$mins, 1), " min)")

fitrow <- function(nm) { f <- fits[[nm]]$fit; data.table(model = nm, converged = !sf2a_failed(f), improper = sf2a_improper(f),
                                                          as.data.table(sf2a_fit_stats(f))) }
fwrite(rbindlist(lapply(c("G", "G_noW", paste0("esem_k", 2:4)), fitrow), fill = TRUE), out12("method_models_fit.csv"))

PG <- g_parts(fits$G$fit)
PG[, `:=`(source = def$source[item], hypothesized = def$factor[item], keying = ifelse(item %in% rv, "reverse", "forward"),
          names_mode = item %in% def$names_mode, text = def$text[item])]
PG0 <- g_parts(fits$G_noW$fit)
PG[, G_without_W := round(PG0$G[match(item, PG0$item)], 3)]
fwrite(PG[order(-G), lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)], out12("G_loadings.csv"))
fwrite(rbind(sf2a_ri(fits$G$fit)[, factor := "RI, G model"], sf2a_ri(fits$G_noW$fit)[, factor := "RI, G without W"]),
       out12("G_agreement.csv"))

esem_tab <- rbindlist(lapply(2:4, function(k) {
  f <- fits[[paste0("esem_k", k)]]$fit
  if (sf2a_failed(f)) return(data.table(k = k, note = sf2a_msg(f)))
  fn <- paste0("F", seq_len(k))
  L <- lavInspect(f, "std")$lambda[items, fn, drop = FALSE]; Kd <- L; Kd[rv, ] <- -Kd[rv, ]
  prim <- apply(abs(Kd), 1, which.max)
  o <- vapply(seq_len(k), function(j) if (sum(Kd[prim == j, j]) < 0) -1 else 1, 1); Kd <- sweep(Kd, 2, o, `*`)
  fwrite(data.table(item = items, source = def$source[items], hypothesized = def$factor[items],
                    keying = ifelse(items %in% rv, "reverse", "forward"), round(Kd, 3),
                    W = round(lavInspect(f, "std")$lambda[items, "W"], 3)), out12(sprintf("esem_k%d_RI_W_loadings.csv", k)))
  load <- Kd[cbind(seq_len(nrow(Kd)), prim)]
  rbindlist(lapply(seq_len(k), function(j) {
    on <- items[prim == j & load >= .40]; comp <- table(factor(def$factor[on], levels = def$factors))
    top <- on[order(-Kd[on, j])][seq_len(min(4, length(on)))]
    data.table(k = k, factor = j, n_ge_40 = length(on), n_reversed = sum(on %in% rv),
               composition = paste(sprintf("%s %d", names(comp), comp)[comp > 0], collapse = ", "),
               top = paste(def$source[top], collapse = " "))
  }))
}), fill = TRUE)
fwrite(esem_tab, out12("esem_RI_W_factors.csv"))

# ---- 3. the strongest scale ----------------------------------------------------------------
partner <- function(i) { p <- def$pairs[def$pairs[, "a"] == i | def$pairs[, "b"] == i, , drop = FALSE]; setdiff(c(p), i) }
greedy_omega <- function(P, n, min_rev = 0, pool = P$item) {
  sel <- character()
  while (length(sel) < n) {
    left <- setdiff(pool, c(sel, unlist(lapply(sel, partner))))
    need <- min_rev - sum(sel %in% rv)
    if (need > 0 && n - length(sel) <= need) left <- intersect(left, rv)
    if (!length(left)) break
    sc <- vapply(left, function(i) set_omega(P, c(sel, i))[["omega_content"]], 0)
    sel <- c(sel, names(which.max(sc)))
  }
  sel
}
alpha_keyed <- function(rows, set) psych::alpha(A$items_keyed[rows, set], warnings = FALSE)$total$raw_alpha
acq <- rowMeans(X, na.rm = TRUE)       # raw mean over all 58: a rough agreement index (39 forward, 19 reversed)
g_score <- tryCatch({ s <- lavPredict(fits$G$fit)[, c("G", "RI")]; s }, error = function(e) NULL)
grp <- A$meta$practice_group
dval <- function(y) { a <- y[grp %in% "current"]; b <- y[grp %in% "never"]
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2) }
cands <- list()
for (n in c(6, 8, 10, 12)) for (mr in unique(c(0, n %/% 3, n %/% 2))) {
  set <- greedy_omega(PG, n, mr)
  # validation: choose on half a, evaluate content omega on half b's fit and alpha on half b
  cv <- rbindlist(lapply(seq_along(halves), function(h) {
    Pa <- g_parts(fits[[paste0("G_half", h, "a")]]$fit); Pb <- g_parts(fits[[paste0("G_half", h, "b")]]$fit)
    s <- greedy_omega(Pa, n, mr); rows_b <- setdiff(seq_len(nrow(X)), halves[[h]])
    data.table(omega_content_b = set_omega(Pb, s)[["omega_content"]], alpha_b = alpha_keyed(rows_b, s),
               overlap = length(intersect(s, set)) / n)
  }))
  om <- set_omega(PG, set); sc <- rowMeans(A$items_keyed[, set], na.rm = TRUE)
  cands[[length(cands) + 1]] <- data.table(
    n = n, min_reversed = mr, n_reversed = sum(set %in% rv), n_names_mode = sum(set %in% def$names_mode),
    omega_content = round(om[["omega_content"]], 3), share_agreement = round(om[["share_agreement"]], 3),
    share_wording = round(om[["share_wording"]], 3), alpha = round(alpha_keyed(seq_len(nrow(X)), set), 3),
    cv_omega_content = round(mean(cv$omega_content_b), 3), cv_alpha = round(mean(cv$alpha_b), 3),
    cv_overlap = round(mean(cv$overlap), 2),
    r_rough_agreement = round(cor(sc, acq, use = "pair"), 3),
    r_RI_score = if (!is.null(g_score)) round(cor(sc, g_score[, "RI"], use = "pair"), 3) else NA_real_,
    d_practice = round(dval(sc), 3), r_ius = round(cor(sc, A$meta$ius_total, use = "pair"), 3),
    statements = paste(def$source[set], collapse = " "), ids = paste(set, collapse = " "))
}
# Reference sets: the registered proposal's broad factor, and the trial's statements in the scale.
ref <- list(`registered proposal F3 (6)` = c("sf2a_40", "sf2a_24", "sf2a_23", "sf2a_04", "sf2a_26", "sf2a_05"),
            `registered proposal F1 (5)` = c("sf2a_02", "sf2a_37", "sf2a_09", "sf2a_15", "sf2a_30"),
            `trial statements, reworded or Study 1 wording (10)` = c("sf2a_01", "sf2a_02", "sf2a_03", "sf2a_04", "sf2a_05",
                                                                     "sf2a_23", "sf2a_24", "sf2a_25", "sf2a_26", "sf2a_27"))
for (nm in names(ref)) { set <- ref[[nm]]; om <- set_omega(PG, set); sc <- rowMeans(A$items_keyed[, set], na.rm = TRUE)
  cands[[length(cands) + 1]] <- data.table(n = length(set), min_reversed = NA, n_reversed = sum(set %in% rv),
    n_names_mode = sum(set %in% def$names_mode), omega_content = round(om[["omega_content"]], 3),
    share_agreement = round(om[["share_agreement"]], 3), share_wording = round(om[["share_wording"]], 3),
    alpha = round(alpha_keyed(seq_len(nrow(X)), set), 3), r_rough_agreement = round(cor(sc, acq, use = "pair"), 3),
    r_RI_score = if (!is.null(g_score)) round(cor(sc, g_score[, "RI"], use = "pair"), 3) else NA_real_,
    d_practice = round(dval(sc), 3), r_ius = round(cor(sc, A$meta$ius_total, use = "pair"), 3),
    statements = paste(def$source[set], collapse = " "), ids = paste(set, collapse = " "), reference = nm) }
best <- rbindlist(cands, fill = TRUE)
fwrite(best, out12("strongest_scale.csv"))
saveRDS(list(G = fits$G$fit, G_noW = fits$G_noW$fit, esem = fits[paste0("esem_k", 2:4)]), sf2a_derived("followup12"))

print(plain[, .(k, factor, n_ge_40, n_reversed, composition, top)])
print(fread(out12("method_models_fit.csv")))
print(esem_tab)
print(head(PG[order(-G), .(item, source, keying, G = round(G, 2), RI = round(RI, 2), W = round(W, 2), G_without_W)], 20))
print(best[, !c("statements", "ids")])
