# Model helpers shared by the Study 2 scripts, so every refit (subsets, bootstrap
# resamples, the holdout half) is the identical model.
#
# Registered settings applied in every fit (Section 5, "General settings"):
#   items ordinal, WLSMV, theta parameterization, missing = "pairwise" (set
#   explicitly: Study 1 left lavaan's listwise default and lost 12 cases);
#   the random intercept (RI): loadings fixed at 1 on every UNRECODED item in the
#   model, free variance, uncorrelated with the content factors.
# conditional.x = FALSE everywhere: with ordinal indicators lavaan's default
# (conditional.x = TRUE whenever an observed exogenous covariate is present) gave
# wrong ESEM regression estimates in testing (signs flipped, spurious effects on
# factors with no true effect); FALSE recovered the generating values.

sf2_member <- function(items, factor_of, fn) {
  m <- outer(factor_of[items], fn, `==`); dimnames(m) <- list(items, fn); m
}

# Partially specified target (Section 5, H1): an item's own factor is free, every
# other loading is targeted to 0.
sf2_target <- function(member) {
  list(target = matrix(0, nrow(member), ncol(member), dimnames = dimnames(member)),
       mask = ifelse(member, 0, 1))
}

sf2_ri_syntax <- function(items, fn) {
  c(paste0("RI =~ ", paste0("1*", items, collapse = " + ")),
    paste0("RI ~~ ", paste0("0*", fn, collapse = " + ")), "RI ~~ RI")
}

# ESEM: correlated content factors in one efa block (target rotation), optionally the RI.
sf2_esem_syntax <- function(items, fn, ri = TRUE) {
  s <- paste0(paste0('efa("content")*', fn, collapse = " + "), " =~ ", paste(items, collapse = " + "))
  if (ri) s <- c(s, sf2_ri_syntax(items, fn))
  s
}

# Simple-structure CFA: each item on its own factor only (Model A); factor variances 1.
sf2_cfa_syntax <- function(member, ri = TRUE) {
  items <- rownames(member); fn <- colnames(member)
  s <- vapply(fn, function(f) { it <- items[member[, f]]
    paste0(f, " =~ NA*", it[1], if (length(it) > 1) paste0(" + ", paste(it[-1], collapse = " + ")) else "") }, "")
  s <- c(unname(s), paste0(fn, " ~~ 1*", fn))
  if (ri) s <- c(s, sf2_ri_syntax(items, fn))
  s
}

# Validity (and other observed) variables as single-indicator latents L_<v>, so they
# can correlate with rotated factors; ri = "zero" keeps them uncorrelated with the RI
# (the registered primary analysis), "free" estimates that correlation.
sf2_observed_syntax <- function(vars, fn, ri = c("zero", "free", "none")) {
  ri <- match.arg(ri)
  unlist(lapply(vars, function(v) c(
    paste0("L_", v, " =~ 1*", v), paste0(v, " ~~ 0*", v),
    paste0("L_", v, " ~~ ", paste(fn, collapse = " + ")),
    switch(ri, zero = paste0("L_", v, " ~~ 0*RI"), free = paste0("L_", v, " ~~ RI"), none = NULL))))
}

# The one fitting call. fast = TRUE skips standard errors and the test statistic
# (bootstrap resamples need point estimates only).
#
# ESEM estimation (target given) is run from two starting points and the solution with
# the lower minimum of the fit function is kept: lavaan's default starts, and the
# estimates of the same model with the efa block replaced by its simple-structure CFA
# (each item on its targeted factor(s)), fitted to the same data. Found in the dry run:
# on a 511-case half, the 39-item ESEM with the random intercept converged from the
# default starts to a local optimum (fit function .66 vs .27, CFI .94 vs .99, RI
# variance -0.11) in which one content factor had absorbed the agreement factor; from
# the CFA starts it reached the solution the full sample gives. The rotation's 30
# random starts do not help there: they act after estimation, on the rotation only.
sf2_sem <- function(model, data, ordered, target = NULL, orthogonal = FALSE, fast = FALSE, group = NULL, ...) {
  model <- paste(model, collapse = "\n")
  args <- list(data = data, ordered = ordered, estimator = "WLSMV", parameterization = "theta",
               missing = "pairwise", conditional.x = FALSE, ...)
  if (fast) args <- c(args, list(se = "none", test = "none"))
  if (!is.null(group)) args$group <- group
  if (is.null(target)) return(suppressWarnings(do.call(lavaan::sem, c(list(model = model), args))))
  rot <- list(rotation = "pst", rotation.args = list(target = target$target, target_mask = target$mask,
                                                     orthogonal = orthogonal, rstarts = SF2$rstarts))
  fits <- list()
  cfa <- tryCatch(suppressWarnings(do.call(lavaan::sem, c(list(model = sf2_cfa_start_syntax(model, target, orthogonal)),
                                                         modifyList(args, list(se = "none", test = "none"))))), error = function(e) e)
  if (inherits(cfa, "lavaan"))
    fits$cfa_start <- tryCatch(suppressWarnings(do.call(lavaan::sem, c(list(model = model, start = cfa), args, rot))), error = function(e) e)
  # Default starts, capped at 300 iterations when the CFA start succeeded (converging
  # fits take under 150): on bootstrap resamples they sometimes ran ~10,000 iterations
  # (minutes) into a worse optimum.
  ctl <- if (inherits(fits$cfa_start, "lavaan")) list(control = list(iter.max = 300)) else list()
  fits$default <- tryCatch(suppressWarnings(do.call(lavaan::sem, c(list(model = model), args, rot, ctl))), error = function(e) e)
  ok <- Filter(function(f) inherits(f, "lavaan") && lavInspect(f, "converged"), fits)
  if (!length(ok)) { f <- fits$default; if (inherits(f, "error")) stop(f); return(f) }
  fx <- vapply(ok, function(f) lavInspect(f, "optim")$fx, 0)
  ok[[which.min(fx)]]
}

# The efa line of an ESEM syntax replaced by its simple-structure CFA: each factor
# loads the items its target leaves free (mask 0), variances fixed at 1.
sf2_cfa_start_syntax <- function(model, target, orthogonal = FALSE) {
  lines <- strsplit(model, "\n", fixed = TRUE)[[1]]
  fn <- colnames(target$mask); items <- rownames(target$mask)
  cfa <- unlist(lapply(fn, function(f) { it <- items[target$mask[, f] == 0]
    c(paste0(f, " =~ NA*", it[1], if (length(it) > 1) paste0(" + ", paste(it[-1], collapse = " + ")) else ""), paste0(f, " ~~ 1*", f)) }))
  if (orthogonal && length(fn) > 1) for (a in 1:(length(fn) - 1)) for (b in (a + 1):length(fn)) cfa <- c(cfa, paste0(fn[a], " ~~ 0*", fn[b]))
  paste(c(cfa, lines[!grepl("^efa\\(", lines)]), collapse = "\n")
}

sf2_try <- function(expr) tryCatch(expr, error = function(e) e)
sf2_failed <- function(x) inherits(x, "error") || (inherits(x, "lavaan") && !lavInspect(x, "converged"))
sf2_msg <- function(x) if (inherits(x, "error")) trimws(gsub("[[:space:]]+", " ", conditionMessage(x))) else if (sf2_failed(x)) "did not converge" else ""

# Standardized (std.all) loadings, items x factors.
sf2_loadings <- function(fit, fn, items = NULL) {
  L <- lavInspect(fit, "std")$lambda
  if (is.null(items)) items <- rownames(L)
  L[items, fn, drop = FALSE]
}

# Keyed direction (reverse-keyed rows sign-flipped) and orientation: a factor's sign
# is arbitrary (lavaan reflects each rotated factor so its raw loadings sum positive,
# which turns Letting be, mostly reversed items, into "needing to know"). Each factor
# is oriented so that the KEYED loadings of the items targeted to it sum positive.
sf2_orient <- function(L, member, reversed) {
  K <- L; rv <- rownames(K) %in% reversed; K[rv, ] <- -K[rv, ]
  o <- vapply(colnames(member), function(f) { s <- sum(K[rownames(member)[member[, f]], f]); if (s < 0) -1 else 1 }, 1)
  list(keyed = sweep(K, 2, o, `*`), raw = sweep(L, 2, o, `*`), orient = o)
}

# Model-based ordinal omega of a factor from the fitted model's standardized loadings
# (latent-response scale), content loading only: omega = (sum l)^2 / ((sum l)^2 +
# sum(1 - l^2)), l = each item's keyed loading on the factor. The agreement factor
# and any cross-loadings count as error. PREREG: this is Study 1's formula
# (12_access_valuing.R), with signed keyed loadings instead of absolute values (an
# item loading against its keying lowers omega rather than raising it).
sf2_omega <- function(l) { l <- l[!is.na(l)]; sum(l)^2 / (sum(l)^2 + sum(1 - l^2)) }

# Per-factor rule (Section 5, H1): confirmed when at least four of its five items
# (n - 1 of n) load >= .40 on it in their keyed direction and none of its items loads
# more strongly (in absolute value) on another factor. Scored for individuals when
# confirmed and omega >= .70.
sf2_factor_rule <- function(Kor, member) {
  rbindlist(lapply(colnames(member), function(f) {
    its <- rownames(member)[member[, f]]
    own <- Kor[its, f]
    other <- if (ncol(Kor) > 1) apply(abs(Kor[its, setdiff(colnames(Kor), f), drop = FALSE]), 1, max) else rep(0, length(its))
    n40 <- sum(own >= .40); stronger <- sum(other > abs(own)); om <- sf2_omega(own)
    conf <- n40 >= length(its) - 1 && stronger == 0
    data.table(factor = f, n_items = length(its), n_keyed_ge_40 = n40, n_stronger_elsewhere = stronger,
               items_stronger_elsewhere = paste(its[other > abs(own)], collapse = " "),
               min_keyed_loading = round(min(own), 3), confirmed = conf, omega = round(om, 3),
               scored_for_individuals = conf && om >= .70)
  }))
}

SF2_FIT_INDICES <- c("ntotal", "npar", "chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled",
                     "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr", "cfi.robust", "tli.robust",
                     "rmsea.robust", "rmsea.ci.upper.robust")
sf2_fit_stats <- function(fit) {
  if (sf2_failed(fit)) return(setNames(as.list(rep(NA_real_, length(SF2_FIT_INDICES))), SF2_FIT_INDICES))
  fm <- suppressWarnings(fitMeasures(fit, SF2_FIT_INDICES))
  as.list(round(unclass(fm), 4))
}

sf2_ri <- function(fit) {
  pe <- parameterEstimates(fit)
  r <- pe[pe$lhs == "RI" & pe$op == "~~" & pe$rhs == "RI", ]
  if (!nrow(r)) return(NULL)
  std <- lavInspect(fit, "std")$lambda[, "RI"]
  list(variance = r$est[1], se = r$se[1], z = r$z[1], p = r$pvalue[1],
       mean_std_loading = mean(std[std != 0]), mean_pct_item_variance = 100 * mean(std[std != 0]^2))
}

# H1 decision (Section 5). PREREG: "partial support" = not fully supported, scaled CFI
# >= .90 and RMSEA <= .08 (the text says "CFI .90-.95 with RMSEA <= .08"; a model with
# CFI >= .95 that misses another full criterion is read as partial, not as no support).
# The random-intercept condition is part of full support only; its p is lavaan's
# two-sided Wald p.
sf2_h1_decision <- function(fs, ri) {
  ri_ok <- !is.null(ri) && isTRUE(ri$variance > 0 && ri$p < .05)
  ok <- function(x) isTRUE(!is.na(x) && x)
  full <- ok(fs$cfi.scaled >= .95) && ok(fs$rmsea.scaled <= .06) && ok(fs$rmsea.ci.upper.scaled <= .08) &&
    ok(fs$srmr <= .08) && ri_ok
  partial <- !full && ok(fs$cfi.scaled >= .90) && ok(fs$rmsea.scaled <= .08)
  if (full) "supported" else if (partial) "partially supported" else "not supported"
}

# Observed scores: mean of available items after reverse-scoring ((min + max) - x),
# when at least 80% of the scale's items are answered; "sum" scales are the mean x k
# (a prorated sum). PREREG: the 80% rule (Section 5, "Missing data") is applied to
# every observed score, including the validity scales.
sf2_score <- function(X, items, reversed = character(), range = c(1, 6), method = "mean",
                      min_prop = SF2$min_prop_answered) {
  M <- X[, items, drop = FALSE]; storage.mode(M) <- "double"
  rv <- intersect(items, reversed)
  if (length(rv)) M[, rv] <- (range[1] + range[2]) - M[, rv]
  n_ok <- rowSums(!is.na(M)); need <- ceiling(min_prop * length(items) - 1e-9)
  s <- rowMeans(M, na.rm = TRUE); s[n_ok < need] <- NA
  if (method == "sum") s <- s * length(items)
  s
}

# Sparse categories (Section 5, step 4): a response category with fewer than 5 answers
# in ANY group is merged with its neighbour toward the scale's centre, for that item,
# in all groups; repeated until no category is sparse (or two remain). Categories are
# then renumbered 1..m. The most extreme sparse category is merged first.
sf2_merge_sparse <- function(X, group = NULL, min_n = SF2$min_category_n, centre = 3.5) {
  if (is.null(group)) group <- rep(1L, nrow(X))
  group <- factor(group); log <- list()
  for (j in colnames(X)) {
    x <- X[, j]
    repeat {
      cats <- sort(unique(x[!is.na(x)]))
      if (length(cats) <= 2) break
      tab <- table(factor(x, levels = cats), group)
      sparse <- cats[apply(tab, 1, min) < min_n]
      if (!length(sparse)) break
      c0 <- sparse[which.max(abs(sparse - centre))]
      into <- if (c0 < centre) min(cats[cats > c0]) else max(cats[cats < c0])
      log[[length(log) + 1]] <- data.table(item = j, category = c0, merged_into = into, n_min_group = min(tab[as.character(c0), ]))
      x[!is.na(x) & x == c0] <- into
    }
    X[, j] <- match(x, sort(unique(x[!is.na(x)])))
  }
  list(X = X, log = rbindlist(log))
}

# ---- invariance (Section 5, step 4; Wu & Estabrook 2016, theta parameterization) ----
#
# pattern: items x factors, "" = no loading, "free", or a number (a loading fixed at
# that value in every group: the ESEM-within-CFA referents). The four nested models:
#   configural  thresholds group-specific; latent-response intercepts 0 and residual
#               variances 1 in every group; factor means 0, variances 1
#   thresholds  thresholds equal; intercepts and residual variances free in groups 2..G
#   loadings    + loadings equal; factor variances free in groups 2..G
#   intercepts  + intercepts equal (0); factor means free in groups 2..G
# The RI keeps loadings 1 in every group, its variance free in every group, and is
# uncorrelated with the content factors; PREREG: its mean is freed with the content
# factor means at the intercepts step (the text does not say). free = items released
# from a step's constraint (partial invariance); an item whose thresholds are released
# keeps intercept 0 and residual variance 1 in every group (not identified otherwise).
SF2_INV_STEPS <- c("configural", "thresholds", "loadings", "intercepts")
sf2_inv_syntax <- function(pattern, ncat, G, step, free = list(), ri = TRUE) {
  items <- rownames(pattern); fn <- colnames(pattern); lev <- match(step, SF2_INV_STEPS)
  gl <- function(x) paste0("c(", paste(x, collapse = ", "), ")")
  grp <- function(lab) paste0(lab, "_g", seq_len(G))
  s <- character()
  for (f in fn) {
    terms <- character()
    for (i in items) {
      p <- pattern[i, f]
      if (p == "") next
      if (p == "free") {
        eq <- lev >= 3 && !(i %in% free$loadings)
        terms <- c(terms, paste0(gl(if (eq) rep(paste0("l_", i, "_", f), G) else grp(paste0("l_", i, "_", f))), "*", i))
      } else terms <- c(terms, paste0(p, "*", i))
    }
    s <- c(s, paste0(f, " =~ ", paste(terms, collapse = " + ")))
  }
  if (ri) s <- c(s, paste0("RI =~ ", paste0("1*", items, collapse = " + ")))
  for (i in items) {
    thr_eq <- lev >= 2 && !(i %in% free$thresholds)
    for (k in seq_len(ncat[[i]] - 1))
      s <- c(s, paste0(i, " | ", gl(if (thr_eq) rep(paste0("t_", i, "_", k), G) else grp(paste0("t_", i, "_", k))), "*t", k))
    can_free <- thr_eq && ncat[[i]] >= 3
    int_eq <- lev >= 4 && !(i %in% free$intercepts)
    s <- c(s, paste0(i, " ~ ", gl(if (can_free && !int_eq) c("0", rep("NA", G - 1)) else rep("0", G)), "*1"),
           paste0(i, " ~~ ", gl(if (can_free) c("1", rep("NA", G - 1)) else rep("1", G)), "*", i))
  }
  for (f in fn) s <- c(s, paste0(f, " ~~ ", gl(if (lev >= 3) c("1", rep("NA", G - 1)) else rep("1", G)), "*", f),
                       paste0(f, " ~ ", gl(if (lev >= 4) c("0", rep("NA", G - 1)) else rep("0", G)), "*1"))
  if (length(fn) > 1) for (a in 1:(length(fn) - 1)) for (b in (a + 1):length(fn)) s <- c(s, paste0(fn[a], " ~~ ", fn[b]))
  if (ri) s <- c(s, paste0("RI ~~ ", gl(rep("NA", G)), "*RI"),
                 paste0("RI ~ ", gl(if (lev >= 4) c("0", rep("NA", G - 1)) else rep("0", G)), "*1"),
                 paste0("RI ~~ ", gl(rep("0", G)), "*", fn))
  paste(s, collapse = "\n")
}

# Longitudinal version (baseline vs follow-up, one group): items named <item>_t1 and
# <item>_t2; each occasion has its own factors (<f>_t1, <f>_t2) and its own RI; the
# same item's residuals are correlated across occasions; occasion 1 is the reference.
sf2_long_syntax <- function(member, ncat, step, free = list()) {
  items <- rownames(member); fn <- colnames(member); lev <- match(step, SF2_INV_STEPS)
  s <- character()
  for (tt in 1:2) {
    sfx <- paste0("_t", tt)
    for (f in fn) {
      its <- items[member[, f]]
      lab <- vapply(its, function(i) if (lev >= 3 && !(i %in% free$loadings)) paste0("l_", i) else paste0("l_", i, sfx), "")
      s <- c(s, paste0(f, sfx, " =~ ", paste0(lab, "*", its, sfx, collapse = " + ")))
    }
    s <- c(s, paste0("RI", sfx, " =~ ", paste0("1*", items, sfx, collapse = " + ")), paste0("RI", sfx, " ~~ RI", sfx),
           paste0("RI", sfx, " ~~ ", paste0("0*", fn, rep(c("_t1", "_t2"), each = length(fn)), collapse = " + ")))
    for (i in items) {
      thr_eq <- lev >= 2 && !(i %in% free$thresholds)
      for (k in seq_len(ncat[[i]] - 1))
        s <- c(s, paste0(i, sfx, " | ", if (thr_eq) paste0("t_", i, "_", k) else paste0("t_", i, "_", k, sfx), "*t", k))
      can_free <- thr_eq && ncat[[i]] >= 3 && tt == 2
      int_eq <- lev >= 4 && !(i %in% free$intercepts)
      s <- c(s, paste0(i, sfx, " ~ ", if (can_free && !int_eq) "NA" else "0", "*1"),
             paste0(i, sfx, " ~~ ", if (can_free) "NA" else "1", "*", i, sfx))
    }
    for (f in fn) s <- c(s, paste0(f, sfx, " ~~ ", if (tt == 2 && lev >= 3) "NA" else "1", "*", f, sfx),
                         paste0(f, sfx, " ~ ", if (tt == 2 && lev >= 4) "NA" else "0", "*1"))
    s <- c(s, paste0("RI", sfx, " ~ ", if (tt == 2 && lev >= 4) "NA" else "0", "*1"))
  }
  s <- c(s, "RI_t1 ~~ RI_t2", paste0(items, "_t1 ~~ ", items, "_t2"))
  lv <- paste0(fn, rep(c("_t1", "_t2"), each = length(fn)))
  for (a in 1:(length(lv) - 1)) for (b in (a + 1):length(lv)) s <- c(s, paste0(lv[a], " ~~ ", lv[b]))
  paste(s, collapse = "\n")
}

# Runs the four steps, comparing each with the previous accepted model; a step is
# non-invariant when delta scaled CFI < -.010 AND delta scaled RMSEA > .015 (PREREG:
# "with" read as both, Chen 2007). Then partial invariance item by item: every item
# still constrained at that step is released in turn; the release with the best CFI
# is kept; repeated until the step passes or SF2$partial_max_free items are released.
# build(step, free) returns syntax; fit(syntax, ...) returns a lavaan fit.
# Every model is fitted from lavaan's default starts and from the previous accepted
# model's estimates (for the configural model, from start0 if given: the CFA's
# configural fit, for the ESEM-within-CFA models), keeping the lower fit-function
# minimum. In the dry run the ESEM-within-CFA models, whose loading matrix is almost
# free beside a random intercept, otherwise stopped at local optima or failed to
# converge (a spurious loadings-step misfit of delta CFI -.035 that vanished when any
# one item was released).
sf2_inv_best <- function(fit, syn, starts) {
  # Every fit capped at 3,000 iterations (converging fits here took under 200; the ones
  # that ran to lavaan's 10,000 never converged, and took minutes each).
  starts <- Filter(function(s0) inherits(s0, "lavaan"), starts)
  ctl <- list(iter.max = 3000)
  cands <- lapply(starts, function(s0) sf2_try(fit(syn, start = s0, control = ctl)))
  cands <- c(cands, list(sf2_try(fit(syn, control = ctl))))
  ok <- Filter(function(f) inherits(f, "lavaan") && lavInspect(f, "converged"), cands)
  if (!length(ok)) return(cands[[1]])
  ok[[which.min(vapply(ok, function(f) lavInspect(f, "optim")$fx, 0))]]
}
sf2_inv_sequence <- function(build, fit, items, label, start0 = NULL) {
  res <- list(); free <- list(thresholds = character(), loadings = character(), intercepts = character())
  prev <- NULL; prev_fit <- start0; fits <- list(); stop_here <- FALSE
  row <- function(f, st, released) {
    fs <- sf2_fit_stats(f)
    data.table(analysis = label, step = st, released = paste(released, collapse = " "), converged = !sf2_failed(f),
               error = sf2_msg(f), chisq = fs$chisq.scaled, df = fs$df, cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled,
               rmsea_upper = fs$rmsea.ci.upper.scaled, srmr = fs$srmr,
               d_cfi = if (!is.null(prev)) round(fs$cfi.scaled - prev$cfi, 4) else NA_real_,
               d_rmsea = if (!is.null(prev)) round(fs$rmsea.scaled - prev$rmsea, 4) else NA_real_)
  }
  for (st in SF2_INV_STEPS) {
    if (stop_here) break
    f <- sf2_inv_best(fit, build(st, free), list(prev_fit))
    r <- row(f, st, character())
    if (!r$converged && st == "configural") { res[[length(res) + 1]] <- r; stop_here <- TRUE; next }
    # A constrained step that does not converge is treated like a non-invariant one:
    # partial invariance is examined item by item.
    r[, noninvariant := !converged | (!is.na(d_cfi) & d_cfi < -.010 & d_rmsea > .015)]
    res[[length(res) + 1]] <- r
    if (isTRUE(r$noninvariant)) {
      key <- st; released <- character()
      repeat {
        if (length(released) >= SF2$partial_max_free) break
        cand <- setdiff(items, c(released, unlist(free)))
        tries <- sf2_par_lapply(cand, function(i) { fr <- free; fr[[key]] <- c(fr[[key]], released, i)
          f2 <- sf2_inv_best(fitf, buildf(stp, fr), list(pf)); if (sf2_failed(f2)) NA_real_ else fitMeasures(f2, "cfi.scaled")[[1]] },
          data = list(free = free, key = key, released = released, fitf = fit, buildf = build, stp = st, pf = prev_fit))
        cf <- unlist(tries)
        if (all(is.na(cf))) break
        released <- c(released, cand[which.max(cf)])
        fr <- free; fr[[key]] <- c(fr[[key]], released)
        f <- sf2_inv_best(fit, build(st, fr), list(prev_fit))
        r <- row(f, st, released); r[, noninvariant := !converged | (!is.na(d_cfi) & d_cfi < -.010 & d_rmsea > .015)]
        r[, step := paste0(st, " (partial)")]
        res[[length(res) + 1]] <- r
        if (!isTRUE(r$noninvariant)) break
      }
      free[[key]] <- c(free[[key]], released)
      if (isTRUE(r$noninvariant) || !r$converged) stop_here <- TRUE   # no acceptable partial model: stop the sequence
    }
    fits[[st]] <- f; prev_fit <- f
    fs <- sf2_fit_stats(f); prev <- list(cfi = fs$cfi.scaled, rmsea = fs$rmsea.scaled)
  }
  list(table = rbindlist(res, fill = TRUE), fits = fits, free = free, completed = !stop_here)
}

# Latent means and variances of the content factors in groups 2..G of a fitted
# multi-group intercepts-step model (group 1 = 0 / 1), oriented by the content items.
sf2_latent_means <- function(f, fn, orient = setNames(rep(1, length(fn)), fn)) {
  pe <- parameterEstimates(f)
  m <- pe[pe$op == "~1" & pe$lhs %in% fn & pe$group > 1, c("lhs", "group", "est", "se", "pvalue")]
  m$est <- m$est * orient[m$lhs]
  as.data.table(m)
}

# ---- known groups (Section 5, H6) -------------------------------------------------------
# lm(score ~ group + covariates); d = the group coefficient over the pooled SD of the
# score (within the two groups); its 95% CI is the coefficient's CI over that SD.
sf2_kg <- function(y, g, covars = NULL) {
  d <- data.table(y = y, g = as.integer(g)); if (!is.null(covars)) d <- cbind(d, covars)
  d <- d[complete.cases(d)]
  cv <- setdiff(names(d), c("y", "g"))
  cv <- cv[vapply(cv, function(v) length(unique(d[[v]])) > 1, TRUE)]
  m <- lm(reformulate(c("g", cv), "y"), data = d)
  s <- summary(m)$coefficients["g", ]; ci <- confint(m)["g", ]
  n1 <- sum(d$g == 1); n0 <- sum(d$g == 0)
  sp <- sqrt(((n1 - 1) * var(d$y[d$g == 1]) + (n0 - 1) * var(d$y[d$g == 0])) / (n1 + n0 - 2))
  data.table(n_group = n1, n_reference = n0, mean_group = mean(d$y[d$g == 1]), mean_reference = mean(d$y[d$g == 0]),
             coef = s[["Estimate"]], se = s[["Std. Error"]], p = s[["Pr(>|t|)"]], pooled_sd = sp,
             d = s[["Estimate"]] / sp, d_lower = ci[[1]] / sp, d_upper = ci[[2]] / sp,
             covariates = paste(cv, collapse = " + "))
}

# A function whose environment holds exactly the objects named, so it can be sent to
# parallel workers whatever environment the calling script runs in (Rscript runs a
# script in the global environment, which is never serialized).
sf2_closure <- function(fun, ...) { environment(fun) <- list2env(list(...), parent = globalenv()); fun }

# ---- parallel map --------------------------------------------------------------------------
# FUN runs on PSOCK workers that have config.R and this file sourced; `data` becomes
# FUN's environment, so nothing else from the calling script is serialized.
sf2_par_lapply <- function(X, FUN, data = list(), cores = SF2$cores) {
  environment(FUN) <- list2env(data, parent = globalenv())
  if (cores <= 1 || length(X) < 2) return(lapply(X, FUN))
  cl <- makeCluster(min(cores, length(X)))
  on.exit(stopCluster(cl))
  clusterCall(cl, function(wd, env) { setwd(wd); do.call(Sys.setenv, as.list(env)); NULL }, getwd(),
              c(SF2_MODE = SF2$mode, SF2_DATA_DIR = SF2$data_dir, SF2_RSTARTS = SF2$rstarts, SF2_CORES = "1"))
  clusterEvalQ(cl, { suppressMessages({ source("study2/analysis/config.R"); source("study2/analysis/sf2_models.R") }); NULL })
  parLapplyLB(cl, X, FUN)
}

# Bootstrap index sets, drawn in the calling process with a recorded seed so the
# result does not depend on how workers are scheduled.
sf2_boot_indices <- function(n, B, seed, strata = NULL) {
  set.seed(seed)
  lapply(seq_len(B), function(b) {
    if (is.null(strata)) return(sample.int(n, n, replace = TRUE))
    unlist(lapply(split(seq_len(n), strata), function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]), use.names = FALSE)
  })
}

sf2_r_ci <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 4) return(list(n = sum(ok), r = NA_real_, lower = NA_real_, upper = NA_real_, p = NA_real_))
  t <- cor.test(x[ok], y[ok])
  list(n = sum(ok), r = unname(t$estimate), lower = t$conf.int[1], upper = t$conf.int[2], p = t$p.value)
}

r3 <- function(x) round(x, 3)

# ---- multi-group invariance runs (06, 07) ---------------------------------------------------------
sf2_pattern_cfa <- function(its, factor_of, fn) {
  m <- sf2_member(its, factor_of, fn); P <- ifelse(m, "free", ""); dimnames(P) <- dimnames(m); P
}
# X: item matrix (rows = participants); group: one value per row, NA = left out.
sf2_run_mg <- function(X, its, P, group, label, start0 = NULL) {
  keep <- !is.na(group)
  mg <- sf2_merge_sparse(X[keep, its, drop = FALSE], group[keep])
  ncat <- apply(mg$X, 2, function(x) length(unique(x[!is.na(x)])))
  d <- as.data.frame(mg$X); d$.group <- as.character(group[keep]); G <- length(unique(d$.group))
  d <- d[order(d$.group), ]   # lavaan takes groups in order of appearance: the alphabetically first is group 1, the reference
  build <- sf2_closure(function(step, free) sf2_inv_syntax(P, ncat, G, step, free, ri = TRUE), P = P, ncat = ncat, G = G)
  # auto.fix.first = FALSE: the syntax frees every loading (labels alone do not stop sem()
  # fixing each factor's first loading at 1); factor variances identify the scale.
  fit <- sf2_closure(function(syn, ...) sf2_sem(syn, d, its, group = ".group", auto.fix.first = FALSE, ...), d = d, its = its)
  sf2_log("invariance: ", label, " (", paste(names(table(d$.group)), table(d$.group), sep = " = ", collapse = ", "), ")")
  r <- sf2_inv_sequence(build, fit, its, label, start0)
  r$merges <- if (nrow(mg$log)) cbind(analysis = label, mg$log) else NULL
  r$groups <- as.list(table(d$.group)); r$pattern <- P
  r
}
