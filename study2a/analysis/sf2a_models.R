# Model helpers shared by every Study 2a script, so a refit on a subset is the
# identical model. Adapted from Study 2's sf2_models.R (the same lavaan settings
# and the two-start ESEM fit, which the Study 2 dry run showed is needed).
#
# Settings in every fit (Section 5, "Statistical models" and "Missing data"):
#   statements ordinal, WLSMV, theta parameterization, missing = "pairwise" (set
#   explicitly: Study 1 left lavaan's listwise default and lost 12 cases);
#   the random intercept (RI): loadings fixed at 1 on every UNRECODED statement,
#   free variance, uncorrelated with the content factors.

# Target membership for k factors (Section 5, step 4): for five factors, the
# hypothesized membership; for any other number, each statement's primary factor
# in the k-factor geomin EFA from 04.
sf2a_membership <- function(def, k) {
  if (k == length(def$factors)) {
    fn <- sf2a_fname(def$factors)
    m <- outer(def$factor[def$items], def$factors, `==`)
  } else {
    E <- fread(sf2a_out(sprintf("04_efa_k%02d_loadings.csv", k)))
    L <- as.matrix(E[, -1]); rownames(L) <- E$item; L <- L[def$items, , drop = FALSE]
    fn <- paste0("F", seq_len(k)); m <- t(apply(abs(L), 1, function(r) seq_len(k) == which.max(r)))
  }
  dimnames(m) <- list(def$items, fn); m
}
sf2a_fname <- function(f) gsub("[^A-Za-z]", "", f)    # "Letting be" -> "Lettingbe"

sf2a_target <- function(member) {
  list(target = matrix(0, nrow(member), ncol(member), dimnames = dimnames(member)),
       mask = ifelse(member, 0, 1))
}

# ESEM syntax: correlated content factors in one efa block, plus optionally the
# RI and a wording factor W (free loadings on the reverse-keyed statements,
# variance 1, uncorrelated with the content factors and the RI; Section 5, step 5).
sf2a_esem_syntax <- function(items, fn, ri = TRUE, wording = NULL) {
  s <- paste0(paste0('efa("content")*', fn, collapse = " + "), " =~ ", paste(items, collapse = " + "))
  if (ri) s <- c(s, paste0("RI =~ ", paste0("1*", items, collapse = " + ")),
                 paste0("RI ~~ ", paste0("0*", fn, collapse = " + ")), "RI ~~ RI")
  if (length(wording)) {
    s <- c(s, paste0("W =~ NA*", wording[1], " + ", paste(wording[-1], collapse = " + ")), "W ~~ 1*W",
           paste0("W ~~ ", paste0("0*", c(fn, if (ri) "RI"), collapse = " + ")))
  }
  s
}

# The efa line replaced by its simple-structure CFA (each factor loads the
# statements its target leaves free), for the second set of starting values.
sf2a_cfa_start_syntax <- function(model, target) {
  lines <- strsplit(model, "\n", fixed = TRUE)[[1]]
  fn <- colnames(target$mask); items <- rownames(target$mask)
  cfa <- unlist(lapply(fn, function(f) { it <- items[target$mask[, f] == 0]
    c(paste0(f, " =~ NA*", it[1], if (length(it) > 1) paste0(" + ", paste(it[-1], collapse = " + ")) else ""),
      paste0(f, " ~~ 1*", f)) }))
  paste(c(cfa, lines[!grepl("^efa\\(", lines)]), collapse = "\n")
}

# The one fitting call. An ESEM (target given) is run from lavaan's default
# starts and from the estimates of its simple-structure CFA, and the converged
# solution with the lower fit-function minimum is kept (Study 2's dry run: from
# default starts a 39-item ESEM with the RI reached a local optimum in which a
# content factor absorbed the agreement factor).
sf2a_sem <- function(model, data, ordered, target = NULL, fast = FALSE, ...) {
  model <- paste(model, collapse = "\n")
  args <- list(data = as.data.frame(data), ordered = ordered, estimator = "WLSMV", parameterization = "theta",
               missing = "pairwise", ...)
  if (fast) args <- c(args, list(se = "none", test = "none"))
  if (is.null(target)) return(suppressWarnings(do.call(lavaan::sem, c(list(model = model), args))))
  rot <- list(rotation = "pst", rotation.args = list(target = target$target, target_mask = target$mask,
                                                     orthogonal = FALSE, rstarts = SF2A$rstarts))
  fits <- list()
  cfa <- tryCatch(suppressWarnings(do.call(lavaan::sem, c(list(model = sf2a_cfa_start_syntax(model, target)),
                                                         modifyList(args, list(se = "none", test = "none"))))),
                  error = function(e) e)
  if (inherits(cfa, "lavaan"))
    fits$cfa_start <- tryCatch(suppressWarnings(do.call(lavaan::sem, c(list(model = model, start = cfa), args, rot))),
                               error = function(e) e)
  fits$default <- tryCatch(suppressWarnings(do.call(lavaan::sem, c(list(model = model), args, rot))), error = function(e) e)
  ok <- Filter(function(f) inherits(f, "lavaan") && lavInspect(f, "converged"), fits)
  if (!length(ok)) { f <- fits$default; if (inherits(f, "error")) stop(f); return(f) }
  ok[[which.min(vapply(ok, function(f) lavInspect(f, "optim")$fx, 0))]]
}

sf2a_fit_esem <- function(X, member, ri = TRUE, wording = NULL, fast = FALSE) {
  items <- rownames(member)
  sf2a_sem(sf2a_esem_syntax(items, colnames(member), ri, wording), X[, items, drop = FALSE], ordered = items,
           target = sf2a_target(member), fast = fast)
}

sf2a_failed <- function(x) inherits(x, "error") || (inherits(x, "lavaan") && !lavInspect(x, "converged"))
sf2a_msg <- function(x) if (inherits(x, "error")) trimws(gsub("[[:space:]]+", " ", conditionMessage(x))) else if (sf2a_failed(x)) "did not converge" else ""

# An improper solution (Section 5, step 5): a negative variance among the latent
# variables, or a standardized loading or factor correlation beyond 1.
sf2a_improper <- function(fit) {
  if (sf2a_failed(fit)) return(TRUE)
  pe <- parameterEstimates(fit, standardized = TRUE)
  lv <- lavNames(fit, "lv")
  own <- pe$op == "~~" & pe$lhs == pe$rhs
  any(own & pe$lhs %in% lv & pe$est < 0) ||
    any(pe$op %in% c("=~", "~~") & pe$lhs %in% lv & !own & abs(pe$std.all) > 1, na.rm = TRUE)
}

# Standardized loadings, statements x factors.
sf2a_loadings <- function(fit, fn, items = NULL) {
  L <- lavInspect(fit, "std")$lambda
  if (is.null(items)) items <- rownames(L)
  L[items, fn, drop = FALSE]
}

# Keyed direction (reverse-keyed rows sign-flipped) and orientation: a factor's
# sign is arbitrary (Study 1's Letting be came back as "needing to know"). Each
# factor is oriented so the KEYED loadings of the statements targeted to it sum
# positive.
sf2a_orient <- function(L, member, reversed) {
  K <- L; rv <- rownames(K) %in% reversed; K[rv, ] <- -K[rv, ]
  o <- vapply(colnames(member), function(f) { s <- sum(K[rownames(member)[member[, f]], f]); if (s < 0) -1 else 1 }, 1)
  list(keyed = sweep(K, 2, o, `*`), raw = sweep(L, 2, o, `*`), orient = o)
}

# Model-based ordinal omega of a factor (Section 5, step 7) from keyed standardized
# loadings: (sum l)^2 / ((sum l)^2 + sum(1 - l^2)). Study 2's formula: agreement and
# cross-loadings count as error.
sf2a_omega <- function(l) { l <- l[!is.na(l)]; sum(l)^2 / (sum(l)^2 + sum(1 - l^2)) }

SF2A_FIT <- c("ntotal", "npar", "chisq.scaled", "df", "cfi.scaled", "tli.scaled", "rmsea.scaled",
              "rmsea.ci.upper.scaled", "srmr")
sf2a_fit_stats <- function(fit) {
  if (sf2a_failed(fit)) return(setNames(as.list(rep(NA_real_, length(SF2A_FIT))), SF2A_FIT))
  as.list(round(unclass(suppressWarnings(fitMeasures(fit, SF2A_FIT))), 4))
}

sf2a_ri <- function(fit, lv = "RI") {
  pe <- parameterEstimates(fit)
  r <- pe[pe$lhs == lv & pe$op == "~~" & pe$rhs == lv, ]
  if (!nrow(r)) return(NULL)
  std <- lavInspect(fit, "std")$lambda[, lv]
  data.table(variance = round(r$est[1], 4), se = round(r$se[1], 4), p = signif(r$pvalue[1], 3),
             mean_std_loading = round(mean(std[std != 0]), 3), mean_pct_item_variance = round(100 * mean(std[std != 0]^2), 1))
}

# Tucker congruence between two loading matrices on shared statements, each
# reference factor matched to its best-congruent partner (sign ignored).
sf2a_congruence <- function(L_ref, L_new) {
  items <- intersect(rownames(L_ref), rownames(L_new))
  C <- psych::factor.congruence(L_ref[items, , drop = FALSE], L_new[items, , drop = FALSE])
  data.table(factor = colnames(L_ref), matched = colnames(L_new)[apply(abs(C), 1, which.max)],
             congruence = round(apply(abs(C), 1, max), 3))
}

# Which hypothesized factor a fitted factor is: the one most of its statements
# (primary keyed loading >= .40 on it) were written for. PREREG: the text names
# factors by their hypothesized content without saying how a factor of a k != 5
# solution is identified; this is the reading used.
sf2a_label <- function(Kor, def) {
  prim <- colnames(Kor)[apply(abs(Kor), 1, which.max)]
  vapply(colnames(Kor), function(f) {
    on <- rownames(Kor)[prim == f & abs(Kor[, f]) >= SF2A$load_min]
    if (!length(on)) return(NA_character_)
    tb <- sort(table(def$factor[on]), decreasing = TRUE); names(tb)[1]
  }, "")
}

# Observed score: mean of available statements after recoding (7 - x).
sf2a_score <- function(X, items, reversed = character(), top = 7) {
  M <- X[, items, drop = FALSE]; storage.mode(M) <- "double"
  rv <- intersect(items, reversed); if (length(rv)) M[, rv] <- top - M[, rv]
  rowMeans(M, na.rm = TRUE)
}
