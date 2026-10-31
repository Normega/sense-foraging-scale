# 20: exploratory, not preregistered (Norm, 2026-10-31): how do the final 8 and its mode-free
# variants behave as a one-factor CFA (what 2b would fit)? In-sample: the statements were
# chosen on these data, so fit is optimistic. Real mode only, after 19.
#   one factor; + the random intercept (RI); + RI and a residual covariance between the two
#   statements that name a mode; standardized loadings and the largest modification indices

source("study2a/analysis/config.R")
source("study2a/analysis/sf2a_models.R")
stopifnot(SF2A$mode == "real")
A <- readRDS(sf2a_derived("sf2a"))
def <- A$definition; X <- A$items_raw; rv <- def$reversed
out20 <- function(f) sf2a_out(paste0("20_", f))
id <- function(src) def$items[match(src, def$source[def$items])]
V <- fread(sf2a_out("19_mode_free.csv"))
sets <- list(`final 8` = V[set == "final 8 (17)", ids], `sf6_44 -> N6` = V[set == "final 8, sf6_44 -> N6", ids],
             `mode-free (N6 + sfc_31)` = V[set == "final 8, both -> N6 + sfc_31", ids])
sets <- lapply(sets, function(s) strsplit(s, " ")[[1]])
model <- function(set, ri, rescov) {
  m <- c(paste0("S =~ NA*", paste(set, collapse = " + ")), "S ~~ 1*S")
  if (ri) m <- c(m, paste0("RI =~ ", paste0("1*", set, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI")
  mm <- intersect(set, def$names_mode)
  if (rescov && length(mm) == 2) m <- c(m, paste(mm[1], "~~", mm[2]))
  m
}
fs <- c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
rows <- list(); loads <- list()
for (nm in names(sets)) for (v in list(c(F, F), c(T, F), c(T, T))) {
  if (v[2] && length(intersect(sets[[nm]], def$names_mode)) < 2) next
  lab <- paste0(nm, if (v[1]) " + RI" else "", if (v[2]) " + mode residual" else "")
  f <- sf2a_sem(model(sets[[nm]], v[1], v[2]), X[, sets[[nm]]], ordered = sets[[nm]])
  mi <- tryCatch({ x <- modindices(f, sort. = TRUE, maximum.number = 3); paste(sprintf("%s%s%s %.0f", def$source[x$lhs] %||% x$lhs, x$op,
                                                                                   def$source[x$rhs] %||% x$rhs, x$mi), collapse = "; ") },
                 error = function(e) NA_character_)
  rc <- if (v[2]) { pe <- parameterEstimates(f, standardized = TRUE); pe[pe$op == "~~" & pe$lhs != pe$rhs & pe$lhs %in% def$names_mode, "std.all"][1] } else NA
  rows[[lab]] <- data.table(model = lab, t(round(fitMeasures(f, fs), 3)), mode_residual_r = round(rc, 3), top_mi = mi)
  l <- lavInspect(f, "std")$lambda[sets[[nm]], "S"]; l <- l * ifelse(names(l) %in% rv, -1, 1); if (sum(l) < 0) l <- -l
  loads[[lab]] <- data.table(model = lab, source = def$source[names(l)], keyed_loading = round(l, 2))
}
fit <- rbindlist(rows); fwrite(fit, out20("fit.csv"))
ld <- dcast(rbindlist(loads), source ~ model, value.var = "keyed_loading"); fwrite(ld, out20("loadings.csv"))
print(fit); print(ld)
