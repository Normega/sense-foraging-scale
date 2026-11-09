# 07: H6, known groups (Section 5, model 5). Current practitioners against people who
# have never practised regularly, on the observed 8: Cohen's d (pooled SD) with 95% CI
# and a two-sided Student t-test (PREREG: Student's, so the test matches the pooled-SD d;
# d's CI is the coefficient's CI over the pooled SD). Supported when the difference favours practitioners,
# p < .05 and d >= 0.50.
# Sensitivity: adjusting for spirituality and religiosity; dedicated practitioners
# against never; the latent difference with agreement modelled (a MIMIC: the H1 model
# with S regressed on the group, S's residual variance 1, so the coefficient is the
# difference in within-group latent SDs); without flagged participants. The 10 is
# reported beside the 8.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
A <- readRDS(sf2b_derived("sf2b")); D <- A$defs; M <- A$meta; S <- A$scores

g_cur <- fifelse(M$practice_group == "current", 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
g_ded <- fifelse(M$dedicated, 1L, fifelse(M$practice_group == "never", 0L, NA_integer_))
rows <- list()
add <- function(label, form, g, covars = NULL, keep = rep(TRUE, nrow(M))) {
  y <- S[[form]]; y[!keep] <- NA
  rows[[length(rows) + 1L]] <<- cbind(analysis = label, form = form, sf2b_kg(y, g, covars))
}
for (form in c("sf8", "sf10")) {
  add("H6: current vs never", form, g_cur)
  add("adjusted for spirituality and religiosity", form, g_cur, M[, .(spiritual, religious)])
  add("dedicated vs never", form, g_ded)
}
add("without flagged participants", "sf8", g_cur, keep = !M$flagged)
kg <- rbindlist(rows)
kg[, supported := analysis == "H6: current vs never" & coef > 0 & p < .05 & d >= SF2B$h6_d]

# Latent difference (MIMIC).
mimic <- function(items) {
  d <- as.data.frame(A$scale_raw[, items]); d$grp <- g_cur; d <- d[!is.na(d$grp), ]
  syn <- c(paste0("S =~ NA*", paste(items, collapse = " + ")), "S ~~ 1*S",
           paste0("RI =~ ", paste0("1*", items, collapse = " + ")), "RI ~~ RI", "S ~~ 0*RI", "S ~ b*grp", "RI ~ 0*grp")
  f <- sf2b_try(sf2b_sem(syn, d, ordered = items))
  if (sf2b_failed(f)) return(data.table(n = nrow(d), latent_d = NA_real_, error = sf2b_msg(f)))
  pe <- parameterEstimates(f); b <- pe[pe$label == "b", ]
  data.table(n = nrow(d), latent_d = b$est, lower = b$ci.lower, upper = b$ci.upper, p = b$pvalue, error = "")
}
lat <- rbind(cbind(form = "sf8", mimic(D$sf8)), cbind(form = "sf10", mimic(D$sf10)), fill = TRUE)
num <- function(DT) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, signif, 4), .SDcols = cols]; DT }
fwrite(num(copy(kg)), sf2b_out("07_known_groups.csv")); fwrite(num(copy(lat)), sf2b_out("07_latent_difference.csv"))
h <- kg[analysis == "H6: current vs never" & form == "sf8"]
write_json(list(H6 = isTRUE(h$supported), d = r3(h$d), d_ci = r3(c(h$d_lower, h$d_upper)), p = signif(h$p, 3),
                n_current = h$n_group, n_never = h$n_reference), sf2b_out("07_h6_decision.json"), auto_unbox = TRUE, pretty = TRUE)
print(kg[, .(analysis, form, n_group, n_reference, d = r3(d), d_lower = r3(d_lower), p = signif(p, 3), supported)]); print(lat)
