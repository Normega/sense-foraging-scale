# 05: H4, criterion validity (Section 5, model 3), on the observed 8 (and the 10, reported).
#   (a) Spearman rho of the 8 with past-week days of deliberately pausing (crit_pause_days)
#   (b) Spearman rho of the 8 with the engaged pause (0/1), among those who reached the
#       offer; and, reported beside it, a logistic regression of the engaged pause on
#       the 8 with device class as a covariate (PREREG: desktop, phone or tablet, as
#       the pause recorded it)
# Each is supported when rho > 0 with p < .05 (two-sided) after Holm correction across
# (a) and (b). Also without flagged participants.

source("study2b/analysis/config.R")
source("study2b/analysis/sf2b_models.R")
A <- readRDS(sf2b_derived("sf2b")); M <- A$meta; S <- A$scores
stopifnot(identical(M$participant, S$participant))

h4 <- function(form, keep = rep(TRUE, nrow(S))) {
  a <- sf2b_r_ci(S[[form]][keep], S$crit_pause_days[keep], "spearman")
  b <- sf2b_r_ci(S[[form]][keep], S$engaged_pause[keep], "spearman")
  holm <- p.adjust(c(a$p, b$p), "holm")
  d <- data.table(y = S$engaged_pause, x = as.numeric(scale(S[[form]])), device = M$pause_device)[keep][complete.cases(y, x, device)]
  d[, device := droplevels(device)]
  lg <- glm(if (nlevels(d$device) > 1) y ~ x + device else y ~ x, family = binomial, data = d)
  co <- summary(lg)$coefficients["x", ]
  list(a = c(a, p_holm = holm[1], supported = isTRUE(a$r > 0 && holm[1] < .05)),
       b = c(b, p_holm = holm[2], supported = isTRUE(b$r > 0 && holm[2] < .05)),
       b_logistic = list(n = nrow(d), odds_ratio_per_sd = exp(co[["Estimate"]]),
                         lower = exp(co[["Estimate"]] - qnorm(.975) * co[["Std. Error"]]),
                         upper = exp(co[["Estimate"]] + qnorm(.975) * co[["Std. Error"]]), p = co[["Pr(>|z|)"]],
                         covariate = if (nlevels(d$device) > 1) "device class" else "none (one device class)"))
}
out <- list(sf8 = h4("sf8"), sf10 = h4("sf10"), sf8_unflagged = h4("sf8", !M$flagged),
            pause_counts = list(reached_offer = sum(!is.na(S$engaged_pause)), took = sum(M$pause_taken, na.rm = TRUE),
                                engaged = sum(S$engaged_pause, na.rm = TRUE)))
rnd <- function(x) if (is.list(x)) lapply(x, rnd) else if (is.numeric(x)) signif(x, 4) else x
write_json(rnd(out), sf2b_out("05_h4_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
sf2b_log("H4a rho ", r3(out$sf8$a$r), " (Holm p ", signif(out$sf8$a$p_holm, 3), "): ", out$sf8$a$supported,
         "; H4b rho ", r3(out$sf8$b$r), " (Holm p ", signif(out$sf8$b$p_holm, 3), "): ", out$sf8$b$supported)
