# 09: every confirmatory decision on one line (09_decisions.csv), and the numbers behind
# them with the sample and the exploratory results (09_report_study2b.json, aggregates only).

source("study2b/analysis/config.R")
rd <- function(f) fromJSON(sf2b_out(f), simplifyVector = FALSE)
h1 <- rd("03_h1_decision.json"); v <- rd("04_validity_decisions.json"); h4 <- rd("05_h4_decision.json")
h5 <- rd("06_h5_decision.json"); h6 <- rd("07_h6_decision.json")
yn <- function(x) if (isTRUE(x)) "supported" else "not supported"
f <- function(x, d = 3) if (is.null(x)) "NA" else formatC(as.numeric(x), format = "f", digits = d)
ci <- function(l, u) sprintf("[%s, %s]", f(l), f(u))
h1s <- h1$h1
rows <- list(
  list("H1", "one factor + agreement fits the 8", yn(h1s$supported),
       sprintf("CFI %s; largest residual %s (freeing it: omega %s, loadings %s); close-fit p %s; RMSEA %s %s",
               f(h1s$cfi), f(h1s$max_residual), f(h1s$freeing_it_changes_omega_by), f(h1s$freeing_it_changes_a_loading_by_up_to),
               f(h1s$closefit_p), f(h1s$rmsea), ci(h1s$rmsea_ci[[1]], h1s$rmsea_ci[[2]]))))
for (x in v$sf8$h2) rows[[length(rows) + 1]] <- list("H2", paste("r with", x$score, "> .30 (lower bound)"), yn(x$h2), sprintf("r %s %s", f(x$r), ci(x$lower, x$upper)))
rows[[length(rows) + 1]] <- list("H2", sprintf("at least %d of the 4 lower bounds > .30", SF2B$h2_needed), yn(v$sf8$H2),
                                 paste(v$sf8$H2_correlations_holding, "of 4"))
rows[[length(rows) + 1]] <- list("H3a", "each H2 upper bound < .80", yn(v$sf8$H3a), paste(vapply(v$sf8$h2, function(x) sprintf("%s (%s)", f(x$upper), x$rc_class), ""), collapse = ", "))
for (x in v$sf8$h3b) rows[[length(rows) + 1]] <- list(paste0("H3b-", x$contrast), paste(x$first, ">", x$second), yn(x$holds),
                                                      sprintf("difference %s %s", f(x$difference), ci(x$lower, x$upper)))
rows[[length(rows) + 1]] <- list("H3b", "at least 4 of the 5 contrasts", yn(v$sf8$H3b), paste(v$sf8$H3b_contrasts_holding, "of 5"))
for (x in v$sf8$h3c) rows[[length(rows) + 1]] <- list("H3c", paste(x$score, "within +/-.30"), yn(x$holds), sprintf("r %s %s", f(x$r), ci(x$lower, x$upper)))
rows[[length(rows) + 1]] <- list("H3d", "multiple R with the four neighbours < .80 (upper bound)", yn(v$sf8$H3d),
                                 sprintf("R %s %s, %s resamples; Ronkko & Cho class: %s", f(v$sf8$h3d$R), ci(v$sf8$h3d$lower, v$sf8$h3d$upper), v$sf8$h3d$resamples, v$sf8$h3d$rc_class))
for (x in v$sf8$h3e) rows[[length(rows) + 1]] <- list("H3e", sprintf("%s within (%s, %s)", x$score, f(x$bound_lower, 2), f(x$bound_upper, 2)), yn(x$holds),
                                                      sprintf("r %s %s", f(x$r), ci(x$lower, x$upper)))
rows[[length(rows) + 1]] <- list("H4a", "rho with past-week pause days > 0 (Holm)", yn(h4$sf8$a$supported), sprintf("rho %s, Holm p %s, n %s", f(h4$sf8$a$r), signif(h4$sf8$a$p_holm, 3), h4$sf8$a$n))
rows[[length(rows) + 1]] <- list("H4b", "rho with the engaged pause > 0 (Holm)", yn(h4$sf8$b$supported), sprintf("rho %s, Holm p %s, n %s", f(h4$sf8$b$r), signif(h4$sf8$b$p_holm, 3), h4$sf8$b$n))
h5a <- fread(sf2b_out("06_h5a.csv"))[form == "sf8"]; h5b <- fread(sf2b_out("06_h5b.csv"))[form == "sf8"]
for (k in seq_len(nrow(h5a))) rows[[length(rows) + 1]] <- list("H5a", paste("delta R2 >= .02 for", h5a$outcome[k]), yn(h5a$supported[k]),
                                                                sprintf("delta R2 %s, Holm p %s, n %s", f(h5a$delta_r2[k]), signif(h5a$p_holm[k], 3), h5a$n[k]))
rows[[length(rows) + 1]] <- list("H5a", "at least one outcome", yn(h5$H5a), "")
rows[[length(rows) + 1]] <- list("H5b", "practice frequency beyond the covariates", yn(h5$H5b), sprintf("OR per SD %s, LR p %s, n %s", f(h5b$odds_ratio_per_sd, 2), signif(h5b$p, 3), h5b$n))
rows[[length(rows) + 1]] <- list("H6", "current vs never, d >= .50", yn(h6$H6), sprintf("d %s %s, p %s (%s vs %s)", f(h6$d, 2), ci(h6$d_ci[[1]], h6$d_ci[[2]]), h6$p, h6$n_current, h6$n_never))
dec <- rbindlist(lapply(rows, function(r) data.table(hypothesis = r[[1]], test = r[[2]], result = r[[3]], estimate = r[[4]])))
fwrite(dec, sf2b_out("09_decisions.csv"))

report <- list(mode = SF2B$mode, written = format(Sys.time(), "%Y-%m-%d %H:%M"),
               software = list(R = R.version.string, lavaan = as.character(packageVersion("lavaan")),
                               psych = as.character(packageVersion("psych")), MASS = as.character(packageVersion("MASS"))),
               flow = fread(sf2b_out("01_flow.csv")), sample = rd("02_sample.json"), h1 = h1, validity = v, h4 = h4, h5 = h5, h6 = h6,
               exploratory = rd("08_exploratory.json"))
write_json(report, sf2b_out("09_report_study2b.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(dec, nrows = 100)
