# Runs the Study 2a pipeline in order. From the repo root:
#   SF2A_MODE=permuted  Rscript study2a/analysis/run_all.R   (default)
#   SF2A_MODE=simulated ...   then 99_check_recovery.R runs too
#   SF2A_MODE=real      ...
# SF2A_K sets the number of content factors in the primary model (default 5),
# chosen from 04's output by decision K.

source("study2a/analysis/config.R")
steps <- c("01_prepare", "02_items", "03_vignettes", "04_factor_number", "05_primary_model", "06_wording",
           "07_ius", "08_reliability", "09_selection", "10_sensitivity", "11_exploratory")
from <- Sys.getenv("SF2A_FROM", "")
if (nzchar(from)) steps <- steps[which(steps == from):length(steps)]
if (SF2A$mode == "simulated" && SF2A$k_primary == 5) steps <- c(steps, "99_check_recovery")
for (s in steps) {
  sf2a_log("==== ", s)
  t0 <- Sys.time()
  source(sprintf("study2a/analysis/%s.R", s), local = new.env())
  sf2a_log(s, " took ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
sf2a_log("done; outputs in ", sf2a_out(""))
