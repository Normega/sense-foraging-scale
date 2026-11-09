# Runs the Study 2b pipeline in order. From the repo root:
#   SF2B_MODE=simulated Rscript study2b/analysis/run_all.R   00_simulate first, 99_check_recovery last
#   SF2B_MODE=permuted  Rscript study2b/analysis/run_all.R   the real extraction, answers shuffled within item (default)
#   SF2B_MODE=real      Rscript study2b/analysis/run_all.R   only after the preregistration is filed
# SF2B_BOOT sets the H3d resamples (registered 2,000), SF2B_CORES the workers,
# SF2B_DFI_REPS the model-specific-cutoff replications (200), SF2B_FROM=<step> starts part-way.

source("study2b/analysis/config.R")
steps <- c("01_prepare", "02_descriptives", "03_structure", "04_validity", "05_criterion", "06_incremental",
           "07_known_groups", "08_exploratory", "09_report")
if (SF2B$mode == "simulated") steps <- c("00_simulate", steps, "99_check_recovery")
from <- Sys.getenv("SF2B_FROM", "")
if (nzchar(from)) steps <- steps[which(steps == from):length(steps)]
times <- list()
for (s in steps) {
  sf2b_log("==== ", s)
  t0 <- Sys.time()
  source(sprintf("study2b/analysis/%s.R", s), local = new.env())
  times[[s]] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
}
fwrite(data.table(step = names(times), minutes = unlist(times), resamples_h3d = SF2B$boot_h3d, cores = SF2B$cores),
       sf2b_out(if (nzchar(from)) "run_times_partial.csv" else "run_times.csv"))
sf2b_log("done in ", sum(unlist(times)), " minutes; outputs in ", sf2b_out(""))
