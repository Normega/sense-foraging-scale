# Runs the Study 2 pipeline in order. From the repo root:
#   SF2_MODE=simulated Rscript study2/analysis/run_all.R   00_simulate first, 99_check_recovery last
#   SF2_MODE=permuted  Rscript study2/analysis/run_all.R   the real extraction, answers shuffled within item (default)
#   SF2_MODE=real      Rscript study2/analysis/run_all.R   only after the preregistration is filed
# Resample counts: SF2_BOOT_H3 (default 2000), SF2_BOOT_STEP7 (default 200); cores: SF2_CORES.

source("study2/analysis/config.R")
steps <- c("01_prepare", "02_descriptives", "03_structure", "04_validity", "05_incremental", "06_invariance",
           "07_known_groups", "08_stability", "09_candidates", "10_revised_form", "11_forage", "12_av14",
           "13_exploratory", "14_report_data")
if (SF2$mode == "simulated") steps <- c("00_simulate", steps, "99_check_recovery")
times <- list()
for (s in steps) {
  sf2_log("==== ", s)
  t0 <- Sys.time()
  source(sprintf("study2/analysis/%s.R", s), local = new.env())
  times[[s]] <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
}
fwrite(data.table(step = names(times), minutes = unlist(times), resamples_h3 = SF2$boot_h3, resamples_step7 = SF2$boot_step7,
                  cores = SF2$cores), sf2_out("run_times.csv"))
sf2_log("done in ", sum(unlist(times)), " minutes; outputs in ", sf2_out(""))
