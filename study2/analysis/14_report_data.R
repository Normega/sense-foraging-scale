# 14: every reported number in one JSON file (aggregates only: no participant-level
# rows, no participant codes), plus one line per confirmatory decision.

source("study2/analysis/config.R")
source("study2/analysis/sf2_models.R")
A <- readRDS(sf2_derived("sf2"))
rj <- function(f) if (file.exists(sf2_out(f))) fromJSON(sf2_out(f)) else NULL
rc <- function(f) if (file.exists(sf2_out(f))) fread(sf2_out(f)) else NULL
h1 <- rj("03_h1_decision.json"); v <- rj("04_validity_decisions.json"); h5 <- rj("05_h5_decision.json")
h6 <- rj("07_h6_decision.json"); h7 <- rj("08_h7_decision.json")

yn <- function(x) if (isTRUE(x)) "supported" else "not supported"
dec <- rbindlist(list(
  list("H1", "structure (ESEM + random intercept)", h1$h1),
  list("H1", paste0("factors confirmed: ", paste(h1$factors$factor[h1$factors$confirmed], collapse = ", ")),
       paste0("scored for individuals: ", paste(h1$factors$factor[h1$factors$scored_for_individuals], collapse = ", "))),
  list("H1 (secondary)", "Letting be separable from reversed wording", yn(h1$wording$letting_be_separable))),
  use.names = FALSE)
setnames(dec, c("hypothesis", "test", "result"))
if (!is.null(v)) dec <- rbind(dec,
  data.table(hypothesis = "H2", test = paste(v$h2$factor, "with", v$h2$score), result = sapply(v$h2$h2_supported, yn)),
  data.table(hypothesis = "H3a", test = paste(v$h2$factor, "with", v$h2$score, "< .80"), result = sapply(v$h2$h3a_supported, yn)),
  data.table(hypothesis = "H3b", test = "Valuing: SBI > BFI-2 Open-Mindedness", result = yn(v$h3b$supported)),
  data.table(hypothesis = "H3c", test = "no factor beyond |r| = .30 with BIDR SDE or IM", result = yn(v$h3c_supported)),
  data.table(hypothesis = "H3d", test = paste(v$h3d$factor, "not redundant with its neighbours"), result = sapply(v$h3d$supported, yn)),
  data.table(hypothesis = "H4", test = paste(v$h4$factor, "with", v$h4$score), result = sapply(v$h4$supported, yn)))
if (!is.null(h5)) dec <- rbind(dec,
  data.table(hypothesis = "H5a", test = paste0(h5$h5a$outcome, ": Access + Valuing beyond FFMQ-15 and MAIA-2, dR2 >= .02 (dR2 = ",
                                               round(h5$h5a$delta_r2, 3), ")"), result = sapply(h5$h5a$supported, yn)),
  data.table(hypothesis = "H5a", test = "at least one of the three outcomes", result = yn(h5$h5a_supported)),
  data.table(hypothesis = "H5b", test = "Access + Valuing predicts practice frequency beyond covariates (ordinal, LR test)",
             result = yn(h5$h5b_supported)))
if (!is.null(h6)) dec <- rbind(dec, data.table(hypothesis = "H6", test = paste(h6$h6$score, "dedicated > never, d >= .50"), result = sapply(h6$h6$supported, yn)))
if (!is.null(h7)) dec <- rbind(dec, data.table(hypothesis = "H7", test = paste(h7$h7$score, "ICC(2,1) >= .70"), result = sapply(h7$h7$supported, yn)))
fwrite(dec, sf2_out("14_decisions.csv"))

ver <- function(p) as.character(packageVersion(p))
out <- list(
  generated = format(Sys.time(), "%Y-%m-%d %H:%M %Z"), mode = SF2$mode, extraction = A$extraction,
  software = list(R = R.version.string, lavaan = ver("lavaan"), semTools = ver("semTools"), psych = ver("psych"), irr = ver("irr"),
                  data.table = ver("data.table"), GPArotation = ver("GPArotation")),
  settings = list(seeds = SF2[grep("^seed", names(SF2))], resamples_h3 = SF2$boot_h3, resamples_step7 = SF2$boot_step7,
                  rstarts = SF2$rstarts),
  decisions = dec,
  flow = A$flow, retest_flow = A$retest_flow, flags = rc("01_flags.csv"), submissions = rc("01_submissions.csv"),
  sample = rj("02_sample.json"), missing_by_measure = rc("02_missing_by_measure.csv"), internal_consistency = rc("02_internal_consistency.csv"),
  items = rc("02_item_descriptives.csv")[, !"text"],
  structure = list(decision = h1, fit = rc("03_fit.csv"), loadings = rc("03_h1_loadings.csv"),
                   factor_correlations = rc("03_h1_factor_correlations.csv"), wording_loadings = rc("03_wording_loadings.csv"),
                   model_a_modification_indices = rc("03_modelA_modification_indices.csv")),
  validity = list(decisions = v, latent = rc("04_latent_validity.csv"), observed = rc("04_observed_validity.csv"),
                  random_intercept = rc("04_random_intercept_validity.csv")),
  incremental = h5, invariance = rj("06_invariance_summary.json"), invariance_steps = rc("06_invariance.csv"),
  latent_means = rc("06_latent_means.csv"),
  known_groups = list(decision = h6, all = rc("07_known_groups.csv"), frequency_trend = rc("07_frequency_trend.csv"),
                      practice_invariance = rc("07_practice_invariance.csv"), practice_latent_means = rc("07_practice_latent_means.csv")),
  stability = list(decision = h7, icc = rc("08_stability.csv"), responders = rc("08_responders.csv")),
  step7 = rj("09_step7.json"), step7_candidates = rc("09_candidates.csv"), step7_bootstrap = rc("09_bootstrap_rates.csv"),
  step8 = rj("10_revised_form.json"), step9_forage = rj("11_forage.json"), step10_av14 = rj("12_av14.json"),
  exploratory = rj("13_exploratory.json"), following_week = rc("13_following_week.csv"))
write_json(out, sf2_out("14_report_study2.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", dataframe = "rows")
print(dec, nrows = 100)
sf2_log("report data: ", sf2_out("14_report_study2.json"))
