# Follow-up invitations (Section 2, "Randomization"): a simple random sample of
# eligible completers, drawn by script with a recorded seed (SF2$seed_retest).
# Not part of run_all.R: run once, about 14 days after the main sessions, on a fresh
# extraction (real mode):
#   SF2_MODE=real Rscript study2/analysis/draw_retest_sample.R [n = 400]
# PREREG: "eligible completers" = reached the end of the main session (the validity
# block), not a staff test, and not excluded for failing both attention checks.
# Writes raw/retest_invited_<date>.csv: the Study 2 codes of those invited. Turning
# them into the Prolific allowlist needs the identity table (lab only); 08 reads this
# file to compare responders with invited non-responders.

source("study2/analysis/config.R")
args <- commandArgs(trailingOnly = TRUE)
n_invite <- if (length(args)) as.integer(args[1]) else 400L   # Section 2: about 400 invited for about 300 completions
stopifnot(SF2$mode != "permuted")
source("study2/analysis/01_prepare.R", local = new.env())   # the registered exclusions, on the latest extraction
A <- readRDS(sf2_derived("sf2"))
eligible <- A$meta[analysable == TRUE & reached_end_validity == TRUE, participant]
set.seed(SF2$seed_retest)
invited <- sort(sample(eligible, min(n_invite, length(eligible))))
f <- sf2_path("raw", sprintf("retest_invited_%s.csv", format(Sys.Date(), "%Y%m%d")))
fwrite(data.table(participant = invited), f)
sf2_log(length(invited), " of ", length(eligible), " eligible completers invited (seed ", SF2$seed_retest, ") -> ", f)
