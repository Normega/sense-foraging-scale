# 99 (simulated mode only): does the pipeline recover what 00_simulate.R planted?
# Pass or fail per check, in 99_recovery.csv.
#   flow      every planted exclusion, flag and repeat submission handled as registered
#   H1        content and agreement loadings against the population; the rule passes
#             (the model is true)
#   validity  each latent correlation's 95% CI covers its population value (expect
#             about 95%) and its sample value (the correlation of the generated latent
#             scores in this sample; any shortfall there is the code's, not chance's)
#   H3b       each contrast's CI covers the population difference
#   invariance  the planted loading difference by gender is found (released at the loadings step)

source("study2b/analysis/config.R")
stopifnot(SF2B$mode == "simulated")
tr <- readRDS(sf2b_derived("simulation_truth")); A <- readRDS(sf2b_derived("sf2b")); D <- A$defs
st <- readRDS(sf2b_derived("structure")); va <- readRDS(sf2b_derived("validity"))
M <- A$meta_all; pl <- tr$planted
chk <- list()
add <- function(check, pass, detail) chk[[length(chk) + 1L]] <<- data.table(check = check, pass = isTRUE(pass), detail = detail)

add("enrollments", nrow(M) == pl$enrolled, sprintf("%d of %d", nrow(M), pl$enrolled))
add("staff tests excluded", sum(M$excl_test) == pl$tests, sprintf("%d of %d", sum(M$excl_test), pl$tests))
add("no scale excluded", sum(M$excl_no_scale) == pl$no_scale, sprintf("%d of %d", sum(M$excl_no_scale), pl$no_scale))
add("two or more failed checks excluded", sum(M$excl_attention) == pl$two_checks, sprintf("%d of %d", sum(M$excl_attention), pl$two_checks))
add("one failed check flagged", sum(M$flag_one_check & M$analysable) == pl$one_check, sprintf("%d of %d", sum(M$flag_one_check & M$analysable), pl$one_check))
add("analysable", sum(M$analysable) == pl$analysable, sprintf("%d of %d", sum(M$analysable), pl$analysable))
flagged_set <- function(col, codes) all(M[participant %in% codes, get(col)])
add("speeders flagged", flagged_set("flag_speeder", pl$speeders), paste(sum(M$flag_speeder), "flagged"))
add("longstrings flagged", flagged_set("flag_longstring", pl$longstrings), paste(sum(M$flag_longstring), "flagged"))
add("many 'Prefer not to answer' flagged", flagged_set("flag_pna", pl$many_pna), paste(sum(M$flag_pna), "flagged"))
sub <- A$submissions[slug == "sf2b-scale"]
add("flagged resubmissions seen and dropped", sum(sub$flagged_resubmission) == length(pl$resubmission), paste(sum(sub$flagged_resubmission), "flagged"))
# The re-entry carried different answers 15 minutes later: the first submission must be the one used.
ext <- fromJSON(sf2b_raw_file(), simplifyVector = FALSE)$rows
first_ok <- all(vapply(pl$reentry, function(p) {
  rs <- Filter(function(r) identical(r$participant, p) && identical(r$slug, "sf2b-scale"), ext)
  rs <- rs[order(vapply(rs, `[[`, "", "received_at"))]
  first <- vapply(D$items, function(i) { v <- rs[[1]]$responses[[i]]; if (is.numeric(v)) v else NA_real_ }, 0)
  isTRUE(all.equal(unname(first), unname(A$scale_raw[p, D$items])))
}, TRUE))
add("a later re-entry does not replace the first submission", first_ok, paste(length(pl$reentry), "re-entries"))

# H1
# The item with the planted gender difference has no single population loading in the
# pooled sample (it is a mixture), so it is left out of the content-loading check.
pop <- tr$population; L <- st$h1$sf8$loadings
ix <- match(L$item, pop$items); inv_item <- L$item == tr$noninvariance$item
rmse <- sqrt(mean((L$content - pop$loading_S[ix])[!inv_item]^2)); rmse_ri <- sqrt(mean((L$agreement - pop$loading_RI[ix])^2))
add("H1 content loadings recovered (RMSE < .05)", rmse < .05, sprintf("RMSE %.3f over the 7 statements without the planted difference", rmse))
add("H1 agreement loadings recovered (RMSE < .05)", rmse_ri < .05, sprintf("RMSE %.3f", rmse_ri))
add("H1 rule passes (model true)", st$h1$sf8$supported, sprintf("CFI %.3f, max residual %.3f, close-fit p %.3f", st$h1$sf8$stats$cfi.scaled, st$h1$sf8$max_residual, st$h1$sf8$closefit_p))

# Validity
lat <- va$latent[form == "sf8" & sample == "all" & score %in% names(tr$rho)]
lat[, `:=`(truth = tr$rho[score], sample_truth = tr$sample_r[score])]
lat[, `:=`(covers = truth >= ci.lower & truth <= ci.upper, covers_sample = sample_truth >= ci.lower & sample_truth <= ci.upper)]
fwrite(lat[, .(score, est, ci.lower, ci.upper, truth, sample_truth = round(sample_truth, 3), covers, covers_sample)], sf2b_out("99_validity_recovery.csv"))
single <- c("lifesat", "engaged_pause", "crit_pause_days", "crit_notice_first", "crit_outdoors")
multi <- lat[!(score %in% single)]
add("latent validity CIs cover the population value (>= 90%)", mean(multi$covers) >= .90, sprintf("%d of %d", sum(multi$covers), nrow(multi)))
add("latent validity CIs cover the sample value (>= 90%)", mean(multi$covers_sample) >= .90, sprintf("%d of %d", sum(multi$covers_sample), nrow(multi)))
add("every validity factor oriented as keyed (no model refused)", all(vapply(va$results, `[[`, TRUE, "ok")),
    paste(sum(!vapply(va$results, `[[`, TRUE, "ok")), "models failed or refused"))
h3b <- va$decisions$sf8$h3b
h3b[, truth := tr$rho[first] - tr$rho[second]]
h3b_cov <- h3b[!is.na(truth), truth >= lower & truth <= upper]
add("H3b contrast CIs cover the population difference", mean(h3b_cov) >= .75, sprintf("%d of %d (contrast 3's total has no exact population value)", sum(h3b_cov), length(h3b_cov)))

# Invariance (gender): the planted item should be the one released at the loadings step.
inv <- fread(sf2b_out("08_invariance.csv"))[analysis == "gender (man / woman)"]
rel <- inv[grepl("^loadings", step) & released != "", released]
add("planted gender loading difference found", length(rel) > 0 && grepl(tr$noninvariance$item, tail(rel, 1)),
    if (length(rel)) paste("released:", tail(rel, 1)) else "no loadings-step release")

res <- rbindlist(chk)
fwrite(res, sf2b_out("99_recovery.csv"))
print(res, nrows = 100)
sf2b_log(sum(res$pass), " of ", nrow(res), " recovery checks pass")
