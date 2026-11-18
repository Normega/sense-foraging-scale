# 33: the 8 by country (six listed countries vs other), and race/ethnicity with country
# held constant. EXPLORATORY, NOT REGISTERED.
#
# Not part of the filed registration (osf.io/d4a5u) and not in its list of exploratory
# analyses. Specified on 18 November 2026 (placement card S3-30), after the registered
# analyses and 10-32 had been seen.
#
# Country is the background question's option (01: M$country2): one of the six listed
# countries (au, ca, ie, nz, uk, us) 173, other 545, missing 2. "Other" is not broken
# down: that needs the free-text answers or Prolific's field (decision X-2b, open), and
# neither is read here. Six countries (173) is below SF2B$country_min (200), which is why
# the registered country comparison did not run; this comparison is under-powered.
#   invariance   24's procedure unchanged (H8's steps, then 24's intercepts step and
#                latent means), groups six countries (reference) vs other.
#   race given country   MIMIC on the H1 model (WLSMV, theta; S's residual variance 1;
#                the agreement factor regressed too): S ~ not_white, then S ~ not_white +
#                other_country, and the same for Black vs White; the coefficient of race
#                is in residual-SD units of S (07's MIMIC scale). Cases with race answered
#                and a country.
# Race/ethnicity as 22-24 read it (option ids only; free text never read).
# Writes aggregates only:
#   output/real/33_invariance.csv, 33_latent_means.csv, 33_race_given_country.csv,
#   33_crosstab.csv, 33_exploratory_country.json
# Run from the repo root:  SF2B_MODE=real Rscript study2b/analysis/33_exploratory_country.R

# Set Up ---------
## Load libraries ---------
packages <- c("jsonlite", "data.table", "lavaan")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings: 24's data and functions, nothing fitted ---------
ex24 <- parse("study2b/analysis/24_exploratory_invariance_race.R", keep.source = FALSE)
run_at <- which(vapply(ex24, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("INV")), TRUE))
stopifnot(length(run_at) == 1)
for (i in seq_len(run_at - 1)) eval(ex24[[i]], envir = globalenv())
stopifnot(exists("run_means"), exists("race_me"), !exists("INV"))
T0 <- Sys.time()
RACE <- GROUPS

cg <- M$country2
GROUPS <- list(country = list(label = "country: six listed countries vs other",
                              g = fifelse(cg == "six_countries", "1_six", fifelse(cg == "other", "2_other", NA_character_))))
n_c <- table(GROUPS$country$g); GROUPS$country$n <- as.list(n_c); GROUPS$country$below_min <- names(n_c)[n_c < SF2B$country_min]
sf2b_log("country: ", paste(names(n_c), n_c, sep = " = ", collapse = ", "))

## Invariance and latent means ---------
inv <- sf2b_run_mg(X, I8, GROUPS$country$g, GROUPS$country$label)
if (!"noninvariant" %in% names(inv$table)) inv$table[, noninvariant := NA]
steps <- cbind(comparison = "country", inv$table); means <- data.table()
last <- tail(inv$table, 1)
if (grepl("^loadings", last$step) && !isTRUE(last$noninvariant)) {
  mm <- run_means("country", inv); steps <- rbind(steps, cbind(comparison = "country", mm$table), fill = TRUE); means <- mm$means
} else sf2b_log("country: loadings step not held; means not compared")

## Race given country (MIMIC) ---------
other_c <- ifelse(cg == "other", 1L, ifelse(cg == "six_countries", 0L, NA_integer_))
xt <- as.data.table(table(race = RACE$race3$g, country = GROUPS$country$g, useNA = "no"))
mimic <- function(race_key, focal) {
  r <- RACE[[race_key]]$g
  nw <- ifelse(r == focal, 1L, ifelse(r == "1_white", 0L, NA_integer_))
  d <- as.data.frame(X[, I8]); d$race <- nw; d$other_c <- other_c
  d <- d[!is.na(d$race) & !is.na(d$other_c), ]
  out <- list()
  for (adj in c(FALSE, TRUE)) {
    rhs <- if (adj) "race + other_c" else "race"
    syn <- c(sf2b_scale_syntax(I8), paste0("S ~ bS*", sub("race", "race", rhs)), paste0("RI ~ bR*", rhs))
    syn <- sub("S ~ bS\\*race \\+ other_c", "S ~ bS*race + cS*other_c", syn)
    syn <- sub("RI ~ bR\\*race \\+ other_c", "RI ~ bR*race + cR*other_c", syn)
    f <- sf2b_try(sf2b_sem(syn, d, ordered = I8))
    if (sf2b_failed(f)) { sf2b_log(race_key, " adj=", adj, " FAILED: ", sf2b_msg(f)); next }
    pe <- as.data.table(parameterEstimates(f))
    g <- function(l) { x <- pe[label == l]; if (!nrow(x)) c(NA, NA, NA) else c(x$est, x$ci.lower, x$ci.upper) }
    out[[length(out) + 1]] <- data.table(comparison = paste0(focal, " vs 1_white"), adjusted_for_country = adj, n = nrow(d),
      S_race = g("bS")[1], S_race_lo = g("bS")[2], S_race_hi = g("bS")[3],
      S_other_country = g("cS")[1], S_other_country_lo = g("cS")[2], S_other_country_hi = g("cS")[3],
      RI_race = g("bR")[1], RI_other_country = g("cR")[1])
  }
  rbindlist(out)
}
rgc <- rbind(mimic("race2", "2_not_white"), mimic("race_bw", "2_black"))

## Write ---------
num <- function(DT, k = 4) { cols <- names(DT)[vapply(DT, is.numeric, TRUE)]; DT[, (cols) := lapply(.SD, function(x) round(x, k)), .SDcols = cols]; DT }
fwrite(num(copy(steps)), sf2b_out("33_invariance.csv")); fwrite(num(copy(means)), sf2b_out("33_latent_means.csv"))
fwrite(num(copy(rgc)), sf2b_out("33_race_given_country.csv")); fwrite(xt, sf2b_out("33_crosstab.csv"))
write_json(list(
  analysis = "33: the 8 by country (six listed vs other) and race/ethnicity with country held constant. EXPLORATORY, NOT REGISTERED (osf.io/d4a5u does not include it); specified 18 November 2026 after 10-32 were seen (card S3-30). Six countries n = 173, below the registered 200: under-powered. 'Other' not broken down (X-2b open).",
  n = as.list(n_c), released = inv$free, steps = num(copy(steps)), latent_means = num(copy(means)),
  race_given_country = num(copy(rgc)), crosstab = xt,
  minutes = round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1)),
  sf2b_out("33_exploratory_country.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
print(xt)
print(num(copy(steps))[, .(step, released, cfi, rmsea, d_cfi, d_rmsea, noninvariant)])
if (nrow(means)) print(num(copy(means))[, .(model, group, n, S_d, S_d_lo, S_d_hi, RI_d, RI_d_lo, RI_d_hi)])
print(num(copy(rgc), 3))
sf2b_log("33 done in ", round(as.numeric(difftime(Sys.time(), T0, units = "mins")), 1), " min")
