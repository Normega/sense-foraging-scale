# Builds the Study 2b OSF filing package from one commit of this repository.
#
#   Rscript study2/filing/build_filing.R          (from the repo root)
#
# Writes to SF2B_FILING_DIR (default F:/radlab_data/sf_study2b/filing/<commit>/):
#   sf2b_code_<commit>.zip       git archive of the code, materials, power simulations,
#                                decision log and preregistration text at that commit
#   preregistration_2b.html/.pdf the preregistration as rendered from that commit, with the
#                                Code paragraph filled in (commit, archive, SHA-256)
#   osf_fields.txt               the text for each field of OSF's Preregistration template,
#                                as plain text to paste
#   manifest.txt                 commit, build time, every output's SHA-256
#
# The text is read from the commit (git show), not the working tree, so the PDF, the
# fields and the archive always match. The Code paragraph's [AT FILING] markers are filled
# here (a commit cannot name its own hash). Any other [AT FILING] marker left in the text,
# or uncommitted changes to the archived paths, make this a DRAFT build: every output is
# marked DRAFT and named *_DRAFT. File only a build without that mark.

# Set Up ---------
## Load libraries ---------
packages <- c("commonmark", "digest")
new_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(new_packages)) install.packages(new_packages)
for (thispack in packages) {
  library(thispack, character.only = TRUE, quietly = TRUE, verbose = FALSE)
}

## Settings ---------
PREREG <- "study2/preregistration_draft.md"
ARCHIVED <- c("study2b/analysis", "study2b/materials", "study2/materials", "study2/power",
              "study2/decisions_after_2a.md", "study2/prereg_2b_plan.md", PREREG, "study2/filing")
CHROME <- Sys.getenv("SF2B_CHROME", "C:/Program Files/Google/Chrome/Application/chrome.exe")

git <- function(...) {
  out <- system2("git", c(...), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status"))) stop("git ", paste(c(...), collapse = " "), " failed:\n", paste(out, collapse = "\n"))
  out
}
sha256 <- function(f) digest(file = f, algo = "sha256")

commit <- git("rev-parse", "HEAD")
short <- substr(commit, 1, 7)
dirty <- git("status", "--porcelain", "--", ARCHIVED)
out_dir <- file.path(Sys.getenv("SF2B_FILING_DIR", "F:/radlab_data/sf_study2b/filing"), short)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## The code archive ---------
zip_name <- sprintf("sf2b_code_%s.zip", short)
zip_path <- file.path(out_dir, zip_name)
git("archive", "--format=zip", paste0("--prefix=sense-foraging-scale-", short, "/"), "-o", zip_path, commit, "--", ARCHIVED)
zip_sha <- sha256(zip_path)

## The text, from the commit ---------
md <- paste(git("show", paste0(commit, ":", PREREG)), collapse = "\n")
md <- sub("sf2b_code_[AT FILING].zip", zip_name, md, fixed = TRUE)
md <- sub("at commit [AT FILING]", paste("at commit", commit), md, fixed = TRUE)
md <- sub("SHA-256 [AT FILING]", paste("SHA-256", zip_sha), md, fixed = TRUE)
open_markers <- lengths(regmatches(md, gregexpr("[AT FILING]", md, fixed = TRUE)))
draft <- open_markers > 0 || length(dirty) > 0
tag <- if (draft) "_DRAFT" else ""
if (draft) {
  why <- c(if (open_markers > 0) sprintf("%d [AT FILING] entries still open in the text", open_markers),
           if (length(dirty)) sprintf("uncommitted changes under the archived paths:\n  %s", paste(dirty, collapse = "\n  ")))
  message("DRAFT build: ", paste(why, collapse = "; "))
}

## HTML and PDF ---------
body <- markdown_html(md, extensions = TRUE)
css <- "
@page { size: Letter; margin: 18mm 16mm; }
body { font-family: 'Segoe UI', Arial, sans-serif; font-size: 10.5pt; line-height: 1.45; color: #111; max-width: 760px; margin: 0 auto; }
h1 { font-size: 17pt; line-height: 1.25; } h2 { font-size: 13.5pt; margin-top: 1.6em; border-bottom: 1px solid #999; }
h3 { font-size: 11.5pt; margin-top: 1.2em; }
table { border-collapse: collapse; margin: .6em 0; font-size: 9.5pt; } th, td { border: 1px solid #aaa; padding: 3px 6px; vertical-align: top; text-align: left; }
th { background: #eee; } code { font-family: Consolas, monospace; font-size: 9pt; }
.draft { position: fixed; top: 40%; left: 0; right: 0; text-align: center; font-size: 80pt; color: rgba(200,0,0,.12); transform: rotate(-30deg); pointer-events: none; }
.stamp { font-size: 8.5pt; color: #555; border-top: 1px solid #ccc; margin-top: 2em; padding-top: .4em; }
"
stamp <- sprintf("<p class='stamp'>Built %s from commit %s%s.</p>", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), commit,
                 if (draft) " &mdash; DRAFT, not for filing" else "")
html <- paste0("<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Sense Foraging Scale Study 2b preregistration</title><style>",
               css, "</style></head><body>", if (draft) "<div class='draft'>DRAFT</div>", body, stamp, "</body></html>")
html_path <- file.path(out_dir, paste0("preregistration_2b", tag, ".html"))
writeLines(html, html_path, useBytes = TRUE)
pdf_path <- file.path(out_dir, paste0("preregistration_2b", tag, ".pdf"))
if (file.exists(pdf_path)) file.remove(pdf_path)
system2(CHROME, c("--headless=new", "--disable-gpu", "--mute-audio", "--no-pdf-header-footer",
                  shQuote(paste0("--print-to-pdf=", normalizePath(pdf_path, mustWork = FALSE))),
                  shQuote(paste0("file:///", normalizePath(html_path, winslash = "/")))), stdout = FALSE, stderr = FALSE)
if (!file.exists(pdf_path)) stop("Chrome did not write the PDF")

## OSF fields, plain text ---------
lines <- strsplit(md, "\n", fixed = TRUE)[[1]]
heads <- grep("^#{2,3} ", lines)
sec <- lapply(seq_along(heads), function(k) {
  end <- if (k < length(heads)) heads[k + 1] - 1 else length(lines)
  list(head = sub("^#+ ", "", lines[heads[k]]), body = if (end > heads[k]) lines[(heads[k] + 1):end] else character())
})
names(sec) <- vapply(sec, `[[`, "", "head")
get <- function(hs) unlist(lapply(seq_along(hs), function(k) { h <- hs[k]
  if (is.null(sec[[h]])) stop("no section: ", h)
  c(if (k > 1) c("", toupper(h)), sec[[h]]$body) }))   # later sections keep their heading
plain <- function(x) {
  x <- x[!grepl("^\\|[-| ]+\\|$", x)]                 # table separator rows
  x <- gsub("\\*\\*([^*]+)\\*\\*", "\\1", x)           # bold
  x <- gsub("(^|[^*])\\*([^* ][^*]*)\\*", "\\1\\2", x) # italics
  x <- gsub("`([^`]+)`", "\\1", x)
  x <- sub("^\\| ?(.*?) ?\\|$", "\\1", x, perl = TRUE)  # table rows as  a | b | c
  paste(trimws(paste(x, collapse = "\n"), "both"), collapse = "")
}
FIELDS <- list(
  "Title" = "Title", "Description" = "Description", "Hypotheses" = "Hypotheses",
  "Study type" = "Study type", "Blinding" = "Blinding", "Study design" = "Study design", "Randomization" = "Randomization",
  "Existing data (choose: Registration prior to analysis of the data) and Explanation of existing data" = "Existing data",
  "Data collection procedures" = "Data collection procedures", "Sample size" = "Sample size",
  "Sample size rationale" = "Sample size rationale", "Stopping rule" = "Stopping rule",
  "Manipulated variables" = "Manipulated variables", "Measured variables" = "Measured variables", "Indices" = "Indices",
  "Statistical models" = c("Statistical models", "Predicted validity matrix", "Specifications fixed in the code"),
  "Transformations" = "Transformations", "Inference criteria" = "Inference criteria", "Data exclusion" = "Data exclusion",
  "Missing data" = "Missing data", "Exploratory analysis" = "Exploratory analysis",
  "Other" = c("6. Other", "References", "Appendix A. The scale, in presented order"))
txt <- c(sprintf("OSF Preregistration template, Study 2b. Commit %s.%s", commit, if (draft) " DRAFT, NOT FOR FILING." else ""),
         "Contributors: add on OSF (Norman Farb, University of Toronto; co-investigators per REB Protocol #00051180).",
         sprintf("Attach: %s (SHA-256 %s) and %s.", zip_name, zip_sha, basename(pdf_path)), "")
for (f in names(FIELDS)) txt <- c(txt, strrep("=", 78), paste0("FIELD: ", f), strrep("=", 78), plain(get(FIELDS[[f]])), "")
osf_path <- file.path(out_dir, paste0("osf_fields", tag, ".txt"))
writeLines(enc2utf8(txt), osf_path, useBytes = TRUE)

## Manifest ---------
outs <- c(zip_path, html_path, pdf_path, osf_path)
man <- c(sprintf("commit %s", commit), sprintf("built %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
         sprintf("status %s", if (draft) "DRAFT" else "final"), if (draft) paste("draft because:", paste(why, collapse = "; ")), "",
         sprintf("%s  %s", vapply(outs, sha256, ""), basename(outs)))
writeLines(man, file.path(out_dir, paste0("manifest", tag, ".txt")))
cat(paste(man, collapse = "\n"), "\n\nWritten to", out_dir, "\n")
