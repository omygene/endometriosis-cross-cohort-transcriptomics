## ===========================================================================
## RUN_journal_plates.R -- one-shot runner (v1.0, 2026-09-30)
##   Step 1: re-render FigS1 with the corrected y-label (full certified
##           re-derivation from the GEO cache, exactly as the original).
##   Step 2: assemble ALL journal plates into
##           "FIGURES PLATES FOR JOURNAL SUBMISSION".
## Usage (project root):
##   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts/RUN_journal_plates.R
## ===========================================================================
set.seed(42)

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
say("== STEP 1/2: FigS1 exploratory decidualization (corrected y-label) ==")
source("scripts/make_figS1_exploratory_drugs_FIXED.R")

say("== STEP 2/2: journal plates ==")
source("scripts/make_journal_plates.R")

say("== ALL DONE -- see folder 'FIGURES PLATES FOR JOURNAL SUBMISSION' ==")
