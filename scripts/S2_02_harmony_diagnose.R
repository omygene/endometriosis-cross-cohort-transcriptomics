# ============================================================================
# S2_02_harmony_diagnose.R — isolated measurement of the EXACT registered
# Harmony call at full scale (233,459 cells), using the real checkpoint saved
# by S2_02_integrate.R v2.4 (data/processed/tmp_S2_02/checkpoint_pre_harmony.rds).
# Purpose (R9): settle by DIRECT MEASUREMENT whether harmony 2.0.5 hangs at
# full scale (as the v2.2 run suggested) or is merely slow on single-threaded
# Windows BLAS. Diagnostic only: writes NOTHING except the log; produces no
# project numbers (R1); changes no parameter (the call below is verbatim the
# registered one, including project.dim = FALSE from Addendum S2_02-B.5).
# RUN ONLY IF the main run was confirmed dead (0% CPU) and killed — never in
# parallel with it (RAM contention would corrupt the measurement).
# ============================================================================
set.seed(42)
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}
suppressPackageStartupMessages(library(Seurat))

chk <- "data/processed/tmp_S2_02/checkpoint_pre_harmony.rds"
if (!file.exists(chk))
  stop("checkpoint not found: ", chk,
       " — the main v2.4 run must reach 'checkpoint saved: pre_harmony' first")

merged <- readRDS(chk)
say("checkpoint loaded: ", ncol(merged), " cells | reductions: ",
    paste(Reductions(merged), collapse = ", "))
gc(verbose = FALSE)

say("starting the EXACT registered Harmony call (theta=c(2,2), ",
    "covariates sample+dataset, dims 1:30, project.dim=FALSE)")
t0 <- Sys.time()
merged <- harmony::RunHarmony(merged, group.by.vars = c("sample", "dataset"),
                              theta = c(2, 2), reduction = "pca",
                              dims.use = 1:30, reduction.save = "harmony",
                              project.dim = FALSE)
say("Harmony COMPLETED in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
gc(verbose = FALSE)
say("verdict: harmony 2.0.5 works at full scale — any main-run stall is a ",
    "speed issue, not a hang. DIAGNOSE DONE.")
