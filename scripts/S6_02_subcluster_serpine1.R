#!/usr/bin/env Rscript
# ============================================================================
# S6_02_subcluster_serpine1.R -- EXPLORATORY dilution-hypothesis test
#   (S6b): SERPINE1 + Tier-1 signature across the 18 endothelial subclusters
#   of GSE213216.
#
# Opened by PI written directive 2026-09-23 ("whatever is scientifically
# defensible and serves the paper -- execute it"). SCOPE IS EXPLORATORY:
# any signal here is hypothesis-generating ONLY and can never upgrade the
# lineage-level S6_01 verdict (P1: NOT CONFIRMED stands forever). This is
# printed in every output (R1/R6).
#
# Rationale (documented before seeing any result): the internal niche finding
# (GSE214411) was SUBCLUSTER-specific (Endothelial A/B out of 5), while S6_01
# tested the whole endothelial lineage. Subcluster aggregation can dilute or
# reverse a subcluster-restricted signal. S6b tests exactly that, on the
# dedicated endothelial object (18 subclusters, 23,226 cells -- measured in
# S6_00 to be the identical cell set as the auxiliary endothelial lineage).
#
# Charter: R1 zero fabrication | R2 registry md5 re-verified | R4 locked
# inputs read-only | R5 BH-FDR<0.05 within THIS exploratory family is the
# reporting gate | R6 negatives/NOT_TESTABLE reported | R8 PI runs it.
# ============================================================================

say  <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
fail <- function(...) { message(paste0(...)); quit(save = "no", status = 1) }
need <- function(p) if (!requireNamespace(p, quietly = TRUE))
  fail("FATAL: package '", p, "' not available -- install it and re-run (R8)")
need("Matrix"); need("SeuratObject")

ROOT <- getwd()
RES  <- file.path(ROOT, "results", "S6_meta")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

MIN_CELLS <- 20L
MIN_ARM   <- 4L

say("== S6_02 (S6b): EXPLORATORY endothelial-subcluster dilution test ==")
say("SCOPE: hypothesis-generating only -- cannot upgrade S6_01 P1 verdict")

# --- 1. registry re-verification (R2/R4) -------------------------------------
reg_path <- file.path(ROOT, "scripts", "S6_input_registry.tsv")
if (!file.exists(reg_path))
  fail("FATAL: scripts/S6_input_registry.tsv missing -- run S6_00 first (R2)")
reg <- readLines(reg_path, warn = FALSE)
ln  <- reg[grepl("endothelial_clusters.shared.rds", reg, fixed = TRUE)]
if (length(ln) != 1L)
  fail("FATAL: registry entry for endothelial_clusters.shared.rds not unique/missing")
reg_md5 <- tolower(regmatches(ln, regexpr("[0-9a-fA-F]{32}", ln)))
rds <- list.files(file.path(ROOT, "data"),
                  pattern = "^endothelial_clusters\\.shared\\.rds$",
                  recursive = TRUE, full.names = TRUE)
if (length(rds) != 1L)
  fail("FATAL: expected exactly 1 endothelial_clusters.shared.rds, found ",
       length(rds))
fresh_md5 <- tolower(unname(tools::md5sum(rds)))
if (fresh_md5 != reg_md5)
  fail("FATAL: endothelial object CHANGED since S6_00 audit -- halt (R2/R4)")
say("registry verified: endothelial object unchanged | md5 ", fresh_md5)

# --- 2. locked Tier-1 core (read-only) ---------------------------------------
t1_path <- file.path(ROOT, "results", "S4_meta", "S4_meta_Tier1_core.csv")
if (!file.exists(t1_path)) fail("FATAL: locked Tier-1 core missing: ", t1_path)
t1 <- read.csv(t1_path, stringsAsFactors = FALSE)
if (!all(c("gene", "M") %in% names(t1)))
  fail("FATAL: Tier-1 core lacks gene/M columns (R4)")
UP   <- t1$gene[t1$M > 0]; DOWN <- t1$gene[t1$M < 0]
if (length(UP) != 121L || length(DOWN) != 364L)
  fail("FATAL: Tier-1 split ", length(UP), "/", length(DOWN),
       " -- locked registry says 121/364 (R4)")
say("locked Tier-1 core verified: 121 up / 364 down")

# --- 3. load + locked mapping (from S6_00 measurements) -----------------------
say("loading ", rds, " (23,226 cells) ...")
obj <- readRDS(rds)
if (!is(obj, "Seurat")) fail("FATAL: object class is ", class(obj)[1])
meta <- obj[[]]
for (cl in c("Patient.No.", "Major.Class", "seurat_clusters", "DF.classifications"))
  if (!cl %in% names(meta))
    fail("FATAL: expected column '", cl, "' missing -- audit mapping violated (R9)")
cts <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
                error = function(e)
                  tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts"),
                           error = function(e2) NULL))
if (is.null(cts)) fail("FATAL: RNA counts layer unavailable (R9)")
rm(obj); invisible(gc(verbose = FALSE))

common <- intersect(colnames(cts), rownames(meta))
if (!length(common)) fail("FATAL: zero cell overlap counts/metadata (R9)")
meta <- meta[common, , drop = FALSE]
cts  <- cts[, common, drop = FALSE]
g <- rownames(cts)
if (!("SERPINE1" %in% g)) fail("FATAL: SERPINE1 absent from object (R9)")
if (sum(c(UP, DOWN) %in% g) != 463L)
  fail("FATAL: Tier-1 coverage ", sum(c(UP, DOWN) %in% g),
       "/485 -- audit measured 463 (R9)")

CLASS_MAP <- c("Endometriosis"                = "Peritoneal lesion",
               "Extra-ovarian endometriosis"  = "Peritoneal lesion",
               "Endometrioma"                 = "Endometrioma",
               "Eutopic Endometrium"          = "Eutopic Endometrium",
               "Unaffected ovary"             = "Unaffected ovary",
               "No endometriosis detected"    = "Disease-free peritoneum")
keep <- meta$DF.classifications == "Singlet" &
        !is.na(meta$Patient.No.) & !is.na(meta$seurat_clusters) &
        !is.na(meta$Major.Class) & meta$Major.Class %in% names(CLASS_MAP)
meta <- meta[keep, , drop = FALSE]
cts  <- cts[, keep, drop = FALSE]
meta$tissue <- unname(CLASS_MAP[as.character(meta$Major.Class)])
if (any(is.na(meta$tissue))) fail("FATAL: unmapped Major.Class value (R9)")
say("cells kept (singlets): ", nrow(meta))

# --- 4. pseudobulk per patient x endothelial subcluster -----------------------
meta$grp <- paste(meta$Patient.No., meta$seurat_clusters, meta$tissue, sep = "||")
G <- Matrix::sparse.model.matrix(~ 0 + grp, data = data.frame(grp = meta$grp))
colnames(G) <- sub("^grp", "", colnames(G))
nc <- Matrix::colSums(G)
G  <- G[, nc >= MIN_CELLS, drop = FALSE]
if (!ncol(G)) fail("FATAL: zero pseudobulk groups survive (R9)")
pb   <- cts %*% G
rm(cts, G); invisible(gc(verbose = FALSE))
lib  <- Matrix::colSums(pb)
lcpm <- log1p(t(t(as.matrix(pb)) / lib) * 1e6)
rm(pb); invisible(gc(verbose = FALSE))
parts <- do.call(rbind, strsplit(colnames(lcpm), "||", fixed = TRUE))
pm <- data.frame(patient = parts[, 1], cluster = parts[, 2],
                 tissue = parts[, 3], stringsAsFactors = FALSE)
say("pseudobulks (>=", MIN_CELLS, " cells): ", nrow(pm),
    " | subclusters present: ", length(unique(pm$cluster)))

# --- 5. endpoints per cluster --------------------------------------------------
CASE_T <- c("Peritoneal lesion", "Endometrioma"); CTRL_T <- "Eutopic Endometrium"
run_test <- function(x_case, x_ctrl) {
  if (length(x_case) < MIN_ARM || length(x_ctrl) < MIN_ARM)
    return(list(testable = FALSE, n_case = length(x_case), n_ctrl = length(x_ctrl),
                p = NA_real_, delta = NA_real_))
  wt <- suppressWarnings(wilcox.test(x_case, x_ctrl, exact = FALSE))
  list(testable = TRUE, n_case = length(x_case), n_ctrl = length(x_ctrl),
       p = wt$p.value, delta = median(x_case) - median(x_ctrl))
}
out <- list(); k <- 0L
for (cl in sort(unique(pm$cluster))) {
  idx <- which(pm$cluster == cl)
  z   <- t(scale(t(lcpm[, idx, drop = FALSE]))); z[is.na(z)] <- 0
  sig <- colMeans(z[UP[UP %in% rownames(z)], , drop = FALSE]) -
         colMeans(z[DOWN[DOWN %in% rownames(z)], , drop = FALSE])
  serp <- lcpm["SERPINE1", idx]
  d  <- pm[idx, , drop = FALSE]
  for (ep in c("SERPINE1_logCPM", "tier1_signature")) {
    v <- if (ep == "SERPINE1_logCPM") serp else sig
    r <- run_test(v[d$tissue %in% CASE_T], v[d$tissue %in% CTRL_T])
    k <- k + 1L
    out[[k]] <- data.frame(cluster = cl, endpoint = ep,
      n_case = r$n_case, n_ctrl = r$n_ctrl, delta_median = r$delta, p = r$p,
      status = ifelse(r$testable, "TESTED", "NOT_TESTABLE(<4/arm)"),
      stringsAsFactors = FALSE)
  }
}
E <- do.call(rbind, out)
if (!any(E$status == "TESTED"))
  fail("FATAL: no testable subcluster at all -- stage meaningless (R9)")
E$FDR <- NA_real_
ok <- E$status == "TESTED"
E$FDR[ok] <- p.adjust(E$p[ok], method = "BH")   # exploratory family only
say("tests: ", nrow(E), " | tested ", sum(ok), " | NOT_TESTABLE ", sum(!ok))

# --- 6. machine-printed exploratory verdict -----------------------------------
sig_hits <- E[ok & !is.na(E$FDR) & E$FDR < 0.05, ]
has_serp <- nrow(sig_hits) > 0L &&
            any(sig_hits$endpoint == "SERPINE1_logCPM" & sig_hits$delta_median > 0)
if (has_serp) {
  verdict <- paste0("EXPLORATORY SIGNAL PRESENT (hypothesis-generating ONLY -- ",
                    "P1 verdict of S6_01 is unchanged: NOT CONFIRMED at lineage level)")
} else {
  verdict <- paste0("NO exploratory subcluster signal for SERPINE1 at FDR<0.05 -- ",
                    "the dilution hypothesis is NOT supported; P1 NOT CONFIRMED stands")
}

# --- 7. outputs (new files only -- R4) ------------------------------------------
f_csv <- file.path(RES, "S6_02_exploratory_endothelial_subclusters.csv")
f_md  <- file.path(RES, "S6_02_exploratory_summary.md")
f_man <- file.path(RES, "S6_02_manifest.tsv")
write.csv(E, f_csv, row.names = FALSE)
md <- c(
  "# S6_02 (S6b) -- EXPLORATORY endothelial-subcluster dilution test",
  "",
  "**SCOPE: hypothesis-generating ONLY. Nothing here upgrades the locked",
  "S6_01 verdict (P1: NOT CONFIRMED). Opened by PI written directive 2026-09-23.**",
  "",
  paste0("Run: ", Sys.time(), " | input md5 verified: ", fresh_md5),
  paste0("Pseudobulks: ", nrow(pm), " | subclusters tested: ",
         sum(E$endpoint == "SERPINE1_logCPM" & E$status == "TESTED"),
         " | exploratory BH family over ", sum(ok), " tests"),
  "",
  paste0("## Verdict (machine-printed, locked rule): ", verdict),
  "",
  "## FDR < 0.05 rows (exploratory family)",
  "",
  if (nrow(sig_hits)) paste0("- ", sig_hits$endpoint, " | cluster ",
                             sig_hits$cluster, " | delta=",
                             signif(sig_hits$delta_median, 3),
                             " | FDR=", signif(sig_hits$FDR, 3),
                             " | n=", sig_hits$n_case, "v", sig_hits$n_ctrl)
  else "- none",
  "",
  "All rows incl. negatives/NOT_TESTABLE: S6_02_exploratory_endothelial_subclusters.csv",
  "",
  "S6_02 DONE -- send both outputs back to the session."
)
writeLines(md, f_md)
man <- data.frame(file = basename(c(f_csv, f_md)),
                  md5  = unname(tools::md5sum(c(f_csv, f_md))))
write.table(man, f_man, sep = "\t", row.names = FALSE, quote = FALSE)
say("wrote ", f_csv); say("wrote ", f_md); say("wrote ", f_man)
say("verdict: ", verdict)
say("S6_02 DONE -- send S6_02_exploratory_summary.md + the csv back to the session")
