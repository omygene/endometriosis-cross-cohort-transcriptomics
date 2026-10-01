# ============================================================================
# S6_00_inspect.R -- GSE213216 structural audit (Stage 6, approved by PI
#   in writing 2026-09-22: "الخطوة القادمة التحقق من خلال GSE213216").
# PURE MEASUREMENT (R9): nothing is assumed about the six shared Seurat
#   objects; every fact used later by S6_01 (assay names, metadata columns,
#   class labels, lineage labels, sample IDs, gene ID style) is measured
#   here and written to results/S6_meta/. S6_01 design mapping is locked
#   ONLY after this audit is reviewed -- same pattern as the S2 audits.
# Charter rules honored:
#   R1 zero fabrication: no numbers are hard-coded; everything is measured.
#   R2 actual files verified, not chat text: fail-loud on any missing input.
#   R4 no results file is ever modified; this script writes ONLY under
#      results/S6_meta/ (+ first-run input registry under scripts/).
#   R5 no new thresholds here (audit only); S6_01 inherits the locked gates.
# Writes: results/S6_meta/S6_00_structure_audit.md
#         results/S6_meta/S6_00_meta_values.csv
#         scripts/S6_input_registry.tsv (first run: loud create; later: verify)
# Usage: Rscript scripts/S6_00_inspect.R    (from project root)
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")   # CFG + find_senmayo_gmt() + read_gmt()
source("scripts/S4_common.R")   # CFG4 + say()

OUT <- "results/S6_meta"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REGISTRY <- "scripts/S6_input_registry.tsv"

## ---- expected objects (names measured in S2 audit: results/S2_qc_GSE213216_audit.csv)
EXPECTED <- c("auxiliary.seurat.shared.rds",
              "endothelial_clusters.shared.rds",
              "EnEpi_cells.shared.rds",
              "EnS_cells.shared.rds",
              "epithelial.annotated.shared.rds",
              "mesenchymal.annotated.shared.rds")

say("== S6_00: GSE213216 structural audit (pure measurement) ==")

## ---- locate the objects under data/ (never assume a path) ------------------
all_rds <- list.files("data", pattern = "\\.rds$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
hits <- lapply(EXPECTED, function(nm) all_rds[basename(all_rds) == nm])
names(hits) <- EXPECTED
missing <- EXPECTED[vapply(hits, length, integer(1)) == 0]
if (length(missing))
  stop("FATAL: expected GSE213216 object(s) not found under data/ (R2):\n  ",
       paste(missing, collapse = "\n  "),
       "\nFound .rds files:\n  ", paste(all_rds, collapse = "\n  "), call. = FALSE)
dup <- EXPECTED[vapply(hits, length, integer(1)) > 1]
if (length(dup))
  stop("FATAL: duplicate copies of object(s) under data/ -- resolve first (R2):\n  ",
       paste(dup, collapse = "\n  "), call. = FALSE)
paths <- vapply(hits, `[[`, character(1), 1)
say("objects located:")
for (nm in EXPECTED) say("  ", nm, " -> ", paths[[nm]])

## ---- input registry (loud create / verify / halt -- same pattern as S5) ----
md5_of <- function(p) unname(tools::md5sum(p))
cur <- data.frame(file = unname(paths), md5 = vapply(paths, md5_of, character(1)),
                  stringsAsFactors = FALSE)
if (file.exists(REGISTRY)) {
  prev <- read.csv(REGISTRY, stringsAsFactors = FALSE, check.names = FALSE)
  m <- merge(prev, cur, by = "file", suffixes = c("_prev", "_cur"), all = TRUE)
  bad <- is.na(m$md5_prev) | is.na(m$md5_cur) | m$md5_prev != m$md5_cur
  if (any(bad))
    stop("FATAL: S6 input registry mismatch (R5). Changed/new files:\n  ",
         paste(m$file[bad], collapse = "\n  "), call. = FALSE)
  say("input registry verified (", nrow(cur), " files)")
} else {
  write.csv(cur, REGISTRY, row.names = FALSE, quote = TRUE)
  say("input registry created loud: ", REGISTRY, " (", nrow(cur), " files)")
}

## ---- locked reference gene sets (read-only) --------------------------------
AXES <- list(
  cytotoxicity         = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GNLY"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"),
  senescence_dormancy  = c("CDKN1A", "CDKN2A", "GADD45A", "SERPINE1"),
  proliferation        = c("MKI67", "TOP2A", "PCNA"),
  stromal_ecm          = c("VIM", "COL1A1", "COL3A1", "FN1", "ACTA2")
)
AXIS_GENES <- unique(unlist(AXES))
if (length(AXIS_GENES) != 24L)
  stop("FATAL: frozen axis registry drifted (expected 24 genes, got ",
       length(AXIS_GENES), ") -- Amendment required (R5)", call. = FALSE)

T1_CSV <- "results/S4_meta/S4_meta_Tier1_core.csv"
if (!file.exists(T1_CSV))
  stop("FATAL: locked Tier-1 core missing: ", T1_CSV, " (R2)", call. = FALSE)
t1 <- read.csv(T1_CSV, stringsAsFactors = FALSE)
if (!"gene" %in% names(t1))
  stop("FATAL: ", T1_CSV, " lacks the locked 'gene' column (R2)", call. = FALSE)
T1_GENES <- t1$gene
say("locked references: 24 axis genes | Tier-1 core = ", length(T1_GENES), " genes")

GMT <- find_senmayo_gmt()
SENMAYO_GENES <- character(0)
if (!is.na(GMT)) {
  gsets <- read_gmt(GMT)
  SENMAYO_GENES <- unique(unlist(lapply(gsets, `[[`, "genes")))
  say("SenMayo GMT reused (never recomputed): ", GMT, " | ", length(SENMAYO_GENES),
      " genes in ", length(gsets), " set(s)")
} else {
  say("NOTE: SenMayo GMT not found under reference/ or data/ -- reported, not fatal at audit stage")
}

## ---- per-object audit -------------------------------------------------------
have_seurat <- requireNamespace("SeuratObject", quietly = TRUE)
say("SeuratObject available: ", have_seurat)

CAND_RX <- paste0("class|type|sample|patient|donor|subject|lineage|cell|",
                  "cluster|site|tissue|group|condit|origin|diagnos|status|",
                  "lesion|endo|arm|batch")

audit_md <- c("# S6_00 -- GSE213216 structural audit (pure measurement)",
              "",
              paste0("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
              paste0("SeuratObject available: ", have_seurat),
              "")
meta_rows <- list()

for (nm in EXPECTED) {
  p <- paths[[nm]]
  say("auditing ", nm, " ...")
  obj <- readRDS(p)
  cls <- class(obj)[1]
  ## dims + assays: Seurat-aware, with plain-attribute fallback (no assumptions)
  if (have_seurat) {
    dm <- dim(obj)
    assays <- tryCatch(SeuratObject::Assays(obj), error = function(e) NA_character_)
    md <- as.data.frame(obj[[]])
    genes <- tryCatch(rownames(obj[["RNA"]]), error = function(e) rownames(obj))
  } else {
    assays_slot <- attr(obj, "assays")
    assays <- if (is.null(assays_slot)) NA_character_ else names(assays_slot)
    md <- attr(obj, "meta.data")
    if (is.null(md)) md <- data.frame(row.names = attr(obj, "cell.names"))
    genes <- if (!is.null(assays_slot) && "RNA" %in% names(assays_slot))
      rownames(attr(assays_slot[["RNA"]], "counts")) else character(0)
    dm <- c(length(genes), nrow(md))
  }
  n_genes <- dm[1]; n_cells <- dm[2]
  say("  class=", cls, " | genes=", n_genes, " | cells=", n_cells,
      " | assays=", paste(assays, collapse = ","))

  ov_axis <- sum(AXIS_GENES %in% genes)
  ov_t1   <- sum(T1_GENES %in% genes)
  ov_sm   <- if (length(SENMAYO_GENES)) sum(SENMAYO_GENES %in% genes) else NA_integer_
  id_style <- if (length(genes) && any(grepl("^ENSG", genes))) "ENSEMBL-like"
              else if (length(genes) && any(grepl("^[A-Z0-9]+$", head(genes, 500)))) "symbol-like"
              else "unknown"

  audit_md <- c(audit_md,
    paste0("## ", nm),
    paste0("- path: `", p, "`"),
    paste0("- class: ", cls, " | genes: ", n_genes, " | cells: ", n_cells),
    paste0("- assays: ", paste(assays, collapse = ", ")),
    paste0("- gene ID style: ", id_style),
    paste0("- axis genes present: ", ov_axis, "/24",
           " | Tier-1 core present: ", ov_t1, "/", length(T1_GENES),
           if (is.na(ov_sm)) " | SenMayo: GMT not found"
           else paste0(" | SenMayo present: ", ov_sm, "/", length(SENMAYO_GENES))),
    paste0("- metadata columns (", ncol(md), "): ",
           paste(names(md), collapse = ", ")),
    "")

  ## candidate classification columns: dump measured values (top 40 by count)
  cand <- names(md)[grepl(CAND_RX, names(md), ignore.case = TRUE)]
  cand <- cand[vapply(md[cand], function(x) is.character(x) || is.factor(x) ||
                        is.logical(x), logical(1))]
  for (cc in cand) {
    tb <- sort(table(md[[cc]], useNA = "ifany"), decreasing = TRUE)
    tb <- head(tb, 40)
    meta_rows[[length(meta_rows) + 1L]] <- data.frame(
      object = nm, column = cc, n_unique = length(unique(md[[cc]])),
      n_cells = n_cells,
      value_counts = paste(sprintf("%s=%d", names(tb), as.integer(tb)),
                           collapse = "; "),
      stringsAsFactors = FALSE)
    audit_md <- c(audit_md,
      paste0("  - `", cc, "` (", length(unique(md[[cc]])), " unique): ",
             paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")))
  }
  audit_md <- c(audit_md, "")
  rm(obj); invisible(gc())
}

if (!length(meta_rows))
  stop("FATAL: no categorical candidate metadata columns measured in any object ",
       "-- cannot design S6_01 contrast mapping without a measurement (R9)",
       call. = FALSE)

meta_df <- do.call(rbind, meta_rows)
write.csv(meta_df, file.path(OUT, "S6_00_meta_values.csv"), row.names = FALSE)

writeLines(c(audit_md,
  "## Next step (locked sequence)",
  "",
  "Send this audit + S6_00_meta_values.csv back to the session. The S6_01",
  "contrast/lineage/sample mapping is locked ONLY from these measured values",
  "(R2/R9), then S6_01_validate.R is finalized for the PI to run (R8)."),
  file.path(OUT, "S6_00_structure_audit.md"))

say("wrote ", file.path(OUT, "S6_00_structure_audit.md"))
say("wrote ", file.path(OUT, "S6_00_meta_values.csv"), " (", nrow(meta_df), " rows)")
say("S6_00 DONE -- send both outputs to the session before S6_01")
