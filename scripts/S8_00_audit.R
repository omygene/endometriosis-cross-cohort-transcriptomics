## ===========================================================================
## S8_00_audit.R -- المرحلة الثامنة: تدقيق LINCS/iLINCS قبل التحقق العكسي
## v1.0 (2026-09-23)
## ---------------------------------------------------------------------------
## السؤال: هل مساحة L1000/iLINCS تسمح باختبار عكسي محكم لإشارة الـ66 جين؟
##   (1) إعادة اشتقاق النواة الصلبة من ملفات S4 المقفلة (R2 -- ليست مكتوبة يدويًا)
##   (2) تغطية الجينات المقفلة في مساحة L1000 (12,328 جينًا / 978 landmark)
##   (3) توافر التوقيعات المطلوبة في iLINCS (KD لـSERPINE1 / estradiol L1000)
##   (4) حكم GO/NO-GO مطبوع آليًا بمعايير مقفلة مسبقًا (في وثيقة التصميم)
##
## الميثاق: R1 بلا اختلاق | R2 تحقق من الملفات الفعلية | R4 قراءة فقط للمقفلة
##          R6 إبلاغ صريح عن أي فشل | R7 مخرجات + manifest | R8 الـPI تشغّل
##
## التشغيل (على جهاز الـPI):
##   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts/S8_00_audit.R
## الحزم المطلوبة: httr, jsonlite (تُثبَّت تلقائيًا إن غابت)
## ===========================================================================

pkgs <- c("httr", "jsonlite")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("installing missing package: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({library(httr); library(jsonlite)})

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

## ---------------------------------------------------------------------------
## 0) الإعدادات المقفلة
## ---------------------------------------------------------------------------
API   <- "https://www.ilincs.org/api"
OUT   <- "results/S8_perturb"
CACHE <- file.path(OUT, "cache_geneinfo")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

GENE_INFO_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/",
  "GSE92742_Broad_LINCS_gene_info.txt.gz")
GENE_INFO_LOCAL <- file.path(CACHE, "GSE92742_Broad_LINCS_gene_info.txt.gz")

## سجل md5 للمدخلات المقفلة (مقيس من أرشيف نتائج الـPI المقفل)
REG <- c(
  "results/S4_meta/S4_meta_Tier1_core.csv"        = "3afbc9d85bec0b4409137cb73cc854bb",
  "results/S4_meta/S4_sens_LOCO_B-G120103_all.csv" = "c5d6f8fc305da678006a71ba5900cb68",
  "results/S4_meta/S4_sens_LOCO_B-G25628_all.csv"  = "b4b62c1f39d5657b4cb0635b69669c86",
  "results/S4_meta/S4_sens_LOCO_B-G51981_all.csv"  = "326f4ad7f7df0631321819d7c7ed1f2e",
  "results/S4_meta/S4_sens_LOCO_B-G6364_all.csv"   = "488f98e302ad3fe88719d7633211c3d0",
  "results/S4_meta/S4_sens_LOCO_B-G7305_all.csv"   = "42c60e220e917cd5e1044ae0be8b020c",
  "results/S4_meta/S4_sens_LOCO_B-G7307_all.csv"   = "14ec9b59df4dfed07721b3e987a10cd1",
  "results/S4_meta/S4_sens_sansG7307_all.csv"      = "14ec9b59df4dfed07721b3e987a10cd1")

## المعايير المقفلة مسبقًا (وثيقة تصميم S8، قسم 6 -- GO/NO-GO)
CRIT <- list(min_kd_lines   = 5L,   ## G1: خطوط خلية مستقلة لـSERPINE1-KD
             min_t1_landmark = 40L, ## G2: جينات Tier-1 قابلة للقياس في L1000 landmark
             min_estradiol_lincs = 3L) ## G3: تواقيع estradiol في مكتبة L1000 الحديثة

md5_file <- function(p) {
  unname(tools::md5sum(p))
}

## ---------------------------------------------------------------------------
## 1) التحقق من سجل المدخلات (R2)
## ---------------------------------------------------------------------------
say("S8_00: verifying locked input registry")
for (f in names(REG)) {
  if (!file.exists(f))
    stop("FATAL: locked input missing: ", f, call. = FALSE)
  m <- md5_file(f)
  if (m != REG[[f]])
    stop("FATAL: md5 mismatch for ", f, "\n  expected ", REG[[f]], "\n  got      ", m,
         call. = FALSE)
}
say("  registry OK (", length(REG), " locked files )")

## ---------------------------------------------------------------------------
## 2) إعادة اشتقاق النواة الصلبة من الملفات المقفلة (ليست ثابتة مكتوبة)
##    القاعدة المقفلة: tier == Tier1 في الأساس + كل سيناريوهات LOCO الستة + sans-G7307
## ---------------------------------------------------------------------------
say("S8_00: re-deriving the robust core from locked S4 files")
f_base <- "results/S4_meta/S4_meta_Tier1_core.csv"
base <- read.csv(f_base, stringsAsFactors = FALSE)
stopifnot(all(c("gene", "tier") %in% names(base)))
base_t1 <- base$gene[base$tier == "Tier1"]

scen_files <- c(
  "results/S4_meta/S4_sens_LOCO_B-G120103_all.csv",
  "results/S4_meta/S4_sens_LOCO_B-G25628_all.csv",
  "results/S4_meta/S4_sens_LOCO_B-G51981_all.csv",
  "results/S4_meta/S4_sens_LOCO_B-G6364_all.csv",
  "results/S4_meta/S4_sens_LOCO_B-G7305_all.csv",
  "results/S4_meta/S4_sens_LOCO_B-G7307_all.csv",
  "results/S4_meta/S4_sens_sansG7307_all.csv")

core <- base_t1
for (f in scen_files) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  stopifnot(all(c("gene", "tier") %in% names(d)))
  core <- intersect(core, d$gene[d$tier == "Tier1"])
}
core_M <- base$M[match(core, base$gene)]
n_up   <- sum(core_M > 0); n_dn <- sum(core_M < 0)

## بوابة اتساق داخلية: النتيجة المعروفة من الملفات المقفلة = 66 (14 فوق / 52 تحت)
if (length(core) != 66L || n_up != 14L || n_dn != 52L)
  stop("FATAL: re-derived robust core = ", length(core), " (", n_up, " up / ", n_dn,
       " down); locked value is 66 (14/52). Input files may have changed (R2).",
       call. = FALSE)
say("  robust core re-derived: ", length(core), " genes (", n_up, " up / ", n_dn,
    " down) -- matches locked value")

## ---------------------------------------------------------------------------
## 3) مساحة جينات L1000: تنزيل ملف gene_info (صغير ~200KB) وتغطية الجينات المقفلة
## ---------------------------------------------------------------------------
say("S8_00: downloading L1000 gene info (GSE92742, ~200KB)")
if (!file.exists(GENE_INFO_LOCAL)) {
  ok <- tryCatch({
    download.file(GENE_INFO_URL, GENE_INFO_LOCAL, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(GENE_INFO_LOCAL))
    stop("FATAL: could not download L1000 gene info from GEO ftp (network?).",
         " Manual fix: place GSE92742_Broad_LINCS_gene_info.txt.gz in ", CACHE, call. = FALSE)
}
gi <- read.delim(gzfile(GENE_INFO_LOCAL), stringsAsFactors = FALSE)
stopifnot(all(c("pr_gene_id", "pr_gene_symbol", "pr_is_lm") %in% names(gi)))
n_lm  <- sum(gi$pr_is_lm == 1)
say("  L1000 space: ", nrow(gi), " genes | ", n_lm, " landmarks")

lm_set    <- gi$pr_gene_symbol[gi$pr_is_lm == 1]
lincs_set <- gi$pr_gene_symbol

core_in_lincs <- core %in% lincs_set
core_in_lm    <- core %in% lm_set
t1 <- base[base$tier == "Tier1", ]
t1_up_lm <- sum(t1$M > 0 & t1$gene %in% lm_set)
t1_dn_lm <- sum(t1$M < 0 & t1$gene %in% lm_set)
axes <- list(
  cytotoxicity         = c("CD8A","CD8B","GZMB","PRF1","NKG7","GNLY"),
  antigen_presentation = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2"),
  senescence_dormancy  = c("CDKN1A","CDKN2A","GADD45A","SERPINE1"),
  proliferation        = c("MKI67","TOP2A","PCNA"),
  stromal_ecm          = c("VIM","COL1A1","COL3A1","FN1","ACTA2"))
axis_cov <- data.frame(
  axis  = names(axes),
  genes = sapply(axes, length),
  in_lincs = sapply(axes, function(gs) sum(gs %in% lincs_set)),
  in_landmark = sapply(axes, function(gs) sum(gs %in% lm_set)),
  stringsAsFactors = FALSE)
serp_in_lm <- "SERPINE1" %in% lm_set

cov_tab <- data.frame(
  measure = c("robust_core_in_L1000", "robust_core_landmark",
              "tier1_up_landmark", "tier1_down_landmark",
              "tier1_total_landmark", "SERPINE1_landmark"),
  n = c(sum(core_in_lincs), sum(core_in_lm),
        t1_up_lm, t1_dn_lm, t1_up_lm + t1_dn_lm, as.integer(serp_in_lm)),
  stringsAsFactors = FALSE)
write.csv(cov_tab, file.path(OUT, "S8_00_coverage.csv"), row.names = FALSE)
write.csv(axis_cov,  file.path(OUT, "S8_00_axis_coverage.csv"), row.names = FALSE)
write.csv(data.frame(gene = core, M = base$M[match(core, base$gene)],
                     in_L1000 = core_in_lincs, is_landmark = core_in_lm),
          file.path(OUT, "S8_00_core_gene_coverage.csv"), row.names = FALSE)
say("  coverage: core ", sum(core_in_lincs), "/66 in L1000 | ", sum(core_in_lm),
    "/66 landmark | Tier-1 landmark ", t1_up_lm + t1_dn_lm, " (", t1_up_lm, " up / ",
    t1_dn_lm, " down)")

## ---------------------------------------------------------------------------
## 4) توافر التوقيعات في iLINCS (استعلامات حقيقية، فشل صريح عند سقوط الشبكة)
## ---------------------------------------------------------------------------
say("S8_00: querying iLINCS API for perturbation availability")
ilincs_get <- function(path, filter = NULL) {
  q <- list()
  if (!is.null(filter)) q$filter <- filter
  r <- GET(paste0(API, path), query = q,
           user_agent("endometriosis-S8-audit/1.0"), timeout(120))
  if (http_error(r))
    stop("FATAL: iLINCS API error ", status_code(r), " for ", path, call. = FALSE)
  txt <- content(r, as = "text", encoding = "UTF-8")
  fromJSON(txt, simplifyVector = TRUE)
}
jl <- function(x) toJSON(x, auto_unbox = TRUE)

avail <- list()

## 4a) SERPINE1 knockdown (كل الأنواع، توثيق pert_type)
d <- ilincs_get("/SignatureMeta", jl(list(where = list(compound = "SERPINE1"))))
if (length(d) == 0L) d <- data.frame()
kd <- data.frame(
  signatureid = if (length(d)) d$signatureid else character(),
  pert_type   = if (length(d)) ifelse(is.null(d$pert_type), NA, d$pert_type) else character(),
  cellline    = if (length(d)) d$cellline else character(),
  libraryid   = if (length(d)) d$libraryid else character(),
  time        = if (length(d)) ifelse(is.null(d$time), NA, d$time) else character(),
  stringsAsFactors = FALSE)
kd_cgs <- kd[kd$pert_type == "trt_sh.cgs" & kd$libraryid == "LIB_6", ]
avail$SERPINE1_KD <- kd
say("  SERPINE1: ", nrow(kd), " signatures | ", nrow(kd_cgs), " CGS-KD across ",
    length(unique(kd_cgs$cellline)), " cell lines")

## 4b) estradiol: مكتبة L1000 الحديثة (LIB_5) + مكتبة CMAP القديمة (LIB_2)
f_est <- '{"where":{"and":[{"compound":"estradiol"},{"lincspertid":{"neq":null}}]},"limit":10000}'
est_l5 <- ilincs_get("/SignatureMeta", f_est)
if (length(est_l5) == 0L) est_l5 <- data.frame()
est_l5 <- data.frame(
  signatureid = if (nrow(est_l5)) est_l5$signatureid else character(),
  cellline    = if (nrow(est_l5)) est_l5$cellline else character(),
  time        = if (nrow(est_l5)) ifelse(is.null(est_l5$time), NA, est_l5$time) else character(),
  libraryid   = if (nrow(est_l5)) est_l5$libraryid else character(),
  stringsAsFactors = FALSE)
est_l5 <- est_l5[est_l5$libraryid == "LIB_5", ]
avail$estradiol_L1000 <- est_l5
say("  estradiol LIB_5 (L1000): ", nrow(est_l5), " signatures across ",
    length(unique(est_l5$cellline)), " cell lines")

## 4c) مرشحون سُجِّلوا مسبقًا ومتوقع غيابهم (توثيق صريح بدل الصمت R6)
cands <- c("tiplaxtinin", "TM5441", "PAI-039", "forskolin")
cand_rows <- lapply(cands, function(nm) {
  d <- tryCatch(ilincs_get("/SignatureMeta",
      jl(list(where = list(compound = nm), limit = 1000))),
    error = function(e) NULL)
  n <- if (is.null(d)) NA_integer_ else length(d)
  data.frame(candidate = nm, n_signatures = n, stringsAsFactors = FALSE)
})
avail$missing_candidates <- do.call(rbind, cand_rows)
write.csv(avail$SERPINE1_KD,      file.path(OUT, "S8_00_avail_serpine1_kd.csv"), row.names = FALSE)
write.csv(avail$estradiol_L1000,  file.path(OUT, "S8_00_avail_estradiol_L1000.csv"), row.names = FALSE)
write.csv(avail$missing_candidates, file.path(OUT, "S8_00_avail_missing_candidates.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5) حكم GO/NO-GO (معايير مقفلة من وثيقة التصميم)
## ---------------------------------------------------------------------------
G1 <- length(unique(kd_cgs$cellline)) >= CRIT$min_kd_lines
G2 <- (t1_up_lm + t1_dn_lm)           >= CRIT$min_t1_landmark
G3 <- nrow(est_l5)                    >= CRIT$min_estradiol_lincs
GO <- G1 && G2 && G3

if (GO) {
  verdict <- "GO -- S8_01 reverse-confirmation may proceed (all locked criteria met)"
} else {
  verdict <- paste0("NO-GO -- S8_01 blocked; failed criteria: ",
                    paste(c("G1_kd_lines"[!G1], "G2_t1_landmark"[!G2],
                            "G3_estradiol"[!G3]), collapse = ", "))
}

summary_md <- c(
  "# S8_00 -- LINCS/iLINCS audit summary (machine-printed)",
  "",
  paste0("Date: ", Sys.time()),
  "",
  "## Robust core re-derived from locked files",
  paste0("- robust core: **", length(core), "** genes (", n_up, " up / ", n_dn,
         " down) -- Tier-1 retained under all 6 LOCO + sans-G7307"),
  "- internal consistency gate vs locked value (66 = 14/52): PASS",
  "",
  "## L1000 coverage (GSE92742 gene info)",
  paste0("- L1000 space: ", nrow(gi), " genes (", n_lm, " landmarks)"),
  paste0("- robust core in L1000: ", sum(core_in_lincs), "/66; landmark: ",
         sum(core_in_lm), "/66 (", paste(core[core_in_lm], collapse = ", "), ")"),
  paste0("- Tier-1 in landmark space: ", t1_up_lm + t1_dn_lm,
         " (", t1_up_lm, " up / ", t1_dn_lm, " down) [criterion G2 >= ",
         CRIT$min_t1_landmark, ": ", if (G2) "PASS" else "FAIL", "]"),
  paste0("- SERPINE1 in landmark space: ", serp_in_lm),
  "",
  "## iLINCS availability",
  paste0("- SERPINE1-KD (CGS, LIB_6): ", nrow(kd_cgs), " signatures / ",
         length(unique(kd_cgs$cellline)), " cell lines [criterion G1 >= ",
         CRIT$min_kd_lines, ": ", if (G1) "PASS" else "FAIL", "]"),
  paste0("- estradiol (LIB_5, L1000): ", nrow(est_l5), " signatures [criterion G3 >= ",
         CRIT$min_estradiol_lincs, ": ", if (G3) "PASS" else "FAIL", "]"),
  "- pre-registered absent candidates (documented, not silently skipped):",
  paste0("  - ", avail$missing_candidates$candidate, ": n=",
         avail$missing_candidates$n_signatures, collapse = "\n"),
  "",
  "## VERDICT",
  paste0("**", verdict, "**"),
  "",
  "*Every number in this file was measured at run time from the locked inputs and",
  " the live iLINCS API. Nothing was copied from chat text (R1/R2).*")
writeLines(summary_md, file.path(OUT, "S8_00_audit_summary.md"))

## ---------------------------------------------------------------------------
## 6) manifest
## ---------------------------------------------------------------------------
outs <- c("S8_00_audit_summary.md", "S8_00_coverage.csv",
          "S8_00_axis_coverage.csv", "S8_00_core_gene_coverage.csv",
          "S8_00_avail_serpine1_kd.csv", "S8_00_avail_estradiol_L1000.csv",
          "S8_00_avail_missing_candidates.csv")
man <- data.frame(
  file = outs,
  md5 = vapply(file.path(OUT, outs), md5_file, character(1)),
  stringsAsFactors = FALSE)
write.table(man, file.path(OUT, "S8_00_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

say("S8_00 DONE")
say("  verdict: ", verdict)
say("  outputs in ", OUT, "/ (+ manifest)")
