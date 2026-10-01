## ===========================================================================
## S8_01_reverse_confirm.R -- التحقق العكسي: هل إسكات SERPINE1 أو تحفيز
##   الإستروجين يعكسان إشارة الآفة المقفلة؟ (مساحة L1000 landmark عبر iLINCS)
## v1.0 (2026-09-23)
## ---------------------------------------------------------------------------
## يعمل فقط بعد S8_00 بحكم GO. الفرضيات مقفلة مسبقًا (وثيقة تصميم S8):
##   H1/P1a: تواقيع SERPINE1-KD معكوسة الاتجاه مع Tier-1 signature للآفة (r<0)
##   H2/P1b: SERPINE1-KD يرفع جينات proliferation (TOP2A/PCNA)
##   H3/P2a: تواقيع estradiol معكوسة الاتجاه مع Tier-1 signature للآفة (r<0)
##   H4/P2b: estradiol يرفع جينات proliferation (TOP2A/PCNA)
## بوابة الادعاء الوحيدة: BH-FDR < 0.05 عبر العائلة الأربعة + اتجاه مطابق للتسجيل المسبق.
##
## الميثاق: R1/R2/R4/R6/R7/R8 كما في S8_00.
## التشغيل:
##   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts/S8_01_reverse_confirm.R
## ===========================================================================

pkgs <- c("httr", "jsonlite")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({library(httr); library(jsonlite)})

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

API <- "https://www.ilincs.org/api"
OUT <- "results/S8_perturb"
CACHE <- file.path(OUT, "cache_sig")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

## حد أدنى لجينات Tier-1 landmark الواجب حضورها داخل كل توقيع محمّل
MIN_OVERLAP <- 40L
## حد أدنى لتواقيع estradiol داخل خط خلية واحد لقبول تجميعه
MIN_PER_LINE <- 3L
POLIF_GENES <- c("TOP2A", "PCNA")

## ---------------------------------------------------------------------------
## 0) سجل المدخلات المقفلة (نفس سجل S8_00)
## ---------------------------------------------------------------------------
REG <- c(
  "results/S4_meta/S4_meta_Tier1_core.csv"        = "3afbc9d85bec0b4409137cb73cc854bb",
  "results/S4_meta/S4_sens_LOCO_B-G120103_all.csv" = "c5d6f8fc305da678006a71ba5900cb68",
  "results/S4_meta/S4_sens_LOCO_B-G25628_all.csv"  = "b4b62c1f39d5657b4cb0635b69669c86",
  "results/S4_meta/S4_sens_LOCO_B-G51981_all.csv"  = "326f4ad7f7df0631321819d7c7ed1f2e",
  "results/S4_meta/S4_sens_LOCO_B-G6364_all.csv"   = "488f98e302ad3fe88719d7633211c3d0",
  "results/S4_meta/S4_sens_LOCO_B-G7305_all.csv"   = "42c60e220e917cd5e1044ae0be8b020c",
  "results/S4_meta/S4_sens_LOCO_B-G7307_all.csv"   = "14ec9b59df4dfed07721b3e987a10cd1",
  "results/S4_meta/S4_sens_sansG7307_all.csv"      = "14ec9b59df4dfed07721b3e987a10cd1")

say("S8_01: verifying locked input registry")
for (f in names(REG)) {
  if (!file.exists(f)) stop("FATAL: locked input missing: ", f, call. = FALSE)
  m <- unname(tools::md5sum(f))
  if (m != REG[[f]])
    stop("FATAL: md5 mismatch for ", f, call. = FALSE)
}
say("  registry OK")

## ---------------------------------------------------------------------------
## 1) مخرجات S8_00 (يجب أن تكون GO) + مساحة الجينات
## ---------------------------------------------------------------------------
man <- read.delim(file.path(OUT, "S8_00_manifest.tsv"), stringsAsFactors = FALSE)
summ <- readLines(file.path(OUT, "S8_00_audit_summary.md"))
if (!any(grepl("GO --", summ)))
  stop("FATAL: S8_00 verdict is not GO. S8_01 must not run (locked rule).",
       call. = FALSE)
say("  S8_00 verdict: GO -- proceeding")

base <- read.csv("results/S4_meta/S4_meta_Tier1_core.csv", stringsAsFactors = FALSE)
t1 <- base[base$tier == "Tier1", ]

gi <- read.delim(gzfile(file.path(OUT, "cache_geneinfo",
                                  "GSE92742_Broad_LINCS_gene_info.txt.gz")),
                 stringsAsFactors = FALSE)
lm_set <- gi$pr_gene_symbol[gi$pr_is_lm == 1]
t1_lm <- t1[t1$gene %in% lm_set, ]
stopifnot(nrow(t1_lm) >= 40L)
core_cov <- read.csv(file.path(OUT, "S8_00_core_gene_coverage.csv"),
                     stringsAsFactors = FALSE)
core8 <- core_cov$gene[core_cov$is_landmark == TRUE]
say("  Tier-1 landmark genes: ", nrow(t1_lm), " | core landmark genes: ",
    length(core8))

## ---------------------------------------------------------------------------
## 2) تحميل توقيع واحد من iLINCS (مع cache على القرص؛ فشل صريح)
## ---------------------------------------------------------------------------
ilincs_post_json <- function(url, body_list, retries = 3L) {
  for (i in seq_len(retries)) {
    r <- tryCatch(POST(url, body = body_list, encode = "form",
                       user_agent("endometriosis-S8/1.0"), timeout(180)),
                  error = function(e) NULL)
    if (!is.null(r) && !http_error(r)) {
      txt <- content(r, as = "text", encoding = "UTF-8")
      out <- tryCatch(fromJSON(txt, simplifyVector = TRUE), error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    say("    retry ", i, " for ", url)
    Sys.sleep(2 * i)
  }
  stop("FATAL: iLINCS call failed after ", retries, " retries: ", url, call. = FALSE)
}

fetch_signature <- function(sig_id) {
  cf <- file.path(CACHE, paste0(sig_id, ".json"))
  if (file.exists(cf)) return(fromJSON(cf, simplifyVector = TRUE)$data$signature)
  out <- ilincs_post_json(paste0(API, "/ilincsR/downloadSignature"),
                          list(sigID = sig_id, noOfTopGenes = "Inf"))
  sig <- out$data$signature
  if (is.null(sig) || nrow(sig) == 0L)
    stop("FATAL: empty signature returned for ", sig_id, call. = FALSE)
  writeLines(toJSON(out), cf)
  Sys.sleep(0.3)
  sig
}

## مقاييس توقيع واحد مقابل الآفة المقفلة
score_signature <- function(sig, label, arm) {
  v <- sig$Value_LogDiffExp
  names(v) <- sig$Name_GeneSymbol
  ov  <- intersect(t1_lm$gene, names(v))
  if (length(ov) < MIN_OVERLAP) {
    say("    WARNING: ", label, " overlap = ", length(ov), " (< ", MIN_OVERLAP,
        ") -- scored NA (reported, not dropped silently)")
    return(data.frame(arm = arm, signatureid = label, n_overlap = length(ov),
                      r_tier1 = NA_real_, core8_mean = NA_real_,
                      prolif_mean = NA_real_, stringsAsFactors = FALSE))
  }
  r  <- cor(t1_lm$M[match(ov, t1_lm$gene)], v[ov], method = "pearson")
  c8 <- intersect(core8, names(v))
  p8 <- intersect(POLIF_GENES, names(v))
  data.frame(arm = arm, signatureid = label, n_overlap = length(ov),
             r_tier1 = r,
             core8_mean = if (length(c8)) mean(v[c8]) else NA_real_,
             prolif_mean = if (length(p8)) mean(v[p8]) else NA_real_,
             stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## 3) الذراع الأول: SERPINE1 knockdown (7 خطوط خلية)
## ---------------------------------------------------------------------------
say("S8_01: downloading SERPINE1-KD signatures (7 cell lines)")
kd_avail <- read.csv(file.path(OUT, "S8_00_avail_serpine1_kd.csv"),
                     stringsAsFactors = FALSE)
kd_scores <- do.call(rbind, lapply(kd_avail$signatureid, function(sid) {
  say("  KD ", sid, " (", kd_avail$cellline[kd_avail$signatureid == sid], ")")
  score_signature(fetch_signature(sid), sid, "SERPINE1_KD")
}))
kd_scores$cellline <- kd_avail$cellline[match(kd_scores$signatureid,
                                              kd_avail$signatureid)]
write.csv(kd_scores, file.path(OUT, "S8_01_kd_scores.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 4) الذراع الثاني: estradiol LIB_5 -- تجميع median لكل خط خلية (المقياس=خط الخلية)
## ---------------------------------------------------------------------------
say("S8_01: downloading estradiol L1000 signatures and aggregating per cell line")
est_avail <- read.csv(file.path(OUT, "S8_00_avail_estradiol_L1000.csv"),
                      stringsAsFactors = FALSE)
tab <- table(est_avail$cellline)
use_lines <- names(tab)[tab >= MIN_PER_LINE]
say("  estradiol cell lines with >= ", MIN_PER_LINE, " signatures: ",
    length(use_lines), " of ", length(tab))
est_scores <- do.call(rbind, lapply(use_lines, function(cl) {
  sids <- est_avail$signatureid[est_avail$cellline == cl]
  per <- lapply(sids, function(sid) {
    sig <- fetch_signature(sid)
    v <- sig$Value_LogDiffExp; names(v) <- sig$Name_GeneSymbol
    v
  })
  all_genes <- unique(unlist(lapply(per, names)))
  med <- vapply(all_genes, function(g) {
    xs <- vapply(per, function(v) if (g %in% names(v)) v[[g]] else NA_real_, numeric(1))
    median(xs, na.rm = TRUE)
  }, numeric(1))
  ov <- intersect(t1_lm$gene, names(med))
  r <- if (length(ov) >= MIN_OVERLAP)
    cor(t1_lm$M[match(ov, t1_lm$gene)], med[ov], method = "pearson") else NA_real_
  c8 <- intersect(core8, names(med)); p8 <- intersect(POLIF_GENES, names(med))
  data.frame(arm = "estradiol", signatureid = cl, n_overlap = length(ov),
             r_tier1 = r,
             core8_mean = if (length(c8)) mean(med[c8]) else NA_real_,
             prolif_mean = if (length(p8)) mean(med[p8]) else NA_real_,
             stringsAsFactors = FALSE)
}))
write.csv(est_scores, file.path(OUT, "S8_01_estradiol_scores.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5) الاختبارات المقفلة: عائلة واحدة من 4 فرضيات، Wilcoxon موقّع، BH-FDR
## ---------------------------------------------------------------------------
say("S8_01: running the pre-registered hypothesis family")
wilcox_vs0 <- function(x, alternative) {
  x <- x[!is.na(x)]
  if (length(x) < 4L) return(list(p = NA_real_, n = length(x)))
  p <- suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative)$p.value)
  list(p = p, n = length(x))
}

t1a <- wilcox_vs0(kd_scores$r_tier1,   "less")
t1b <- wilcox_vs0(kd_scores$prolif_mean, "greater")
t2a <- wilcox_vs0(est_scores$r_tier1,  "less")
t2b <- wilcox_vs0(est_scores$prolif_mean, "greater")

pvals <- c(P1a = t1a$p, P1b = t1b$p, P2a = t2a$p, P2b = t2b$p)
fdr <- p.adjust(pvals, method = "BH")

hyp <- data.frame(
  hypothesis = names(pvals),
  test = c("SERPINE1-KD r_tier1 < 0 (7 cell lines)",
           "SERPINE1-KD prolif(TOP2A,PCNA) > 0",
           "estradiol per-line r_tier1 < 0 (cell-line medians)",
           "estradiol prolif > 0"),
  n = c(t1a$n, t1b$n, t2a$n, t2b$n),
  p_value = unname(pvals), BH_FDR = unname(fdr),
  stringsAsFactors = FALSE)
hyp$verdict <- ifelse(!is.na(hyp$BH_FDR) & hyp$BH_FDR < 0.05,
                      "CONFIRMED (FDR<0.05, pre-registered direction)",
                      "NOT CONFIRMED")
write.csv(hyp, file.path(OUT, "S8_01_hypotheses.csv"), row.names = FALSE)

med_kd_r   <- median(kd_scores$r_tier1, na.rm = TRUE)
med_est_r  <- median(est_scores$r_tier1, na.rm = TRUE)
med_kd_pr  <- median(kd_scores$prolif_mean, na.rm = TRUE)
med_est_pr <- median(est_scores$prolif_mean, na.rm = TRUE)

if (hyp$verdict[1] == "CONFIRMED (FDR<0.05, pre-registered direction)") {
  story <- "P1 CONFIRMED: SERPINE1 knockdown reverses the locked lesion Tier-1 signature (external, pharmacologically-directed evidence)"
} else {
  story <- "P1 NOT CONFIRMED: SERPINE1 knockdown does not significantly reverse the lesion signature -- reported as such"
}

if (hyp$verdict[3] == "CONFIRMED (FDR<0.05, pre-registered direction)") {
  story2 <- "P2 CONFIRMED: estradiol stimulation reverses the locked lesion Tier-1 signature"
} else {
  story2 <- "P2 NOT CONFIRMED: estradiol reversal not significant -- reported as such"
}

summary_md <- c(
  "# S8_01 -- Reverse-confirmation summary (machine-printed)",
  "",
  paste0("Date: ", Sys.time()),
  "",
  "## Pre-registered hypothesis family (one BH family, sole gate FDR < 0.05)",
  "",
  "| H | test | n | p | BH-FDR | verdict |",
  "|---|------|---|---|--------|---------|",
  paste0("| ", hyp$hypothesis, " | ", hyp$test, " | ", hyp$n, " | ",
         signif(hyp$p_value, 3), " | ", signif(hyp$BH_FDR, 3), " | ",
         hyp$verdict, " |"),
  "",
  "## Median effects",
  paste0("- SERPINE1-KD: median r_tier1 = ", round(med_kd_r, 3),
         " | median prolif = ", round(med_kd_pr, 3)),
  paste0("- estradiol:   median r_tier1 = ", round(med_est_r, 3),
         " | median prolif = ", round(med_est_pr, 3)),
  "",
  "## Story verdicts (machine-printed)",
  paste0("- **", story, "**"),
  paste0("- **", story2, "**"),
  "",
  "*Negative r_tier1 = the perturbation moves landmark-measurable Tier-1 genes in the",
  " direction OPPOSITE to the lesion (reversal). Proliferation readout uses the two",
  " proliferation-axis genes measurable in L1000 landmark space (TOP2A, PCNA; MKI67",
  " is absent from the 978). The 66-gene robust core is represented by its 8 landmark",
  " genes in this space (coverage measured in S8_00). Nothing here is a claim beyond",
  " the locked FDR gate.*")
writeLines(summary_md, file.path(OUT, "S8_01_summary.md"))

## ---------------------------------------------------------------------------
## 6) manifest
## ---------------------------------------------------------------------------
outs <- c("S8_01_kd_scores.csv", "S8_01_estradiol_scores.csv",
          "S8_01_hypotheses.csv", "S8_01_summary.md")
man2 <- data.frame(
  file = outs,
  md5 = vapply(file.path(OUT, outs), function(p) unname(tools::md5sum(p)),
               character(1)),
  stringsAsFactors = FALSE)
write.table(man2, file.path(OUT, "S8_01_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

say("S8_01 DONE")
say("  ", story)
say("  ", story2)
say("  outputs + manifest in ", OUT)
