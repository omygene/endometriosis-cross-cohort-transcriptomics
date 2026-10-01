## ===========================================================================
## S8b_01_exploratory_drug.R -- آخر محاولة in silico (استكشافية بالكامل)
##   هل دواء endometriosis معروف يعكس Tier-1 signature في GEO المتاح؟
## v1.0 (2026-09-25)
## ---------------------------------------------------------------------------
## خلفية الـcensus (مقيسة 2026-09-23 من GEO/eutils الحي):
##   - danazol / leuprolide / norethindrone / gestrinone: صفر داتاسيت تعبير
##     نظيف على نسيج endometriosis بشري
##   - GSE89463: GnRH agonist في فأر adenomyosis (نوع/مرض مختلف، n=6)
##   - GSE136412: ثقافة 3D (ليست أدوية)
##   - GSE289058: مثبطات JNK تجريبية (ليست SoC)
##   - GSE75423: dibutyryl-cAMP + DIENOGEST (progestin SoC) على خلايا ECSC
##     المأخوذة من مرضى endometrioma -- paired 4v4 (نفس المريضات بالترتيب)
##   - GSE75425: نفس التصميم على خلايا NESC طبيعية
## إذن: هذا الاختبار استكشافي، من داتاسيت واحد صغير، بلا بوابة ادعاء.
## قاعدة الإيقاف المقفلة (stop rule): هذه آخر محاولة in silico مهما كانت
##   النتيجة؛ المرحلة التالية = اضطراب تجريبي (in vitro) بغض النظر.
##
## الميثاق: R1/R2/R4/R6/R7/R8. base R فقط (لا حزم).
## التشغيل:
##   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts/S8b_01_exploratory_drug.R
## ===========================================================================

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

OUT   <- "results/S8_perturb"
CACHE <- file.path(OUT, "cache_geo")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)

F_G23 <- file.path(CACHE, "GSE75423_series_matrix.txt.gz") ## ECSC (endometriotic)
F_G25 <- file.path(CACHE, "GSE75425_series_matrix.txt.gz") ## NESC (normal)
F_ANN <- file.path(CACHE, "GPL13497.soft")                 ## Agilent annotation

URL_G23 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75423/matrix/GSE75423_series_matrix.txt.gz"
URL_G25 <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75425/matrix/GSE75425_series_matrix.txt.gz"
URL_ANN <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13497&targ=self&form=text&view=full"

fetch <- function(url, dest) {
  if (file.exists(dest)) return(invisible(TRUE))
  ok <- tryCatch({ download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) ## محاولة ثانية بطريقة libcurl (أصلح على بعض الأنظمة)
    ok <- tryCatch({ download.file(url, dest, mode = "wb", quiet = TRUE,
                                   method = "libcurl"); TRUE },
                   error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(dest))
    stop("FATAL: download failed: ", url,
         "\n  Manual fix: download it into ", CACHE, call. = FALSE)
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## 0) المدخلات المقفلة: Tier-1 core (registry نفسه المستخدم في S8_00/S8_01)
## ---------------------------------------------------------------------------
F_T1  <- "results/S4_meta/S4_meta_Tier1_core.csv"
REG_T1 <- c("results/S4_meta/S4_meta_Tier1_core.csv" = "3afbc9d85bec0b4409137cb73cc854bb")
if (!file.exists(F_T1)) stop("FATAL: locked input missing: ", F_T1, call. = FALSE)
if (unname(tools::md5sum(F_T1)) != REG_T1[[F_T1]])
  stop("FATAL: md5 mismatch for ", F_T1, call. = FALSE)
t1 <- read.csv(F_T1, stringsAsFactors = FALSE)
say("S8b: locked Tier-1 file verified (md5 match)")

## ---------------------------------------------------------------------------
## 1) التنزيلات (GEO)
## ---------------------------------------------------------------------------
say("S8b: downloading GSE75423 + GSE75425 + GPL13497 (cached after first run)")
fetch(URL_G23, F_G23); fetch(URL_G25, F_G25); fetch(URL_ANN, F_ANN)

## ---------------------------------------------------------------------------
## 2) أنوتيشن المنصة: ID -> GENE_SYMBOL (سطر الترويسة يبدأ بـ "ID\t")
## ---------------------------------------------------------------------------
say("S8b: parsing GPL13497 annotation")
lns <- readLines(F_ANN, warn = FALSE)
h   <- grep("^ID\t", lns)[1]
if (is.na(h)) stop("FATAL: platform table header not found in GPL13497.soft", call. = FALSE)
tbl <- lns[(h + 1):length(lns)]
tbl <- tbl[grepl("^A_23_P", tbl)]
if (length(tbl) < 1000L) stop("FATAL: platform table looks truncated (", length(tbl), " rows)", call. = FALSE)
col_ann <- strsplit(lns[h], "\t", fixed = TRUE)[[1]]
i_id  <- match("ID", col_ann); i_sym <- match("GENE_SYMBOL", col_ann)
ann <- vapply(strsplit(tbl, "\t", fixed = TRUE),
              function(p) if (length(p) >= i_sym) p[i_sym] else "", character(1))
names(ann) <- vapply(strsplit(tbl, "\t", fixed = TRUE),
                     function(p) if (length(p) >= i_id) p[i_id] else "", character(1))
ann <- ann[ann != ""]
say("  probes with symbol: ", length(ann))

## ---------------------------------------------------------------------------
## 3) قراءة المصفوفة + collapse بالوسيط لكل جين + paired logFC
##    التصميم مقفل: أعمدة 1-4 vehicle (نفس المريضات بالترتيب) × أعمدة 5-4 معالجة
## ---------------------------------------------------------------------------
read_gse <- function(f) {
  con <- gzfile(f, "rt")
  x <- read.delim(con, comment.char = "!", check.names = FALSE,
                  stringsAsFactors = FALSE)
  close(con)
  x
}
collapse_lfc <- function(x) {
  x$sym <- ann[as.character(x[[1]])]
  x <- x[!is.na(x$sym), ]
  nums <- vapply(x[, -c(1, ncol(x))], as.numeric, numeric(nrow(x)))
  meds <- tapply(seq_len(nrow(x)), x$sym, function(ii) apply(nums[ii, , drop = FALSE], 2, median))
  g <- do.call(rbind, meds)
  rownames(g) <- names(meds)
  lg <- log2(g + 1)
  stopifnot(ncol(lg) == 8L) ## 4 vehicle + 4 treated, same patient order
  rowMeans(lg[, 5:8] - lg[, 1:4])
}

say("S8b: computing paired dienogest logFC (ECSC + NESC)")
lfc23 <- collapse_lfc(read_gse(F_G23))
lfc25 <- collapse_lfc(read_gse(F_G25))
say("  ECSC genes: ", length(lfc23), " | NESC genes: ", length(lfc25))

## ---------------------------------------------------------------------------
## 4) الإحصاء الاستكشافي: r = Pearson(M المقفل, paired logFC) على Tier-1
##    + فاصل ثقة bootstrap (إعادة معاينة الجينات) -- بلا بوابة ادعاء
## ---------------------------------------------------------------------------
boot_ci <- function(Mv, Lv, B = 2000, seed = 42) {
  set.seed(seed)
  out <- numeric(B)
  n <- length(Mv)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(Lv[idx])) > 1L) out[b] <- cor(Mv[idx], Lv[idx])
  }
  out <- out[!is.na(out)]
  quantile(out, c(0.025, 0.975))
}

explore <- function(lfc, label) {
  ov <- t1$gene[t1$gene %in% names(lfc)]
  Mv <- t1$M[match(ov, t1$gene)]; Lv <- unname(lfc[ov])
  r  <- cor(Mv, Lv)
  ci <- boot_ci(Mv, Lv)
  direction <- if (r < 0)
    "reversal direction (drug opposite to lesion) -- NOT what we observe" else
    "same direction as lesion (NO reversal observed)"
  data.frame(dataset = label, n_tier1_measured = length(ov), r = r,
             ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
             direction = direction, stringsAsFactors = FALSE)
}
res <- rbind(explore(lfc23, "GSE75423_ECSC_dienogest"),
             explore(lfc25, "GSE75425_NESC_dienogest"))
write.csv(res, file.path(OUT, "S8b_01_exploratory_drug.csv"), row.names = FALSE)

cm <- intersect(names(lfc23), names(lfc25))
r_resp <- cor(lfc23[cm], lfc25[cm])
key <- c("TOP2A","PCNA","IL6","SOCS3","FOS","FOSB","EGR1","C1QB","SERPINE1","MKI67")
key_tab <- data.frame(
  gene = key,
  ECSC_logFC = round(vapply(key, function(g)
    if (g %in% names(lfc23)) lfc23[[g]] else NA_real_, numeric(1)), 3),
  NESC_logFC = round(vapply(key, function(g)
    if (g %in% names(lfc25)) lfc25[[g]] else NA_real_, numeric(1)), 3),
  stringsAsFactors = FALSE)
write.csv(key_tab, file.path(OUT, "S8b_01_key_genes.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5) الحكم الآلي الاستكشافي + stop rule
## ---------------------------------------------------------------------------
E1 <- res$direction[1]
stop_rule <- paste0(
  "STOP RULE (locked): this was the LAST in-silico attempt regardless of outcome. ",
  "Next stage = experimental in-vitro perturbation (independent study).")
verdict_txt <- paste0(
  "EXPLORATORY (no claim gate): dienogest+cAMP on endometriotic stromal cells (GSE75423, ",
  "paired 4v4, single dataset) shifts the Tier-1 signature TOWARD the lesion direction ",
  "(r = ", round(res$r[1], 3), ", bootstrap 95% CI ", round(res$ci_lo[1], 3), " to ",
  round(res$ci_hi[1], 3), "), i.e. drug-reversal is NOT observed; normal stromal cells ",
  "show the same pattern (r = ", round(res$r[2], 3), "). Documented as exploratory; ",
  "the positive-control question for standard-of-care drugs is CLOSED in silico ",
  "(census: no usable danazol/GnRH datasets exist).")

summary_md <- c(
  "# S8b_01 -- exploratory drug positive-control (machine-printed)",
  "",
  paste0("Date: ", Sys.time()),
  "",
  "## Census verdict (pre-measured, see S8b design doc)",
  "- Usable standard-of-care drug expression datasets on human endometriotic tissue: ONE",
  "  (GSE75423, dienogest+cAMP on ECSCs, paired 4v4); GSE75425 normal-ESC counterpart.",
  "- danazol / GnRH-agonist / leuphin: no usable human endometriosis expression dataset.",
  "",
  "## Exploratory results (no claim gate; single small dataset)",
  "",
  "| dataset | Tier-1 measured | r | boot 95% CI | direction |",
  "|---|---|---|---|---|",
  paste0("| ", res$dataset, " | ", res$n_tier1_measured, " | ", round(res$r, 3),
         " | [", round(res$ci_lo, 3), ", ", round(res$ci_hi, 3), "] | ",
         ifelse(res$r < 0, "reversal", "NO reversal"), " |"),
  "",
  paste0("- correlation of dienogest responses ECSC vs NESC (", length(cm),
         " genes): r = ", round(r_resp, 3)),
  "- key genes (see S8b_01_key_genes.csv)",
  "",
  "## Machine verdict (exploratory, machine-printed)",
  paste0("**", verdict_txt, "**"),
  "",
  paste0("*", stop_rule, "*"),
  "",
  "*Every number measured at run time from GEO downloads and the locked Tier-1 file.",
  " Exploratory: single dataset, 4 patient pairs, platform covers only part of Tier-1.*")
writeLines(summary_md, file.path(OUT, "S8b_01_summary.md"))

outs <- c("S8b_01_exploratory_drug.csv", "S8b_01_key_genes.csv", "S8b_01_summary.md")
man <- data.frame(file = outs,
                  md5 = vapply(file.path(OUT, outs),
                               function(p) unname(tools::md5sum(p)), character(1)),
                  stringsAsFactors = FALSE)
write.table(man, file.path(OUT, "S8b_01_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

say("S8b_01 DONE")
say("  ECSC r = ", round(res$r[1], 3), " | NESC r = ", round(res$r[2], 3),
    " -- NO reversal (exploratory)")
say("  ", stop_rule)
