# تصميم مرحلة S3 المُفصل — المقارنات الأولية (مجمّد، بانتظار اعتمادك الكتابي قبل أي رن)

تاريخ التجميد: 2026-09-18 · الحاكم: الصلابة العلمية أمام المراجع + الدفاع المنهجي + عدم حذف أي معلومة ذات قيمة للقصة (الأقل دفاعًا = حساسيات لا أساس؛ لا إقصاء إلا بسبب موثق).

## 1. المعايير الموقّعة (معتمدة، حاكمة على كل مخرجات S3)

| # | المعيار | التنفيذ في السكريبتات |
|---|---------|------------------------|
| ① | العتبة: BH-FDR < 0.05 + طبقة تأويلية \|log2FC\| ≥ 0.5 بدون إقصاء + توافق الاتجاه مُبلَّغ | `CFG$fdr_gate` / `CFG$lfc_tier` في كل المخرجات؛ لا سطر يستبعد بحجم الأثر؛ S3_06 يبلغ التوافق ولا يفلتر به |
| ② | scRNA: pseudobulk per sample أساس + Wilcoxon per lineage حساسية | S3_03 (أساس) / S3_04 (حساسية فقط، كل صف `sensitivity=TRUE` + تنبيه pseudoreplication) |
| ③ | ssGSEA للـ bulk + z-mean للـ scRNA | S3_02 (GSVA ssgseaParam، موثّق بالسورس) / S3_05 (تعريف z-mean مجمّد أدناه) |
| ④ | الدمج في S4: random-effects أساسي + fixed حساسية + I²/τ² | ليست جزءًا من S3 — S3 يخرج جداول effect sizes وسكورات جاهزة لها فقط |
| ⑤ | S5–S7 تُسجَّل عند فتحها | لا شيء هنا |
| ⑥ | المجلة مؤجلة؛ Main ≤ 6 أشكال | لا أشكال في S3 أصلًا |
| ⑦ | powercfg بلا تغيير حتى ظهور رن > ساعتين | أقصى تقدير لأي رنة في S3 ≪ ساعتين (جدول R10) — القرار متسق |

## 2. جدول الـ contrasts المجمّد (18 contrast: 7 أساسية + 4 ثانوية + 4 حساسية + 3 استكشافية)

الملف الآلي: `scripts/S3_contrasts_frozen.tsv` (MD5 بالأسفل). القراءة البشرية:

### أساسية (primary)
| id | cohort | المقارنة | التصميم | ملاحظات |
|----|--------|----------|---------|---------|
| B-G7305 | GSE7305 | Endometriosis vs Normal endometrium | غير مقترن | 10 vs 10 |
| B-G6364 | GSE6364 | Endometriosis eutopic vs Normal eutopic | غير مقترن + covariate phase (3 مستويات مجمّدة من عمود detail) | 21 vs 16 |
| B-G25628 | GSE25628 | كل المريضات (Eutopic + Ectopic) vs Normal endometrium | غير مقترن فقط | 16 vs 6 |
| B-G7307 | GSE7307 | Endometriosis vs Normal endometrium | غير مقترن | 18 vs 23؛ حساسية إقصاء الكوهورت كاملة في S3_06 (تحفظ E2) |
| B-G51981 | GSE51981 | كل الـ Endometriosis (77 بما فيها 2 severity-NA — لا سبب موثق للإقصاء فتُبقى) vs Non-EMS no pelvic pathology | غير مقترن | 77 vs 34؛ أنظف مجموعة ضابط |
| B-G11691 | GSE11691 | Ectopic lesion vs Eutopic (patient) | مقترن 9 سيدات، duplicateCorrelation | ملف الاقتران المجمّد مع تدقيق ذاتي 9×2 |
| B-G120103 | GSE120103 | (2A+2B) vs Fertile control | غير مقترن | 13 vs 9 بعد الإقصاءات الموثقة (2B=4) |

### ثانوية (secondary، نتائج منشورة لكنها ليست الأساس)
- B-G51981-MILD / B-G51981-SEV: طبقات الشدة ضد الضابط النظيف
- B-G120103-2A / B-G120103-2B (n=4، report-only صراحةً في الملاحظات)

### حساسية (sensitivity — تُنشر بجانب الأساس)
- B-G51981-S1: باستبعاد GSM1256659_134 (المُعلَّم QC) · B-G51981-S2: باستبعاد الزوج الحدّي _103+_134 · B-G51981-S3: بضابط مجمّع (no pathology + pelvic pathology)
- B-G120103-IC: (2A+2B) vs **Infertile control** — يحفظ قصة العقم في القصة البيولوجية (لا حذر لمعلومة ذات قيمة)

### استكشافية (exploratory — hypothesis-generating، صراحةً ليست نتائج)
- B-G25628-EXPL: Ectopic vs Eutopic غير مقترن · B-G120103-FCIC: Fertile vs Infertile control (أثر العقم نفسه)

### scRNA (من S3_03/05، بذراعي N1 المجمدين Amendment S2-B/B3)
- A (أساس): pseudobulk لكل dataset × lineage (+ ساق ALL بلا تقسيم) — N1 **مشمول**
- B: نفسه مع **استبعاد N1** — والذراعان تُنشران
- A-HCOV / A-NOTX (GSE179640 فقط): covariate الهرمونات إلزامي + ذراع باستبعاد المعالجات — معلّقة على تأكيد عمود الهرمونات من الـ census

## 3. تعريفات مجمّدة أخرى
- **z-mean (S3_05):** لكل خلية: متوسط z-scores الجين-لجين (محسوبة على log1p عبر كل خلايا الأطلس) لجينات SenMayo الموجودة في الكائن. تلخيص لكل dataset × lineage × sample، ثم مقارنة sample-level بـ limma.
- **الـ collapse للـ ssGSEA (S3_02):** لكل رمز جيني، يُحتفظ بالـ probe صاحب أعلى متوسط تعبير؛ تُصبح rownames هي الرموز (شرط مطابقة GSVA). خرائط Affy عبر hgu133plus2.db / hgu133a.db / hgu133a2.db، وAgilent عبر عمود الرموز في genes (الـ census يحدده، وغياب أي منها = STOP بانتظار موافقتك R5).
- **فلتر pseudobulk:** ≥10 خلايا لكل (lineage × sample)؛ فحص CPM>1 في ≥ نصف العينات قبل voom.
- **GSE11691 pairing:** 9 سيدات مثبتات في `scripts/S3_GSE11691_pairs.tsv` مع تدقيق ذاتي 9×2=18 (أي مخالفة = FATAL، وأي تصحيح = Amendment لا تعديل صامت).
- **GSVA:** مركبة الـ ssGSEA هي GSVA (معيار قابل للاستشهاد)؛ الفرع الرئيسي `ssgseaParam`+`gsva` (GSVA ≥ 1.50، تم التحقق من السورس 2.6.6)، وفرع legacy للإصدارات الأقدم. غياب GSVA = STOP بانتظار موافقة كتابية — **لا** يُستبدل بمُنفِّذ يدوي بلا Amendment.

## 4. سير التنفيذ (لا أي رن قبل اعتمادك)

```
scripts\run_S3.bat
  → S3_smoke.R        إلزامي أولًا (B.4): 8 أرجل اصطناعية، أي FAIL يوقف المرحلة
  → S3_00_preflight.R census قراءة-فقط: يطبع قيم الإعدادات المطلوب تعبئتها
  → [تعبئة scripts/S3_common.R من لوج الـ census — خطوة مراجعة منك]
  → S3_01 (cohort×7 + collect) → S3_02 → S3_03 → S3_04 → S3_05 → S3_06
```

أوامر يدوية بديلة (paste-ready، من جذر المشروع):
```
"C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_smoke.R > logs\S3_smoke.txt 2>&1
"C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_00_preflight.R > logs\S3_00_preflight.txt 2>&1
"C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_01_bulk_contrasts.R GSE51981 > logs\S3_01_GSE51981.txt 2>&1
```

## 5. تقديرات R10 (16GB، نهاري — القرار ⑦ متسق: لا شيء هنا يقترب من ساعتين)

| رنة | الذروة التقديرية | الزمن التقديري |
|-----|------------------|----------------|
| S3_smoke | < 1 GB | < دقيقة |
| S3_00 preflight | < 2 GB (تحميل الأطلس للـ census) | < 3 دقائق |
| S3_01 لكل cohort | < 1.5 GB (GSE51981 أكبر) | ثوانٍ–دقائق لكل cohort؛ الكل < 10 دقائق |
| S3_02 ssGSEA | < 3 GB | < 10 دقائق |
| S3_03 pseudobulk | ≈ 4–5 GB (أكبر رنة) | 10–15 دقيقة |
| S3_04 Wilcoxon (حساسية) | ≈ 5–6 GB | 20–40 دقيقة |
| S3_05 z-mean | ≈ 4–5 GB | < 10 دقائق |
| S3_06 collect | < 1 GB | < دقيقتين |

## 6. بنود مفتوحة تُغلق من لوج S3_00 (تعبئة 5 أسطر في `S3_common.R` + أي موافقات)

1. `n1_sample` — معرف N1 في GSE214411 (الـ census يطبع كل قيم `N1` تلقائيًا).
2. `hormonal_col` + `hormonal_treated_lbl` — عمود الهرمونات وقيمة "معالج" (للذراعين A-HCOV / A-NOTX).
3. `group_col` + `ems_labels` + `ctrl_labels` — عمود الشرط وقيم EMS/Control.
4. `lineage_col` — عمود الأنماط الثمانية من S2_02.
5. حضور GSVA و`.db` packages والـ GMT — إن غاب أي منها: STOP بانتظار موافقتك الكتابية (R5).
6. تأكيدك البصري لملف اقتران GSE11691 ضد characteristics المصفوفة (أي اختلاف = Amendment).

## 7. جرد الملفات و MD5

| ملف | MD5 |
|-----|-----|
| S3_contrasts_frozen.tsv | c23feac0880fa6f832ee21a1b0c647f2 |
| S3_GSE11691_pairs.tsv | 4ff386f6633385fde3f9bec3118aac45 |
| S3_common.R | 0da7bb5c64347f740bda1ebd9ef220f7 |
| S3_00_preflight.R | 000a3685eff654fa006277750ea3e693 |
| S3_01_bulk_contrasts.R | 5c646432616d38837bd2bcfad380626f |
| S3_02_ssgsea_bulk.R | a4ccd8c1e816cd14f813153c4bc2c34c |
| S3_03_scrna_pseudobulk.R | 899126ae24d51747bdac2f123390d080 |
| S3_04_scrna_wilcoxon_sens.R | 1126e6eaf8c5c494c69fe64bd8968f03 |
| S3_05_scrna_zmean.R | f59c5fade1cbd366dcf138f25265c61d |
| S3_06_collect.R | a7c89523a41cd05863e8c02881bd582c |
| S3_smoke.R | 0432d2d663f2beda9f470f01bc18fb91 |
| run_S3.bat | c211655f0f5052822c3d757c228b084e |

ملاحظة حوكمة: سكريبتات scRNA تسحب الطبقات عبر `get_rna()` (JoinLayers عند الحاجة — Seurat 5)؛ S3_03 يشتغل per-dataset (تصميم أكثر دفاعًا، والدمج بين الدراسات يحدث في S4 بالعشوائيات المقفلة ④). لا يوجد أي إقصاء عينة في S3 إلا المُعلَّمات الموثقة (GSM1256659_134 في ذراع S1 فقط، والخمس GE1-v5 في GSE120103 المسبقة التوثيق).
