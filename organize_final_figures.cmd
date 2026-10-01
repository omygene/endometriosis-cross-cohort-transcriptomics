@echo off
setlocal

set "P=C:\endometriosis_immunoediting"
set "OUT=%P%\figures\FINAL_MANUSCRIPT"

if not exist "%OUT%" mkdir "%OUT%"

copy /Y "%P%\figures\Fig-1.jpg" "%OUT%\Fig1_Study_design.jpg"
copy /Y "%P%\figures\fig2_tier1_core.jpg" "%OUT%\Fig2_Gene_meta_analysis.jpg"
copy /Y "%P%\figures\Fig3d_robust_core_summary.jpg" "%OUT%\Fig3_Robust_core.jpg"
copy /Y "%P%\figures\fig3_sensitivity.jpg" "%OUT%\FigS1_Sensitivity_analysis.jpg"
copy /Y "%P%\figures\fig5_pathway_meta_AB.jpg" "%OUT%\Fig4A_Pathway_meta.jpg"
copy /Y "%P%\figures\fig5_pathway_meta_CD.jpg" "%OUT%\Fig4B_Paired_concordance.jpg"
copy /Y "%P%\figures\fig4_localization.jpg" "%OUT%\Fig5_scRNA_localization.jpg"
copy /Y "%P%\figures\Fig6a_validation_design.jpg" "%OUT%\Fig6A_Validation_design.jpg"
copy /Y "%P%\figures\fig6_GSE213216_validation.jpg" "%OUT%\Fig6B_External_validation.jpg"
copy /Y "%P%\figures\FigS7_test_matrix.jpg" "%OUT%\FigS2_Validation_test_matrix.jpg"
copy /Y "%P%\figures\Fig7a_lincs_coverage.jpg" "%OUT%\FigS3_LINCS_coverage.jpg"
copy /Y "%P%\figures\fig7_s8_reverse_confirmation.jpg" "%OUT%\FigS4_LINCS_reversal_negative.jpg"
copy /Y "%P%\figures\Fig7d_hypotheses.jpg" "%OUT%\FigS5_Pre_registered_hypotheses.jpg"
copy /Y "%P%\figures\figS1_drugs_exploratory.jpg" "%OUT%\FigS6_Exploratory_decidualization.jpg"

echo.
echo DONE
echo Final figures are here:
echo %OUT%
start "" "%OUT%"
pause