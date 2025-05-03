# About

This repository hosts a series of scripts that analyze and visualize data for the Beyond Auxotrophy project at the Bertrand Lab of marine microbial proteomics and trace nutrient biogeochemistry. For this experiment, a diatom (*F. cylindrus*) were cultured with and without the addition of vitamin B<sub>12</sub>. Through targeted mass spectrometry, we investigate the intracellular concentrations of a suite of metabolites when the facultative cobalamin consumer is grown with and without cobalamin

# Folder Directory

## 0_TSQ_Frag_Metab_BMIS

This script in this folder, [TSQ_Frag_Metab_BMIS.R](0_TSQ_Frag_Metab_BMIS/TSQ_Frag_Metab_BMIS) uses the quality control (QC) samples to choose a best matched internal standard (BMIS). When normalized to it's BMIS, a metabolite's coefficient of variation in the QC samples injected throughout the course of the run is decreased by at least 30%.


### Inputs:

-   Raw skyline data for calibration samples [ER3_149_Frag_Catalina_01_Cal03022022_output.csv](TSQ_Frag_Metab_RAW/ER3_149_Frag_Catalina_01_Cal03022022_output.csv)

### Outputs

- [A list of BMIS per analyte](0_TSQ_Frag_Metab_BMIS/QC_BMIS_results_sum.csv) and information about decreases in CV and final normalized peaks if a BMIS was chosen

## 1_TSQ_Frag_Metab_CalCurveQuant

[TSQ_Frag_Metab_CalCurveQuant.R](1_TSQ_Frag_Metab_CalCurveQuant/TSQ_Frag_Metab_CalCurveQuant.R) performs calibrations and normalizes to heavy standards for absolute quantification. It computes limits of detection and quantification for each compound. Finally, it filters for only metabolites which can be quantified (either relative quantification in the form of peak per cell or absolute quantification in moles per cell).

### Inputs:

-   [TSQ_Metab_CalSum.csv](1_TSQ_Frag_Metab_CalCurveQuant/TSQ_Metab_CalSum.csv), which contains information about whether we did a calibration curve for each molecule and the initial spike amount for each metabolite
-   [spike_amount_calcurve.csv](1_TSQ_Frag_Metab_CalCurveQuant/spike_amount_calcurve.csv), a lookup file that relates sample names to amount spiked. Here, the cal_amount column contains 1 or 5, and denotes how much of each sample was spiked. (if 1, amounts in sample name are correct. ex: p5fmol = 0.5 fmol spike. if 5, 5x that. ex: p5fmol = 2.5 fmol).
-   [B_sample_info_loading_220116.csv](1_TSQ_Frag_Metab_CalCurveQuant/B_sample_info_loading_220116.csv), which contains information about the number of cells loaded on column and the amount of carbon per cell ()

### Outputs

-   [LODQ_export.csv](1_TSQ_Frag_Metab_CalCurveQuant/LODQ_export.csv), a list of compounds and their respective LOD/LOQ's after normalization
-   [rsq_metab.csv](1_TSQ_Frag_Metab_CalCurveQuant/rsq_metab.csv), a list of all of the compounds, if they were normalized, and their R-squared value when readings were compared to the BMIS if r-squared is less than .5, no calibration is applied
-   [recommended_norm_df.csv](1_TSQ_Frag_Metab_CalCurveQuant/recommended_norm_df.csv), a list of recommended normalization by analyte and the type of quantification we can obtain based on LODs/LOQs

## 2_TSQ_Frag_Metab_CobalaminStackedBar
[TSQ_Frag_Metab_CobalaminStackedBar.R](2_TSQ_Frag_Metab_CobalaminStackedBar/TSQ_Frag_Metab_CobalaminStackedBar.R) 

## TSQ_Frag_Metab_RAW




