# About

This repository hosts a series of scripts that analyze and visualize data for the Beyond Auxotrophy project at the Bertrand Lab of marine microbial proteomics and trace nutrient biogeochemistry. For this experiment, a diatom (F. cylindrus) and bacterial isolate (Pibocella sp.) were cultured with and without the addition of vitamin B12.

# Folder Directory

## 0_TSQ_Frag_Metab_BMIS

This script uses the quality control (QC) samples to choose a best matched internal standard (BMIS). When normalized to it's BMIS, a metabolite's coefficient of variation in the QC samples injected throughout the course of the run is decreased by at least 30%.

### Inputs:

-   Raw skyline data (ER3_149_Frag_Catalina_01_Cal03022022_output.csv)

### Outputs

-   Lists of BMIS per metabolite, with information about decreases in CV and final normalized peaks if a BMIS was chosen (QC_Norm_Export_Sum.csv and QC_Norm_Export_Sum.csv)

## 1_TSQ_Frag_Metab_CalCurveQuant

This script performs calibrations and normalizes to heavy standards for absolute quantification. It computes limits of detection and quantification for each compound. Finally, it filters for only metabolites which can be quantified (either relative quantification in the form of peak per cell or absolute quantitication in moles per cell).

### Inputs:

-   A list of molecules, which tells whether we did a calibration curve for them and the initial spike amount for each metabolite (TSQ_Metab_CalSum.csv)
-   A lookup file that relates sample names to amount spiked (spike_amount_calcurve.csv). Here, the cal_amount column contains 1 or 5, and denotes how much of each sample was spiked. (if 1, amounts in sample name are correct. ex: p5fmol = 0.5 fmol spike. if 5, 5x that. ex: p5fmol = 2.5 fmol ).
-   A file that contains information about the number of cells loaded on column and the amount of carbon per cell (B_sample_info_loading_220116.csv)

### Outputs

-   A list of compounds and their respective LOD/LOQ's after normalization (LODQ_export.csv)
-   A list of all of the compounds, if they were normalized, and their R-squared value when readings were compared to the BMIS if r-squared is less than .5, no calibration is applied (rsq_metab.csv)
