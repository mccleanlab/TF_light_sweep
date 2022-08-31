# Instructions
This repository contains the code to make plots for the paper "Transcription Factor Localization Dynamics and DNA Binding Drive Distinct Promoter Interpretations"

Due to file size restrictions:
-1. Raw microscopy images (.ND2) are available only upon request
-2. Single-cell fluorescence measurements extracted from the .ND2 images are available: INSERT MENDELEY DATA URL WHEN UPLOAD STOPS FAILING
-3. Gene expression model solutions are available here: INSERT MENDELEY DATA URL WHEN UPLOAD STOPS FAILING

If you have the raw microscopy images and want to calculate everything from scratch, run the scripts as follows:
-1. Run calc_1_run_cellfinder_in_folder_light_sweep.m to segment and measure yeast from .ND2 files (requires github.com/mccleanlab/cell_finder)
-2. Run calc_2_collect_measurements_light_sweep.m to collect/organize single-cell measurements from many experiments
-3. Run calc_3_analyze_collected_measurements_light_sweep.m to process measurements (will export simpler .mat files to work with)
-4. Run calc_4_run_promoter_model_LHS.m to calculate solutions to gene expression model (slow, depending on number of guesses)

If you just want to make plots from paper
-1. Download single-cell fluorescence measurements (collected_measurements.zip) and gene expression model solutions (model_solutions.zip) and unzip into this folder
-2. Run calc_3_analyze_collected_measurements_light_sweep.m to process measurements (will export simpler .mat files to work with)
-3. Run plot_figures_*.m (requires https://www.mathworks.com/matlabcentral/fileexchange/54465-gramm-complete-data-visualization-toolbox-ggplot2-r-like)

