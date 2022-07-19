# BHT_Unct
Uncertainty based model to correct Bottom Hole Temperatures of deep wells to static formation temperatures with conventional analytical and numerical correction methods.

# About
This Python tool is designed to estimate the formation temperature from thermally disturbed Bottom Hole Temperatures (BHTs). BHT data are usually of poor quality, therefore it is important to take into account the uncertainty of the required input parameters to the different models that can be qpplied for correction. This tool allows for specifying the quality of input data (uncertainty range), performing a Saltelli sampling (Saltellis extension of Sobol) and running one of six conventional BHT correction schemes.
The output is provided as density distribution of the result space. p10, p50, p90 or modal value can be used to describe the corrected temperature value with its uncertainty, dependent on the input ranges.


# Licence
BHT_Unct is distributed under the GNU GENERAL PUBLIC LICENSE v3.

# important files
BHT_Unct.py         main file. User input is requested here
bht_functions.py    includes BHT correction methods
Iterativ_Goetzl.py  includes the iterative solved Linearization method. It is highly recommended to use multiprocessing if Linearization method is used.
plotting_stuff.py   includes functions for plotting and saving plots

# Usage
User input is requested in the file BHT_Unct.py.
- Give a name for the well. This will appear in the final result plot
- Choose a correction method, e.g. Brennand Method (BM), Forward Modelling (FM), Linearization Method (LM). One BHT correction scheme (1BHTM) will be automatically       chosen if only one BHT is given later.
- Specify the depth of the BHT measurements (m TVD). This will appear in the final plot. If the chosen model is Brennand 'BM' the circulation time is estimated from     the defined depth.
- Insert the borehole radius
- Specify the measured BHTs (°C) in BHT_list and the respective shut-in times (h) in the list t_list
- Specify the mud temperature Tmud (°C) if known. When this value is unknown it will be varied between 20 and 80°C
- Specify the quality of input data BHT and shut-in time. 
    if quality_BHT = "low":  uncertainty is +/- 3%
    if quality_BHT = "med":  uncertainty is +/- 2%
    if quality_BHT = "high": uncertainty is +/- 1%
    if quality_ts = "low":   uncertainty is +/- 1h
    if quality_ts = "med":   uncertainty is +/- 30min
    if quality_ts = "high":  uncertainty is +/- 15min    
- Specify the number of samples for Saltelli sampling in the list "sample_variation"

After the script is run, the Sobol results (first order index, second order index and total order index) will be saved into the folder "SobolInformationen"
Under "plots" the main result (density and box plot) will be saved including expected value (median, p50), worst-case (p10) and best-case prediction (p90) as well as the ranges p50-p10 and p90-p50 that can be used to describe the uncertainty of the expected value.
Under "sobol_plots", the user can check the Sobol plots for convergence of the respective model. If the total order index does still change with increasing number of samples, the user should extend the list "sample_variation" and rerun the model.
