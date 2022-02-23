# BHT_uncert
Uncertainty based model to correct Bottom Hole Temperatures of deep wells to static formation temperatures with conventional analytical and numerical correction methods.

# About
This Python tool is designed to correct thermally disturbed Bottom Hole Temperatures (BHTs) to the true formation temperature. BHT data are usually of poor quality, therefore it is important to take into account the uncertainty of the required input parameters to the different models. This tool allows for specifying the quality of input data (uncertainty range), performs a Saltelli sampling (Saltellis extension of Sobol) and runs one of six conventional BHT correction schemes.
The output is provided as density distribution of the result space. p10, p50, p90 or modal value can be used to describe the corrected temperature value with its uncertainty, dependent on the input ranges.

For at least two BHT available, the user can switch between Horner method, Lachenbruch&Brewer, Brennand method, Forward modelling or linearized model.
For only one BHT available, a 1BHT correction scheme after Zekiri () is applied.

# Licence
BHT_uncert is distributed under the GNU GENERAL PUBLIC LICENSE v3.

# Usage
User input is requested in the file BHT_Correction.py. 
- Specify the number of BHTs as n
- Specify the model to use if n ≥ 2 
    'HM' for Horner Plot method
    'LBM' for Lachenbruch & Brewer method
    'BM' for Brennand Plot method
    'FM' for Forward modelling
    'LM' for Linearized method

- Specify the quality of input data. This refers to the shut-in time, as this is the most sensitive parameter.
- Specify the depth of the BHT measurements (TVD)
- Specify the borehole radius (radius)
- Specify all BHTs in the list (BHT_list)
- Specify all BHT related shut-in times in the list (t_list)
