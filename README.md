# BHT_Unct
Uncertainty based model to correct Bottom Hole Temperatures of deep wells to static formation temperatures with conventional analytical and numerical correction methods.

# About
This Python tool is designed to estimate the formation temperature from thermally disturbed Bottom Hole Temperatures (BHTs). BHT data are usually of poor quality, therefore it is important to take into account the uncertainty of the required input parameters to the different models that can be qpplied for correction. This tool allows for specifying the quality of input data (uncertainty range), performing a Saltelli sampling (Saltellis extension of Sobol) and running one of six conventional BHT correction schemes.
The output is provided as density distribution of the result space. p10, p50, p90 or modal value can be used to describe the corrected temperature value with its uncertainty, dependent on the input ranges.


# Licence
BHT_Unct is distributed under the GNU GENERAL PUBLIC LICENSE v3.

# Usage
User input is requested in the file BHT_Correction.py. 
- Give a name for the well. This will appear in the final result plot
- Specify the number of BHTs as n
- if n = 1 1BHTM is automatically chosen as model
- Specify the model to use if n ≥ 2 
    'HM' for Horner Plot method
    'LBM' for Lachenbruch & Brewer method
    'BM' for Brennand Plot method
    'FM' for Forward modelling
    'LM' for Linearized method

- Specify the quality of input data. This refers to the shut-in time, as this is the most sensitive parameter. The uncertainty can be edited for each defined quality
- Specify the depth of the BHT measurements (TVD). This will appear in the final plot. If the chosen model is Horner 'HM' or Brennand 'BM' the circulation time is estimated from the defined depth.
- Specify the borehole radius 'radius'
- Specify all BHTs in the list 'BHT_list'
- Specify all BHT related shut-in times in the list 't_list'
