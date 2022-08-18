# Application of Partially Linear Fréchet Single Index model to NHANES data on wearable devices.

1. First create a folder in your computer and keep it as the working directory for your RStudio operations. In that folder, download and save the following files from this repository, in no particular order. 

  - cart2polar.R
    [Description: It computes the the polar coordinates of the corresponding cartesian coordinates for the dimensions 2, 3, 4, 5]
  - polar2cart.R

    [Description: This is a function to compute the p-dimensional cartesian coordinates from the given polar coordinates in (p-1) dimension and radius.
    Inputs: 1) eta  - a (p-1) dimensional vector
            2) r    - radius of the polar coordinates (set default, r=1)
    Output:         - cartesian coordinates in p-dimension.]
    
  - wn_cost.R
    [Description: This file computes the Wn function in the equation 9 in the main document]
  - survey2wassersteinmodel.R
  - cudratico.R
  - adj_fr_r2.R
  - PLFSI_model.R
  - spline_variation.R
  - AnalysisAge30to50_with_BMXWAIST_HEI.Rmd
  - AnalysisAge30to50_with_BMXWAIST_HEI_wTAC.Rmd



