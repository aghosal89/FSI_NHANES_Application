# Application of Partially Linear Fréchet Single Index model to NHANES data on wearable devices.

1. First create a folder in your computer and keep it as the working directory for your RStudio operations. In that folder, download and save the following files from this repository, in no particular order. 

 - cart2polar.R
    
            Description: This function computes polar coordinates from the cartesian coordinates in p-dimensional euclidean 
            space for the dimensions 2, 3, 4, 5.
            
    Input:
    
            cartesian coordinates in p-dimensional space.
             
    Outputs:  
    
            1) r   : the radius of the polar coordinates.
            2) eta : the vector of polar coordinates with length (p-1).

 - polar2cart.R

    Description: 
    
            This is a function to compute the p-dimensional cartesian coordinates from the given polar coordinates in 
            (p-1) dimension and radius r. Computes for p=2,3,4,5. This performs the reverse operation of the 'cart2polar' 
            function above.

    Inputs: 
    
            1) eta : a (p-1) dimensional vector
            2) r   : radius of the polar coordinates (set default, r=1)
    
    
    Output:         
            
            cartesian coordinates in p-dimensional space.
    
  - wn_cost.R

    Description: 
    
            This is the cost function for estimating the index parameter obtained from equation 9 in the main document.

    Inputs:
    
            1. et          : index parameter in polar coordinates
            2. datosfda    : an nxm matrix as the response whose each row represents an observation, each column 
                             represents a quantile, t from [0,1].
            3. datosx      : a data frame with n rows whose columns include the covariates for the model as well as the
                             variables for survey design.
            5. linear_vars : the names of the p covariates for the linear part of model included in datosx.
            6. si_vars     : the names of the q covariates for the SI part of model included in datosx.
            7. tt          : the equidistant grid on [0,1] of length m.
            8. sp          : the degree of polynomial considered for spline regression.
            9. dfs         : degrees of freedom as an alternative to specifying the knots.

  
     Output:
  
            the mean square prediction error as a function of the Index parameter.
            
  - survey2wassersteinmodel.R
  - cudratico.R
  - adj_fr_r2.R
  - PLFSI_model.R
  - spline_variation.R
  - AnalysisAge30to50_with_BMXWAIST_HEI.Rmd
  - AnalysisAge30to50_with_BMXWAIST_HEI_wTAC.Rmd



