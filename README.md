# Application of Partially Linear Fréchet Single Index model to NHANES data on wearable devices.

## Step 1:

First create a folder in your computer and keep it as the working directory for your RStudio operations. In that folder, download and save the following files from this repository, in no particular order. 

 - ### cart2polar.R
    
            Description: This function computes polar coordinates from the cartesian coordinates in p-dimensional euclidean 
            space for the dimensions p = 2, 3, 4, 5.
            
    Input:
    
            cartesian coordinates in p-dimensional space.
             
    Outputs:  
    
            1) r   : the radius of the polar coordinates.
            2) eta : the vector of polar coordinates with length (p-1).

 - ### polar2cart.R

    Description: 
    
            This is a function to compute the p-dimensional cartesian coordinates from the given polar coordinates in 
            (p-1) dimension and radius r. Computes for p=2,3,4,5. This performs the inverse operation of the 'cart2polar' 
            function above.

    Inputs: 
    
            1) eta : a (p-1) dimensional vector
            2) r   : radius of the polar coordinates (set default, r=1)
    
    
    Output:         
            
            cartesian coordinates in p-dimensional space.
    
  - ### cudratico.R
    
    Description: 
    
            This function projects the model predictions into the space of quantiles, the L^2 - Wasserstein space by solving 
            a quadratic programing problem. 
    
    Inputs:
    
            1. prediciones  - model predictions of the GLM.
            2. cotainferior - lower limit of the quantiles.
            3. cotasuperior - upper limit of the quantiles.
    
    Output:
    
            a matrix of same dimension as the input matrix prediciones, of which the rows are projected onto the space of quantiles.
    
  - ### survey2wassersteinmodel.R

    Description:
    
           This function computes the performance metrics, estimates the beta functional parameters, residuals of the model utilizing the 
           survey design proposed by NHANES. 
           
    Inputs: 
    
           1) formula           : a character representing the formula for regression
           2) data_analysis_svy : a survey GLM object 
           3) objetofda         : the functional data representation of the response quantiles

    Outputs:
   
           1) r^2           : the coefficient of determination (R-squared) of the fitted model.
           2) betaj         : functional beta coefficients of the fitted regression models. 
           3) beta_lcl      : pointwise lower 95% Confidence Interval of the beta parameter
           4) beta_ucl      : pointwise upper 95% Confidence Interval of the beta parameter
           5) predicciones  : the prediction of regression model.
           6) residuos      : the residuals of the fitted model obtained for each point on the grid [0,1]
           7) r2vec         : a vector of length same as number of columns in predicciones 
                              above, containing the multiple R-squared for each percentile.
           
    
  - ### wn_cost.R

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
            7. formula_lv  : a character of length=1, the formula of the covariates in the linear part.
            8. tt          : the equidistant grid on [0,1] of length m.
            9. sp          : the degree of polynomial considered for spline regression.
           10. dfs         : degrees of freedom as an alternative to specifying the knots.

  
    Output:
  
            the mean square prediction error as a function of the Index parameter.
            
  - ### adj_fr_r2.R

    Description: 
   
    Inputs:  
   
           1. qin   :  an nxm matrix whose rows are the quantiles corresponding to a distribution
           2. qpred :  an nxm matrix whose rows are predicted quantiles for the Y above. 
           3. tt    :  an equidistant grid of points on [0,1] of length m for quantiles.
           4. q     :  number of covariates in the model.
   
   
    Output 
    
         A list contanining Fréchet R-squared (Frechet_R2), Adjusted Fréchet R-squared (Adj_Frechet_R2).
        
  - ### PLFSI_model.R

    Description: 
    
         This function fits the partially-linear Fréchet Single Index model to the distributional responses.
         
    Sourcing: 
    
         The files 'wn_cost.R', 'polar2cart.R', 'cart2polar.R', 'survey2wassersteinmodel.R', 'cuadratico.R' 
         have to be sourced prior to running the codes in this script.
         
    Inputs:
    
         1. tt          : length m grid spanning [0, 1], used as grid for quantile functions.
         2. datosfda    : nxm matrix of response quantile functions on grid tt.
         3. si_vars     : a p-vector of variables' names to be considered in the Single Index part.
         4. linear_vars : a q-vector of variables' names to be considered in the linear part.
         5. formula_lv  : a character of length=1, the formula of the covariates in the linear part
         6. nsp         : integer giving the number of starting points in each dimension to be used by optim. A 
                          lattice of points will be created by constructing an equally spaced grid for each of 
                          the (p - 1) hyperspherical coordinates used to represent theta in the optimization. 
                          Default is 3.
         7. L           : a list of integers specifying which starting points to use. If L = 0 (default), all of 
                          the starting points in the lattice will be utilized. Otherwise, L of these will be 
                          chosen by row number. If L = -1, the user will have to input a matrix whose rows are 
                          the starting points.
         8. etaStart    : a matrix with (p-1) columns each row indicating a unique starting value used in optimization 
                          for estimating theta. This is input only if L=-1 
         9. datosx      : the dataset of n whose columns include the covariates, survey variables of the model.
        10. sp          : order of spline.
        11. dfs         : degrees of freedem of the spline.

    
    Outputs: A List with the following elements:
    
         1. thetaHat : length p vector giving the estimated coefficient
         2. fnvalue  : achieved minimum value of the criterion function for estimating theta
         3. etaStart : matrix with (p - 1) columns, each row indicating a unique starting value
                       used in optimization for estimating theta
         4. optInf   : list containing information about optimization routine for each
                       starting point
                       
                       
  - ### Plots.R
  - ### AnalysisAge30to50_with_BMXWAIST_HEI.Rmd
  - ### AnalysisAge30to50_with_BMXWAIST_HEI_wTAC.Rmd



