# Functions used for model

First create a folder in your computer and keep it as the working directory for your RStudio operations. In that folder, download and save the files from this repository. The following are the descriptions for running the main model Partially Linear Fréchet Single Index Model (PL-FSI) and performing other tasks.

### cart2polar.R

##### Description: 

 This function computes polar coordinates from the cartesian coordinates in p-dimensional Euclidean space for the dimensions p = 2, 3, 4, 5.

##### Input: 

 - cartesian coordinates in p-dimensional space.

##### Outputs:

  - r      :  the radius of the polar coordinates.
  - eta    :  the vector of polar coordinates with length (p-1).


### polar2cart.R

##### Description: 

This is a function to compute the p-dimensional cartesian coordinates from the given polar coordinates in (p-1) dimension and radius r. Computes for p=2,3,4,5. This performs the inverse operation of the 'cart2polar' function above.

Inputs:
•	eta : a (p-1) dimensional vector.
•	r      : radius of the polar coordinates (set default, r=1).

Output: cartesian coordinates in p-dimensional space.

3.	cudratico.R

Description:     This function projects the model predictions into the space of quantiles, the L^2 - Wasserstein space by solving 
    a quadratic programing problem.

Inputs:         
•	prediciones    : model predictions of the GLM.
•	cotainferior    : lower limit of the quantiles.
•	cotasuperior : upper limit of the quantiles.

Output: a matrix of same dimension as the input matrix ‘prediciones’, of which the rows are projected onto the space of quantiles.
4.	survey2wassersteinmodel_2.R
Description: This function computes the performance metrics, estimates the beta functional parameters, residuals of the model utilizing the survey design proposed by NHANES.
Inputs:
•	formula                        : a character representing the formula for regression.
•	data_analysis_svy  : a survey GLM object. 
•	objetofda                    : the functional data representation of the response quantiles.
•	xout		    : a matrix of order n1 x p whose columns are the same as those of the input covariate matrix in ‘data_analysis_svy’ above, the number of rows n1 being the observation held out, to perform prediction.
Outputs: 
•	r^2                               : the coefficient of determination (R-squared) of the fitted model.
•	betaj                           : functional beta coefficients of the fitted regression models. 
•	Predicciones_In    : the prediction of regression model before projection into the 2-Wasserstein space of distributions.
•	Projection_In          : the projections into the 2-Wasserstein space of distributions.
•	Predictions_Out   : the out-of-sample prediction of regression model after projection to 2-Wasserstein space.
•	residuos                   : the residuals of the fitted model obtained for each point on the grid [0,1].
•	r2vec                          : a vector of length same as number of columns in ‘predicciones’ above, containing the multiple R-squared for each percentile.

5.	wn_cost.R
Description: This is the cost function for estimating the index parameter obtained from equation 9 in the main document.
Inputs: 
•	et             : index parameter in polar coordinates
•	datosfda     : an nxm matrix as the response whose each row represents an observation, each column represents a quantile, t from [0,1].
•	datosx          : a data frame with n rows whose columns include the covariates for the model as well as the variables for survey design.
•	linear_vars : the names of the p covariates for the linear part of model included in datosx.
•	si_vars          : the names of the q covariates for the SI part of model included in datosx.
•	formula_lv  : a character of length=1, the formula of the covariates in the linear part.
•	tt                      : the equidistant grid on [0,1] of length m.
•	sp                    : the degree of polynomial considered for spline regression.
•	dfs                   : degrees of freedom as an alternative to specifying the knots.
Output: the mean square prediction error as a function of the Index parameter.
6.	adj_fr_r2.R
Description: This function is to compute the Fréchet R-squared and the adjusted Fréchet R-squared from the models.
Inputs:
•	qin                             :  an nxm matrix whose rows are the quantiles corresponding to a distribution
•	qpred                       :  an nxm matrix whose rows are predicted quantiles for the Y above. 
•	tt                                 :  an equidistant grid of points on [0,1] of length m for quantiles.
•	q                                 :  number of covariates in the model.
•	Survey_weights   :  Survey weights used for the individuals in NHANES study.
Outputs: 
•	Fréchet R-squared (Frechet_R2).
•	Adjusted Fréchet R-squared (Adj_Frechet_R2).

7.	PLFSI_model.R
Description: This function fits the PL-FSI model to the distributional responses.
Sourcing: The files 
•	wn_cost.R
•	polar2cart.R
•	cart2polar.R
•	survey2wassersteinmodel_2.R
•	cuadratico.R
must be sourced prior to running the codes in this script. The descriptions of these functions are as mentioned above.
Inputs:           
•	tt                     : length m grid spanning [0, 1], used as grid for quantile functions.
•	datosfda     : nxm matrix of response quantile functions on grid tt.
•	si_vars         : a p-vector of variables' names to be considered in the Single Index part.
•	linear_vars : a q-vector of variables' names to be considered in the linear part.
•	formula_lv  : a character of length=1, the formula of the covariates in the linear part
•	nsp                 : integer giving the number of starting points in each dimension to be used by optim. A lattice of points will be created by constructing an equally spaced grid for each of the (p - 1) hyperspherical coordinates used to represent theta in the optimization. Default is 3.
•	 L                     : a list of integers specifying which starting points to use. If L = 0 (default), all of the starting points in the lattice will be utilized. Otherwise, L of these will be chosen by row number. If L = -1, the user will have to input a matrix whose rows are the starting points.
•	etaStart        : a matrix with (p-1) columns each row indicating a unique starting value used in optimization for estimating theta. This is input only if L=-1. 
•	datosx           : the dataset of n whose columns include the covariates, survey variables of the model.
•	sp                    : order of spline.
•	dfs                   : degrees of freedem of the spline.
Outputs:
•	thetaHat  : length p vector giving the estimated coefficient
•	fnvalue     : achieved minimum value of the criterion function for estimating theta
•	etaStart    : matrix with (p - 1) columns, each row indicating a unique starting value used in optimization for estimating theta
•	optInf         : list containing information about optimization routine for each starting point.


PLFSI regression model

To run the main PLFSI model run the “AnalysisAge20to80_PLFSI_noTAC.R”, which computes the model with the participants in the age range 20 - 80 years and BMI range 18.5 - 40. For the hybrid nature of the model, we considered the covariates Age and BMI in the non-linear part, and the covariates Sex, Ethnicity, HEI in the linear part, also considered the interaction between the covariates Sex and Ethnicity.

Sources: 

•	adj_fr_r2.R
•	PLFSI_model.R
•	survey2wassersteinmodel_2.R
•	wn_cost.R

Input:   

       The dataset ‘datosalex(1).csv’ file can be obtained from the link.

Outputs:

•	‘Theta_Hat.csv’: the .csv file contaning the 2 x 1 estimated index parameter for the Single Index part of the model.
•	‘Output_Age20to80_noTAC_residuals.csv’: The model residuals obtained for every order of quantile t.  
•	‘Predictions_before_projection.csv’: The model predictions obtained for every order of quantile t. 
•	‘Predictions_after_projection.csv’: The prediction obtained in ‘Predictions_before_projection.csv’ and projected into the 2-Wasserstein space.
•	‘Output_Age20to80_noTAC_betas.csv’: Estimated beta coefficients for every order of quantile t. 
•	‘Output_Age20to80_noTAC_rnames.csv’: The estimated parameters in the model.
•	‘Output_Age20to80_noTAC_betas_UCL.csv’: pointwise 95% Upper Confidence limits of estimated parameters.
•	‘Output_Age20to80_noTAC_betas_LCL.csv’: pointwise 95% Lower Confidence limits of estimated parameters.
•	‘Output_Age20to80_noTAC_betaeffects.csv’: pointwise estimate of the effects.
•	‘Output_Age20to80_noTAC_betaeffects_UCL.csv’: pointwise 95% Upper Confidence limits of estimated effects.
•	‘Output_Age20to80_noTAC_betaeffects_LCL.csv’: pointwise 95% Lower Confidence limits of estimated effects.


Comparison of the Out-Sample performances of the PLFSI, PLF and GF regression models

To compare the PLF (Partially Linear Fréchet), GF (Global Fréchet) and the PLFSI (Partially Linear Fréchet Single Index) models w.r.t. the out-sample performances, run the following file:

“Combined_PLFSI_PLF_GF_performance_comparison.R”

First, we split the data into in-sample and out-sample parts so that 462 observations are held out chosen randomly without replacement, while the rest 4154 are considered for the in-sample model fitting of the data. We created 40 such splits, while fitting each of the three models PLF, PLFSI, GF in the training split and predicting the activity distributions for the held-out data. The Random split of the data is computed in the file:

“Sample_data.R”

which outputs file “Sample_data.csv” containing a 40 x 462 matrix whose entries represent the observation number that is to be held out into the out-sample for the corresponding split. 

Sources: 

•	adj_fr_r2.R
•	PLFSI_model.R
•	survey2wassersteinmodel_2.R
•	wn_cost.R

Inputs:

•	The ‘Sample_data.csv’ can be found at the link.
•	The dataset ‘datosalex(1).csv’ file can be obtained from the link mentioned above.

Outputs:

•	“Theta_hat_PLFSI.csv”: Contains a matrix of order 40 x 2, whose rows are the estimated Single-Index parameters for each of the 40 data splits.

•	“Theta_convergence.csv”: Contains a matrix of order 40 x 4, whose rows are the numbers indicating convergence and successful estimation of the SI parameter for each of the 4 starting values. Here ‘0’ would mean convergence of the estimation algorithm.

•	“MSPEs_PLFSI.csv”: Contains a vector of length 40 that contains the MSPE of out-sample prediction of the PLFSI model for each of the 40 splits.

•	“MSPEs_GF.csv”: Contains a vector of length 40 of the MSPE of the out-sample prediction of the Global Fréchet regression model for each of the 40 splits.

•	“MSPE_PLF.csv”: Contains a vector of 40 MSPE values corresponding to the out-sample prediction of the PLF model corresponding to the 40 splits.




Evaluation of the PLFSI model through Bootstrap

We performed 100 bootstrap simulations of the PLFSI model to re-estimate the model parameters and their linear combinations. To perform this exercise run the file 

“PLFSI_Bootstrap.R”.

The bootstrap survey weight multipliers are created in the script:

“Survey_data_preparation.R”,

which takes as input the “datosalex(1).csv” (link above) and as output provides “Boot_survey_data.csv” available in the link, which contains the survey weight multipliers for each observation in the study.

For the “PLFSI_bootstrap.R” script,

Sources: 

•	adj_fr_r2.R
•	PLFSI_model.R
•	survey2wassersteinmodel_2.R
•	wn_cost.R

Input:   

•	The dataset “datosalex(1).csv” (link above).
•	“Boot_survey_data.csv” (link above).
•	“Output_Age20to80_noTAC_betaeffects.csv” file as output from the file “AnalysisAge20to80_with_HEI_noTAC.Rmd”. 

Outputs:

•	“Bootstrap_lower95_ConfidenceBound.csv”: the 95% lower confidence bounds for the estimated model parameter effects. 

•	“Bootstrap_upper95_ConfidenceBound.csv”: the 95% upper confidence bounds for the estimated model parameter effects.





Plots and figures in the document

To run all the plots and figures in the document run the
