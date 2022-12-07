
## This is the cost function for estimating the index parameter

## Inputs: 

## 1) et            - index parameter in polar coordinates
## 2) datosfda      - an nxm matrix as the response whose each row represents an observation, 
##                    each column represents a quantile, t from [0,1].
## 3) linear_vars   - the names of the p covariates for the linear part of model included in datosx 
## 4) si_vars       - the names of the q covariates for the SI part of model included in datosx
## 5) formula_lv    - a character of length=1, the formula of the covariates in the linear part
## 5) tt            - the equidistant grid on [0,1] of length m.
## 6) sp            - the degree of polynomial considered for spline regression.
## 7) dfs           - degrees of freedom as an alternative to specifying the knots. 
## 8) datosx        - a data frame with n rows whose columns include the covariates for the model
##                    as well as for survey design.

## Output: the mean square prediction error. 

wn_cost <- function(et, datosfda, linear_vars, si_vars, formula_lv, tt, sp, dfs, datosx) {
  
  # read the necessary libraries:
  library("survey")
  library("splines")
  library("fda.usc")
  library("fdadensity")
  
  source("polar2cart.R", local = knitr::knit_global())
  th<- as.matrix(polar2cart(et, 1))
  bs_fun<- as.matrix(datosx[,c(si_vars)])%*% th
  
  bs_data<- data.frame(bs(x=bs_fun, degree=sp, df=dfs))
  colnames(bs_data)<- paste("BS", seq(1:ncol(bs_data)), sep ="")
  datosx<- data.frame(cbind.data.frame(datosx, bs_data))
  
  dat<- cbind.data.frame(datosx, datosfda)
  data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                          weights= ~survey_wt, data = dat, nest = TRUE)
  
  objetofda= fdata(datosfda, argvals = tt)
  formula<- paste(formula_lv, paste(colnames(bs_data),collapse="+"),sep="+")
  source("survey2wassersteinmodel.R", local =knitr::knit_global())
  res=  survey2wassersteinmodel(formula=formula, data_analysis_svy=data_analysis_svy, objetofda = objetofda)
  Yhat<- res$predicciones
  
  #get array of Wasserstein distances between the response and its prediction
  RSS <- sapply(1:nrow(datosfda), function(i) fdadensity:::trapzRcpp(X = tt, Y = (datosfda[i, ] - Yhat[i, ])^2))
  
  w<- data_analysis_svy$variables$survey_wt
  return(sum(w*RSS/sum(w)))
}
