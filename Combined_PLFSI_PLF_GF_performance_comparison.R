
#########################################################
# Comparison of the PLF, PL-FSI, GF regression models
#########################################################

# Set your R directory
setwd("~/Documents/2024_09_21_Aritra_codeForMainAnalysis")

### Read the libraries 
library("optimParallel")
library("tidyverse")
library("survey")
library("fda.usc")
library("sandwich")
library("compiler")
library("sampling")
library("statip")
library("fdadensity")
library("Rcpp")
library("osqp")
library("Matrix")
library("splines")
library("statip")

# Read the necessary functions to obtain the results
source("adj_fr_r2.R")
source("PLFSI_model.R")
source("survey2wassersteinmodel_2.R")
source("wn_cost.R")
source("polar2cart.R")
source("cuadratico.R")

# read the dataset
datos= read.csv("datosalex(1).csv")

# column indices for distributional representations
indices= round(seq(188,687, length=500))
colnames(datos)[indices]= paste("X", 1:500, sep="")
datos= datos[, c(1:187, indices)]
tt= seq(0,1, length=500)

# obtain integrated R-squared for the quantile range [0, 0.97].
m97 <- floor(which(tt>=0.97))[1]
tt97 <- tt[1:m97]

# subset data for age 20-80 years
datos <- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  

# subset data for BMI in the range 18.5 - 40
datos <- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  

# create the response data as functional data
datosfda <- datos[,188:687] 
datosfda <- as.matrix(datosfda)

datosfda[which(datosfda>280)] <- datosfda[which(datosfda>280)-1]+0.5
datos[,188:687] <- datosfda

# read the training-testing splits for the response and covariate data
Sample_data <- read.csv("Sample_data.csv", header = TRUE)[,-1]

si_vars <- c("BMXBMI", "RIDAGEYR.y")  # covariates for the Single Index part
linear_vars<- c("RIAGENDR", "RIDRETH3", "HEI")  # covariates for Linear part
linear_vars_numerical <- "HEI"  # numeric covariates for Linear part
linear_vars_categorical <- c("RIAGENDR", "RIDRETH3")  # categorical covariates for Linear part

# the survey variables
datosx <- cbind.data.frame(scale(datos[,si_vars]), datos[,linear_vars_categorical], 
                           HEI=scale(datos[,linear_vars_numerical]), 
                           survey_wt= datos[,"wtmec4yr_adj_norm"], 
                           survey_id=datos[,"SDMVPSU"], 
                           survey_strata= datos[,"SDMVSTRA"])

datosx$RIAGENDR <- as.factor(datosx$RIAGENDR)
datosx$RIDRETH3 <- as.factor(datosx$RIDRETH3)
datosx$survey_strata<- as.factor(datosx$survey_strata)
datosx$survey_id<- as.factor(datosx$survey_id)

## compute the regression formula for the Global Frechet model
formula_lv <- paste(c(paste(linear_vars_categorical, collapse="*"),
                      linear_vars_numerical), collapse = "+")
formula_gfr <- paste(formula_lv, paste(paste(si_vars),collapse="+"),sep="+")


######
# Out-of-sample performance evaluation for the PL-FSI model
######

# spline parameters
sp<- 4
dfs<- 9

# Traditional R-squared for each data-split
R2_plfsi <- rep(NA, len= nrow(Sample_data))

# Frechet R-squared for each data-split
Frechet_R2_plfsi <- rep(NA, len= nrow(Sample_data))  

# Adjusted Frechet R-squared for each data split
Adj_Frechet_R2_plfsi <- rep(NA, len= nrow(Sample_data))

# MSPE for each data split to test the Out-of-sample prediction
MSPEs_plfsi <- rep(NA, len= nrow(Sample_data))

# store the estimated single index parameter for each data split
Theta_hat <- matrix(NA, nrow = nrow(Sample_data), ncol = 2)

# store the convergence record for each estimation process
theta_convergence<- matrix(NA, nrow =nrow(Sample_data), ncol=4)

# Run the PL-FSI model to run for each data split considered
for (i1 in 1:nrow(Sample_data)) {
  
  print(i1)
  s1<- as.integer(Sample_data[i1,])
  # create the in-sample data for regression 
  datosx_in<- datosx[-s1, ] # prediction
  datosfda_in<- datosfda[-s1, ] # response
  
  # create the out-of-sample data for regression
  yout<- datosfda[s1, ] # response 
  datosx_out<- datosx[s1,] # predictor variables
  
  # Estimate the model parameters
  lmod<- PLFSI_model(si_vars = si_vars, linear_vars = linear_vars, 
                     datosfda = datosfda_in, tt =tt, nsp = 4, L = 0, 
                     etaStart = NULL, sp=4, dfs=9, formula_lv=formula_lv, 
                     datosx= datosx_in)
  
  # report the Frechet R^2
  Frechet_R2_plfsi[i1] <- lmod$Frechet_R2
  
  # report the Adjusted Frechet R^2
  Adj_Frechet_R2_plfsi[i1] <- lmod$Adj_Frechet_R2
  
  # report the R^2
  R2_plfsi[i1]<- lmod$R2
  
  # estimate of index parameter theta
  Theta_hat[i1,]<- lmod$thetaHat
  
  # To check if all points converged, 0 indicates convergence
  theta_convergence[i1,] <- lmod$converge
  
  # obtain the residuals with patient id
  res_id<- cbind.data.frame(SEQN =datos$SEQN[-s1], lmod$residuals)
  colnames(res_id)<- c("SEQN", paste("X", c(1:ncol(lmod$residuals)), sep = ""))
  
  # obtain the global Frechet model
  objetofda_in<- fdata(datosfda_in, argvals = tt)
  
  bs_fun <- as.matrix(datosx[,si_vars])%*% Theta_hat[i1,]
  
  bs_data <- data.frame(bs(x=bs_fun, degree=sp, df=dfs))
  colnames(bs_data) <- paste("BS", seq(1:ncol(bs_data)), sep = "")
  bs_data_in <- bs_data[-s1,]
  bs_data_out <- bs_data[s1,]
  data_temp <- cbind.data.frame(datosx_in, bs_data_in)
  
  xout <- cbind(datosx_out, bs_data_out)
  
  ## modifying the inputs of the algorithm
  
  formula<- paste(formula_lv, paste(colnames(bs_data),collapse="+"), sep="+")
  data_temp<- cbind.data.frame(data_temp, datosfda_in)
  data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                                 weights= ~survey_wt, data = data_temp, nest = TRUE)
  
  res<- survey2wassersteinmodel_2(formula=formula, data_analysis_svy, 
                                  objetofda=objetofda_in,
                                  xout = xout)
  
  MSPEs_plfsi[i1]<- mean(sapply(1: nrow(res$predicciones_Out), function(i) {
    fdadensity:::trapzRcpp(X=tt97, Y=(res$predicciones_Out[i,1:m97] -
                                        yout[i,1:m97])^2) 
  } ))
}

write.csv(Theta_hat, "Theta_hat.csv")
write.csv(theta_convergence, "Theta_convergence.csv")
write.csv(R2_plfsi, "Rsquared_PLFSI.csv")
write.csv(Frechet_R2_plfsi, "Frechet_Rsquared_PLFSI.csv")
write.csv(Adj_Frechet_R2_plfsi, "Adj_FRsquared_PLFSI.csv")
write.csv(MSPEs_plfsi, "MSPEs_PLFSI.csv")


################################################################
# Run the out-of-sample performance for the PLF regression model
################################################################

# for the first variable in the single index part of the model 
bs_data1<- data.frame(bs(x=datosx[,si_vars[1]], degree=sp, df=dfs))
colnames(bs_data1) <- paste("BS", "_1_", seq(1:ncol(bs_data1)), sep ="")

# for the second variable in the single index part of the model
bs_data2<- data.frame(bs(x=datosx[,si_vars[2]], degree=sp, df=dfs))
colnames(bs_data2) <- paste("BS", "_2_", seq(1:ncol(bs_data2)), sep = "")

# combine the data for further modeling exercise
datosx_plf<- cbind.data.frame(datosx, bs_data1, bs_data2)

# the main computation part
R2<- rep(NA, len= nrow(Sample_data))
Frechet_R2<- rep(NA, len= nrow(Sample_data))
Adj_Frechet_R2<- rep(NA, len= nrow(Sample_data))
MSPEs_plf<- rep(NA, len= nrow(Sample_data))

for (i2 in 1:nrow(Sample_data)) {
  
  sdat<- as.numeric(Sample_data[i2,])
  yout <- datosfda[sdat,]
  
  datosfda_temp <- datosfda[-sdat, ]
  objetofda_temp <- fdata(datosfda_temp,argvals = tt)
  
  # obtain the global Frechet model
  
  ## compute the regression formula for the Global Frechet model
  formula_lv <- paste(c(paste(linear_vars_categorical, collapse="*"),
                        linear_vars_numerical), 
                        collapse = "+")
  
  xout_temp <- datosx_plf[sdat,]
  datosx_temp <- datosx_plf[-sdat,]
  
  si_varnames <- c(colnames(bs_data1), colnames(bs_data2))
  
  formula_gfr <- paste(formula_lv, paste(si_varnames, collapse ="+"), sep="+")
  data_temp <- cbind.data.frame(datosx_temp, datosfda_temp)
  data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                                 weights= ~survey_wt, data = data_temp, 
                                 nest = TRUE)
  
  # run the Global Frechet model with the survey weights
  res <- survey2wassersteinmodel_2(formula=formula_gfr, data_analysis_svy, 
                                   objetofda=objetofda_temp, xout=xout_temp)   
  
  # obtain R-squared over all the quantiles 
  R2[i2] <- res$r2
  
  # Frechet R-squared of the GF model
  frsq<- adj_fr_r2(qin=datosfda_temp, qpred=res$predicciones_I, tt=tt, 
                   q=nrow(res$betaj), survey_weights=datosx_temp$survey_wt)
  
  Frechet_R2[i2] <- frsq$Frechet_R2
  
  # Adjusted Frechet R-squared of GF model
  Adj_Frechet_R2[i2] <- frsq$Adj_Frechet_R2
  
  ## Out-of-sample performance 
  
  # Mean Square Prediction Error (MSPE) for the out-of-sample prediction
  MSPEs_plf[i2]<- mean(sapply(1: nrow(res$predicciones_Out), function(i) {
    fdadensity:::trapzRcpp(X=tt97, Y=(res$predicciones_Out[i,1:m97] -
                                        yout[i,1:m97])^2) 
  } ))
  
}

write.csv(MSPEs, "MSPEs_PLF.csv")
write.csv(Frechet_R2, "Frechet_R2_PLF.csv")
write.csv(Adj_Frechet_R2, "Adj_Frechet_R2_PLF.csv")
write.csv(R2, "Traditional_R2_PLF.csv")

######
# Out-of-sample performance evaluation for the GF model
######

# Traditional R-squared for each data-split
R2_gf<- rep(NA, len= nrow(Sample_data))

# Frechet R-squared for each data-split
Frechet_R2_gf<- rep(NA, len= nrow(Sample_data))  

# Adjusted Frechet R-squared for each data split
Adj_Frechet_R2_gf<- rep(NA, len= nrow(Sample_data))

# MSPE for each data split to test the Out-of-sample prediction
MSPEs_gf<- rep(NA, len= nrow(Sample_data))

si_vars <- c("BMXBMI","RIDAGEYR.y")  # covariates for Single Index part
linear_vars<- c("RIAGENDR","RIDRETH3","HEI")  # covariates for Linear part
linear_vars_numerical <- "HEI"  # numeric covariates for Linear part
linear_vars_categorical <- c("RIAGENDR","RIDRETH3")  # categorical covariates for Linear part

# the survey variables
datosx <- cbind.data.frame(scale(datos[,si_vars]), datos[,linear_vars_categorical], 
                           HEI=scale(datos[,linear_vars_numerical]), 
                           survey_wt= datos[,"wtmec4yr_adj_norm"], 
                           survey_id=datos[,"SDMVPSU"], 
                           survey_strata= datos[,"SDMVSTRA"])

datosx$RIAGENDR <- as.factor(datosx$RIAGENDR)
datosx$RIDRETH3 <- as.factor(datosx$RIDRETH3)
datosx$survey_strata<- as.factor(datosx$survey_strata)
datosx$survey_id<- as.factor(datosx$survey_id)

# In-sample and out-of sample performance of the 
# Global Fréchet regression model 

for (i3 in 1:nrow(Sample_data)) {
  
  s1<- as.integer(Sample_data[i3,])
  datosx_in<- datosx[-s1, ]
  datosfda_in<- datosfda[-s1, ]
  
  yout<- datosfda[s1, ]
  datosx_out<- datosx[s1,]
  
  # obtain the global Frechet model
  objetofda_in <- fdata(datosfda_in, argvals = tt)
  
  ## compute the regression formula for the Global Frechet model
  formula_lv <- paste(c(paste(linear_vars_categorical, collapse="*"),
                        linear_vars_numerical), collapse = "+")
  formula_gfr <- paste(formula_lv, paste(paste(si_vars),collapse="+"),sep="+")
  data_temp <- cbind.data.frame(datosx_in, datosfda_in)
  
  data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                                 weights= ~survey_wt, data = data_temp, nest = TRUE)
  
  # run the Global Frechet model with the survey weights
  res <- survey2wassersteinmodel_2(formula=formula_gfr, data_analysis_svy, 
                                   objetofda=objetofda_in,
                                   xout= datosx_out)   
  
  # obtain R-squared over all the quantiles 
  R2_gf[i3] <- res$r2
  
  # Frechet R-squared of the GF model for the corresponding data split
  frsq<- adj_fr_r2(qin=datosfda_in, qpred=res$predicciones_I, tt=tt, q=nrow(res$betaj), 
                   survey_weights=datosx_in$survey_wt)
  
  Frechet_R2_gf[i3] <- frsq$Frechet_R2
  
  # Adjusted Frechet R-squared of GF model for the corresponding data split
  Adj_Frechet_R2_gf[i3] <- frsq$Adj_Frechet_R2
  
  MSPEs_gf[i3]<- mean(sapply(1: nrow(res$predicciones_Out), function(i) {
    fdadensity:::trapzRcpp(X=tt97, Y=(res$predicciones_Out[i,1:m97] -
                                        yout[i,1:m97])^2) 
  } ))
  
}

write.csv(R2_gf, "Rsquared_GF.csv")
write.csv(MSPEs_gf, "MSPEs_GF.csv")
write.csv(Frechet_R2_gf, "Frechet_Rsquared_GF.csv")
write.csv(Adj_Frechet_R2_gf, "Adj_Frechet_Rsquared_GF.csv")

