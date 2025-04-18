
########
# Analysis of the Physical activity data for age range 20-80 years with 
# covariates BMXBMI, RIDAGEYR.y in Single Index; RIAGENDR, RIDRETH3, HEI in linear part
########

# Set your R directory


### Read the libraries
library("optimParallel")
library("tidyverse")
library("survey")
library("fda.usc")
library("sandwich")
library("compiler")
library("sampling")
library("osqp")
library('Rcpp')
library('numbers')
library('splines')
library('fdadensity')
library('Matrix')

# source the necessary functions
source("adj_fr_r2.R")
source("PLFSI_model.R")
source("survey2wassersteinmodel_2.R")
source("wn_cost.R")
source("cuadratico.R")
source("polar2cart.R")

# read the dataset
datos= read.csv("datosalex(1).csv")

# column indices for distributional representations
indices= round(seq(188,687, length=500))
colnames(datos)[indices]= paste("X", 1:500, sep="")
datos= datos[, c(1:187, indices)]

# subset data for age 20-80 years
datos<- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  

# subset data for BMI in the range 18.5 - 40
datos<- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  

# create the response data as functional data
datosfda= datos[,188:687] 
datosfda= as.matrix(datosfda)

datosfda[which(datosfda>280)]= datosfda[which(datosfda>280)-1]+0.5
datos[,188:687]= datosfda

#relevel(data_analysis_svy$variables$RIDRETH3,3) 

# read the sampling data for training-testing split

si_vars <- c("BMXBMI","RIDAGEYR.y")  # covariates for Single Index part
linear_vars<- c("RIAGENDR","RIDRETH3","HEI")  # covariates for Linear part
linear_vars_numerical <- "HEI"  # numeric covariates for Linear part
linear_vars_categorical <- c("RIAGENDR","RIDRETH3")  # categorical covariates for Linear part

# the survey variables
datosx <- cbind.data.frame(scale(datos[,si_vars]), datos[,linear_vars_categorical], HEI=scale(datos[,linear_vars_numerical]), survey_wt= datos[,"wtmec4yr_adj_norm"], survey_id=datos[,"SDMVPSU"], survey_strata= datos[,"SDMVSTRA"])

datosx$RIAGENDR<- as.factor(datosx$RIAGENDR)
datosx$RIDRETH3<- as.factor(datosx$RIDRETH3)
datosx$survey_strara<- as.factor(datosx$survey_strata)
datosx$survey_id<- as.factor(datosx$survey_id)

# obtain the global Frechet model
tt= seq(0,1,length=500)
objetofda= fdata(datosfda,argvals = tt)

## compute the regression formula for the Global Frechet model
formula_lv<- paste(c(paste(linear_vars_categorical, collapse="*"),linear_vars_numerical), collapse = "+")
formula_gfr<- paste(formula_lv, paste(paste(si_vars),collapse="+"),sep="+")
data_temp <- cbind.data.frame(datosx, datosfda)
data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                               weights= ~survey_wt, data = data_temp, nest = TRUE)

# run the Global Frechet model with the survey weights
res<- survey2wassersteinmodel_2(formula=formula_gfr, data_analysis_svy, objetofda=objetofda)   

# obtain R-squared over all the quantiles 
res$r2

# obtain integrated R-squared across the orders of quantile (tt) in the range [0, 0.97].
m97 <- floor(which(tt>=0.97))[1]
tt97 <- tt[1:m97]
fdadensity:::trapzRcpp(X=tt97, Y=res$R2_vector[1:m97])

# Frechet R-squared of the GF model
frsq<- adj_fr_r2(qin=datosfda, qpred=res$predicciones, tt=tt, q=nrow(res$betaj), survey_weights=datosx$survey_wt)
frsq$Frechet_R2

# Adjusted Frechet R-squared of GF model
frsq$Adj_Frechet_R2

# run the Partially Linear Frechet Single-Index regression model
lmod<- PLFSI_model(si_vars = si_vars, linear_vars = linear_vars, datosfda = datosfda, tt =tt, 
                   nsp = 4, L = 0, etaStart = NULL, sp=4, dfs=9, formula_lv=formula_lv, 
                   datosx= datosx)

# report the Frechet R^2
lmod$Frechet_R2

# report the Adjusted Frechet R^2
lmod$Adj_Frechet_R2

# report the R^2
lmod$R2

# obtain integrated R-squared for the quantile range [0, 0.97].
fdadensity:::trapzRcpp(X=tt97, Y=lmod$R2_vector[1:m97])

# estimate of index parameter theta
lmod$thetaHat
write.csv(lmod$thetaHat, "Theta_Hat.csv")

# To check if all points converged, 0 indicates convergence
lmod$converge

# obtain the residuals with patient id
res_id<- cbind.data.frame(SEQN =datos$SEQN, lmod$residuals)
colnames(res_id)<- c("SEQN", paste("X", c(1:ncol(lmod$residuals)), sep = ""))

# save the residuals for further analysis
write.csv(res_id, "Output_Age20to80_noTAC_residuals.csv") 
write.csv(lmod$predictions, "Output_Age20to80_noTAC_predictions.csv")

# save the beta estimates for plotting later
betas  <- lmod$betaj
rnames <- rownames(betas)

write.csv(betas, "Output_Age20to80_noTAC_betas.csv" )
write.csv(rnames, "Output_Age20to80_noTAC_rnames.csv")

## In-sample performance evaluation

# using the estimated index parameter from the prior chunk

# the survey variables
datosx <- cbind.data.frame(scale(datos[,si_vars]), datos[,linear_vars_categorical], HEI=scale(datos[,linear_vars_numerical]), survey_wt= datos[,"wtmec4yr_adj_norm"], survey_id=datos[,"SDMVPSU"], survey_strata= datos[,"SDMVSTRA"])

datosx$RIAGENDR<- as.factor(datosx$RIAGENDR)
datosx$RIDRETH3<- as.factor(datosx$RIDRETH3)
datosx$survey_strara<- as.factor(datosx$survey_strata)
datosx$survey_id<- as.factor(datosx$survey_id)

bs_fun <- as.matrix(datosx[,si_vars])%*%(lmod$thetaHat)

# obtain the global Frechet model

tt= seq(0,1,length=500)

sp<- 4  # order of B-spline bases
dfs<- 9 # degrees of freedom 

bs_data<- data.frame(bs(x=bs_fun, degree=sp, df=dfs))
colnames(bs_data) <- paste("BS", seq(1:ncol(bs_data)), sep = "")
data_temp<- cbind.data.frame(datosx, bs_data)

## modifying the inputs of the algorithm

formula<- paste(formula_lv, paste(colnames(bs_data),collapse="+"),sep="+")
data_temp<- cbind.data.frame(data_temp, datosfda)
data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                               weights= ~survey_wt, data = data_temp, nest = TRUE)

res<- survey2wassersteinmodel_2(formula=formula, data_analysis_svy, objetofda=objetofda)

# Save the predictions before and after projection to the the 2-Wasserstein space
# before projection
write.csv(res$predicciones_In, "Predictions_before_projection.csv") 

# after projection
write.csv(res$projection_In, "Predictions_after_projection.csv") 

q_n<- nrow(res$betaj) 

# obtain R-squared over all the quantiles 
res$r2

# obtain integrated R-squared across the orders of quantile (tt) in the range [0, 0.97].
fdadensity:::trapzRcpp(X=tt97, Y=res$R2_vector[1:m97])

# Frechet R-squared of the GF model
frsq<- adj_fr_r2(qin= as.matrix(datosfda), qpred=res$predicciones, tt=tt, q=q_n, 
                 survey_weights= datosx$survey_wt)
frsq$Frechet_R2

# Adjusted Frechet R-squared of GF model
frsq$Adj_Frechet_R2


### Computation of the 95% Confidence Intervals for the effects

formulaaux= paste("X", 1, sep="")
formulaaux= paste(formulaaux, "~", sep="")

formulafinal= paste(formulaaux, formula,sep="")
formulafinal= as.formula(formulafinal)  

m2 = svyglm(formulafinal, design=data_analysis_svy, family=stats::gaussian())

# the following function creates a vector of 0s and 1s with 1s specified at the 

vec1s<- function(vec_pos, lenv) {
  vec_temp<- rep(0, lenv)
  vec_temp[vec_pos]<- 1
  return(vec_temp)
}

l2<- length(m2$coefficients)

# define the matrix for added effects 
C= matrix(c(vec1s(c(1),l2),             # male, mexican american
            vec1s(c(1,2),l2),           # female, mexican american
            vec1s(c(1,3),l2),           # male, other hispanic
            vec1s(c(1,2,3,18),l2),      # female, other hispanic 
            vec1s(c(1,4),l2),           # male, non-hispanic white 
            vec1s(c(1,2,4,19),l2),      # female, non-hispanic white
            vec1s(c(1,5),l2),           # male, non-hispanic black
            vec1s(c(1,2,5,20),l2),      # female, non-hispanic black
            vec1s(c(1,6),l2),           # male, non-hispanic asian
            vec1s(c(1,2,6,21),l2),      # female, non-hispanic asian
            vec1s(c(1,2,7,22),l2),      # male, other races incl. multi-racial
            vec1s(c(1,7),l2),           # female, other races incl. multi-racial
            vec1s(c(8),l2)),            # Healthy Eating Index
          byrow =TRUE, ncol=l2)

betaj= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
betaj_ucl= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
betaj_lcl= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))

beta_effects= matrix(0, nrow=dim(C)[1], ncol=ncol(objetofda$data))
betaeffects_ucl= matrix(0, nrow=dim(C)[1], ncol=ncol(objetofda$data))
betaeffects_lcl= matrix(0, nrow=dim(C)[1], ncol=ncol(objetofda$data))

t3= proc.time()

# looping over the order of quantile tt in [0,1].
for(i in 1:ncol(objetofda$data)) {
  
  formulaaux= paste("X", i, sep="")
  formulaaux= paste(formulaaux,"~", sep="")
  
  formulafinal= paste(formulaaux, formula, sep="")
  formulafinal= as.formula(formulafinal)
  m = svyglm( formulafinal, design=data_analysis_svy, family=stats::gaussian())
  vcovm <- survey:::vcov.svyglm(m)
  meffect <- matrix( as.vector(m$coefficients), ncol=1)
  
  betaj[,i]<- meffect
  betaj_lcl[,i]<- confint(m)[,1]
  betaj_ucl[,i]<- confint(m)[,2]
  
  beta_effects[,i]<- as.vector(C%*%(meffect))
  betaeffects_ucl[,i]= as.vector(C %*% (meffect)) + qt(0.975, m$df.residual) * sqrt( diag( C%*% vcovm %*% t(C)))
  betaeffects_lcl[,i]= as.vector(C %*% (meffect)) - qt(0.975, m$df.residual) * sqrt( diag( C%*% vcovm %*% t(C)))
  
}

proc.time()-t3
write.csv(betaj, "Output_Age20to80_noTAC_betas.csv")
write.csv(betaj_ucl, "Output_Age20to80_noTAC_betas_UCL.csv" )
write.csv(betaj_lcl, "Output_Age20to80_noTAC_betas_LCL.csv" )

write.csv(beta_effects, "Output_Age20to80_noTAC_betaeffects.csv")
write.csv(betaeffects_ucl, "Output_Age20to80_noTAC_betaeffects_UCL.csv" )
write.csv(betaeffects_lcl, "Output_Age20to80_noTAC_betaeffects_LCL.csv" )



