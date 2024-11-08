
################################################################################
# This R script performs the Bootstrap of the parameters estimated in the PL-FSI 
# model run in the age range 20-80 years with BMI, AGE in the single index; 
# GENDER, Ethnicity, HEI in the linear part.

# Set your R directory
setwd("~/Documents/2024_09_21_Aritra_codeForMainAnalysis")

# Read the libraries
library("optimParallel")
library("tidyverse")
library("survey")
library("fda.usc")
library("sandwich")
library("compiler")
library("sampling")
library('numbers')

# source the necessary files from the R directory read above
source("adj_fr_r2.R")
source("PLFSI_model.R")
source("survey2wassersteinmodel_2.R")
source("wn_cost.R")

# read the dataset contining the response distributions and the covariate data
datos <- read.csv("datosalex(1).csv")

# column indices for distributional representations
indices <- round(seq(188,687, length=500))
colnames(datos)[indices] <- paste("X", 1:500, sep="")
datos <- datos[, c(1:187, indices)]

# subset data for age 20-80 years
datos <- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  

# subset data for BMI in the range 18.5 - 40
datos<- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  

# convert the response data to functional data
datosfda <- datos[,188:687] 
datosfda <- as.matrix(datosfda)

datosfda[which(datosfda>280)] <- datosfda[which(datosfda>280)-1]+0.5
datos[,188:687] <- datosfda

#relevel(data_analysis_svy$variables$RIDRETH3,3) 

# read the sampling data for training-testing split
si_vars <- c("BMXBMI","RIDAGEYR.y")  # covariates for Single Index part
linear_vars <- c("RIAGENDR","RIDRETH3","HEI")  # covariates for Linear part
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
datosx$survey_strata <- as.factor(datosx$survey_strata)
datosx$survey_id <- as.factor(datosx$survey_id)

# obtain the global Frechet model
tt= seq(0,1,length=500)

# obtain the quantile (tt) in the range [0, 0.97].
m97 <- floor(which(tt>=0.97))[1]
tt97 <- tt[1:m97]
objetofda= fdata(datosfda,argvals = tt)

## compute the regression formula for the Global Frechet model
formula_lv<- paste(c(paste(linear_vars_categorical, collapse="*"),
                     linear_vars_numerical), 
                   collapse = "+")

###############################################
# Bootstrapping of the PL-FSI model parameters
###############################################

# read the dataset containing the multiplicative factors for the NHANES survey weights
boot_survey_data <- read.csv("Boot_survey_data.csv", header = TRUE)[,-1]

bn<- dim(boot_survey_data)[2]  # number of bootstrap samples

# estimate of the single index for each bootstrap sample
thetahat_boot <- matrix(nrow = ncol(boot_survey_data), ncol=2)
frechet_r2_boot <- rep(NA, ncol(boot_survey_data))
adjfrechet_r2_boot <- rep(NA, ncol(boot_survey_data))

betas_bootstrap <- NA
R2_integral_boot <- rep(NA, ncol(boot_survey_data))
thetaconverge_boot <- matrix(NA, nrow = ncol(boot_survey_data), ncol= 4)
Predictions_boot <- NA

#data_temp <- cbind.data.frame(datosx, datosfda)

for (i in 1:bn) {
  
  print(i)
  datosx_new <- datosx
  datosx_new$survey_wt <- datosx$survey_wt*boot_survey_data[,i]
  
  # run the Partially Linear Fréchet Single-Index regression model
  lmod <- PLFSI_model(si_vars = si_vars, linear_vars = linear_vars, 
                     datosfda = datosfda, tt = tt, nsp = 4, L = 0, 
                     etaStart = NULL, sp=4, dfs=9, formula_lv = formula_lv, 
                     datosx = datosx_new)
  
  # estimate of index parameter theta
  thetahat_boot[i,] <- lmod$thetaHat
  
  # To check if all points converged, 0 indicates convergence
  thetaconverge_boot[i,] <- lmod$converge
   
  # save the beta estimates for plotting later
  betas_bootstrap <- rbind.data.frame(betas_bootstrap, lmod$betaj)
  
  # obtain integrated R-squared for the quantile range [0, 0.97].
  R2_integral_boot[i] <- fdadensity:::trapzRcpp(X=tt97, Y=lmod$R2_vector[1:m97])
  
  # Adjusted Frechet R-squared over all the quantiles 
  adjfrechet_r2_boot[i]<- lmod$Adj_Frechet_R2
  
  # Frechet R-squared over all the quantiles
  frechet_r2_boot[i] <- lmod$Frechet_R2

  # model predictionas collated together to form a combined dataset
  Predictions_boot<- rbind(Predictions_boot, lmod$predictions)

}

#Predictions_boot <- Predictions_boot[-1,]
betas_bootstrap <- betas_bootstrap[-1,]

# save the beta estimtes, Frechet R-squared and the Adjusted Frechet R-squared
# from every bootstrap sample.

write.csv(betas_bootstrap, "Beta_estimates_bootstrap_all.csv")
write.csv(frechet_r2_boot, "Frechet_R2_boot.csv")
write.csv(adjfrechet_r2_boot, "AdjFrechet_R2_boot.csv")

#########################################################
# Computation of the bootstrap C.I. for the PL-FSI model
#########################################################

# read the true beta estimates and their effects from the PL-FSI model run
betas <- read.csv("Output_Age20to80_noTAC_betas.csv", header=T)[,-1]
beta_effects <- read.csv("Output_Age20to80_noTAC_betaeffects.csv", header = TRUE)[,-1]

beffects_boot_u95CI <- matrix(NA, nrow=dim(beta_effects)[1], ncol = dim(beta_effects)[2])
beffects_boot_l95CI <- matrix(NA, nrow=dim(beta_effects)[1], ncol = dim(beta_effects)[2])

q95_boot <- vector(length = dim(beta_effects)[1])

vec1s<- function(vec_pos, lenv) {
  vec_temp<- rep(0, lenv)
  
  for (i in 1:length(vec_pos)) {
    vec_temp[abs(vec_pos[i])] <- sign(vec_pos[i])
  }
  
  return(vec_temp)
}

l2<- dim(betas)[1]  # number of parameters in the full model

C <- matrix(c(vec1s(c(2),l2),     # male,female difference, mexican american
            vec1s(c(2,3),l2),   # male,female difference, other hispanic
            vec1s(c(2,4),l2),   # male,female difference, non-hispanic white
            vec1s(c(2,5),l2),   # male,female difference, non-hispanic black 
            vec1s(c(2,6),l2),   # male,female difference, non-hispanic asian 
            vec1s(c(2,7),l2),   # male,female difference, ORIMR
            
            vec1s(c(3),l2),     # male difference, MA vs OH 
            vec1s(c(4),l2),     # male difference, MA vs NHW
            vec1s(c(5),l2),     # male difference, MA vs HNB
            vec1s(c(6),l2),     # male difference, MA vs NHA 
            vec1s(c(7),l2),     # male difference, MA vs ORIMR
            vec1s(c(3,-4),l2),  # male difference, OH vs NHW
            vec1s(c(3,-5),l2),  # male difference, OH vs NHB
            vec1s(c(3,-6),l2),  # male difference, OH vs NHA
            vec1s(c(3,-7),l2),  # male difference, OH vs ORIMR
            vec1s(c(4,-5),l2),  # male difference, NHW vs NHB
            vec1s(c(4,-6),l2),  # male difference, NHW vs NHA
            vec1s(c(4,-7),l2),  # male difference, NHW vs ORIMR
            vec1s(c(5,-6),l2),  # male difference, NHB vs NHA
            vec1s(c(5,-7),l2),  # male difference, NHB vs ORIMR
            vec1s(c(6,-7),l2),  # male difference, NHA vs ORIMR
            
            vec1s(c(3,18),l2),  # female difference, MA vs OH
            vec1s(c(4,19),l2),  # female difference, MA vs NHW
            vec1s(c(5,20),l2),  # female difference, MA vs NHB
            vec1s(c(6,21), l2), # female difference, MA vs NHA
            vec1s(c(7,22),l2),  # female difference, MA vs ORIMR
            vec1s(c(3,18,-4,-19), l2),  # female difference, OH vs NHW
            vec1s(c(3,18,-5,-20), l2),  # female difference, OH vs NHB
            vec1s(c(3,18,-6,-21), l2),  # female difference, OH vs NHA
            vec1s(c(3,18,-7,-22), l2),  # female difference, OH vs ORIMR
            vec1s(c(4,19,-5,-20), l2),  # female difference, NHW vs NHB
            vec1s(c(4,19,-6,-21), l2),  # female difference, NHW vs NHA
            vec1s(c(4,19,-7,-22), l2),  # female difference, NHW vs ORIMR
            vec1s(c(5,20,-6,-21), l2),  # female difference, NHB vs NHA
            vec1s(c(5,20,-7,-22), l2),  # female difference, NHB vs ORIMR
            vec1s(c(6,21,-7,-22), l2), # female difference, NHA vs ORIMR
            vec1s(c(8), l2)),   # Healthy Eating Index
          
          byrow = TRUE, ncol=l2)

beffects_bootstrap<- NA

for (j1 in 1:bn) {
  
  boot_temp <- as.matrix(betas_bootstrap[seq(22*(j1-1)+1, 22*j1), ])
  boot_temp_effects<- C %*% boot_temp
  
  beffects_bootstrap <- rbind(beffects_bootstrap, boot_temp_effects)
  
}
  
beffects_bootstrap <- beffects_bootstrap[-1,]

for (j2 in 1:nrow(C)) {
  
  beffects_temp <- beffects_bootstrap[seq(j2, nrow(beffects_bootstrap), by=37), ]
  beffects_temp2 <- t(apply(beffects_temp, 1, function(x) {matrix(as.numeric(x) 
                     - as.numeric(beta_effects[j2,]), ncol=1) } ))
  
  # take the standard deviation across the simulations
  boot_stdev <- apply(beffects_temp2, 2, sd)   
  
  beffects_temp3 <- t(apply(beffects_temp2, 1, function(x) {
    max(as.numeric(na.omit(as.numeric(abs(x))/ as.numeric(boot_stdev))))
  }) )
  
  q95_boot[j2] <- quantile(beffects_temp3, prob=0.95)
  
  beffects_boot_u95CI[j2,] <- as.numeric(beta_effects[j2,] + 
                                           q95_boot[j2] * boot_stdev)
  beffects_boot_l95CI[j2,] <- as.numeric(beta_effects[j2,] - 
                                           q95_boot[j2] * boot_stdev)

}

# write the Bootstrap CIs for further plotting and exploratory analyses
write.csv(beffects_boot_l95CI, "Bootstrap_lower95_ConfidenceBound.csv")
write.csv(beffects_boot_u95CI, "Bootstrap_upper95_ConfidenceBound.csv")
