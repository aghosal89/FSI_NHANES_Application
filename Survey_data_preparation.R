
# read the library for survey bootstrapping
library('surveybootstrap')

# Set your R directory


# read the dataset
datos <- read.csv("datosalex(1).csv")

# column indices for distributional representations
indices <- round(seq(188,687, length=500))
colnames(datos)[indices] <- paste("X", 1:500, sep="")
datos <- datos[, c(1:187, indices)]

# subset data for age 20-80 years
datos <- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  

# subset data for BMI in the range 18.5 - 40
datos <- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  

# create the response data as functional data
datosfda <- datos[,188:687]
datosfda <- as.matrix(datosfda)

# transforming the data for further analysis
datosfda[which(datosfda>280)] <- datosfda[which(datosfda>280)-1]+0.5
datos[,188:687] <- datosfda

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

data_temp <- cbind.data.frame(datosx, datosfda)
data_analysis_svy <- svydesign(id= ~survey_id, strata = ~survey_strata, 
                               weights= ~survey_wt, data = data_temp, 
                               nest = TRUE)

datosx$RIAGENDR <- as.factor(datosx$RIAGENDR)
datosx$RIDRETH3 <- as.factor(datosx$RIDRETH3)
datosx$survey_strara <- as.factor(datosx$survey_strata)
datosx$survey_id <- as.factor(datosx$survey_id)

# obtain the global Frechet model
tt= seq(0, 1, length=500)

# obtain the quantile (tt) in the range [0, 0.97]
m97 <- floor(which(tt>=0.97))[1]
tt97 <- tt[1:m97]

###################################################################################
# Bootstrapping of the PL-FSI model parameters to understand the model performance
###################################################################################

set.seed(1723)

# create the bootstrapped weights for each sample
boot_surveys <- rescaled.bootstrap.sample(survey.data = data_analysis_svy$variables,
                                          survey.design = ~survey_strata,
                                          num.reps = 100)

boot_survey_data <- matrix(NA, nrow =nrow(datosx), ncol = length(boot_surveys))

for (i in 1:length(boot_surveys)) {
  boot_survey_data[,i] <- boot_surveys[[i]]$weight.scale
}

boot_survey_data <- as.data.frame(boot_survey_data)
colnames(boot_survey_data) <- paste("Survey", 1:length(boot_surveys), sep = "_")

write.csv(boot_survey_data, "Boot_survey_data.csv")
