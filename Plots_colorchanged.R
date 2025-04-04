
################################################################################
# This script is used to create the plots in the figures 1 - 8 in the document
################################################################################

# Set working directory for R
setwd("~/Downloads/FSI_NHANES_Application-main")

# read the R libraries
library("latex2exp")
library("survey")
library("fda.usc")
library("viridis")
library("ggplot2")
library("energy")

library('numbers')
library('splines')
library('fdadensity')
library('osqp')
library('Rcpp')
library('Matrix')

# source the functions
source("PLFSI_model.R")
source("cuadratico.R")
source("adj_fr_r2.R")
source("survey2wassersteinmodel_2.R")
source("wn_cost.R")
source("polar2cart.R")

# read the NHANES dataset
datos= read.csv("datosalex(1).csv")
indices= round(seq(188,687, length=500))

colnames(datos)[indices]= paste("X", 1:500, sep="")
datos= datos[, c(1:187, indices)]

# subset data for age 30-80 years
datos <- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  

# subset data for BMI in the range 18.5 - 40
datos <- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  
#datos <- datos[-c(734, 1216, 2094, 2699, 2958, 3092, 3933, 3963, 4232), ]

# create the response data as functional data
datosfda <- datos[,188:687]
datosfda <- as.matrix(datosfda)
datosfda[which(datosfda>280)] <- datosfda[which(datosfda>280)-1]+0.5
datos[,188:687] <- datosfda

tt <- seq(0,1,length=500)
objetofda <- fdata(datosfda,argvals = tt)

#relevel(data_analysis_svy$variables$RIDRETH3,3) 

si_vars <- c("BMXBMI","RIDAGEYR.y")  # covariates for Single Index part
linear_vars<- c("RIAGENDR","RIDRETH3","HEI")  # covariates for Linear part
linear_vars_numerical <- "HEI"  # numeric covariates for Linear part
linear_vars_categorical <- c("RIAGENDR","RIDRETH3")  # categorical covariates for Linear part
#datosx$Gender<- as.factor(datosx$RIAGENDR)

# the survey variables
datosx <- cbind.data.frame(scale(datos[,si_vars]), 
                           datos[,linear_vars_categorical], 
                           HEI=scale(datos[,linear_vars_numerical]), 
                           survey_wt= datos[,"wtmec4yr_adj_norm"],
                           survey_id=datos[,"SDMVPSU"], 
                           survey_strata= datos[,"SDMVSTRA"])

datosx$RIAGENDR <- as.factor(datosx$RIAGENDR)

datosx$RIDRETH3 <- as.factor(datosx$RIDRETH3)
datosx$survey_strara <- as.factor(datosx$survey_strata)
datosx$survey_id <- as.factor(datosx$survey_id)
#datosx$RIAGENDR <- relevel(datosx$RIAGENDR, ref = 2)  # Apply relevel function
# obtain the global Frechet model

################################################################################
# Codes to create figure: Raw physical activity data from wearable accelerometer 
# device, obtained from the NHANES
####

# read the data frame contaning the raw accelerometer readings 
datosacel<- read.csv("datosaritra.csv")
# read data frame containing the quantile distributions of individual physical 
# activity representations
obs= c()

for(i in 1:8){
  obs= c(obs,as.numeric(unlist(datosacel[i,])))
}

# create the top panel (a) of the figure in the .eps format
setEPS()
postscript(file = "grafico_rawactivity.eps", width = 8, height = 4, paper = 'special')

t= seq(0,8,length= length(obs))
plot(t, obs, bty="n", xlab="days",ylab="Activity Count", 
     main="Raw physical activity data", pch=16, cex=0.25)
abline(v=0)
abline(v=1)
abline(v=2)
abline(v=3)
abline(v=4)
abline(v=5)
abline(v=6)
abline(v=7)
abline(v=8)

grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

dev.off()


#####################################################################################
# Codes to create figure: Physical activity representation for 1 and all participants
####

# create the panel (c) of the figure in the .eps format
postscript(file = "grafico_quantiles_all.eps", width = 4, height = 4, paper = 'special')

grid= seq(0,1,length=500)

# the quantile distributions of all the participants' physical activity representation 
plot(grid, datosfda[1,], col= 1, xlab = "t", ylim= c(0, 285),
     ylab="Empirical Quantiles", main="Physical activity representations", type = "l")
for(i in 2:nrow(datosfda)) {
  lines(grid, datosfda[i,], type = "l", col= "grey")
}


grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

dev.off()


# create the panel (b) of the figure in the .eps format
postscript(file = "grafico_quantiles_1.eps", width = 4, height = 4, paper = 'special')

# the quantile distribution of the chosen indivdual's physical activity representation 
plot(grid, quantile(obs, probs= grid,na.rm= TRUE), ylab = "Empirical Quantile",xlab="t", 
     main="Physical activity representation", type ="l", ylim= c(0, 285), lwd=2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

dev.off()


###############################################################################
# codes to create the figure: Out-of-sample evaluation, performance comparison 
# of PLF, PLFSI, GF
####

# first the adjusted Frechet R-squared
mspe_plfsi <- read.csv("MSPE_PLFSI.csv", header = TRUE)[,-1]
mspe_plf <- read.csv("MSPEs_PLF.csv", header = TRUE)[,-1]
mspe_gf <- read.csv("MSPEs_GF.csv", header =TRUE)[,-1]

mspe_diff <- data.frame("MSPEdiff" = c(mspe_gf-mspe_plfsi, mspe_plf-mspe_plfsi),
                        "Difference"= rep(c("GF-PLFSI","PLF-PLFSI"), 
                                          each= length(mspe_gf)))

dat_mspe<- cbind.data.frame("MSPE"= c(mspe_gf, mspe_plf, mspe_plfsi),
                            model= rep(c("GF","PLF","PLFSI"), each=length(mspe_gf)))

setEPS()
postscript(file = "Rplot_MSPE_variation_PLF_PLFSI_GF.eps", width = 5, 
           height = 4, paper = 'special')

par(mfrow=c(1,1))
boxplot(MSPE ~ model, data=dat_mspe, ylab ="Model MSPE", pch=16, cex=0.7)
dev.off()

postscript(file = "Rplot_MSPEdiff_PLF_PLFSI_GF.eps", width = 5, 
           height = 4, paper = 'special')

boxplot(MSPEdiff ~ Difference, data=mspe_diff, ylab ="MSPE difference",
        xlab= "Comparison pair", pch=16, cex=0.7)
abline(h=0)

dev.off()


#####################################################################
# codes to create the figure: HEI 95% pointwise confidence intervals
####

# create an equidistant grid on the interval [0,1]
tt<- seq(0,1,length=ncol(datosfda))

# consider the quantiles upto order t=0.97
qrnt <- floor(which(tt>=0.97))[1]
tt1<- tt[1:qrnt]

# read the functional parameter estimates for the covariate HEI
b_effects_hei <- read.csv("Output_Age20to80_noTAC_betaeffects_HEI.csv")[,-1]

# read the pointwise 95% lower confidence limit of the functional parameter estimates
# for the covariate HEI
b_effects_l95_hei<- read.csv("Output_Age20to80_noTAC_betaeffects_LCL_HEI.csv")[,-1]

# read the pointwise 95% upper confidence limit of the functional parameter estimates
# for the covariate HEI
b_effects_u95_hei<- read.csv("Output_Age20to80_noTAC_betaeffects_UCL_HEI.csv")[,-1]

# consider the estimated functional parameters for the covariate HEI and their 
# functional lower and upper 95% Confidence limits upto t=0.97.
b1_eff_hei<- b_effects_hei[,c(1:qrnt)]
b_l95_1eff_hei<- b_effects_l95_hei[,c(1:qrnt)]
b_u95_1eff_hei<- b_effects_u95_hei[,c(1:qrnt)]

plotdf_hei<- data.frame(Estimate= as.numeric(t(b1_eff_hei)), 
                        u95CI= as.numeric(t(b_u95_1eff_hei)), 
                        l95CI= as.numeric(t(b_l95_1eff_hei)),
                        t=tt1)

plotdf_hei<- data.frame(plotdf_hei, Coefficient=rep(c("Healthy Eating Index"), 
                                                    each= length(tt1)),
                        Ord = rep(c(1:nrow(b1_eff_hei)), each= length(tt1)))

gen_hei<- c(min(b_l95_1eff_hei), max(b_u95_1eff_hei))

blank_data_hei<- data.frame(Coefficient= rep(c("Healthy Eating Index"), each= 2), 
                            t = 0, Coefficients = gen_hei, 
                            CI= rep(c("Lower 95%","Upper 95%"), nrow(b1_eff_hei)),
                            Ord=1)

ggplot(plotdf_hei, aes(x = t, y = Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI),
              fill = "#999999", alpha = 0.6) + 
  geom_line(size = .6, color = "black") +
  labs(x = "t", title = "Healthy Eating Index (HEI)") +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_hei, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_HEI_2.png", 
       width = 5.5, height = 3.5) 


################################################################################
# codes to create the figure: pointwise 95% Confidence Intervals for Female-Male 
# difference in physical activity for the various ethnicities
####

# read the functional effect differences between Males and Females 
b_effects_mf <- read.csv("Output_Age20to80_noTAC_betaeffects_MF.csv")[,-1]

# read the pointwise 95% lower confidence limit of the functional effect differences 
# between Males and Females
b_effects_l95_mf<- read.csv("Output_Age20to80_noTAC_betaeffects_LCL_MF.csv")[,-1]

# read the pointwise 95% upper confidence limit of the functional effect differences 
# between Males and Females
b_effects_u95_mf<- read.csv("Output_Age20to80_noTAC_betaeffects_UCL_MF.csv")[,-1]

# consider the estimated effect differences and their lower and upper 95% Confidence
# Intervals between males and females upto t=0.97.
b1_eff_mf <- b_effects_mf[,c(1:qrnt)]
b_l95_1eff_mf<- b_effects_l95_mf[,c(1:qrnt)]
b_u95_1eff_mf<- b_effects_u95_mf[,c(1:qrnt)]

# create the data frame that gathers the difference of intercepts and their pointwise 
# 95% Confidence Intervals in a format helpful in creating the plots in figure 3. The 
# goal is to understand the disparity in physical activity levels between Males and 
# Females of different ethnicities.

plotdf_mf <- data.frame(Estimate=NaN, u95CI=NaN, l95CI=NaN, t=NaN)

for (i in 1:nrow(b1_eff_mf)) {
  plotdf_temp<- data.frame(Estimate= as.numeric(t(b1_eff_mf[i,])), 
                           u95CI= as.numeric(t(b_u95_1eff_mf[i,])), 
                           l95CI= as.numeric(t(b_l95_1eff_mf[i,])),
                           t=tt1)
  
  plotdf_mf<- rbind(plotdf_mf, plotdf_temp)
}

plotdf_mf<- plotdf_mf[-1,]

plotdf_mf <- data.frame(plotdf_mf, 
                  Coefficient=rep(c("Mexican American",
                                    "Other Hispanic", 
                                    "Non-Hispanic White",
                                    "Non-Hispanic Black",
                                    "Non-Hispanic Asian",
                                    "Other Races-IMR"), 
         each= length(tt1)), Ord = rep(c(1:nrow(b1_eff_mf)), each= length(tt1)))

# For better visual comparison, the y-axes were fixed, so that the disparity among 
# different ethnicities were studied. Here read from data the upper and lower limits
# of the y-axis and create an auxiliary dataset to fix the order of the panels to 
# appear in the figure.

gen_eth_mf<- rep(c(min(b_l95_1eff_mf),max(b_u95_1eff_mf)), nrow(b1_eff_mf)) 

blank_data_mf <- data.frame(Coefficient=rep(c("Mexican American",
                                              "Other Hispanic", 
                                              "Non-Hispanic White",
                                              "Non-Hispanic Black",
                                              "Non-Hispanic Asian",
                                              "Other Races-IMR"), each=2),
                            t = 0, Coefficients = gen_eth_mf,
                            CI= rep(c("Lower 95%","Upper 95%"), 
                                    nrow(b1_eff_mf)), Ord=1)

# create the main figure and save it in the .png format
ggplot(data = plotdf_mf, aes(x = t, y = Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI), fill = "#999999",
              alpha = 0.6) + 
  geom_line(size = .6, color = "black") +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(vars(Coefficient), nrow= 5) +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_mf, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_MF_2.png", 
       width = 6.5, height = 5) 


######################################################################################
# codes to create the figure: pointwise 95% Confidence Intervals for physical activity
# of females of each ethnicity.
###############################

# read the functional effect differences between Females of each pair of ethnicities 
b_effects_f <- read.csv("Output_Age20to80_noTAC_betaeffects_F.csv")[,-1]

# read the pointwise 95% lower confidence limits of the functional effect differences 
# between Females of each pair of ethnicities
b_effects_l95_f<- read.csv("Output_Age20to80_noTAC_betaeffects_LCL_F.csv")[,-1]

# read the pointwise 95% upper confidence limits of the functional effect differences 
# between Males of each pair of ethnicities
b_effects_u95_f<- read.csv("Output_Age20to80_noTAC_betaeffects_UCL_F.csv")[,-1]

# consider the estimated effect differences and their lower and upper 95% Confidence
# Intervals between females of different pairs of ethnicities upto t=0.97.
b1_eff_f<- b_effects_f[,c(1:qrnt)]
b_l95_1eff_f<- b_effects_l95_f[,c(1:qrnt)]
b_u95_1eff_f<- b_effects_u95_f[,c(1:qrnt)]

plotdf_f <- data.frame(Estimate=NaN, u95CI=NaN, l95CI=NaN, t=NaN)

for (i in 1:nrow(b1_eff_f)) {
  plotdf_temp<- data.frame(Estimate= as.numeric(t(b1_eff_f[i,])), 
                           u95CI= as.numeric(t(b_u95_1eff_f[i,])), 
                           l95CI= as.numeric(t(b_l95_1eff_f[i,])),
                           t=tt1)
  
  plotdf_f<- rbind(plotdf_f, plotdf_temp)
}

plotdf_f<- plotdf_f[-1,]


plotdf_f<-    data.frame(plotdf_f, 
                         Coefficient= rep(c("OH - MA",
                                            "NHW - MA", 
                                            "NHB - MA",
                                            "NHA - MA",
                                            "ORIMR - MA",
                                            "OH - NHW",
                                            "OH - NHB",
                                            "OH - NHA", 
                                            "OH - ORIMR",
                                            "NHW - NHB",
                                            "NHW - NHA",
                                            "NHW - ORIMR",
                                            "NHB - NHA",
                                            "NHB - ORIMR",
                                            "NHA - ORIMR"),
                                          each= length(tt1)), Ord = rep(c(1:nrow(b1_eff_f)),each= length(tt1)))

# create a data frame which contains the nature of the differences considered in the 
# order which appears in the data frame 'plotdf_all'. The abbreviations for the ethnicities
# OH   : Other Hispanic
# MA   : Mexican American
# NHW  : Non-Hispanic White
# NHB  : Non-Hispanic Black
# NHA  : Non-Hispanic Asian
# ORIMR: Other Races Including Multi-Racial.


gen_eth_f<- rep(c(min(b_l95_1eff_f),max(b_u95_1eff_f)), nrow(b1_eff_f))

# For better visual comparison, the y-axes were fixed, so that the disparity among 
# different ethnicities were studied. Here read from data the upper and lower limits
# of the y-axis and create an auxiliary dataset to fix the order of the panels to 
# appear in the figure.


blank_data_f<- data.frame(Coefficient=rep(c("OH - MA",
                                            "NHW - MA", 
                                            "NHB - MA",
                                            "NHA - MA",
                                            "ORIMR - MA",
                                            "OH - NHW",
                                            "OH - NHB",
                                            "OH - NHA", 
                                            "OH - ORIMR",
                                            "NHW - NHB",
                                            "NHW - NHA",
                                            "NHW - ORIMR",
                                            "NHB - NHA",
                                            "NHB - ORIMR",
                                            "NHA - ORIMR"), each=2),
                          t = 0, Coefficients = gen_eth_f,
                          CI= rep(c("Lower 95%","Upper 95%"), nrow(b1_eff_f)),
                          Ord=1)


# the main plotting code, here the data frames 'plotdf_all' and the 'blank_data'
# are both used to modify the aesthetics of the plot.

ggplot(data = plotdf_f, aes(x = t, y = Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI), fill = "#999999", alpha = 0.6) + 
  geom_line(size = .5, color = "black") +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(vars(Coefficient), nrow= 5) +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_f, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_F_2.png", 
       width = 8.5, height = 9) 


######################################################################################
# codes to create the figure: pointwise 95% Confidence Intervals for physical activity
# of males of each ethnicity.
####

# read the functional effect differences between Males of each pair of ethnicities 
b_effects_m <- read.csv("Output_Age20to80_noTAC_betaeffects_M.csv")[,-1]

# read the pointwise 95% lower confidence limits of the functional effect differences 
# between Males of each pair of ethnicities
b_effects_l95_m<- read.csv("Output_Age20to80_noTAC_betaeffects_LCL_M.csv")[,-1]

# read the pointwise 95% upper confidence limits of the functional effect differences 
# between Males of each pair of ethnicities
b_effects_u95_m<- read.csv("Output_Age20to80_noTAC_betaeffects_UCL_M.csv")[,-1]

# consider the estimated effect differences and their lower and upper 95% Confidence
# Intervals between males of different pairs of ethnicities upto t=0.97.
b1_eff_m<- b_effects_m[,c(1:qrnt)]
b_l95_1eff_m<- b_effects_l95_m[,c(1:qrnt)]
b_u95_1eff_m<- b_effects_u95_m[,c(1:qrnt)]

plotdf_m <- data.frame(Estimate=NaN, u95CI=NaN, l95CI=NaN, t=NaN)

for(i in 1:nrow(b1_eff_m)) {
  plotdf_temp<- data.frame(Estimate= as.numeric(t(b1_eff_m[i,])), 
                           u95CI= as.numeric(t(b_u95_1eff_m[i,])), 
                           l95CI= as.numeric(t(b_l95_1eff_m[i,])),
                           t=tt1)
  plotdf_m<- rbind(plotdf_m, plotdf_temp)
}

plotdf_m<- plotdf_m[-1,]

plotdf_m <- data.frame(plotdf_m,         
                       Coefficient= rep(c("OH - MA",
                                          "NHW - MA", 
                                          "NHB - MA",
                                          "NHA - MA",
                                          "ORIMR - MA",
                                          "OH - NHW",
                                          "OH - NHB",
                                          "OH - NHA", 
                                          "OH - ORIMR",
                                          "NHW - NHB",
                                          "NHW - NHA",
                                          "NHW - ORIMR",
                                          "NHB - NHA",
                                          "NHB - ORIMR",
                                          "NHA - ORIMR"),
                                        each= length(tt1)),
                       Ord = rep(c(1:nrow(b1_eff_m)), each= length(tt1)))

# For better visual comparison, the y-axes were fixed, so that the disparity among 
# different ethnicities were studied. Here read from data the upper and lower limits
# of the y-axis and create an auxiliary dataset to fix the order of the panels to 
# appear in the figure.

gen_eth_m<- rep(c(min(b_l95_1eff_m),max(b_u95_1eff_m)), nrow(b1_eff_m)) 

blank_data_m<- data.frame(Coefficient= rep(c("OH - MA",
                                             "NHW - MA", 
                                             "NHB - MA",
                                             "NHA - MA",
                                             "ORIMR - MA",
                                             "OH - NHW",
                                             "OH - NHB",
                                             "OH - NHA", 
                                             "OH - ORIMR",
                                             "NHW - NHB",
                                             "NHW - NHA",
                                             "NHW - ORIMR",
                                             "NHB - NHA",
                                             "NHB - ORIMR",
                                             "NHA - ORIMR"), each=2),
                          t = 0, Coefficients = gen_eth_m,
                          CI= rep(c("Lower 95%","Upper 95%"), 
                                  nrow(b1_eff_m)), Ord=1)

# create the main figure and save it in .png format.
ggplot(data = plotdf_m, aes(x = t, y = Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI), fill = "#999999",
              alpha = 0.6) + 
  geom_line(size = .6, color = "black") +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(vars(Coefficient), nrow= 5) +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_m, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_M_2.png", 
       width = 8.5, height = 9) 


##############################################################################
# Codes to create figure: Heatmap of the PLFSI model prediction for ranges of
# Age and BMI of patients
####

# The length of the equidistant grid on the standardized ranges of both Age and BMI
# for our model
m2<- 500

# create the grid on standardized BMI range
x_ <- seq(min(datos$BMXBMI), max(datos$BMXBMI), len=m2)

# create the grid on standardized age range
y_ <- seq(min(datos$RIDAGEYR.y), max(datos$RIDAGEYR.y), len=m2)

# read the functional regression coed=fficients of the partially linear Frechet 
# single index model
b <- read.csv("Output_Age20to80_noTAC_betas.csv", header = TRUE)[,-1]

# consider the variables on a 2-dimensional grid to construct the heatmap plot. 
temp_si_data<- data.frame(BMXBMI= seq(min(datosx$BMXBMI), max(datosx$BMXBMI), len=m2), 
                          RIDAGEYR.y=seq(min(datosx$RIDAGEYR.y),max(datosx$RIDAGEYR.y),len=m2))
tsd_grid<- expand.grid(temp_si_data)
tsd_grid_level<- expand.grid(data.frame(BMI=x_, Age=y_))

# read the estimated index parameter
th <- matrix(read.csv("Theta_Hat.csv")[,-1], ncol =1)

# store the index corresponding to 50,75,90,97th percentiles on the grid tt
q50<- floor(which(tt>=0.50))[1]    
q75<- floor(which(tt>=0.75))[1]   
q90<- floor(which(tt>=0.90))[1]   
q97<- floor(which(tt>=0.97))[1]

# set the order and degrees of the bspline basis
sp <- 4 
dfs <- 9

# get the bspline basis expanstion of single index portion
u_ <- as.matrix(tsd_grid)%*%th
bs_data1 <- data.frame(splines2::dbs(x=u_, degree=sp, df=dfs))
colnames(bs_data1) <- paste("BS", seq(1:ncol(bs_data1)), sep ="")

# add intercept, fix the HEI variable to its median level and add columns of the 
# bspline basis functions to form a single data frame. 
x_males_eth1 <- cbind.data.frame(Intercept=1, HEI=median(datosx$HEI), bs_data1)

# creating the predictions for all quantiles on tt using the functional beta 
# coefficients of the partially linear frechet single index model
yh_males_eth1 <- matrix(NA, nrow = dim(x_males_eth1)[1], ncol = length(tt))

for (r in 1:length(tt)) {
  yh_males_eth1[,r] <- as.matrix(x_males_eth1)%*% matrix(b[,r][c(1,8:17)],ncol=1)
}

# transporting the predictions to elements of L^2-Wasserstein space if necessary
t<- proc.time()
indt <- apply(yh_males_eth1, 1, function(x) {all(diff(x)>= 0)==FALSE})
for(j in 1:length(indt)) {
  if(indt[j]==TRUE) {
    yh_males_eth1[j,] <- cuadratico(matrix(yh_males_eth1[j,], nrow = 1))
  }
}
proc.time() - t

# save the quantiles for future use 
write.csv(yh_males_eth1, "Pedicted_transport_percentiles.csv")
Yh<- as.data.frame(yh_males_eth1[,c(q50, q75, q90, q97)])
colnames(Yh)<- c("Q50","Q75","Q90","Q97")

# create the data frames necessary for creating the plots 
# for the 50th percentile
df_q50<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                       AGE=tsd_grid_level$Age,
                       MIMS=Yh$Q50)

# for the 75th precentile
df_q75<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          MIMS=Yh$Q75)

# for the 90th percentile
df_q90<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          MIMS=Yh$Q90)

# for the 97th percentile
df_q97<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          MIMS=Yh$Q97)

# the main plotting codes for the 50,75,90,97th quantiles 
ggplot(df_q50, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = MIMS)) +
  #scale_color_grey() +
  scale_fill_gradientn(colors = gray.colors(100,
                                            start = 0,
                                            end = .9, rev = TRUE)) +
  theme(axis.text=element_text(size=15),
      axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("50th Quantile")
ggsave(file = "Rplot_quantile50_prediction.png", width = 6.5, height = 5.5)


ggplot(df_q75, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = MIMS)) +
  scale_fill_gradientn(colors = gray.colors(100,
                                            start = 0,
                                            end = .9, rev = TRUE)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("75th Quantile")
ggsave(file = "Rplot_quantile75_prediction.png", width = 6.5, height = 5.5)


ggplot(df_q90, aes(BMI, AGE)) +
  geom_raster(aes(fill = MIMS)) +
  scale_fill_gradientn(colors = gray.colors(100,
                                            start = 0,
                                            end = .9, rev = TRUE)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("90th Quantile")
ggsave(file = "Rplot_quantile90_prediction.png", width = 6.5, height = 5.5)


ggplot(df_q97, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = MIMS)) +
  scale_fill_gradientn(colors = gray.colors(100,
                                            start = 0,
                                            end = .9, rev = TRUE)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("97th Quantile")
ggsave(file = "Rplot_quantile97_prediction.png", width = 6.5, height = 5.5)


################################################################################
# Codes to create figure: Heatmap of the integral of the quantiles of the PLFSI 
# model prediction for ranges of Age and BMI of patients
####

yh_males_eth1_int<- sapply(1:nrow(yh_males_eth1), function(i) {
  fdadensity:::trapzRcpp(X=tt, Y= yh_males_eth1[i,])
})

df_int<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          MIMS=yh_males_eth1_int)

ggplot(df_int, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = MIMS)) +
  scale_fill_gradientn(colors = gray.colors(100,
                       start = 0, end = .9, rev = TRUE)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("Integral of Quantiles of MIMS")
ggsave(file = "Rplot_quantile_integrals.png", width = 6.5, height = 5.5)


######################################################
# Codes to create figure: derivative of the quantiles 
####

# the length of the equidistant grid of the single index consisting of standardized 
# BMI and age for the partially linear Frechet single index model.
m4<- 60

# consider equidistant grid of variables BMXBMI, RIDAGEYR.y over their ranges. 
# For each variables consider of length = m4. 

x_ <- seq(min(datos$BMXBMI), max(datos$BMXBMI), len=m4)
y_ <- seq(min(datos$RIDAGEYR.y), max(datos$RIDAGEYR.y), len=m4)

temp_si_data<- data.frame(BMXBMI= seq(min(datosx$BMXBMI),max(datosx$BMXBMI), len=m4), 
                RIDAGEYR.y=seq(min(datosx$RIDAGEYR.y),max(datosx$RIDAGEYR.y),len=m4))

# consider the variables on a 2 dimensional grid to construct the contour plot. 
tsd_grid<- expand.grid(temp_si_data)
tsd_grid_level<- expand.grid(data.frame(BMI=x_, Age=y_))

# univariate single index obtained 
u_2<- as.matrix(cbind(datosx$BMXBMI, datosx$RIDAGEYR.y))%*%th
u_2 <- sort(u_2)

# compute the bspline basis expansion of the single index 
aspl1.1<- splines2::dbs(u_2, df= dfs, knots = NULL, degree = sp, 
                      intercept = FALSE)

#####################
# The following are the bondary knots and the rest of the knot sequence for the 
# spline basis computation above
bknots<- attributes(aspl1.1)$Boundary.knots
knots <- attributes(aspl1.1)$knots
knots_all<- c(bknots[1], knots, bknots[2])

# compute the first derivative of the basis spline as in equation (4.15) in the 
# paper with respect to the single index 
aspl2<- splines2::dbs(x=u_2, derivs = 1L, df= dfs, knots = knots, 
                      degree = sp, Boundary.knots= bknots,
                      intercept = FALSE)
  
# check index for 16.6667th percentile of the ordered single index
row_ind_fbo<- floor(which(u_2>=knots_all[2]))[1]
  
# check index for 83.3333rd percentile of the ordered single index
row_ind_lbo<- floor(which(u_2>=knots_all[length(knots_all)-1]))[1]

# obtain the derivative of the model predictions for prior to mapping them to 
# L^2 Wasserstein space  
drv2<- matrix(NA, nrow =nrow(aspl2), ncol =length(tt)+1)
for (r in 1:length(tt)) {
  drv2[,r] <- aspl2%*% matrix(b[,r][c(9:17)],ncol=1)
}

# obtain the integral of all derivative quantiles over their range [0,1]
drv2[,ncol(drv2)]<- sapply(1:nrow(drv2[,1:length(tt)]), function(i) {
  fdadensity:::trapzRcpp(X=tt, Y= drv2[i,1:length(tt)]) })

# since we are interested to study the derivative of the bspline basis with 
# respect to the single index, hence we draw vertical lines at all the 
# knot locations. Outside the lower and upper boundary knots, we present the 
# derivative of the bspline bases in dotted lines, while within the boundary knots,
# we present the lines as solid.


df_in <- drv2[c(row_ind_fbo:row_ind_lbo),]
u_2_in<- u_2[c(row_ind_fbo:row_ind_lbo)]
u_2_boundary_l<- u_2[c(1:row_ind_fbo)]
df_boundary_l <- drv2[c(1:row_ind_fbo),]
u_2_boundary_r <- u_2[c(row_ind_lbo:length(u_2))]
df_boundary_r <- drv2[c(row_ind_lbo:length(u_2)),]

png(file = "Rplot_qderiv_SI.png", width = 5, height = 5, units = "in", res =300)

#"#E69F00", "#56B4E9", "#009E73", "#0072B2", #CCBB44

plot(u_2_in, df_in[,q50], type ="l", ylab ="Derivative of Quantiles", 
     xlab="Single Index", ylim=c(-10,25), xlim =c(-2.5,2.7), lwd=1.5, col="black")
abline(h=0, lty="dotted")
lines(u_2_in, df_in[,q75], col="black", lwd=2)
lines(u_2_in, df_in[,q90], col="black", lwd=2.5)
lines(u_2_in, df_in[,q97], col="black", lwd=3)
#lines(u_2_in, df_in[,ncol(drv2)], col="#CCBB44", lwd=3.5)

lines(u_2_boundary_l, df_boundary_l[,q50], lwd=1.5, lty=4, col="grey")
lines(u_2_boundary_l, df_boundary_l[,q75], col="grey", lwd=2, lty =4)
lines(u_2_boundary_l, df_boundary_l[,q90], col="grey", lwd=2.5, lty =4)
lines(u_2_boundary_l, df_boundary_l[,q97], col="grey", lwd=3, lty =4)
#lines(u_2_boundary_l, df_boundary_l[,ncol(drv2)], col="#CCBB44", lwd=3.5, lty =3)

lines(u_2_boundary_r, df_boundary_r[,q50], lwd=1.5, lty=4, col="grey")
lines(u_2_boundary_r, df_boundary_r[,q75], col="grey", lwd=2.5, lty =4)
lines(u_2_boundary_r, df_boundary_r[,q90], col="grey", lwd=3, lty =4)
lines(u_2_boundary_r, df_boundary_r[,q97], col="grey", lwd=3.5, lty =4)
#lines(u_2_boundary_r, df_boundary_r[,ncol(drv2)], col="#CCBB44", lwd=3.5, lty =3)

for (j in 1:length(knots_all)) {
  abline(v=knots_all[j], lty=3)
}

legend(-.5, 25, legend=c("t=0.50","t=0.75","t=0.90","t=0.97"),
       #col=c("#000000", "#D55E00", "#009E73", "#0072B2"), 
       lwd=c(1.5,2,2.5,3), lty=c(1,1,1,1), cex=1.2, text.font=10, bg='lightgrey')

#grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
#     lwd = par("lwd"), equilogs = TRUE)

dev.off()

png(file = "Rplot_histogram_ofSI.png", width = 4.5, height = 4.5, units = "in", 
    res=300)

hist(u_2, freq=FALSE, main="Histogram of Single Index", xlab= "")

dev.off()


#################################################################################
# Codes to create the figure: Analysis of log TAC with respect to three clusters
####

# we read the model residuals from the fitted PLSIFR-model
datosw <- read.csv("Output_Age20to80_noTAC_residuals.csv")

# Compute the TAC variable for each participant
datosw$TAC= apply(datos[,c(188:687)],1,mean) 

x2= kgroups(1/100*datosw[,round(seq(3,qrnt,length=100))],k=3)  # why 40 was deducted?
datosw$grupo= x2$cluster
fda= fdata(datosw[,c(3:502)], argvals = seq(0,1,length=500))
cluster1= fda[x2$cluster==1,1:qrnt]
cluster2= fda[x2$cluster==2,1:qrnt]
cluster3= fda[x2$cluster==3,1:qrnt]

media1= func.mean(cluster1)
media2= func.mean(cluster2)
media3= func.mean(cluster3)

png(file= "clusteringresults.eps", units ="in", height = 6, width = 7, res=300)
par(mfrow= c(2,3))
plot(cluster1[,1:qrnt], col="grey", main="Cluster 1", ylim= c(-30,95), xlim=c(0,1.001),
     ylab="Residuals")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

plot(cluster2[,1:qrnt], col="grey", main="Cluster 2",ylim= c(-30,95), xlim=c(0,1.001),
     ylab="Residuals")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

plot(cluster3[,1:qrnt], col="grey", main="Cluster 3",ylim= c(-30,95), xlim=c(0,1.001)
     , ylab="Residuals")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

boxplot(log(datos$TAC)~datosw$grupo, ylim= c(7,11), ylab= "log TAC", xlab="Cluster",
        pch=16, cex=0.4)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

boxplot(log(datos$LBXSGL_43)~datosw$grupo, ylim= c(3.5,6.8), ylab= "log Glucose, mg/dL", 
        xlab="Cluster", pch=16, cex=0.4)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

#boxplot(datos2$LBXGH_39~datosw$grupo, ylim= c(4,15), ylab= "A1C, %", xlab="Cluster")
boxplot(log(datos$LBXSCR_43)~datosw$grupo, ylim= c(-1,2), ylab= "log Creatinine, mg/dL",
        pch=16, cex=0.4, xlab="Cluster")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

dev.off()


###################################################################
# Codes for the Bootstrap CIs for the parameters and their effects
####

beta_boot_u95 <- read.csv("Bootstrap_upper95_ConfidenceBound.csv")[,-1]
beta_boot_l95 <- read.csv("Bootstrap_lower95_ConfidenceBound.csv")[,-1]
beta_effects <- read.csv("Output_Age20to80_noTAC_betaeffects.csv")[,-1]

# read the pointwise 95% lower/upper bootstrap confidence limit of the functional 
# effect differences between Males and Females
boot_effects_l95_mf<- beta_boot_l95[c(1:6),]
boot_effects_u95_mf<- beta_boot_u95[c(1:6),]

# consider the estimated effect differences and their lower and upper 95% Confidence
# Intervals between males and females upto t=0.97.
boot_l95_1eff_mf<- boot_effects_l95_mf[,c(1:qrnt)]
boot_u95_1eff_mf<- boot_effects_u95_mf[,c(1:qrnt)]

# create the data frame that gathers the difference of intercepts and their pointwise 
# 95% Confidence Intervals in a format helpful in creating the plots in figure 3. The 
# goal is to understand the disparity in physical activity levels between Males and 
# Females of different ethnicities.

plotdf_mf <- data.frame(Estimate=NaN, u95CI=NaN, l95CI=NaN, t=NaN)

for (i in 1:nrow(b1_eff_mf)) {
  plotdf_temp<- data.frame(Estimate= as.numeric(t(b1_eff_mf[i,])), 
                           u95CI= as.numeric(t(boot_u95_1eff_mf[i,])), 
                           l95CI= as.numeric(t(boot_l95_1eff_mf[i,])),
                           t=tt1)
  
  plotdf_mf<- rbind(plotdf_mf, plotdf_temp)
}

plotdf_mf<- plotdf_mf[-1,]

plotdf_mf <- data.frame(plotdf_mf, 
                        Coefficient=rep(c("Mexican American",
                                          "Other Hispanic", 
                                          "Non-Hispanic White",
                                          "Non-Hispanic Black",
                                          "Non-Hispanic Asian",
                                          "Other Races-IMR"), 
                                        each= length(tt1)), 
                        Ord = rep(c(1:nrow(b1_eff_mf)), each= length(tt1)))

# For better visual comparison, the y-axes were fixed, so that the disparity among 
# different ethnicities were studied. Here read from data the upper and lower limits
# of the y-axis and create an auxiliary dataset to fix the order of the panels to 
# appear in the figure.

gen_eth_mf<- rep(c(min(boot_l95_1eff_mf),max(boot_u95_1eff_mf)), nrow(b1_eff_mf)) 

blank_data_mf <- data.frame(Coefficient=rep(c("Mexican American",
                                              "Other Hispanic", 
                                              "Non-Hispanic White",
                                              "Non-Hispanic Black",
                                              "Non-Hispanic Asian",
                                              "Other Races-IMR"), each=2),
                            t = 0, Coefficients = gen_eth_mf,
                            CI= rep(c("Lower 95%","Upper 95%"), 
                                    nrow(b1_eff_mf)), Ord=1)

# create the main figure and save it in the .png format

#####################################################################################
# This is for the male-female differences of physical activities and their pointwise
# 95% Confidence Intervals
####

ggplot(data = plotdf_mf, aes(x = t, y =Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI), fill = "#999999",
              alpha = 0.6) + 
  geom_line(size = .6, color = "black") +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(vars(Coefficient), nrow= 5) +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_mf, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_MF_bootstrap.png", 
       width = 6.5, height = 5) 


# read the pointwise 95% upper/lower bootsrtap confidence limits of the functional effect differences 
# between Males of each pair of ethnicities
boot_effects_u95_m<- beta_boot_u95[c(7:21),]
boot_effects_l95_m<- beta_boot_l95[c(7:21),]

# consider the estimated effect differences and their lower and upper 95% Confidence
# Intervals between males of different pairs of ethnicities upto t=0.97.
boot_l95_1eff_m <- boot_effects_l95_m[,c(1:qrnt)]
boot_u95_1eff_m <- boot_effects_u95_m[,c(1:qrnt)]

plotdf_m <- data.frame(Estimate=NaN, u95CI=NaN, l95CI=NaN, t=NaN)

for(i in 1:nrow(b1_eff_m)) {
  plotdf_temp<- data.frame(Estimate= as.numeric(t(b1_eff_m[i,])), 
                           u95CI= as.numeric(t(boot_u95_1eff_m[i,])), 
                           l95CI= as.numeric(t(boot_l95_1eff_m[i,])),
                           t=tt1)
  plotdf_m<- rbind(plotdf_m, plotdf_temp)
}

plotdf_m<- plotdf_m[-1,]

plotdf_m <- data.frame(plotdf_m,         
                       Coefficient= rep(c("OH - MA",
                                          "NHW - MA", 
                                          "NHB - MA",
                                          "NHA - MA",
                                          "ORIMR - MA",
                                          "OH - NHW",
                                          "OH - NHB",
                                          "OH - NHA", 
                                          "OH - ORIMR",
                                          "NHW - NHB",
                                          "NHW - NHA",
                                          "NHW - ORIMR",
                                          "NHB - NHA",
                                          "NHB - ORIMR",
                                          "NHA - ORIMR"),
                                        each= length(tt1)),
                       Ord = rep(c(1:nrow(b1_eff_m)), each= length(tt1)))

# For better visual comparison, the y-axes were fixed, so that the disparity among 
# different ethnicities were studied. Here read from data the upper and lower limits
# of the y-axis and create an auxiliary dataset to fix the order of the panels to 
# appear in the figure.

gen_eth_m <- rep(c(min(boot_l95_1eff_m), max(boot_u95_1eff_m)), nrow(b1_eff_m)) 

blank_data_m<- data.frame(Coefficient= rep(c("OH - MA",
                                             "NHW - MA", 
                                             "NHB - MA",
                                             "NHA - MA",
                                             "ORIMR - MA",
                                             "OH - NHW",
                                             "OH - NHB",
                                             "OH - NHA", 
                                             "OH - ORIMR",
                                             "NHW - NHB",
                                             "NHW - NHA",
                                             "NHW - ORIMR",
                                             "NHB - NHA",
                                             "NHB - ORIMR",
                                             "NHA - ORIMR"), each=2),
                          t = 0, Coefficients = gen_eth_m,
                          CI= rep(c("Lower 95%","Upper 95%"), 
                                  nrow(b1_eff_m)), Ord=1)

# create the main figure and save it in .png format.
ggplot(data = plotdf_m, aes(x = t, y = Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI), fill = "#999999",
              alpha = 0.6) + 
  geom_line(size = .6, color = "black") +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(vars(Coefficient), nrow= 5) +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_m, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_M_bootstrap.png", 
       width = 8.5, height = 9) 

# read the pointwise 95% lower/upper bootstrap confidence limits of the functional effect differences 
# between Females of each pair of ethnicities
boot_effects_l95_f <- beta_boot_l95[c(22:36), ]
boot_effects_u95_f <- beta_boot_u95[c(22:36), ]

# consider the estimated effect differences and their lower and upper 95% bootstrap CIs
# between females of different pairs of ethnicities upto t=0.97.
boot_l95_1eff_f<- boot_effects_l95_f[, c(1:qrnt)]
boot_u95_1eff_f<- boot_effects_u95_f[, c(1:qrnt)]

plotdf_f <- data.frame(Estimate=NaN, u95CI=NaN, l95CI=NaN, t=NaN)

# read the functional effect differences between Females of each pair of ethnicities 
b_effects_f <- read.csv("Output_Age20to80_noTAC_betaeffects_F.csv")[,-1]
b1_eff_f<- b_effects_f[,c(1:qrnt)]


for (i in 1:nrow(b1_eff_f)) {
  plotdf_temp<- data.frame(Estimate= as.numeric(t(b1_eff_f[i,])), 
                           u95CI= as.numeric(t(boot_u95_1eff_f[i,])), 
                           l95CI= as.numeric(t(boot_l95_1eff_f[i,])),
                           t=tt1)
  
  plotdf_f<- rbind(plotdf_f, plotdf_temp)
}

plotdf_f<- plotdf_f[-1,]

gen_eth_f<- rep(c(min(boot_l95_1eff_f),max(boot_u95_1eff_f)), nrow(b1_eff_f))

plotdf_f <-   data.frame(plotdf_f, 
                         Coefficient= rep(c("OH - MA",
                                            "NHW - MA", 
                                            "NHB - MA",
                                            "NHA - MA",
                                            "ORIMR - MA",
                                            "OH - NHW",
                                            "OH - NHB",
                                            "OH - NHA", 
                                            "OH - ORIMR",
                                            "NHW - NHB",
                                            "NHW - NHA",
                                            "NHW - ORIMR",
                                            "NHB - NHA",
                                            "NHB - ORIMR",
                                            "NHA - ORIMR"),
                                          each= length(tt1)), 
                         Ord = rep(c(1:nrow(b1_eff_f)), each= length(tt1)))

blank_data_f<- data.frame(Coefficient=rep(c("OH - MA",
                                            "NHW - MA", 
                                            "NHB - MA",
                                            "NHA - MA",
                                            "ORIMR - MA",
                                            "OH - NHW",
                                            "OH - NHB",
                                            "OH - NHA", 
                                            "OH - ORIMR",
                                            "NHW - NHB",
                                            "NHW - NHA",
                                            "NHW - ORIMR",
                                            "NHB - NHA",
                                            "NHB - ORIMR",
                                            "NHA - ORIMR"), each=2),
                          t = 0, Coefficients = gen_eth_f,
                          CI= rep(c("Lower 95%","Upper 95%"), nrow(b1_eff_f)),
                          Ord=1)

# plot of the bootstrap confidence intervals for the females patients
ggplot(data = plotdf_f, aes(x = t, y = Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI), fill = "#999999", alpha = 0.6) + 
  geom_line(size = .5, color = "black") +
  guides(fill = FALSE, color = FALSE) +
  facet_wrap(vars(Coefficient), nrow= 5) +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_f, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_F_bootstrap.png", 
       width = 8.5, height = 9) 

# read the pointwise 95% lower/upper bootstrap confidence limit of the functional parameter estimates
# for the covariate HEI
b_effects_l95_hei<- beta_boot_l95[37,]
b_effects_u95_hei<- beta_boot_u95[37,]

# consider the estimated functional parameters for the covariate HEI and their 
# functional lower and upper 95% Confidence limits upto t=0.97.
boot_l95_1eff_hei <- b_effects_l95_hei[,c(1:qrnt)]
boot_u95_1eff_hei <- b_effects_u95_hei[,c(1:qrnt)]

plotdf_hei<- data.frame(Estimate= as.numeric(t(b1_eff_hei)), 
                        u95CI= as.numeric(t(boot_u95_1eff_hei)), 
                        l95CI= as.numeric(t(boot_l95_1eff_hei)),
                        t=tt1)

# create a data frame which contains the nature of the differences considered in the 
# order which appears in the data frame 'plotdf_all'. The abbreviations for the ethnicities
# OH   : Other Hispanic
# MA   : Mexican American
# NHW  : Non-Hispanic White
# NHB  : Non-Hispanic Black
# NHA  : Non-Hispanic Asian
# ORIMR: Other Races Including Multi-Racial.

plotdf_hei<- data.frame(plotdf_hei, Coefficient=rep(c("Healthy Eating Index"), 
                                                    each= length(tt1)),
                        Ord = rep(c(1:nrow(b1_eff_hei)), each= length(tt1)))

gen_hei<- c(min(boot_l95_1eff_hei), max(boot_u95_1eff_hei))

# For better visual comparison, the y-axes were fixed, so that the disparity among 
# different ethnicities were studied. Here read from data the upper and lower limits
# of the y-axis and create an auxiliary dataset to fix the order of the panels to 
# appear in the figure.

blank_data_hei<- data.frame(Coefficient= rep(c("Healthy Eating Index"), each= 2), 
                            t = 0, Coefficients = gen_hei,
                            CI= rep(c("Lower 95%","Upper 95%"), nrow(b1_eff_hei)),
                            Ord=1)

########## HEI bootstrap plot ##############

ggplot(plotdf_hei, aes(x = t, y = Estimate)) +
  geom_ribbon(aes(ymin = l95CI, ymax = u95CI),
              fill = "#999999", alpha = 0.6) + 
  geom_line(size = .6, color = "black") +
  labs(x = "t", title = "Healthy Eating Index (HEI)") +
  theme(axis.text=element_text(size=100),
        axis.title=element_text(size=100)) +
  theme(text=element_text(size=100)) +
  geom_blank(data = blank_data_hei, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimates") + theme_bw() +
  geom_hline(yintercept = 0, col="black", linetype="dotted")

ggsave(file= "Rplot_betafunctional_linear_effects_differences_HEI_bootstrap.png", 
       width = 5.5, height = 3.5) 


#####################################################################################
# Code to create the figure: plot quantiles for disparity between predicted quantiles
# before and after projection into the 2-Wasserstein space
######

# read the predictions before and after projection
predictions<- as.matrix(read.csv("Predictions_before_projection.csv", 
                                 header=TRUE)[,-1])
projections<- as.matrix(read.csv("Predictions_after_projection.csv", 
                                 header =TRUE)[,-1])

# absolute differences between the predictions before and after projection
pred_diff <- abs(predictions - projections)

## summary across the quantiles 
pred_diff_min<- apply(pred_diff, 2, min)
pred_diff_max<- apply(pred_diff, 2, max)
pred_diff_1stq<- apply(pred_diff, 2, quantile, probs=0.25)
pred_diff_med<- apply(pred_diff, 2, median)
pred_diff_3rdq<- apply(pred_diff, 2, quantile, probs=0.75)
dat_combined<- rbind(pred_diff_1stq, pred_diff_med, pred_diff_3rdq, pred_diff_max)
row.names(dat_combined) <- c("25Quantiles", "Medians", "75Quantiles", "Maximums")
write.csv(dat_combined, "PLFSI_Projection_disparity.csv")

setEPS()
postscript(file = "Rplot_PLFSI_projection_disparity.eps", width =7, height = 6, paper = 'special')
par(mfrow= c(1,1))

plot(tt, pred_diff_3rdq, type= "l", lty=2, col="black", ylab="Absolute differences",
     xlab="Order of Quantiles", main="Absolute differences across quantiles", lwd=1.5)
lines(tt, pred_diff_med, type= "l", lty=7, col="black", lwd= 1.5)
lines(tt, pred_diff_1stq, type= "l", lty=3, col="black", lwd=1.5)
legend("topright", legend = c("3rd Quartile", "Median","1st Quartile"), 
       #col = c("red2", "black","blue"), 
       lty = c(2,7,3), lwd=1.5)

dev.off()
 
