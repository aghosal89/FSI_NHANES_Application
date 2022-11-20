
# read the necessary libraries
library("latex2exp")
library("survey")
library("fda.usc")
library("viridis")
library("tikzDevice")
# Set library for R oprations
setwd("~/Documents/FSI/Application /alex/FSI_NHANES/Archive_PLFSI")

###########################################
# Codes to create figure 1 in the document
###########################################

datosacel<- read.csv("datosaritra.csv")
datosfda<- read.csv("Y.csv")[,-1]
obs= c()

for(i in 1:8){
  obs= c(obs,as.numeric(unlist(datosacel[i,])))
}

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



postscript(file = "grafico_quantiles_all.eps", width = 4, height = 4, paper = 'special')

grid= seq(0,1,length=500)

plot(grid, datosfda[1,], xlab = "t", col=1, ylim= c(0, 285),
     ylab="Empirical Quantiles", main="Physical activity representations", type = "l")
for(i in 2:nrow(datosfda)) {
  lines(grid, datosfda[i,], type = "l", col=i)
  
}
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

dev.off()


postscript(file = "grafico_quantiles_1.eps", width = 4, height = 4, paper = 'special')

plot(grid, quantile(obs, probs= grid,na.rm= TRUE), ylab = "Empirical Quantile",xlab="t", 
     main="Physical activity representation", type ="l", ylim= c(0, 285), lwd=2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

dev.off()

###########################################
# Codes to create figure 2 in the document
###########################################

# read the functional beta coefficients
b= read.csv("Output_Age20to80_noTAC_betas3.csv", header = T)[,-1]

# read the functional upper 95% CI of beta coefficients
b_u95= read.csv("Output_Age20to80_noTAC_betas_UCL3.csv", header = T)[,-1]

# read the functional lower 95% CI of beta coefficients
b_l95= read.csv("Output_Age20to80_noTAC_betas_LCL3.csv", header = T)[,-1]

# sequence of quantiles 
tt= seq(0,1,length=ncol(b))


# creating the plots 
qrnt <- floor(which(tt>=0.97))[1]
tt1<- tt[1:qrnt]
b1<- b[,c(1:qrnt)]
b_u95_1<- b_u95[,c(1:qrnt)]
b_l95_1<- b_l95[,c(1:qrnt)]
b1[1,]

plotdf_intercept<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                     ncol=1),t(b1[1,]), t(b_u95_1[1,]), t(b_l95_1[1,]))),
                t=rep(tt1, 4), CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                                                  each= length(tt1)))

plotdf_SexF<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                      ncol=1),t(b1[2,]), t(b_u95_1[2,]), t(b_l95_1[2,]))),
                       t=rep(tt1, 4), 
                       CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                               each= length(tt1)))

plotdf_OH<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                     ncol=1),t(b1[3,]), t(b_u95_1[3,]), t(b_l95_1[3,]))),
                       t=rep(tt1, 4), 
                       CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), each= length(tt1)))

plotdf_NHW<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                      ncol=1),t(b1[4,]), t(b_u95_1[4,]), t(b_l95_1[4,]))),
                        t=rep(tt1, 4), 
                        CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), each= length(tt1)))

plotdf_NHB<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                      ncol=1),t(b1[5,]), t(b_u95_1[5,]), t(b_l95_1[5,]))),
                        t=rep(tt1, 4), 
                        CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), each= length(tt1)))

plotdf_NHA<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                      ncol=1),t(b1[6,]), t(b_u95_1[6,]), t(b_l95_1[6,]))),
                        t=rep(tt1, 4), 
                        CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                        each= length(tt1)))

plotdf_NHM<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                      ncol=1),t(b1[7,]), t(b_u95_1[7,]), t(b_l95_1[7,]))),
                        t=rep(tt1, 4), 
                        CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                        each= length(tt1)))

plotdf_HEI<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                      ncol=1),t(b1[8,]), t(b_u95_1[8,]), t(b_l95_1[8,]))),
                        t=rep(tt1, 4), 
                        CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                        each= length(tt1)))

plotdf<- data.frame(rbind(plotdf_intercept, plotdf_SexF, plotdf_OH, plotdf_NHW,
                          plotdf_NHB, plotdf_NHA, plotdf_HEI), 
                  Coefficient=rep(c("Intercept", "Sex:Female","Other Hispanic",
                    "Non Hispanic White", "Non Hispanic Black",
                   "Non Hispanic Asian", "HEI"), each= nrow(plotdf_intercept)))

postscript(file = "Rplot_betafunctional_linear2.eps", width = 8, height = 8, paper = 'special')

ggplot(data = plotdf, aes(x=t, y=Coefficients, fill=CI)) +
  geom_path(aes(colour=as.factor(CI)), size=.9, alpha=1) +
  labs(colour = "Coefficient") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  facet_wrap(~Coefficient, scales = "free") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11)) +
  theme(text=element_text(size=11)) +
  ylab(label = "Estimated Coefficients") 

dev.off()

###########################################
# Codes to create figure 3 in the document
###########################################

# number of values for each of Age and BMI used to create the grid for figure 3
m2<- 500

# read the covaraiets in the PLFSI model
rnames<- read.csv("Output_Age20to80_noTAC_rnames3.csv", header = T)[,-1]

# read and manipulate the dataset
datos <- read.csv("datosalex(1).csv")
indices <- round(seq(188,687, length=500))

colnames(datos)[indices]= paste("X", 1:500, sep="")
datos= datos[, c(1:187, indices)]
datos<- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  # subset data for age 30-50 years
datos<- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  # subset data for BMI range 18.5-40
datosfda= datos[,188:687]
datosfda= as.matrix(datosfda)
datosfda[which(datosfda>280)]= datosfda[which(datosfda>280)-1]+0.5
datos[,188:687]= datosfda

tt= seq(0,1,length=500)
objetofda= fdata(datosfda,argvals = tt)

si_vars <- c("BMXBMI","RIDAGEYR.y")  # covariates for Single Index part
linear_vars <- c("RIAGENDR","RIDRETH3","HEI")  # covariates for Linear part
# question: why should we choose HEI if not significant. 
# scale 
linear_vars_numerical<- c("HEI")
linear_vars_categorical<- c("RIAGENDR","RIDRETH3")

# the survey variables
datosx <- cbind.data.frame(scale(datos[,si_vars]), HEI=scale(datos[,linear_vars_numerical]),
                           datos[,linear_vars_categorical], 
                           survey_wt= datos[,"wtmec4yr_adj_norm"],
                           survey_id=datos[,"SDMVPSU"], 
                           survey_strata= datos[,"SDMVSTRA"])

datosx$RIAGENDR<- as.factor(datosx$RIAGENDR)
datosx$RIDRETH3<- as.factor(datosx$RIDRETH3)
datosx$survey_strara<- as.factor(datosx$survey_strata)
datosx$survey_id<- as.factor(datosx$survey_id)

# consider equidistant grid of variables BMXBMI, RIDAGEYR.y over their range. 
# For each variables consider of length = 100. 

x_ <- seq(min(datos$BMXBMI), max(datos$BMXBMI), len=m2)
y_ <- seq(min(datos$RIDAGEYR.y), max(datos$RIDAGEYR.y), len=m2)

temp_si_data<- data.frame(BMXBMI= seq(min(datosx$BMXBMI),max(datosx$BMXBMI), len=m2), 
                          RIDAGEYR.y=seq(min(datosx$RIDAGEYR.y),max(datosx$RIDAGEYR.y),len=m2))

# consider the variables on a 2 dimensional grid to construct the contour plot. 
tsd_grid<- expand.grid(temp_si_data)
tsd_grid_level<- expand.grid(data.frame(BMI=x_, Age=y_))
x<- xin<- tsd_grid

# read the estimated index parameter
th<- matrix(read.csv("Theta_Hat.csv")[,-1], ncol =1)
Theta_hat<- th

q50<- floor(which(tt>=0.50))[1]   # check index for 50th percentile 
q75<- floor(which(tt>=0.75))[1]   # check index for 75th percentile
q80<- floor(which(tt>=0.80))[1]   # check index for 80th percentile
q90<- floor(which(tt>=0.90))[1]   # check index for 90th percentile

sp <- 4 
dfs<- 9

u_ <- as.matrix(tsd_grid)%*%th
bs_data1 <- data.frame(splines::bs(x=u_, degree=sp, df=dfs))

colnames(bs_data1)<- paste("BS", seq(1:ncol(bs_data1)), sep ="")

# use survey weights for median
x_males_eth1<- cbind.data.frame(Intercept=1, HEI=median(datosx$HEI), bs_data1)

# removing gender and ethnicity columns since using the baselines, i.e. males 
# for gender and Ethnicity1 for Ethnicity 


# creating the predictions for all quantiles on tt
yh_males_eth1<- matrix(NA, nrow = dim(x_males_eth1)[1], ncol = length(tt))

for (r in 1:length(tt)) {
  yh_males_eth1[,r] <- as.matrix(x_males_eth1)%*% matrix(b[,r][-c(2:7)],ncol=1)
}

source("cuadratico.R", local = knitr::knit_global())
t<- proc.time()
indt <- apply(yh_males_eth1, 1, function(x) {all(diff(x)>= 0)==FALSE})
for(j in 1:length(indt)) {
  if(indt[j]==TRUE){
    yh_males_eth1[j,]<- cuadratico(matrix(yh_males_eth1[j,], nrow = 1))
  }
}
proc.time() - t

# save the quantiles for future use 
Yh<- data.frame(yh_males_eth1[,c(q50, q75, q80, q90)])
colnames(Yh)<- c("Q50","Q75","Q80","Q90")
write.csv(yh_males_eth1, "Pedicted_percentiles.csv")

yh_males_eth1$Integral<- sapply(1:nrow(yh_males_eth1), 
        function(i) fdadensity:::trapzRcpp(X = tt, Y = as.numeric(yh_males_eth1[i,])))
  
Yh<- data.frame(read.csv("Pedicted_percentiles.csv")[,c(q50, q75, q80, q90)])

df_q50<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                       AGE=tsd_grid_level$Age,
                       Quantiles=Yh$Q50)

df_q75<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          Quantiles=Yh$Q75)

df_q80<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          Quantiles=Yh$Q80)

df_q90<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          Quantiles=Yh$Q90)

setEPS()
postscript(file = "Rplot_quantile50_prediction.eps", width = 6.5, height = 5.5, paper = 'special')

ggplot(df_q50, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
      axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("50th Quantile")

dev.off()

postscript(file = "Rplot_quantile75_prediction.eps", width = 6.5, height = 5.5, paper = 'special')

ggplot(df_q75, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("75th Quantile")

dev.off()

postscript(file = "Rplot_quantile80_prediction.eps", width = 6.5, height = 5.5, paper = 'special')

ggplot(df_q80, aes(BMI, AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("80th Quantile")

dev.off()

postscript(file = "Rplot_quantile90_prediction.eps", width = 6.5, height = 5.5, paper = 'special')

ggplot(df_q90, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("90th Quantile")

dev.off()





###########################################
# Codes to create figure 4 in the document
###########################################

# number of values for each of Age and BMI used to create the grid for figure 3
m2<- 60

# read the functional beta coefficients
b= read.csv("Output_Age20to80_noTAC_betas3.csv", header = T)[,-1]
#beta_functional_coeffs<- b

# read the covaraiets in the PLFSI model
rnames<- read.csv("Output_Age20to80_noTAC_rnames3.csv", header = T)[,-1]

# read and manipulate the dataset
datos <- read.csv("datosalex(1).csv")
indices <- round(seq(188,687, length=500))

colnames(datos)[indices] = paste("X", 1:500, sep="")
datos= datos[, c(1:187, indices)]

# subset data for age 20-80 years
datos<- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  

# subset data for BMI 15-40
datos<- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  
datos<- datos[-c(734, 1216, 2094, 2699, 2958, 3092, 3933, 3963, 4232), ]

datosfda= datos[,188:687]
datosfda= as.matrix(datosfda)
datosfda[which(datosfda>280)]= datosfda[which(datosfda>280)-1]+0.5
datos[,188:687]= datosfda

tt= seq(0,1,length=500)
objetofda= fdata(datosfda,argvals = tt)

si_vars <- c("BMXBMI","RIDAGEYR.y")  # covariates for Single Index part
linear_vars <- c("RIAGENDR","RIDRETH3","HEI")  # covariates for Linear part
# question: why should we choose HEI if not significant. 
# scale 
linear_vars_numerical<- c("HEI")
linear_vars_categorical<- c("RIAGENDR","RIDRETH3")

# the survey variables
datosx <- cbind.data.frame(scale(datos[,si_vars]), HEI=scale(datos[,linear_vars_numerical]),
                           datos[,linear_vars_categorical], 
                           survey_wt= datos[,"wtmec4yr_adj_norm"],
                           survey_id=datos[,"SDMVPSU"], 
                           survey_strata= datos[,"SDMVSTRA"])

datosx$RIAGENDR<- as.factor(datosx$RIAGENDR)
datosx$RIDRETH3<- as.factor(datosx$RIDRETH3)
datosx$survey_strara<- as.factor(datosx$survey_strata)
datosx$survey_id<- as.factor(datosx$survey_id)

# consider equidistant grid of variables BMXBMI, RIDAGEYR.y over their range. 
# For each variables consider of length = m2. 


x_ <- seq(min(datos$BMXBMI), max(datos$BMXBMI), len=m2)
y_ <- seq(min(datos$RIDAGEYR.y), max(datos$RIDAGEYR.y), len=m2)

temp_si_data<- data.frame(BMXBMI= seq(min(datosx$BMXBMI),max(datosx$BMXBMI), len=m2), 
                          RIDAGEYR.y=seq(min(datosx$RIDAGEYR.y),max(datosx$RIDAGEYR.y),len=m2))

# consider the variables on a 2 dimensional grid to construct the contour plot. 
tsd_grid<- expand.grid(temp_si_data)
tsd_grid_level<- expand.grid(data.frame(BMI=x_, Age=y_))

# read the estimated index parameter
th<- matrix(read.csv("Theta_Hat.csv")[,-1], ncol =1)
Theta_hat<- th

u_t
u_<- as.matrix(tsd_grid)%*%th
# fix the order of the B-spline basis 
sp <- 4 

# fix the degrees of freedom
dfs<- 9


#####################
# The following are the bondary knots and the rest of the knot sequence for the 
# spline basis computation

aspl<- splines2::dbs(x=u_, derivs = 1L, df= dfs, knots = NULL, degree = sp, 
                     intercept = FALSE)

drv<- matrix(NA, nrow =nrow(aspl), ncol =length(tt))
for (r in 1:length(tt)) {
  drv[,r] <- aspl%*% matrix(b[,r][-c(1:8)],ncol=1)
}

q50<- floor(which(tt>=0.50))[1]   # check index for 50th percentile 
q75<- floor(which(tt>=0.75))[1]   # check index for 75th percentile
q80<- floor(which(tt>=0.80))[1]   # check index for 80th percentile
q90<- floor(which(tt>=0.90))[1]   # check index for 90th percentile


# univariate
u_2<- as.matrix(cbind(datosx$BMXBMI, datosx$RIDAGEYR.y))%*%th
u_2 <- sort(u_2)
#u_2<- seq(min(u_d), max(u_d), len=200)
aspl1.1<- splines::bs(u_2, df= dfs, knots = NULL, degree = sp, 
                      intercept = FALSE)

bknots<- attributes(aspl1.1)$Boundary.knots
knots <- attributes(aspl1.1)$knots
knots_all<- c(bknots[1], knots, bknots[2])

aspl2<- splines2::dbs(x=u_2, derivs = 1L, df= dfs, knots = knots, 
                      degree = sp, Boundary.knots= bknots,
                      intercept = FALSE)


# check index for 15th percentile of the Single Index
row_ind_fbo<- floor(which(u_2>=knots_all[2]))[1]
  
# check index for 83rd percentile of the Single Index
row_ind_lbo<- floor(which(u_2>=knots_all[length(knots_all)-1]))[1]
  
drv2<- matrix(NA, nrow =nrow(aspl2), ncol =length(tt))
for (r in 1:length(tt)) {
  drv2[,r] <- aspl2%*% matrix(b[,r][-c(1:8)],ncol=1)
}

df_in <- drv2[c(row_ind_fbo:row_ind_lbo),]
u_2_in<- u_2[c(row_ind_fbo:row_ind_lbo)]
u_2_boundary_l<- u_2[c(1:row_ind_fbo)]
df_boundary_l <- drv2[c(1:row_ind_fbo),]
u_2_boundary_r <- u_2[c(row_ind_lbo:length(u_2))]
df_boundary_r <- drv2[c(row_ind_lbo:length(u_2)),]


setEPS()
postscript(file = "Rplot_qderiv_SI.eps", width = 5, height = 5, paper = 'special')

plot(u_2_in, df_in[,q50], type ="l", ylab ="Derivative of Quantiles", 
     xlab="Single Index", ylim=c(-10,25), xlim =c(-2.5,2.7), lwd=3)
lines(u_2_in, df_in[,q75], col="green", lwd=3)
lines(u_2_in, df_in[,q80], col="blue", lwd=3)
lines(u_2_in, df_in[,q90], col="yellow2", lwd=3)

lines(u_2_boundary_l, df_boundary_l[,q50], lwd=3, lty=3)
lines(u_2_boundary_l, df_boundary_l[,q75], col="green", lwd=3, lty =3)
lines(u_2_boundary_l, df_boundary_l[,q80], col="blue", lwd=3, lty =3)
lines(u_2_boundary_l, df_boundary_l[,q90], col="yellow2", lwd=3, lty =3)

lines(u_2_boundary_r, df_boundary_r[,q50], lwd=3, lty=3)
lines(u_2_boundary_r, df_boundary_r[,q75], col="green", lwd=3, lty =3)
lines(u_2_boundary_r, df_boundary_r[,q80], col="blue", lwd=3, lty =3)
lines(u_2_boundary_r, df_boundary_r[,q90], col="yellow2", lwd=3, lty =3)

for (j in 1:length(knots_all)) {
  abline(v=knots_all[j], lty=3, col="red")
}

legend(-1, 25, legend=c("t=0.50","t=0.75","t=0.80","t=0.90"),
       col=c("black","green","blue","yellow2"), 
       lwd =c(3,3,3,3),
       lty=c(1,1,1,1), cex=1.2, text.font=10
       #, bg='lightgrey'
)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

dev.off()


postscript(file = "Rplot_histogram_ofSI.eps", width = 4.5, height = 4.5, paper = 'special')

hist(u_d, freq=FALSE, main="Histogram of Single Index", xlab= "")
dev.off()

postscript(file = "Rplot_kdens_SI.eps", width = 4.5, height = 4.5, paper = 'special')

plot(density(u_d), main="")
abline(v=0, lty=3, col="red")
dev.off()



# computation of the integral of the beta coefficients:



