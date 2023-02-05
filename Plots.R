
# read the necessary libraries
library("latex2exp")
library("survey")
library("fda.usc")
library("viridis")
library("ggplot2")
library("energy")

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

# read the estimated main effects, 95% upper and lower confidence limits pointwise
b_effects<- read.csv("Output_Age20to80_noTAC_betaeffects.csv")[,-1]
b_effects_u95<- read.csv("Output_Age20to80_noTAC_betaeffects_ucl.csv")[,-1]
b_effects_l95<- read.csv("Output_Age20to80_noTAC_betaeffects_lcl.csv")[,-1]

# sequence of order of quantiles 
tt= seq(0,1,length=ncol(b_effects))

# creating the plots for order of quantile in the range [0, 0.97].
qrnt <- floor(which(tt>0.97))[1]
tt1<- tt[1:qrnt]

b1_eff<- b_effects[,c(1:qrnt)]
b_u95_1eff<- b_effects_u95[,c(1:qrnt)]
b_l95_1eff<- b_effects_l95[,c(1:qrnt)]


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

############
# Plot for the added effects

plotdf_SexM_E1<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[1,]), t(b_u95_1eff[1,]), t(b_l95_1eff[1,]))),
                  t=rep(tt1, 4), CI= rep(c("0","Estimate","Upper 95%","Lower 95%"),
                  each= length(tt1)))

plotdf_SexF_E1<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[2,]), t(b_u95_1eff[2,]), t(b_l95_1eff[2,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                  each= length(tt1)))

plotdf_SexM_OH<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[3,]), t(b_u95_1eff[3,]), t(b_l95_1eff[3,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                  each= length(tt1)))

plotdf_SexF_OH<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[4,]), t(b_u95_1eff[4,]), t(b_l95_1eff[4,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                  each= length(tt1)))

plotdf_SexM_NHW<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                   ncol=1),t(b1_eff[5,]), t(b_u95_1eff[5,]), t(b_l95_1eff[5,]))),
                   t=rep(tt1, 4), 
                   CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                            each= length(tt1)))

plotdf_SexF_NHW<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[6,]), t(b_u95_1eff[6,]), t(b_l95_1eff[6,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                            each= length(tt1)))

plotdf_SexM_NHB<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[7,]), t(b_u95_1eff[7,]), t(b_l95_1eff[7,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                           each= length(tt1)))

plotdf_SexF_NHB<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[8,]), t(b_u95_1eff[8,]), t(b_l95_1eff[8,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                                    each= length(tt1)))

plotdf_SexM_NHA<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[9,]), t(b_u95_1eff[9,]), t(b_l95_1eff[9,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                                    each= length(tt1)))

plotdf_SexF_NHA<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                  ncol=1),t(b1_eff[10,]), t(b_u95_1eff[10,]), t(b_l95_1eff[10,]))),
                  t=rep(tt1, 4), 
                  CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                                    each= length(tt1)))

plotdf_SexM_ORIMR<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                       ncol=1),t(b1_eff[11,]), t(b_u95_1eff[11,]), t(b_l95_1eff[11,]))),
                       t=rep(tt1, 4), 
                       CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                                    each= length(tt1)))

plotdf_SexF_ORIMR<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                     ncol=1),t(b1_eff[12,]), t(b_u95_1eff[12,]), t(b_l95_1eff[12,]))),
                     t=rep(tt1, 4), 
                     CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                                    each= length(tt1)))

plotdf_HEI<- data.frame(Coefficients= as.vector(rbind(matrix(0, nrow=length(tt1), 
                    ncol=1),t(b1_eff[13,]), t(b_u95_1eff[13,]), t(b_l95_1eff[13,]))),
                            t=rep(tt1, 4), 
                          CI= rep(c("0","Estimate","Upper 95%","Lower 95%"), 
                                    each= length(tt1)))

plotdf <- data.frame(rbind(plotdf_SexM_E1, plotdf_SexF_E1, 
                     plotdf_SexM_OH,  plotdf_SexF_OH,
                     plotdf_SexM_NHW, plotdf_SexF_NHW,
                     plotdf_SexM_NHB, plotdf_SexF_NHB,
                     plotdf_SexM_NHA, plotdf_SexF_NHA,
                     plotdf_SexM_ORIMR, plotdf_SexF_ORIMR,
                     plotdf_HEI), 
                     Coefficient=rep(c("Male,Mexican American",
                     "Female,Mexican American", 
                     "Male,Other Hispanic",
                     "Female,Other Hispanic",
                     "Male,Non-Hispanic White",
                     "Female,Non-Hispanic White", 
                     "Male,Non-Hispanic Black",
                     "Female,Non-Hispanic Black", 
                     "Male,Non Hispanic Asian",
                     "Female,Non Hispanic Asian", 
                     "Male,Other Races(MR)",
                     "Female,Other Races(MR)", 
                     "HEI"), each= nrow(plotdf_HEI)),
                     Ord = rep(c(1:13), each= nrow(plotdf_HEI)))

# minimum/maximum for all the gender, ethnicities
gen_eth_min<- min(b_l95_1eff[1:12,])
gen_eth_max<- max(b_u95_1eff[1:12,])

blank_data <- data.frame(Coefficient=rep(c("Male,Mexican American",
                                           "Female,Mexican American", 
                                           "Male,Other Hispanic",
                                           "Female,Other Hispanic",
                                           "Male,Non-Hispanic White",
                                           "Female,Non-Hispanic White", 
                                           "Male,Non-Hispanic Black",
                                           "Female,Non-Hispanic Black", 
                                           "Male,Non Hispanic Asian",
                                           "Female,Non Hispanic Asian", 
                                           "Male,Other Races(MR)",
                                           "Female,Other Races(MR)", 
                                           "HEI"), each= 2), 
                        t = 0, Coefficients = 
                          c(rep(c(gen_eth_min, gen_eth_max), 12),
                            min(b_l95_1eff[13,]), max(b_u95_1eff[13,])),
                         CI= rep(c("Lower 95%","Upper 95%"), 13),
                         Ord=1)

ggplot(data = plotdf, aes(x=t, y=Coefficients, fill=CI)) +
  geom_path(aes(colour= as.factor(CI)), size=.9, alpha=1) +
  labs(colour = "Coefficient") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  facet_wrap( ~reorder(Coefficient, Ord), scales = "free_y") +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=11)) +
  theme(text=element_text(size=11)) +
  geom_blank(data = blank_data, aes(x = t, y = Coefficients)) +
  ylab(label = "Estimated Coefficients") + theme_bw() 

ggsave(file= "Rplot_betafunctional_linear2_effects.eps", width = 9, height = 10) 


###########################################
# Codes to create figure 3 in the document
###########################################

# number of values for each of Age and BMI used to create the grid for figure 3
m2<- 500

# read the covaraiets in the PLFSI model
rnames<- read.csv("Output_Age20to80_noTAC_rnames.csv", header = TRUE)[,-1]
b<- read.csv("Output_Age20to80_noTAC_betas.csv", header = TRUE)[,-1]

# consider equidistant grid of variables BMXBMI, RIDAGEYR.y over their range. 
# For each variables consider of length = 100. 

x_ <- seq(min(datos$BMXBMI), max(datos$BMXBMI), len=m2)
y_ <- seq(min(datos$RIDAGEYR.y), max(datos$RIDAGEYR.y), len=m2)

temp_si_data<- data.frame(BMXBMI= seq(min(datosx$BMXBMI), max(datosx$BMXBMI), len=m2), 
                          RIDAGEYR.y=seq(min(datosx$RIDAGEYR.y),max(datosx$RIDAGEYR.y),len=m2))

# consider the variables on a 2 dimensional grid to construct the contour plot. 
tsd_grid<- expand.grid(temp_si_data)
tsd_grid_level<- expand.grid(data.frame(BMI=x_, Age=y_))
x<- xin<- tsd_grid

# read the estimated index parameter
th <- matrix(read.csv("Theta_Hat.csv")[,-1], ncol =1)
Theta_hat<- th

q50<- floor(which(tt>=0.50))[1]   # check index for 50th percentile 
q75<- floor(which(tt>=0.75))[1]   # check index for 75th percentile
q90<- floor(which(tt>=0.90))[1]   # check index for 80th percentile
q97<- floor(which(tt>=0.97))[1]   # check index for 90th percentile

sp <- 4 
dfs<- 9

u_ <- as.matrix(tsd_grid)%*%th
bs_data1 <- data.frame(splines2::dbs(x=u_, degree=sp, df=dfs))

colnames(bs_data1)<- paste("BS", seq(1:ncol(bs_data1)), sep ="")

# use survey weights for median
x_males_eth1<- cbind.data.frame(Intercept=1, HEI=median(datosx$HEI), bs_data1)

# removing gender and ethnicity columns since using the baselines, i.e. males 
# for gender and Ethnicity1 for Ethnicity 


# creating the predictions for all quantiles on tt
yh_males_eth1<- matrix(NA, nrow = dim(x_males_eth1)[1], ncol = length(tt))

for (r in 1:length(tt)) {
  yh_males_eth1[,r] <- as.matrix(x_males_eth1)%*% matrix(b[,r][c(1,8:17)],ncol=1)
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
write.csv(yh_males_eth1, "Pedicted_transport_percentiles.csv")
Yh<- as.data.frame(yh_males_eth1[,c(q50, q75, q90, q97)])
colnames(Yh)<- c("Q50","Q75","Q90","Q97")


df_q50<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                       AGE=tsd_grid_level$Age,
                       Quantiles=Yh$Q50)

df_q75<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          Quantiles=Yh$Q75)

df_q90<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          Quantiles=Yh$Q90)

df_q97<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          Quantiles=Yh$Q97)

df_qint<- cbind.data.frame(BMI= tsd_grid_level$BMI,
                           AGE= tsd_grid_level$Age,
                           Integral = yh_mth1_int)

ggplot(df_q50, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
      axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("50th Quantile")
ggsave(file = "Rplot_quantile50_prediction.eps", width = 6.5, height = 5.5)


ggplot(df_q75, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("75th Quantile")
ggsave(file = "Rplot_quantile75_prediction.eps", width = 6.5, height = 5.5)


ggplot(df_q90, aes(BMI, AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("90th Quantile")
ggsave(file = "Rplot_quantile90_prediction.eps", width = 6.5, height = 5.5)


ggplot(df_q97, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = Quantiles)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("97th Quantile")
ggsave(file = "Rplot_quantile97_prediction.eps", width = 6.5, height = 5.5)


########
# Compute the integral of the quantiles

yh_males_eth1_int<- sapply(1:nrow(yh_males_eth1), function(i) {
  fdadensity:::trapzRcpp(X=tt, Y= yh_males_eth1[i,])
})

df_int<- cbind.data.frame(BMI=tsd_grid_level$BMI,
                          AGE=tsd_grid_level$Age,
                          Integrals=yh_males_eth1_int)

ggplot(df_int, aes(x=BMI, y=AGE)) +
  geom_raster(aes(fill = Integrals)) +
  scale_fill_viridis() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15)) +
  theme(text=element_text(size=15)) +
  ggtitle("Integated Quantiles")
ggsave(file = "Rplot_quantile_integrals.eps", width = 6.5, height = 5.5)

###########################################
# Codes to create figure 4 in the document
###########################################

# number of values for each of Age and BMI used to create the grid for figure 3
m2<- 60

# read the covaraiets in the PLFSI model
rnames<- read.csv("Output_Age20to80_noTAC_rnames.csv", header = T)[,-1]

# consider equidistant grid of variables BMXBMI, RIDAGEYR.y over their range. 
# For each variables consider of length = m2. 

x_ <- seq(min(datos$BMXBMI), max(datos$BMXBMI), len=m2)
y_ <- seq(min(datos$RIDAGEYR.y), max(datos$RIDAGEYR.y), len=m2)

temp_si_data<- data.frame(BMXBMI= seq(min(datosx$BMXBMI),max(datosx$BMXBMI), len=m2), 
                          RIDAGEYR.y=seq(min(datosx$RIDAGEYR.y),max(datosx$RIDAGEYR.y),len=m2))

# consider the variables on a 2 dimensional grid to construct the contour plot. 
tsd_grid<- expand.grid(temp_si_data)
tsd_grid_level<- expand.grid(data.frame(BMI=x_, Age=y_))

# univariate
u_2<- as.matrix(cbind(datosx$BMXBMI, datosx$RIDAGEYR.y))%*%th
u_2 <- sort(u_2)

aspl1.1<- splines2::dbs(u_2, df= dfs, knots = NULL, degree = sp, 
                      intercept = FALSE)

#####################
# The following are the bondary knots and the rest of the knot sequence for the 
# spline basis computation

bknots<- attributes(aspl1.1)$Boundary.knots
knots <- attributes(aspl1.1)$knots
knots_all<- c(bknots[1], knots, bknots[2])

aspl2<- splines2::dbs(x=u_2, derivs = 1L, df= dfs, knots = knots, 
                      degree = sp, Boundary.knots= bknots,
                      intercept = FALSE)

# check index for 16.6667th percentile of the Single Index
row_ind_fbo<- floor(which(u_2>=knots_all[2]))[1]
  
# check index for 83.3333rd percentile of the Single Index
row_ind_lbo<- floor(which(u_2>=knots_all[length(knots_all)-1]))[1]
  
drv2<- matrix(NA, nrow =nrow(aspl2), ncol =length(tt))
for (r in 1:length(tt)) {
  drv2[,r] <- aspl2%*% matrix(b[,r][c(9:17)],ncol=1)
}

df_in <- drv2[c(row_ind_fbo:row_ind_lbo),]
u_2_in<- u_2[c(row_ind_fbo:row_ind_lbo)]
u_2_boundary_l<- u_2[c(1:row_ind_fbo)]
df_boundary_l <- drv2[c(1:row_ind_fbo),]
u_2_boundary_r <- u_2[c(row_ind_lbo:length(u_2))]
df_boundary_r <- drv2[c(row_ind_lbo:length(u_2)),]


setEPS()
postscript(file = "Rplot_qderiv_SI.eps", width = 5, height = 5, paper = 'special')

#"#E69F00", "#56B4E9", "#009E73", "#0072B2"

plot(u_2_in, df_in[,q50], type ="l", ylab ="Derivative of Quantiles", 
     xlab="Single Index", ylim=c(-10,25), xlim =c(-2.5,2.7), lwd=1.5, col="#000000")
lines(u_2_in, df_in[,q75], col="#D55E00", lwd=2)
lines(u_2_in, df_in[,q90], col="#009E73", lwd=2.5)
lines(u_2_in, df_in[,q97], col="#0072B2", lwd=3)

lines(u_2_boundary_l, df_boundary_l[,q50], lwd=1.5, lty=3, col="#000000")
lines(u_2_boundary_l, df_boundary_l[,q75], col="#D55E00", lwd=2, lty =3)
lines(u_2_boundary_l, df_boundary_l[,q90], col="#009E73", lwd=2.5, lty =3)
lines(u_2_boundary_l, df_boundary_l[,q97], col="#0072B2", lwd=3, lty =3)

lines(u_2_boundary_r, df_boundary_r[,q50], lwd=1.5, lty=3, col="#000000")
lines(u_2_boundary_r, df_boundary_r[,q75], col="#D55E00", lwd=2, lty =3)
lines(u_2_boundary_r, df_boundary_r[,q90], col="#009E73", lwd=2.5, lty =3)
lines(u_2_boundary_r, df_boundary_r[,q97], col="#0072B2", lwd=3, lty =3)

for (j in 1:length(knots_all)) {
  abline(v=knots_all[j], lty=3, col="red")
}

legend(-1, 25, legend=c("t=0.50","t=0.75","t=0.90","t=0.97"),
       col=c("#000000", "#D55E00", "#009E73", "#0072B2"), 
       lwd=c(1.5,2,2.5,3), lty=c(1,1,1,1), cex=1.2, text.font=10, bg='lightgrey')

#grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
#     lwd = par("lwd"), equilogs = TRUE)

dev.off()


postscript(file = "Rplot_histogram_ofSI.eps", width = 4.5, height = 4.5, paper = 'special')

hist(u_2, freq=FALSE, main="Histogram of Single Index", xlab= "")
dev.off()



# computation of the integral of the beta coefficients:
###########################################
# Codes to create the figure 5 in document
###########################################

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

setEPS()
postscript(file= "clusteringresults.eps", horizontal = FALSE, onefile = FALSE, 
           paper = "special", height = 6, width = 7)
par(mfrow= c(2,3))
plot(cluster1[,1:qrnt], main="Cluster 1", ylim= c(-30,95), xlim=c(0,1.001),
     ylab="Residuals")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

plot(cluster2[,1:qrnt], main="Cluster 2",ylim= c(-30,95), xlim=c(0,1.001),
     ylab="Residuals")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

plot(cluster3[,1:qrnt], main="Cluster 3",ylim= c(-30,95), xlim=c(0,1.001)
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

