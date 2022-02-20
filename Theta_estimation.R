
###### survey estimation of theta #########

library("tidyverse")
library("survey")
library("fda.usc")
library("sandwich")

# Aritra's library for R oprations
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/FSI/Application /alex/FSI_NHANES")

datos= read.csv("datosalex(1).csv") 

colnames(datos)[188:687]= paste("X", 1:500, sep="")

datosfda <- datos[,188:687]
datosfda= as.matrix(datosfda)
datosfda[which(datosfda>280)] <- datosfda[which(datosfda>280)-1]+0.5
datos[,188:687]<- datosfda


ts = seq(0, 1, length=500)
objetofda= fdata(datosfda,argvals = ts)

etiquetaredad= function(x){
  if(x>=18&x<=39){
    return(1)
  }else if(x>=40&x<=60){
    return(2)
  }else if(x>60){
    return(3)
  }else{
    return(0)
  } 
}

edadcategorica= sapply(datos$RIDAGEYR.x, etiquetaredad)
datos$edadcategorica= edadcategorica


data_analysis_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, weights= ~wtmec4yr_adj_norm,
                               data = datos,nest = TRUE)
data_analysis_svy$variables$RIDRETH3= as.factor(data_analysis_svy$variables$RIDRETH3)

## creating a grid of starting values 

numerical<- c("BMXBMI","HEI","RIDAGEYR.y") # numerical covariates
categorical<- c("RIAGENDR","RIDRETH3","edadcategorica") # categorical covariates

# compute dimension of covariates
p <- length(numerical)
nsp <- 3
# specified spacing between staring points in each coordinate
spc <- pi/nsp  
# equally spaced starting points in polar coordinates
f <- lapply(1:(p - 1), function(j) seq(-pi/2 + spc/2, pi/2 - spc/2, by = spc)) 
# create grid of starting values
etaStart <- as.matrix(expand.grid(f))  
#L=0
#if(L != 0) {
#  smp <- sample.int(nrow(etaStart), min(L, nrow(etaStart)))
#  etaStart <- etaStart[smp,]
#}

## modifying the inputs of the algorithm
datos$RIAGENDR<- as.factor(datos$RIAGENDR)
datos$RIDRETH3<- as.factor(datos$RIDRETH3)
datos$edadcategorica<- as.factor(datos$edadcategorica)


## defining the options of optimization as output
optInf <- list()

optim_optns <- list(factr = 1e11, maxit = 100)
WnMin <- rep(NA, nrow(etaStart))
etaMin <- matrix(NA, nrow = nrow(etaStart), ncol = p - 1)
sp<- 4

wn_cost(et=etaStart[6,], datosfda=datosfda, categorical=categorical, 
        numerical=numerical, ts=ts, sp=4, datos =datos, 
        dfs=sp+5)


###### Run main optimization ######
for (k in 1:nrow(etaStart)) {
  print(k)
  WnOpt <- optim(par = etaStart[k, ], fn = wn_cost, method = "L-BFGS-B",
                 lower = -pi/2, upper = pi/2, 
                 control = optim_optns, 
                 datosfda=datosfda, categorical=categorical, 
                 numerical=numerical, ts=ts, sp=4, datos =datos, 
                 dfs=sp+5)
  
  optInf[[k]] <- WnOpt
  
  WnMin[k] <- WnOpt$value
  etaMin[k, ] <- WnOpt$par
  
}



