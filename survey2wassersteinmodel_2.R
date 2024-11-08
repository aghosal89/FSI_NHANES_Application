## This function computes the metrics of model performance utilizing the survey weights

## Inputs:

#  1) formula           - a character representing the formula for regression.
#  2) data_analysis_svy - a survey GLM object. 
#  3) objectofda        - the functional data representation of the response quantiles
#  4) xout              - the covariate values where prediction is required.

## Outputs:

#  1) r^2            - the R^2 of the fitted model
#  2) betaj          - beta coefficients of the fitted regression model
#  3) predicciones   - the prediction of regression model after projection to 
#                      L^2-Wasserstein space.
#  4) residuos       - the residuals of the fitted model. 
#  5) r2vec          - a vector of length same as number of columns in predicciones 
#                      above, containing the multiple R-squared for each percentile.

## R libraries needed to run this function:

# survey

## Other functions from this repository needed to run this function:

# 'cuadratico.R'

survey2wassersteinmodel_2 <- function(formula=NULL, data_analysis_svy=NULL, 
                                    objetofda= NULL, xout=NULL) {
  
  prediciones_in = matrix(0, nrow= dim(data_analysis_svy$variables)[1], 
                       ncol= ncol(objetofda$data))
  
  residuos <- residuos2<- prediciones_in2 <- prediciones_in  # for the in-sample
  
  if(!is.null(xout)) {
  prediciones_out = matrix(0, nrow= dim(xout)[1], ncol= ncol(objetofda$data))
  }
  
  r2vec<- rep(NA, dim(prediciones_in)[2])
  formulaaux= paste("X", 1, sep="")
  formulaaux= paste(formulaaux, "~", sep="")
  
  formulafinal <- paste(formulaaux, formula, sep="")
  formulafinal <- as.formula(formulafinal)  
  
  #xout <- as.data.frame(xout)
  m2 = svyglm(formulafinal, design=data_analysis_svy, family=stats::gaussian())
  #pred<- as.matrix(predict.glm(m2, newdata = data.frame(xout), type="response"))
  
  #print(summary(m2))
  betaj= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
  #betaj_ucl= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
  #betaj_lcl= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
  #pred<- matrix(0, nrow= )
  
  # obtain the survey weights
  w= data_analysis_svy$variables$survey_wt
  w = w/sum(w)
  
  #setup parallel backend to use many processors
  t= proc.time()
  for(i in 1:ncol(objetofda$data)) {
    
    formulaaux= paste("X", i, sep="")
    formulaaux= paste(formulaaux,"~",sep="")
    
    formulafinal= paste(formulaaux, formula,sep="")
    formulafinal= as.formula(formulafinal)
    m=svyglm(formulafinal, design=data_analysis_svy, family=stats::gaussian())
    
    aux= as.numeric(m$residuals)
    betaj[,i] <- as.numeric(m$coefficients)
    
    prediciones_in[,i] <- as.numeric(m$fitted.values)
    
    if(!is.null(xout)) {
      prediciones_out[,i] <- as.matrix(predict.glm(m, 
                                newdata = data.frame(xout), type="response"))
    }
    residuos[,i]= aux
    
    # Compute the Multiple R-squared for each quantile:
    y_mean<- sum(w*objetofda$data[,i])
    TSS<- sum(w*(objetofda$data[,i] - y_mean)^2)
    if(TSS==0) {
      r2vec[i] <- 1
      SSE<- 0
    } else {
      SSE<- sum(w*(aux)^2)
      r2vec[i] <- 1- (SSE/TSS)
    }
  }
  print(t - proc.time())
  
  # check if the predictions are quantiles
  indt_in <- apply(prediciones_in, 1, function(x) {all(diff(x)>= 0)==FALSE})
  if(!is.null(xout)) {
    indt_out <- apply(prediciones_out, 1, function(x) {all(diff(x)>= 0)==FALSE})
  }
  
  for(j in 1:length(indt_in)) {
    if(indt_in[j]==TRUE) {
      prediciones_in2[j,]<- cuadratico(matrix(prediciones_in[j,], nrow = 1))
    } else {
      prediciones_in2[j,]<- prediciones_in[j,]
    }
  }
  
  if(!is.null(xout)) {
    for(j in 1:length(indt_out)) {
      if(indt_out[j]==TRUE) {
        prediciones_out[j,]<- cuadratico(matrix(prediciones_out[j,], nrow = 1))
      } else {
        prediciones_out[j,]<- prediciones_out[j,]
      }
    }
  }
  
  media1 <- apply(w*objetofda$data, 2, sum)
  for (i in 1:nrow(objetofda$data)) {
    residuos2[i,]<- objetofda$data[i,]- media1
  }
  
  disp1= apply(residuos2, 1, function(x){sum(x^2)})
  disp2= apply(residuos, 1, function(x){sum(x^2)})
  
  # The R-squared takes into account the survey weights
  
  r2 <- 1-(sum(disp2*w)/sum(disp1*w))
  
  rownames(betaj) <- names(m$coefficients)
  
  if(!is.null(xout)) {
  return(list("r2"=r2, "betaj"=betaj, "predicciones_In"= prediciones_in,
              "projection_In" = prediciones_in2,
              "predicciones_Out"= prediciones_out,
              "residuos"= residuos, "R2_vector"= r2vec))
  } else {
    return(list("r2"=r2, "betaj"=betaj, "predicciones_In"= prediciones_in,
                "projection_In" = prediciones_in2,
                "residuos"= residuos, "R2_vector"= r2vec))
  }
   
}

