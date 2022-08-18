
## This function computes the metrics of model performance utilizing the survey weights
## Inputs: 1) formula           - a character representing the formula for regression
##         2) data_analysis_svy - a survey GLM object 
##         3) objectofda        - the functional data representation of the response quantiles

## Outputs:1) r^2               - the R^2 of the fitted model
##         2) betaj             - beta coefficients of the fitted regression model
##         3) predicciones      - the prediction of regression model
##         4) residuos          - the residuals of the fitted model


survey2wassersteinmodel <- function(formula=formula, data_analysis_svy=data_analysis_svy, objetofda= objetofda) {
  library("survey")
  
  prediciones = matrix(0, nrow= dim(data_analysis_svy$variables)[1], ncol= ncol(objetofda$data))
  residuos = prediciones
  
  formulaaux= paste("X", 1, sep="")
  formulaaux= paste(formulaaux, "~", sep="")
  
  formulafinal= paste(formulaaux, formula,sep="")
  formulafinal= as.formula(formulafinal)  
  
  m2= svyglm(formulafinal,design=data_analysis_svy, family=stats::gaussian())
  
  betaj= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
  betaj_ucl= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
  betaj_lcl= matrix(0, nrow=length(m2$coefficients), ncol=ncol(objetofda$data))
  
  #setup parallel backend to use many processors
  t= proc.time()
  for(i in 1:ncol(objetofda$data)){
    
    formulaaux= paste("X", i, sep="")
    formulaaux= paste(formulaaux,"~",sep="")
    
    formulafinal= paste(formulaaux, formula,sep="")
    formulafinal= as.formula(formulafinal)
    m=svyglm(formulafinal, design=data_analysis_svy, family=stats::gaussian())
    
    aux= as.numeric(m$residuals)
    betaj[,i]= as.numeric(m$coefficients)
    ci<- confint(m)
    betaj_ucl[,i]= as.numeric(ci[,2])
    betaj_lcl[,i]= as.numeric(ci[,1])
    
    prediciones[,i]= as.numeric(m$fitted.values)
    residuos[,i]= aux
    
  }
  print(t - proc.time())
  
  pesosestandarizados= data_analysis_svy$variables$survey_wt
  pesosestandarizados= pesosestandarizados/sum(pesosestandarizados)
  
  media1= func.mean(objetofda)
  residuos2= objetofda-media1
  
  disp1= apply(residuos2$data, 1, function(x){sum(x^2)})
  disp2= apply(residuos, 1, function(x){sum(x^2)})
  
  # I will must update this with survey weights
  
  r2 <- 1-(sum(disp2*pesosestandarizados)/sum(disp1*pesosestandarizados))
  
  rownames(betaj) <- names(m$coefficients)
  
  return(list("r2"=r2,"betaj"=betaj, "beta_lcl"=betaj_lcl, "beta_ucl"=betaj_ucl, "predicciones"= prediciones,"residuos"= residuos))
  
}
