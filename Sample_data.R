######
# Select the observations for the training/testing splits for the data 
######

# Read the libraries
library('sampling')

# Set your R directory
setwd("~/Documents/2024_09_21_Aritra_codeForMainAnalysis")

# read the dataset
datos= read.csv("datosalex(1).csv")

# subset data for age 20-80 years
datos <- subset(datos, datos$RIDAGEYR.x >= 20 & datos$RIDAGEYR.x <= 80)  

# subset data for BMI in the range 18.5 - 40
datos <- subset(datos, datos$BMXBMI >= 18.5 & datos$BMXBMI <= 40)  

n <- nrow(datos)
ns <- round(n*.10, 0)
seed <- 30190
Sample_data <- NA 
s1 <- 40 # number of training-testing splits

for (i in 1:s1) {
  a<- srswor(ns, n)
  
  set.seed(seed)
  insampl <- which(a==1)
  Sample_data <- rbind(Sample_data, insampl)
  seed <- seed+5
  
}

Sample_data <- Sample_data[-1,]
rownames(Sample_data) <- NULL

write.csv(Sample_data, "Sample_data.csv")
