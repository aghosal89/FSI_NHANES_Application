## This is a function to compute the p-dimensional cartesian coordinates from 
## given polar coordinates in (p-1) dimension and radius.

## Inputs: 

# 1) eta  - a (p-1) dimensional vector.
# 2) r    - radius of the polar coordinates (set default, r=1).

# codes below were written for p= 2,3,4,5; which can be generalized.

polar2cart <- function(eta, r  = 1) {
  # read the length (p-1)
  s <- length(eta)
  # define vector for storing cartesian coordinates
  theta <- vector(length = s+1)
  
  # compute and return cartesian coordinates for a given p.
  if(s==4) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])*cos(eta[4])
    theta[2]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])*sin(eta[4])
    theta[3]<- r*cos(eta[1])*cos(eta[2])*sin(eta[3])
    theta[4]<- r*cos(eta[1])*sin(eta[2])
    theta[5]<- r*sin(eta[1])
    return(theta)
  }
  
  if(s==3) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])
    theta[2]<- r*cos(eta[1])*cos(eta[2])*sin(eta[3])
    theta[3]<- r*cos(eta[1])*sin(eta[2])
    theta[4]<- r*sin(eta[1])
    return(theta)
    
  }
  if(s==2) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])
    theta[2]<- r*cos(eta[1])*sin(eta[2])
    theta[3]<- r*sin(eta[1])
    return(theta)
  }
  
  if(s==1) {
    theta[1]<- r*cos(eta[1])
    theta[2]<- r*sin(eta[1])
    return(theta)
  }
}
