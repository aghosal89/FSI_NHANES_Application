

# This function computes polar coordinates from the cartesian coordinates in
# p-dimensional euclidean space. 
# Input  :  1) x   - is the cartesian coordinate, p-dimensional vector of reals.
# Outputs:  1) r   - the radius of the polar coordinates.
#           2) eta - the vector of polar coordinates with length (p-1).

# This performs the reverse operation of the 'polar2cart' function above. This 
# function is written for p=2,3,4,5; which can be generalized.

cart2polar <- function(x) {
  # read dimension of the cartesian coordinate:
  p =length(x)      
  # define vector for storing polar coordinates of length (p-1)
  eta =  vector(length = p-1)   
  # computing the radius 
  r <- sqrt(sum(x^2))    
  
  # compute the polar coordinates for given p
  if(p==2) {
    eta[1]= atan(x[2]/x[1])
  }
  
  if(p==3) {
    eta[1]= atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[2]= atan(x[2]/x[1])
  }
  
  if(p==4) {
    eta[1] = atan(x[4]/sqrt(sum((c(x[1],x[2],x[3]))^2)))
    eta[2] = atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[3] = atan(x[2]/x[1])
  }
  
  if(p==5) {
    eta[1] = atan(x[5]/sqrt(sum((c(x[1],x[2],x[3],x[4]))^2)))
    eta[2] = atan(x[4]/sqrt(sum((c(x[1],x[2],x[3]))^2)))
    eta[3] = atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[4] = atan(x[2]/x[1])
  }
  return(c(r,eta))
}

