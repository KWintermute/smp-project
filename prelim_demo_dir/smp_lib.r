###
### Copy of the R functions used in the test code
###

firstpasFt <-function(s) {
  tmp1 <-p*ft(s)
  tmp2 <-solve(Id-tmp1)
  tmp3 <-solve(Id*tmp2)
  return(1/s*tmp1%*%tmp2%*%tmp3)
}

firstpasft <-function(s) {
  tmp1 <-p*ft(s)
  tmp2 <-solve(Id-tmp1)
  tmp3 <-solve(Id*tmp2)
  return (tmp1%*%tmp2%*%tmp3)
}

#Definition of function to compute the transient probability matrix in the LT domain transprobt <-function(s) {
transprobt <-function(s) {
  J <-matrix(1,ncol=n,nrow=n)
  tmp1 <-p*ft(s)
  tmp2 <-solve(Id-tmp1)
  return(1/s*tmp2%*%(Id-Id*(tmp1%*%J)))
}

#Definition of function to compute the first passage CDF matrix in the LT domain
SMP.firstpass.cdf <- function(s) {
  tmp1 <- p*ft(s)  
  tmp2 <- solve(Id-tmp1)
  tmp3 <- solve(Id*tmp2)
  return(1/s*tmp1%*%tmp2%*%tmp3)
}

#Definition of function to compute the first passage PDF matrix in the LT domain
SMP.firstpass.pdf <- function(s) {
  tmp1 <- p*ft(s)  
  tmp2 <- solve(Id-tmp1)
  tmp3 <- solve(Id*tmp2)
  return(tmp1%*%tmp2%*%tmp3)
}

#Definition of function to compute the transient probability matrix in the LT domain
SMP.transprob.t <- function(s) {
  J <- matrix(1,ncol=n,nrow=n)
  tmp1 <- p*ft(s)  
  tmp2 <- solve(Id-tmp1) 
  return(1/s*tmp2%*%(Id-Id*(tmp1%*%J))) 
}

#Definition of function to compute the expected visits to a state matrix in the LT domain
SMP.expected.visits.to.state <- function(s) {
  tmp1 <- p*ft(s)  
  tmp2 <- solve(Id-tmp1) 
  return(1/s*(tmp2-Id)) 
}


#Definition of function to compute the expected time in state matrix in the LT domain
SMP.expected.time.in.state <- function(s) {
  J <- matrix(1,ncol=n,nrow=n)
  tmp1 <- p*ft(s)  
  tmp2 <- solve(Id-tmp1) 
  return(1/s^2*tmp2%*%(Id-Id*(tmp1%*%J)))
}


#Definition of function to compute the probability of transistioning to state 0 times matrix in the LT domain
SMP.prob.trans.to.state.0.times <- function(s) {
  J <- matrix(1,ncol=n,nrow=n)
  tmp1 <- p*ft(s)
  tmp2 <- solve(Id-tmp1)
  tmp3 <- solve(Id*tmp2)
  g <- tmp1%*%tmp2%*%tmp3
  return(1/s*(J-g))
}


#Definition of function to compute the probability of transistioning to state 1 or less times matrix in the LT domain
SMP.prob.trans.to.state.1.times <- function(s) {
  J <- matrix(1,ncol=n,nrow=n)
  tmp1 <- p*ft(s)  
  tmp2 <- solve(Id-tmp1)
  tmp3 <- solve(Id*tmp2)
  g <- tmp1%*%tmp2%*%tmp3
  return(1/s*(J-g*(J%*%(Id*g)))) 
}

#Definition of function to compute the probability of transistioning to state 2 or less times matrix in the LT domain
SMP.prob.trans.to.state.2.times <- function(s) {
  J <- matrix(1,ncol=n,nrow=n)
  tmp1 <- p*ft(s)  
  tmp2 <- solve(Id-tmp1)
  tmp3 <- solve(Id*tmp2)
  g <- tmp1%*%tmp2%*%tmp3
  return(1/s*(J-g*(J%*%((Id*g)^2)))) 
}

###########################################################################################


p <- matrix(c(0,0,0.14,0,0.86,0,1,0),nrow=n,ncol=n,byrow=T)

#function that calculates the LT for for the f-tilde matrix
ft <- function(s) {
  #Checking if the LTs were just calculated for this value of s
  if (s == temp_s) {return(temp)}
  
  #Computing the numeric integration for the LT of the 5 Weibull waiting time distributions
  ft1r <- function(t) {dlnorm(t,5.8,.8384)*exp(-Re(s)*t)*cos(Im(s)*t)}
  ft1i <- function(t) {dlnorm(t,5.8,.8384)*exp(-Re(s)*t)*sin(Im(s)*t)}
  temp1 <- integrate(ft1r,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value - 
    integrate(ft1i,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value*1i
  ft2r <- function(t) {dfisk(t,34.231,2.097)*exp(-Re(s)*t)*cos(Im(s)*t)}
  ft2i <- function(t) {dfisk(t,34.231,2.097)*exp(-Re(s)*t)*sin(Im(s)*t)}
  temp2 <- integrate(ft2r,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value - 
    integrate(ft2i,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value*1i
  ft3r <- function(t) {dlnorm(t,5.128,1.766)*exp(-Re(s)*t)*cos(Im(s)*t)}
  ft3i <- function(t) {dlnorm(t,5.128,1.766)*exp(-Re(s)*t)*sin(Im(s)*t)}
  temp3 <- integrate(ft3r,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value - 
    integrate(ft3i,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value*1i
  ft4r <- function(t) {dunif(t,0,10)*exp(-Re(s)*t)*cos(Im(s)*t)}
  ft4i <- function(t) {dunif(t,0,10)*exp(-Re(s)*t)*sin(Im(s)*t)}
  #temp4 <- integrate(ft4r,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value - 
  #          integrate(ft4i,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value*1i
  #ft5r <- function(t) {dchisq(t,2)*exp(-Re(s)*t)*cos(Im(s)*t)}
  #ft5i <- function(t) {dchisq(t,2)*exp(-Re(s)*t)*sin(Im(s)*t)}
  #temp5 <- integrate(ft5r,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value - 
  #          integrate(ft5i,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value*1i
  #ft6r <- function(t) {dweibull(t,1.2,2)*exp(-Re(s)*t)*cos(Im(s)*t)}
  #ft6i <- function(t) {dweibull(t,1.2,2)*exp(-Re(s)*t)*sin(Im(s)*t)}
  #temp6 <- integrate(ft6r,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value - 
  #        integrate(ft6i,0,Inf,subdivisions=10000,rel.tol=1e-10,stop.on.error=FALSE)$value*1i
  
  #Defining a matrix to output from this function
  output <- matrix(
    c(0.000,0.0,
      temp3,0.000,temp2,
      0.000,temp1,0
    ),nrow=n,ncol=n,byrow=T)
  
  #Updating the "global" variable temp_s
  temp_s <<- s
  
  #Updating the "global" matrix variable temp
  temp <<- output
  
  return(output)
}