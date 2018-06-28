#--------------------------------------------------------------------------------------------------------------#
# lognormal project
# Jeroen Elassaiss-Schaap (a, b), Kevin Duisters (c)

# a. Leiden Academic Centre for Drug Research, Leiden University
# b. PD-Value
# c. Mathematical Institute, Leiden University
#--------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------#
# KL functions

# theta: meanlog in lognormal distribution p
# omega: sdlog in lognormal distribution p
# mu: mean in normal distribution q
# sigma: sd in normal distribution q

KL.proxy <- function(theta,omega,mu,sigma,h=0.00001){
  x <- seq(0.001,10,h)
  px <- pmax(dlnorm(x,theta,omega),1e-300)
  qx <- pmax(dnorm(x,mu,sigma),1e-300)
  return(sum(px*log(px/qx))*h)
}

KL.exact <- function(theta,omega,mu,sigma){
  return(-1/2 - theta + log(sigma/omega) + (exp(2*theta+2*omega^2) - 2*mu*exp(theta+omega^2/2)+mu^2  )/(2*sigma^2))
}

#--------------------------------------------------------------------------------------------------------------#
