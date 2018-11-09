#--------------------------------------------------------------------------------------------------------------#
# lognormal project
# Jeroen Elassaiss-Schaap (a, b), Kevin Duisters (c)

# a. Leiden Academic Centre for Drug Research, Leiden University
# b. PD-Value
# c. Mathematical Institute, Leiden University
#--------------------------------------------------------------------------------------------------------------#

# source functions (set working directory to source file location)
#source("functions/KL.R")
# test
#KL.proxy(theta=0.5,omega=0.2,mu=1,sigma=0.3)
#KL.exact(theta=0.5,omega=0.2,mu=1,sigma=0.3)
#--------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------#
# Comparison of three statistics in two interpretations
#--------------------------------------------------------------------------------------------------------------#

# settings
theta <- log(8) # Fix any theta
omegaseq <- seq(0.05,2,by=0.05)
L <- length(omegaseq)
A <- D <- matrix(NA,L,2,dimnames=list(1:L,c("0.10","0.90")))

for(l in 1:L){
  omega <- omegaseq[l]
  
  # Method A (true lognormal interpretation)
  A[l,] <- c(qlnorm(0.10,theta,omega),qlnorm(0.90,theta,omega)) 
  
  # Method D (see paper Table with candidate model interpretations)
  mu <- exp(theta + 0.5*omega^2) # mean
  sigma <- sqrt(exp(2*theta)*(exp(2*omega^2) - exp(omega^2)  ))  # mean match
  D[l,] <- c(max(0,qnorm(0.10,mu,sigma)),qnorm(0.90,mu,sigma))  # cut off at zero (support lognormal)
  
}
par(mfrow=c(1,2))
plot(omegaseq,(D[,1]-A[,1])/A[,1],ylab="rel diff",xlab="omega",main="10",type="l",ylim=c(-1,1))
plot(omegaseq,(D[,2]-A[,2])/A[,2],ylab="rel diff",xlab="omega",main="90",type="l",ylim=c(-1,1))

#--------------------------------------------------------------------------------------------------------------#
# Introducing the statistics
MMR <- function(omega){exp((1/2)*omega^2)}
MDR <- function(omega){(1/omega)*exp((1/2)*(omega^2-1))}
skew <- function(omega){(exp(omega^2)+2)*sqrt(exp(omega^2-1)) }

# first some intuition; this is a reproduction of Jeroen's plot in email Nov 1
par(mfrow=c(1,1))
plot(omegaseq,omegaseq,ylim=c(0,3),type="n")
lines(omegaseq,MDR(omegaseq))
lines(omegaseq,MMR(omegaseq),col="blue",lty=2)
lines(omegaseq,skew(omegaseq),col="magenta",lty=4)

# Now, the question is how well statistics can capture bad relative performance (D-A)/A
# plot along omega values
par(mfrow=c(1,2))
rel10<-(D[,1]-A[,1])/A[,1]
sub10 <- (rel10 < -0.05) & (rel10 > -1) # subinterval to prevent measures from visually 'exploding'
plot(x=rel10[sub10],y=MDR(omegaseq[sub10]),ylab="stat",xlab="rel diff",type="l",ylim=c(0,3),main="10",xlim=c(-1,0))
lines(rel10[sub10],MMR(omegaseq[sub10]),col="blue",lty=2)
lines(rel10[sub10],skew(omegaseq[sub10]),col="magenta",lty=4)

rel90<-(D[,2]-A[,2])/A[,2]
sub90 <- (rel90 > 0.05) & (rel90 < 1) # subinterval to prevent measures from visually 'exploding'
plot(x=rel90[sub90],y=MDR(omegaseq[sub90]),ylab="stat",xlab="rel diff",type="l",ylim=c(0,3),main="90",xlim=c(0,1))
lines(x=rel90[sub90],MMR(omegaseq[sub90]),col="blue",lty=2)
lines(x=rel90[sub90],skew(omegaseq[sub90]),col="magenta",lty=4)

#--------------------------------------------------------------------------------------------------------------#
# Rescaled statistics on subinterval
# this is the figure from my email Nov 8
normf <- function(vec){
  return((max(vec)-vec)/(max(vec)-min(vec)))
}

par(mfrow=c(1,2))
rel10<-(D[,1]-A[,1])/A[,1]
sub10 <- (rel10 < -0.05) & (rel10 > -1)
plot(x=rel10[sub10],y=normf(MDR(omegaseq[sub10])),ylab="Rescaled stat",xlab="rel diff",type="l",ylim=c(0,1),main="10",xlim=c(-1,0))
lines(rel10[sub10],normf(MMR(omegaseq[sub10])),col="blue",lty=2)
lines(rel10[sub10],normf(skew(omegaseq[sub10])),col="magenta",lty=4)

rel90<-(D[,2]-A[,2])/A[,2]
sub90 <- (rel90 > 0.05) & (rel90 < 1)
plot(x=rel90[sub90],y=normf(MDR(omegaseq[sub90])),ylab="Rescaled stat",xlab="rel diff",type="l",ylim=c(0,1),main="90",xlim=c(0,1))
lines(x=rel90[sub90],normf(MMR(omegaseq[sub90])),col="blue",lty=2)
lines(x=rel90[sub90],normf(skew(omegaseq[sub90])),col="magenta",lty=4)
