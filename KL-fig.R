#--------------------------------------------------------------------------------------------------------------#
# lognormal project
# Jeroen Elassaiss-Schaap (a, b), Kevin Duisters (c)

# a. Leiden Academic Centre for Drug Research, Leiden University
# b. PD-Value
# c. Mathematical Institute, Leiden University
#--------------------------------------------------------------------------------------------------------------#
# load libraries
library(plotly)

# source functions (set working directory to source file location)
source("functions/KL.R")

# test
KL.proxy(theta=0.5,omega=0.2,mu=1,sigma=0.3)
KL.exact(theta=0.5,omega=0.2,mu=1,sigma=0.3)


#--------------------------------------------------------------------------------------------------------------#
# visualize

#--------------------------------------------------------------------------------------------------------------#
# An example fitting a normal curve in a lognormal with mode match in densities
par(mfrow=c(1,2))

omega=0.67

for(theta in c(1,10)){
x <- seq(0,2*exp(theta),length.out=1e3)
plot(x,dlnorm(x,theta,omega),ylab="density",type="l",main=paste0("theta = ", theta, ", omega = ",omega))
mu <- exp(theta-omega^2) # mode match
sigma <- sqrt(exp(2*theta)*(exp(2*omega^2) - 2*exp(-0.5*omega^2) + exp(-2*omega^2) )) # KL min for mode match
lines(x,dnorm(x,mu,sigma),col="blue",lty=2,lwd=2)
legend("topright",bty="n",c("lognormal","normal"),col=c(1,"blue"),lty=c(1,2),lwd=c(1,2))
}

#--------------------------------------------------------------------------------------------------------------#

w <- seq(0,2,0.0001)
w[which(abs(exp(0.5*w^2-log(w)-0.5)-1.1)<1e-4)]

# 3D KL
L1 <- L2 <- 1000
thetaseq <- seq(-10,10,length.out=L1)
omegaseq <- seq(0.001,2.5,length.out=L2)
z <- matrix(NA,L1,L2)
plotlist <- vector("list",3)
scenelist <- lapply(1:3,function(i) paste0("scene",i))

par(mfrow=c(1,1))
for(i in 1:3){
  flip <- c("mode","mean","median")[i]
for(r in 1:L1){
  for(c in 1:L2){
    theta <- thetaseq[r]
    omega <- omegaseq[c]
    
    if(flip=="mode"){
      mu <- exp(theta - omega^2) # mode 
      sigma <- sqrt(exp(2*theta)*(exp(2*omega^2) - 2*exp(-0.5*omega^2) + exp(-2*omega^2)   ))  # mode match
    }
    if(flip=="mean"){
      mu <- exp(theta + 0.5*omega^2) # mean
      sigma <- sqrt(exp(2*theta)*(exp(2*omega^2) - exp(omega^2)  ))  # mean match
    }
    if(flip=="median"){
      mu <- exp(theta) # median 
      sigma <- sqrt(exp(2*theta)*(exp(2*omega^2) - 2*exp(omega^2/2) + 1   ))  # median match
    }
    
    z[r,c] <- min(KL.exact(theta,omega,mu,sigma),5)
  }
}

#contour(x=thetaseq,y=omegaseq,z) # alternative
  plotlist[[i]] <- plot_ly(x=~thetaseq,y=~omegaseq,z = ~t(z)) %>% add_surface() %>%layout(
    title = paste("KL divergence with",flip,"match"),
    scene = list(
      xaxis = list(title = "theta"),
      yaxis = list(title = "omega"),
      zaxis = list(title = "KL")
    ))
  
}
plotlist[[1]]
plotlist[[2]]
plotlist[[3]]

#--------------------------------------------------------------------------------------------------------------#