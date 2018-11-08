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
L1 <- 50
L2 <- 200
thetaseq <- sort(c(0,seq(-10,10,length.out=L1)))
theta.zero<-which(thetaseq==0)
omegaseq <- seq(0.001,2.5,length.out=L2)
z <- matrix(NA,L1,L2)
p <- matrix(NA,L1,L2*4)
dim(p)<-c(L1,L2,4)
zlist<-vector("list",3)
plist<-vector("list",3)
plotlist <- vector("list",3)
plotlist2d <- vector("list",3)
scenelist <- lapply(1:3,function(i) paste0("scene",i))

tail.interest<-0.05

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
    zlist[[i]]<-z
    p[r,c,] <- c(mu=mu,sigma=sigma,
                low.over=pnorm(qlnorm(tail.interest,theta,omega),mu,sigma)/tail.interest,
                high.under=(1-pnorm(qlnorm(1-tail.interest,theta,omega),mu,sigma))/tail.interest
    )
    plist[[i]]<- p
    
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
  
  plotlist2d[[i]] <- ggplot(data=data.frame(cbind(omegaseq,z=z[theta.zero,])),aes(omegaseq,z))  + 
                              geom_line() +
                              ggtitle(paste("KL divergence with",flip,"match"))+
                              xlab("Omega")+ylab("KL")+
    geom_hline(yintercept = 3.8/100,linetype=3)
                              
      

}
plotlist2d[[1]]
plotlist2d[[2]]
plotlist2d[[3]]

optima<-c("mode","mean","median")
matplot(omegaseq,cbind(mode=plist[[1]][theta.zero,,3],
                       mean=plist[[2]][theta.zero,,3],
                       median=plist[[3]][theta.zero,,3]
                       )*100-100,
        type="l",xlab='omega',ylab='Percentage error in lower tail',log="xy")
abline(h=25,lty=3)
legend("topleft",col=1:3,lty=1:3,legend=optima)

matplot(omegaseq,cbind(mode=plist[[1]][theta.zero,,4],
                       mean=plist[[2]][theta.zero,,4],
                       median=plist[[3]][theta.zero,,4]
                       ),
type="l",xlab='omega',ylab='Fold error in upper tail',log="xy")
legend("topleft",col=1:3,lty=1:3,legend=optima)

matplot(omegaseq,( cbind(mode=plist[[1]][theta.zero,,4],
                       mean=plist[[2]][theta.zero,,4],
                       median=plist[[3]][theta.zero,,4]
                       ) *100-100 ) %>% abs(),
type="l",xlab='omega',ylab='Fold error in upper tail',log="xy")

matplot(omegaseq,cbind(mode=plist[[1]][theta.zero,,3],
                       mean=plist[[2]][theta.zero,,3],
                       median=plist[[3]][theta.zero,,3]
)*100-100, col=1,
type="l",xlab='omega',ylab='Percentage error in lower and upper tail',log="xy")
abline(h=25,lty=3)
legend("topleft",col=1,lty=1:3,legend=optima)
matplot(omegaseq,( cbind(mode=plist[[1]][theta.zero,,4],
                         mean=plist[[2]][theta.zero,,4],
                         median=plist[[3]][theta.zero,,4]
) *100-100 ) %>% abs(), col=2,
type="l",add=T)

#--------------------------------------------------------------------------------------------------------------#