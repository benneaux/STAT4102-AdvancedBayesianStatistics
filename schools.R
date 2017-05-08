
schoolmodel <- function(){
  for(j in 1:J){
    y[j]~ dnorm(theta[j],tau.y[j])
    theta[j] <- mu.theta + xi*eta[j]
    tau.y[j] <- pow(sigma.y[j],-2)
  }
  xi ~ dnorm(0, tau.xi)
  tau.xi  <- pow(prior.scale,-2)
  for  (j in 1:J){
    eta[j] ~ dnorm(0,tau.eta)
  }
  tau.eta ~ dgamma(0.5,0.5)
  sigma.theta <- abs(xi)/sqrt(tau.eta)
  mu.theta ~ dnorm(0.0,0.0001)
}


write.model(schoolmodel, "schoolmodel.txt")
model.file2 = paste(getwd(),"schoolmodel.txt", sep="/")
J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
prior.scale <- 25
data <- list("J",  "y",  "sigma.y",  "prior.scale")


inits <- function(){list(eta=rnorm(J),mu.theta=rnorm(1),xi=rnorm(1),tau.eta=runif(1))}

WINE="/opt/local/bin/wine"
WINEPATH="/opt/local/bin/winepath"
OpenBUGS.pgm=paste0("/Users/benjamin/Applications/wine/",
                    "drive_c/ProgramFiles/OpenBUGS/OpenBUGS323/OpenBUGS.exe")

parameters <- c("theta", "mu.theta", "sigma.theta")


schools.sim <- bugs(data, 
                    inits, 
                    parameters,  
                    model.file2, 
                    n.chains=3, 
                    n.iter=1000,
                    OpenBUGS.pgm=OpenBUGS.pgm, 
                    WINE=WINE, 
                    WINEPATH=WINEPATH,
                    useWINE=T,
                    codaPkg = T,
                    debug = F)

out.coda <- read.bugs(schools.sim)
library(coda)
library(lattice)
plot(out.coda)
