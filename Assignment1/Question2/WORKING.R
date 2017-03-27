library(tidyverse)

#credible interval for (1D) posterior distribution stored in a numeric
## vector.  assumed unimodal!  pvec (vector of parameter values),
## npost (vector of posterior densities), level, tolerance
ncredint <- function(pvec,npost,level=0.95,tol=0.01,verbose=FALSE) {
  dx = diff(pvec)[1]
  cumdist <- cumsum(npost)*dx
  midpt <- which.min(abs(cumdist-0.5))

  lims <- function(pdens) { ## find lower and upper values for which
                            ## prob dens is closest to target value
    lower <- which.min(abs(npost[1:midpt]-pdens))
    upper <- which.min(abs(npost[(midpt+1):length(npost)]-pdens))+midpt
    c(lower,upper)
  }
  
  limarea <- function(pdens) {
    intlim <- lims(pdens)
    d <- sum(npost[intlim[1]:intlim[2]])*dx
    ##    cat(pdens,intlim,d,"\n")
    d
  }
  ## find credible interval
  v2 <- seq(0,max(npost),by=tol)
  vals <- sapply(v2,limarea)
  w <- which.min(abs(vals-level))
  r = c(pvec[lims(v2[w])],v2[w],limarea(v2[w]))
  names(r) = c("lower","upper","p","area")
  plot(pvec, cumdist, type = 'l')
  abline(h = 0.5)
  abline(v = pvec[midpt])
  abline(v = r["upper"])
  abline(v = r["lower"])
  ## credible intervals; posterior density, area
  if (verbose) return(r) else return(r[1:2])
}



cauchydist <- function(p,x){
  H <- vector(mode = "numeric",length = length(x))
  for(i in 1:length(x)){
    H[i] = (1+(x[i] - p)^2)^(-1)
  }
  return(prod(H))
}

dens <- function(p) vapply(p,cauchydist,0,x=x)
c <- (integrate(dens, -Inf, Inf)$value)^(-1)

x <- c(3,5)
p <- seq(-10,20,1/1000)
d <- integrate(dens, -Inf, Inf)
c <- d$value^(-1)

data <- data.frame(param = p, density = c*dens(p))

ncredint(data$param, data$density, level = 0.95, tol = 0.0001, verbose = TRUE)
