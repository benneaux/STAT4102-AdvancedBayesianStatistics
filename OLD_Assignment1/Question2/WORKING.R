library(tidyverse)



CauchyHPD <- function(x, # vector of samples 
                     p, # vector of possible parameter values
                     alpha=0.95, # HPD interval value
                     tol=0.0001) # level of tolerance for exact HPD interval
  {
  
  cauchydist <- function(x,p){
    H <- vector(mode = "numeric",length = length(x))
    for(i in 1:length(x)){
      H[i] = (1+(x[i] - p)^2)^(-1)
    }
    return(prod(H))
  }
  
  dens <- function(p) vapply(p,cauchydist,0,x=x)
  
  c = (integrate(dens, -Inf, Inf)$value)^(-1)
  
  df = data.frame(vparams = p, postdens = c*dens(p))
  
  
  cumdist = cumsum(df$postdens)*diff(df$vparams)[1]
  post_median = which.min(abs(cumdist-0.5))

  HPDlimits <- function(post_dens) { ## find lower and upper values for which
                               ## prob dens is closest to target value
    lower = which.min(abs(df$postdens[1:post_median]-post_dens))
    upper = which.min(abs(df$postdens[(post_median+1):length(df$postdens)]-post_dens))+post_median
    limits = c(lower,upper)
  }
  
  HPDlimitarea <- function(post_dens) {
    limitints = HPDlimits(post_dens)
    limitarea = sum(df$postdens[limitints[1]:limitints[2]])*diff(df$vparams)[1]
  }
  ## find credible interval
  v2 = seq(0,max(df$postdens),by=tol)
  vals = sapply(v2,HPDlimitarea)
  w = which.min(abs(vals-alpha))
  r = c(df$vparams[HPDlimits(v2[w])])
  names(r) = c("lower","upper")
  
  plot(df$vparams, cumdist, type = 'l')
  abline(h = 0.5)
  abline(v = df$vparams[post_median])
  abline(v = r["upper"])
  abline(v = r["lower"])
  return(r)
}


