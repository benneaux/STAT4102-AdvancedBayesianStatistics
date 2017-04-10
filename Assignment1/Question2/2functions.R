library(tidyverse)

CauchyPercentage <- function(y = NULL,# a value to test Pr[p < y]
                             x, # a vector of samples
                             p, # a vector of possible parameters
                             climlow = -Inf,
                             climhigh = Inf,
                             alpha = 0.025) {
  
  cauchydist <- function(x,p) {
    H = vector(mode = "numeric",length = length(x))
    for(i in 1:length(x)){
      H[i] = (1+(x[i] - p)^2)^(-1)
    }
    return(prod(H))
  }
  dens <- function(p) {
    vapply(p,cauchydist,min(p), x = x)
  }
  cauchydens <- function(x,p, ...){
    c = (integrate(dens, climlow, climhigh)$value)^(-1)
    data.frame(vparams = p, postdens = c*dens(p))
  }
  
  df = cauchydens(x,p)
  
  yind = ifelse(length(which(df$vparams == y))==1,
                which(df$vparams == y),
                max(which(df$vparams < y)))

  c = (integrate(dens, climlow, climhigh)$value)^(-1)
  cumdist = c*integrate(dens,climlow,y)$value
  difference = abs(alpha - cumdist)
  results = c(yind,df[yind,"postdens"],cumdist)
  invisible(results)
  results <<- c(yind,df[yind,"postdens"],cumdist)
  return(difference)
}
CauchyPercentageNO <- function(y = NULL,# a value to test Pr[p < y]
                             x, # a vector of samples
                             p, # a vector of possible parameters
                             climlow = -Inf,
                             climhigh = Inf) {
  
  cauchydist <- function(x,p) {
    H = vector(mode = "numeric",length = length(x))
    for(i in 1:length(x)){
      H[i] = (1+(x[i] - p)^2)^(-1)
    }
    return(prod(H))
  }
  dens <- function(p) {
    vapply(p,cauchydist,min(p), x = x)
  }
  cauchydens <- function(x,p, ...){
    c = (integrate(dens, climlow, climhigh)$value)^(-1)
    data.frame(vparams = p, postdens = c*dens(p))
  }
  
  df = cauchydens(x,p)
  
  yind = ifelse(length(which(df$vparams == y))==1,
                which(df$vparams == y),
                max(which(df$vparams < y)))
  
  c = (integrate(dens, climlow, climhigh)$value)^(-1)
  c*integrate(dens,climlow,y)$value
}

CauchyHPD <- function(x, # vector of samples 
                      p, # vector of possible parameter values
                      alpha = 0.95, # HPD interval value
                      tol = 0.0001, # level of tolerance for exact HPD interval
                      plot = TRUE) { 
  
  cauchydens <- function(x,p){ # returns a df containing dens and param values.
    # Cauchy density described in notes.
    cauchydist <- function(x,p) {
      H = vector(mode = "numeric",length = length(x))
      for(i in 1:length(x)) {
        H[i] = (1+(x[i] - p)^2)^(-1) }
      return(prod(H)) }
    
    dens = function(p) vapply(p,cauchydist,min(p), x = x)
    c = (integrate(dens, -Inf, Inf)$value)^(-1)
    data.frame(vparams = p, postdens = c*dens(p)) }
  
  df = cauchydens(x,p) # Call the above function.
  
  cumdist = cumsum(df$postdens)*diff(df$vparams)[1] # Cumulative distribution
  post_median = which.min(abs(cumdist-0.5)) # The posterior median.
  
  # find lower and upper values for which prob dens is closest to target value
  HPDlimits <- function(post_dens) { 
    lower = which.min(abs(df$postdens[1:post_median]-post_dens))
    upper = which.min(abs(df$postdens[(post_median+1):length(df$postdens)]-post_dens))+post_median
    limits = c(lower,upper) }
  
  HPDlimitarea <- function(post_dens) { # find the area corresponding to that interval
    limitints = HPDlimits(post_dens)
    limitarea = sum(df$postdens[limitints[1]:limitints[2]])*diff(df$vparams)[1]
  }
  
  # Now calculate the HPD using the above functins.
  v2 = seq(0,max(df$postdens),by=tol)
  vals = sapply(v2,HPDlimitarea)
  w = which.min(abs(vals-alpha))
  r = c(df$vparams[HPDlimits(v2[w])])
  names(r) = c("Lower (HPD)","Upper (HPD)")
  if(plot) {
    plot(df$vparams, df$postdens, type = 'l', 
         xlab = expression(theta), ylab = "Density")
    abline(h = df[HPDlimits(v2[w])[1],"postdens"])
    abline(v = r["Upper (HPD)"], col = 'blue')
    abline(v = r["Lower (HPD)"], col = 'blue') }
  
  return(r)
}

