library(tidyverse)

CauchyPercentage <- function(y = NULL,# a value to test Pr[p < y]
                             alpha = 0.025,
                             x, # a vector of samples
                             p # a vector of possible parameters
) {
  
  dens <- function(p) {
    vapply(p,cauchydist,min(p), x = x)
  }
  
  cauchydens <- function(x,p,climlow = -Inf,climhigh = Inf){
    
    cauchydist <- function(x,p) {
      H = vector(mode = "numeric",length = length(x))
      for(i in 1:length(x)){
        H[i] = (1+(x[i] - p)^2)^(-1)
      }
      return(prod(H))
    }

    c = (integrate(dens, climlow, climhigh)$value)^(-1)
    data.frame(vparams = p, postdens = c*dens(p))
  }
  
  df = cauchydens(x,p)
  
  yind = ifelse(length(which(df$vparams == y))==1,
                which(df$vparams == y),
                max(which(df$vparams < y)))
  
  plot(df$vparams, df$postdens, type = 'l')
  abline(v = df$vparams[yind])
  abline(h = df[yind,"postdens"])
  
  cumdist = c*integrate(dens,climlow,y)$value
  difference = abs(alpha - cumdist)
  results = c(yind,df[yind,"postdens"],cumdist)
  invisible(results)
  return(difference)
}

CauchyHPD <- function(x, # vector of samples 
                      p, # vector of possible parameter values
                      alpha = 0.95, # HPD interval value
                      tol = 0.0001) { # level of tolerance for exact HPD interval
  
  cauchydens <- function(x,p){
    
    cauchydist <- function(x,p) {
      H = vector(mode = "numeric",length = length(x))
      for(i in 1:length(x)){
        H[i] = (1+(x[i] - p)^2)^(-1)
      }
      return(prod(H))
    }
    dens = function(p) vapply(p,cauchydist,min(p), x = x)
    c = (integrate(dens, -Inf, Inf)$value)^(-1)
    data.frame(vparams = p, postdens = c*dens(p))
  }
  
  df = cauchydens(x,p)
  
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
  # par(mfrow = c(1,2))
  # plot(df$vparams, cumdist, type = 'l')
  # abline(h = 0.5)
  # abline(v = df$vparams[post_median], col = 'red')
  # abline(v = r["upper"], col = 'blue')
  # abline(v = r["lower"], col = 'blue')
  # plot(df$vparams, df$postdens, type = 'l')
  # abline(v = df$vparams[post_median], col = 'red')
  # abline(h = df[HPDlimits(v2[w])[1],"postdens"])
  # abline(v = r["upper"], col = 'blue')
  # abline(v = r["lower"], col = 'blue')
  return(r)
}




HPD <- function(h, # height
                mode, # how to split the uniroot interval 
                dfunc, 
                pfunc,
                alpha = 0.95,
                plot = TRUE){
  
  lt       = uniroot(f=function(x){ dfunc(x) - h}, 
                     lower = 0, 
                     upper = mode)$root
  ut       = uniroot(f=function(x){ dfunc(x) - h}, 
                     lower = mode,
                     upper = 1)$root
  
  coverage = pfunc(ut) - pfunc(lt)
  
  hpdval   = abs(alpha - coverage)
  
  if (plot) {
    th     = seq(0, 1, length=10000)
    plot(th, dfunc(th),
         t = "l", 
         lty = 1,
         xlab = expression(theta),
         ylab = "posterior Density")
    
    abline(h = h)
    segments(ut, 0, ut, dfunc(ut))
    segments(lt, 0, lt, dfunc(lt))
    title(
      bquote(
        paste("p(",.(round(lt, 5))," < ", theta, " < ",
              .(round(ut,5)), " | " , y, ") = ",
              .(round(coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval)}

y=1; n=100; a=1; b=17; p=0.95
upper <- max(dbeta(seq(0, 1, length=10000),y + a,n-y+b))
interval <- seq(0,upper,by = 0.00001)
dfunc <- function(x) dbeta(x, y + a, n - y + b)
pfunc <-function(x) pbeta(x, y + a, n - y + b)
mode <- (y + a - 1)/(n + a + b - 2)

optimize(HPD,
         interval = interval,
         lower = min(interval),
         tol = .Machine$double.eps,
         mode = mode,
         dfunc = dfunc,
         pfunc = pfunc)
