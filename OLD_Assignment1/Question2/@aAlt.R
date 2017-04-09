library(tidyverse)

CauchyPercentage <- function(y = NULL,# a value to test Pr[p < y]
                             climlow = -Inf,
                             climhigh = Inf)
  {
  cauchydist <- function(p) {
      H = vector(mode = "numeric",length = length(x))
      for(i in 1:length(x)){
        H[i] = (1+(x[i] - p)^2)^(-1)
      }
      return(prod(H))
    }
    dens <- function(p) {
      vapply(p,cauchydist,min(p))
    }
    c = (integrate(dens, climlow, climhigh)$value)^(-1)
    df = data.frame(vparams = p, postdens = c*dens(p))
    
    yind = ifelse(length(which(df$vparams == y))==1,
                  which(df$vparams == y),
                  max(which(df$vparams < y)))
    
    
    cumdist = c*integrate(dens,climlow,y)$value
  return(cumdist)
}
CauchyHPD <- function(p,
                      climlow = -Inf,
                      climhigh = Inf
                      ) { # vector of samples 
    
    cauchydist <- function(p) {
      H = vector(mode = "numeric",length = length(x))
      for(i in 1:length(x)){
        H[i] = (1+(x[i] - p)^2)^(-1)
      }
      return(prod(H))
    }
    dens = function(p) vapply(p,cauchydist,min(p))
    c = (integrate(dens, climlow, climhigh)$value)^(-1)
    df = data.frame(vparams = p, postdens = c*dens(p))
  
  return(df$postdens)
}

optimize(HPD,
         interval = interval,
         lower = min(interval),
         tol = .Machine$double.eps,
         alpha = 0.9,
         mode = mode,
         dfunc = CauchyHPD,
         pfunc = CauchyPercentage,
         th = p,
         one.sided = FALSE,
         plot = TRUE,
         climlow = -Inf,
         climhigh = Inf)

x <- c(3,5)
p <- seq(0,10,by = 1/100)
 x <- -1*c(4.0, 5.5, 7.5, 4.5, 3.0)
 p <- seq(0,3,1/100)
upper <- max(CauchyHPD(p = p, climlow = 0))
mode <- p[which.max(CauchyHPD(p = p,climlow = 0))]
interval <- seq(0,upper,by = 0.00001)

HPD <- function(h, # height
                mode, # how to split the uniroot interval 
                dfunc, 
                pfunc,
                alpha = 0.95,
                plot = TRUE,
                th, #plot range
                one.sided = FALSE,
                climlow = NULL,
                climhigh = NULL,
                ...){
  if(length(climlow)==1|length(climhigh)==1){
  
  if(one.sided){
    ut       = uniroot(f=function(v){ dfunc(v, climlow,climhigh) - h}, 
                       lower = climlow,
                       upper = 100)$root
    coverage = pfunc(ut, climlow,climhigh)
  } else {
  lt       = uniroot(f=function(v){ dfunc(v, climlow,climhigh) - h}, 
                     lower = -100, 
                     upper = mode)$root
  ut       = uniroot(f=function(v){ dfunc(v, climlow,climhigh) - h}, 
                     lower = mode,
                     upper = 100)$root
  coverage = pfunc(ut, climlow,climhigh) - pfunc(lt, climlow,climhigh)
  }

  hpdval   = abs(alpha - coverage)

  if (plot) {
    plot(th, dfunc(th, climlow,climhigh),
         t = "l",
         lty = 1,
         xlab = expression(theta),
         ylab = "posterior Density")

    abline(h = h)
    if(one.sided){
      segments(ut, 0, ut, dfunc(ut,climlow,climhigh))
      title(
        bquote(
          paste("p(", theta, " < ",
                .(round(ut,5)), " | " , y, ") = ",
                .(round(coverage, 5)), ")")))
    } else {
      segments(ut, 0, ut, dfunc(ut, climlow,climhigh))
      segments(lt, 0, lt, dfunc(lt, climlow,climhigh))
      title(
        bquote(
          paste("p(",.(round(lt, 5))," < ", theta, " < ",
                .(round(ut,5)), " | " , y, ") = ",
                .(round(coverage, 5)), ")")))
    }
  }
  } else {if(one.sided){
    ut       = uniroot(f=function(v){ dfunc(v) - h}, 
                       lower = -100,
                       upper = 100)$root
    coverage = pfunc(ut, climlow)
  } else {
    lt       = uniroot(f=function(v){ dfunc(v) - h}, 
                       lower = -100, 
                       upper = mode)$root
    ut       = uniroot(f=function(v){ dfunc(v) - h}, 
                       lower = mode,
                       upper = 100)$root
    coverage = pfunc(ut) - pfunc(lt)
  }
    
    hpdval   = abs(alpha - coverage)
    
    if (plot) {
      plot(th, dfunc(th),
           t = "l",
           lty = 1,
           xlab = expression(theta),
           ylab = "posterior Density")
      
      abline(h = h)
      if(one.sided){
        segments(ut, 0, ut, dfunc(ut))
        title(
          bquote(
            paste("p(", theta, " < ",
                  .(round(ut,5)), " | " , y, ") = ",
                  .(round(coverage, 5)), ")")))
      } else {
        segments(ut, 0, ut, dfunc(ut))
        segments(lt, 0, lt, dfunc(lt))
        title(
          bquote(
            paste("p(",.(round(lt, 5))," < ", theta, " < ",
                  .(round(ut,5)), " | " , y, ") = ",
                  .(round(coverage, 5)), ")")))
      }
    }
  }
  if(one.sided){
    results <<- data.frame(ut,coverage,h)
  } else {
    results <<- data.frame(lt,ut,coverage,h)
  }
  return(hpdval)
  }



y=47; n=100; a=1; b=17; p=0.95
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
         pfunc = pfunc,
         th = seq(0,1,by = 1/1000),
         climlow = NULL,
         climhigh = NULL)
