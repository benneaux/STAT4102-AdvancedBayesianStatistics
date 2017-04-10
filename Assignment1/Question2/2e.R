# Berger 2
require(dplyr)
require(ggplot2)
require(ggthemes)
x <- c(-4.0, -5.5, -7.5, -4.5, -3.0)
p <- seq(0,5,1/1000)

CauchyHPD <- function(x, # vector of samples 
                      p, # vector of possible parameter values
                      alpha = 0.95, # HPD interval value
                      tol = 0.0001,
                      climlow = 0,
                      climhigh = Inf) { # level of tolerance for exact HPD interval
  
    
    cauchydist <- function(x,p) {
      H = vector(mode = "numeric",length = length(x))
      for(i in 1:length(x)){
        H[i] = (1+(x[i] - p)^2)^(-1)
      }
      return(prod(H))
    }
    
  dens = function(p) vapply(p,cauchydist,min(p), x = x)
  c = (integrate(dens, climlow , climhigh)$value)^(-1)
  df = data.frame(vparams = p, postdens = c*dens(p))

  cumdist = cumsum(df$postdens)*diff(df$vparams)[1]
  # cumdist = c*integrate(dens,climlow,y)$value
  difference = abs(alpha - cumdist)
  post_median = which.min(abs(cumdist-0.5))
  
  HPDlimits <- function(post_dens) { ## find lower and upper values for which
    ## prob dens is closest to target value
    lower = 1
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
  r = c(df$vparams[HPDlimits(v2[w])],c*integrate(dens,climlow,df$vparams[HPDlimits(v2[w])][2])$value)
  names(r) = c("lower","upper","dens")

  plot(df$vparams, df$postdens, type = 'l',xlab = expression(theta),
       ylab = "Density")
  abline(h = df[HPDlimits(v2[w])[2],"postdens"])
  abline(v = r["upper"], col = 'blue')
  title(bquote(paste("Cauchy Dens. with ", theta," > 0")))
  
  return(r)
}

CauchyHPD(x=x,p=p, climlow = 0, climhigh = Inf)
