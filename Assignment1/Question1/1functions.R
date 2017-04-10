# HDIofqbeta = function(shape1, 
#                       shape2, 
#                       credint = 0.95,
#                       tol = 1e-8,
#                       plot = FALSE,
#                       one.sided = FALSE) {
#   
#   if(one.sided){  # returns the one-sided interval if required.
#     return(qbeta(credint,
#                  shape1,
#                  shape2))
#   } else {
#     
#     # a function to calculate the interval for specific values of the upper and
#     # lower limits
#     
#     interval = function(lowint,  # for 95%: 0 <= lowint <=0.05
#                         credint, # e.g. 95%
#                         shape1,  # i.e: y - a
#                         shape2   # i.e: n - y + b
#     ){
#       
#       hig = qbeta(credint + lowint, # the upper limit
#                   shape1,
#                   shape2) 
#       
#       low = qbeta(lowint, # the lower limit.
#                   shape1, 
#                   shape2)
#       
#       hig - low # returns the interval.
#     }
#     
#     alpha    = 1 - credint # an upper limit for the interval optimiser.
#     length  = optimize(interval, 
#                        c(0, alpha),
#                        shape1 = shape1,
#                        shape2 = shape2,
#                        credint = 0.95,
#                        tol = 1e-8)
#     
#     HDIinterval = length$minimum # the actual optimised value from 0 to 
#     # the lower limit.
#     
#     lt = qbeta(HDIinterval, shape1, shape2)
#     ut = qbeta(HDIinterval + credint, shape1, shape2)
#     coverage = credint
#     l = HDIinterval
#     
#     
#     if (plot) {
#       th = seq(0, 1, length=10000)
#       plot(th, dbeta(th, shape1, shape2),
#            t="l",xlab=expression(theta),
#            ylab="Posterior Dens.")
#       abline(h = dbeta(results[[1]],shape1,shape2))
#       segments(results[[2]],0,results[[2]],dbeta(results[[2]],shape1,shape2))
#       segments(results[[1]],0,results[[1]],dbeta(results[[1]],shape1,shape2))
#       title(bquote(paste("p(", .(round(results[[1]], 5))," < ", theta, " < ",
#                          .(round(results[[2]],5)), " | " , y, ") = ",
#                          .(round(results[[3]], 5)))))
#     }
#     
#     return(c(qbeta(HDIinterval, shape1, shape2),
#              qbeta(HDIinterval + credint, shape1, shape2)))
#   }
# }
# 
solve.HPD.beta = function(shape1, shape2, credint = 0.95, plot=FALSE, ...){
if(shape1 <= 1){
  lt = 0
  ut = qbeta(credint, shape1, shape2)
  coverage = credint
  results = data_frame("Lower"      = lt,
                       "Upper"    = ut,
                       "Coverage" = coverage,
                       "Height"   = ut)
  if (plot) {
    th = seq(0, 1, length=10000)
    plot(th, dbeta(th, shape1, shape2), t="l",
         xlab=expression(theta), ylab="Posterior Dens.") #Plot the curve
    abline(h=results[[4]]) # plot the optimised 'height'
    segments(results[[2]],0,results[[2]],
             dbeta(results[[2]],shape1,shape2)) # the upper limit
    segments(results[[1]],0,results[[1]],
             dbeta(results[[1]],shape1,shape2)) # the lower limit
    title(bquote(paste("p(", .(round(results[[1]], 9))," < ", theta, " < ",
                       .(round(results[[2]],9)), " | " , y, ") = ",
                       .(round(results[[3]], 9))))) #title
  }
  return(c(results[[1]],results[[2]]))
}
if(shape1 > n){
  lt = qbeta(1-credint, shape1, shape2)
  ut = 1
  coverage = credint
  results = data_frame("Lower"    = lt,
                       "Upper"    = ut,
                       "Coverage" = coverage,
                       "Height"   = lt)
  if (plot) {
    th = seq(0, 1, length=10000)
    plot(th, dbeta(th, shape1, shape2), t="l",
         xlab=expression(theta), ylab="Posterior Dens.") #Plot the curve
    abline(h=results[[4]]) # plot the optimised 'height'
    segments(results[[2]],0,results[[2]],
             dbeta(results[[2]],shape1,shape2)) # the upper limit
    segments(results[[1]],0,results[[1]],
             dbeta(results[[1]],shape1,shape2)) # the lower limit
    title(bquote(paste("p(", .(round(results[[1]], 9))," < ", theta, " < ",
                       .(round(results[[2]],9)), " | " , y, ") = ",
                       .(round(results[[3]], 9))))) #title
  }
  return(c(results[[1]],results[[2]]))
} else {
  hpdfunc <- function(h, shape1, shape2){
    mode = (shape1 - 1)/(shape1 + shape2 - 2)
    lt = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
                 lower=0, upper=mode)$root
    ut = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
                 lower=mode, upper=1)$root
    coverage = pbeta(ut, shape1, shape2) - pbeta(lt, shape1, shape2)
    
    hpdval = abs(credint-coverage)
    
    return(hpdval)
  }
  upper = max(dbeta(seq(0,1, by = 0.001), shape1, shape2)) 
  
  h = optimize(hpdfunc,
               interval = seq(0,upper,by = 0.001),
               lower = 0,
               tol = .Machine$double.eps,
               shape1,
               shape2)
  
  h <- h$minimum
  mode = (shape1 - 1)/(shape1 + shape2 - 2)
  lt = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
               lower=0, upper=mode)$root
  ut = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
               lower=mode, upper=1)$root
  coverage = pbeta(ut, shape1, shape2) - pbeta(lt, shape1, shape2)
  results = data_frame("Lower"    = lt,
                       "Upper"    = ut,
                       "Coverage" = coverage,
                       "Height"   = h)
  if (plot) {
    th = seq(0, 1, length=10000)
    plot(th, dbeta(th, shape1, shape2), t="l",
         xlab=expression(theta), ylab="Posterior Dens.") #Plot the curve
    abline(h=results[[4]]) # plot the optimised 'height'
    segments(results[[2]],0,results[[2]],
             dbeta(results[[2]],shape1,shape2)) # the upper limit
    segments(results[[1]],0,results[[1]],
             dbeta(results[[1]],shape1,shape2)) # the lower limit
    title(bquote(paste("p(", .(round(results[[1]], 9))," < ", theta, " < ",
                       .(round(results[[2]],9)), " | " , y, ") = ",
                       .(round(results[[3]], 9))))) #title
  }
  return(c(results[[1]],results[[2]]))
}}

waldcover <- function(p,x = 0:n) {
  fpx   = dbinom(x, n, p)
  phat  = x / n
  low   = phat - z * sqrt(phat * (1 - phat) / n)
  hig   = phat + z * sqrt(phat * (1 - phat) / n)
  inies = as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
adjwaldcover <- function(p,x = 0:n) {
  fpx   = dbinom(x, n, p)
  nadj  = n + (z^2)
  phat  = (1/nadj)*(x + ((z^2)/2))
  low   = phat - z * sqrt(phat * (1 - phat)/nadj)
  hig   = phat + z * sqrt(phat * (1 - phat)/nadj)
  inies = as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
scorecover <- function(p,x = 0:n) {
  fpx   = dbinom(x, n, p)
  phat  = x/n
  z2    = z*z
  low   = (phat + (z2/2)/n 
           - z * sqrt((phat * (1 - phat) + (z2/4)/n)/n))/(1 + z2/n)
  hig   = (phat + (z2/2)/n 
           + z * sqrt((phat * (1 - phat) + (z2/4)/n)/n))/(1 + z2/n)
  inies = as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
exactcover <- function(p,x = 0:n) {
  fpx   = dbinom(x, n, p)
  low   = qbeta(a/2, x, n - x + 1)
  hig   = qbeta((1-a/2), x + 1, n - x)
  inies = as.numeric(low <= p & p <= hig)
  sum(inies * fpx)
}
jeffreyscover <- function(p,x = 0:n,...) {
  fpx  = dbinom(x,n,p)
    data <-matrix(data = NA, nrow = n+1, ncol = 2)
    for(i in 0:n){
      data[i+1,] = solve.HPD.beta(shape1 = i + 0.5, shape2 = n - i + 0.5,...)
    }
    inies = as.numeric(data[,1] <= p & p <= data[,2])
    sum(inies * fpx)
  }


blcover <- function(p,x = 0:n,...) {
  fpx = dbinom(x,n,p)
    data <-matrix(data = NA, nrow = n+1, ncol = 2)
    for(i in 0:n){
      data[i+1,] = solve.HPD.beta(shape1 = i + 1, shape2 = n - i + 1,...)
    }
    inies = as.numeric(data[,1] <= p & p <= data[,2])
    sum(inies * fpx)
  }

blcover0 <- function(p,x = 0,...) {
  fpx = dbinom(x,n,p)
    hig = HDIofqbeta(shape1 = 1, shape2 = n + 1,one.sided = TRUE,...)
    inies = as.numeric(p <= hig)
    sum(inies * fpx)
}

