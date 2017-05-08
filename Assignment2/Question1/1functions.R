solve.HPD.beta = function(shape1, shape2, credint = 0.95, one.sided = FALSE,...){
if(shape1 <= 1| one.sided == TRUE){
  lt = 0
  ut = qbeta(credint, shape1, shape2)
  coverage = credint
  results = data_frame("Lower"    = lt,
                       "Upper"    = ut,
                       "Coverage" = coverage,
                       "Height"   = ut)

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
  
  h = h$minimum
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

  return(c(results[[1]],results[[2]]))
}}

waldcover <- function(p,x = 0:n) {
  dens   = dbinom(x, n, p)
  phat  = x / n
  low   = phat - z * sqrt(phat * (1 - phat) / n)
  hig   = phat + z * sqrt(phat * (1 - phat) / n)
  intvals = as.numeric(low <= p & p <= hig)
  sum(intvals * dens)
}
adjwaldcover <- function(p,x = 0:n) {
  dens   = dbinom(x, n, p)
  nadj  = n + (z^2)
  phat  = (1/nadj)*(x + ((z^2)/2))
  low   = phat - z * sqrt(phat * (1 - phat)/nadj)
  hig   = phat + z * sqrt(phat * (1 - phat)/nadj)
  intvals = as.numeric(low <= p & p <= hig)
  sum(intvals * dens)
}
scorecover <- function(p,x = 0:n) {
  dens   = dbinom(x, n, p)
  phat  = x/n
  z2    = z*z
  low   = (phat + (z2/2)/n 
           - z * sqrt((phat * (1 - phat) + (z2/4)/n)/n))/(1 + z2/n)
  hig   = (phat + (z2/2)/n 
           + z * sqrt((phat * (1 - phat) + (z2/4)/n)/n))/(1 + z2/n)
  intvals = as.numeric(low <= p & p <= hig)
  sum(intvals * dens)
}
exactcover <- function(p,x = 0:n) {
  dens   = dbinom(x, n, p)
  low   = qbeta(a/2, x, n - x + 1)
  hig   = qbeta((1-a/2), x + 1, n - x)
  intvals = as.numeric(low <= p & p <= hig)
  sum(intvals * dens)
}
jeffreyscover <- function(p,x = 0:n,...) {
  dens  = dbinom(x,n,p)
    data <-matrix(data = NA, nrow = n+1, ncol = 2)
    for(i in 0:n){
      data[i+1,] = solve.HPD.beta(shape1 = i + 0.5, shape2 = n - i + 0.5,...)
    }
    intvals = as.numeric(data[,1] <= p & p <= data[,2])
    sum(intvals * dens)
  }


blcover <- function(p,x = 0:n,...) {
  dens = dbinom(x,n,p)
    data <-matrix(data = NA, nrow = n+1, ncol = 2)
    for(i in 0:n){
      data[i+1,] = solve.HPD.beta(shape1 = i + 1, shape2 = n - i + 1,...)
    }
    intvals = as.numeric(data[,1] <= p & p <= data[,2])
    sum(intvals * dens)
  }

blcover0 <- function(p,x = 0,...) {
  dens = dbinom(x,n,p)
    hig = solve.HPD.beta(shape1 = 1, shape2 = n + 1,...)
    intvals = as.numeric(p <= hig)
    sum(intvals * dens)
}

jeffreyscover0 <- function(p,x = 0,...) {
  dens = dbinom(x,n,p)
  hig = solve.HPD.beta(shape1 = 0.5, shape2 = n + 0.5,...)
  intvals = as.numeric(p <= hig)
  sum(intvals * dens)
}

