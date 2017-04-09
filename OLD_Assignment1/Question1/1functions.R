# solve.HPD.beta2(shape1 = y + a, shape2 = n - a + b)

solve.HPD.beta2<-  function(shape1, shape2, credint = 0.95){
  post_a = shape1
  post_b = shape2
  post_mode = which.max(dbeta(seq(0,1,by = 1/1000000),post_a,post_b))/1000000
  upper <- max(dbeta(seq(0, 1, length=10000), post_a, post_b))
  limits <- function(h,x){
  lt = uniroot(f = function(x){dbeta(x,post_a, post_b) - h}, lower=0, upper=post_mode)$root
  ut = uniroot(f = function(x){dbeta(x,post_a, post_b) - h}, lower = post_mode, upper = 1)$root
  coverage = pbeta(ut, post_a, post_b) - pbeta(lt, post_a, post_b)
  abs(credint-coverage)
  }
  h2 = optimise(limits,seq(0,upper,by = 0.001),tol=1e-10)$minimum
  lt = uniroot(f = function(x){dbeta(x,post_a, post_b) - h2}, lower=0, upper = post_mode)$root
  ut = uniroot(f = function(x){dbeta(x,post_a, post_b) - h2}, lower = post_mode, upper = 1)$root
  return(c(lt,ut))
}

# Coverage Intervals

HDIofqbeta = function(credint=0.95,
                     tol=1e-8,
                     shape1, 
                     shape2, 
                     one.sided = FALSE) {
  
  func <- qbeta
  
  if(one.sided){
    return(func(credint, shape1,shape2))
  } else {
  
  incredint    = 1 - credint
  intervalWidth = function(lowint,
                           func,
                           credint,
                           shape1,
                           shape2) 
  {
    func(credint + lowint, shape1, shape2) - func(lowint, shape1, shape2)
  }
  optInfo        = optimize(intervalWidth, 
                            c(0, incredint) ,
                            func = func,
                            credint = credint, 
                            tol = tol,
                            shape1, shape2)
  HDIlowint    = optInfo$minimum
  
  return(c(func(HDIlowint, shape1, shape2),
           func(credint + HDIlowint, shape1, shape2) ) )
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
    data[i+1,] = HDIofqbeta(shape1 = i + 0.5, shape2 = n - i + 0.5,...)
  }
  inies = as.numeric(data[,1] <= p & p <= data[,2])
  sum(inies * fpx)
}
blcover <- function(p,x = 0:n,...) {
  fpx = dbinom(x,n,p)
  data <-matrix(data = NA, nrow = n+1, ncol = 2)
  for(i in 0:n){
    data[i+1,] = HDIofqbeta(shape1 = i + 1, shape2 = n - i + 1,...)
  }
  inies = as.numeric(data[,1] <= p & p <= data[,2])
  sum(inies * fpx)
}

HDIofICDF = function(ICDF,
                     credMass=0.95,
                     tol=1e-8,
                     ...) {
  
  incredMass    = 1.0 - credMass
  intervalWidth = function(lowTailPr,
                           ICDF,
                           credMass,
                           ...) 
  {
    ICDF(credMass + lowTailPr, ...) - ICDF(lowTailPr, ...)
  }
  optInfo        = optimize(intervalWidth, 
                            c(0, incredMass) ,
                            ICDF = ICDF,
                            credMass = credMass, 
                            tol = tol,
                            ...)
  HDIlowTailPr    = optInfo$minimum
  
  return(c(ICDF(HDIlowTailPr, ...),
           ICDF(credMass + HDIlowTailPr, ...) ) )
}