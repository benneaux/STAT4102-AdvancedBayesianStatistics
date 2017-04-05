
solve.HPD.beta <-  function(h, y, n, a, b, p, plot=T, ...){
  post_a = y + a
  post_b = n - y + b
  if (post_a < 1 | post_b < 1)
    warning("code assumes mode is not at 0 or 1")

  post_mode = (y + a - 1)/(n + a + b - 2)

  lt = uniroot(
    f = function(x){
        dbeta(x,post_a, post_b) - h
        },
    lower=0,
    upper=post_mode)$root

  ut = uniroot(
    f = function(x){
      dbeta(x,post_a, post_b) - h
      },
    lower = post_mode,
    upper = 1)$root

  coverage = pbeta(ut, post_a, post_b) - pbeta(lt, post_a, post_b)
  hpdval = abs(p-coverage)
  if (plot) {
    th = seq(0, 1, length=5000)
    plot(th, dbeta(th, post_a, post_b),
         t="l", lty=1,xlab=expression(theta),
         ylab="posterior Density", ...)
    abline(h=h)
    segments(ut,0,ut,dbeta(ut,post_a,post_b))
    segments(lt,0,lt,dbeta(lt,post_a,post_b))
    title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
                       .(round(ut,5)), " | " , y, ") = ",
                       .(round(coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval)
}

# Coverage Intervals

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
jeffreyscover <- function(p,x = 0:n) {
  fpx  = dbinom(x,n,p)
  data <-matrix(data = NA, nrow = n+1, ncol = 2)
  for(i in 0:n){
    data[i+1,] = HDIofICDF(qbeta,shape1 = i + 0.5, shape2 = n - i + 0.5)
  }
  inies = as.numeric(data[,1] <= p & p <= data[,2])
  sum(inies * fpx)
}
blcover <- function(p,x = 0:n) {
  fpx = dbinom(x,n,p)
  data <-matrix(data = NA, nrow = n+1, ncol = 2)
  for(i in 0:n){
    data[i+1,] = HDIofICDF(qbeta,shape1 = i + 1, shape2 = n - i + 1)
  }
  inies = as.numeric(data[,1] <= p & p <= data[,2])
  sum(inies * fpx)
}
  
  #sum(inies * fpx)
  # x = 1:(n-1)    
  #     low0 = qbeta(0, 1, n + 1)
  #     hig0 = qbeta(1-a, 1, n + 1)
  #     inies0 = as.numeric(low0 <= p & p <= hig0)
  #     
  #     lown = qbeta(a, n + 1, 1)
  #     hign = qbeta(1, n + 1, 1)
  #     iniesn = as.numeric(lown <= p & p <= hign)
  #     
  #     low = qbeta(a/2, x + 1, n - x + 1)
  #     hig = qbeta(1-a/2, x + 1, n - x + 1)
  #     iniesrest = as.numeric(low <= p & p <= hig)
  #     
  #   inies = c(inies0,iniesrest,iniesn)
      #sum(inies * fpx)
#}
