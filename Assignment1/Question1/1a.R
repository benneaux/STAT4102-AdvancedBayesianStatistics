# Beta-Binomial model
options(digits = 8)
solve.HPD.beta = function(shape1, shape2, credint = 0.95, plot=FALSE, ...){
  
  hpdfunc <- function(h, shape1, shape2){
  mode = (shape1 - 1)/(shape1 + shape2 - 2)
  lt = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
               lower=0, upper=mode)$root
  ut = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
               lower=mode, upper=1)$root
  coverage = pbeta(ut, shape1, shape2) - pbeta(lt, shape1, shape2)
  
  hpdval = abs(credint-coverage)
  
  results <<- (c(lt,ut,coverage,h))
  return(hpdval)
  }
  upper = max(dbeta(seq(0,1, by = 0.001), shape1, shape2))
  h = optimize(hpdfunc,
           interval = seq(0,upper,by = 0.001),
           lower = 0,
           tol = .Machine$double.eps,
           shape1,
           shape2) 
  
  if (plot) {
    th = seq(0, 1, length=10000)
    plot(th, dbeta(th, shape1, shape2),
        t="l",xlab=expression(theta),
        ylab="posterior Density")
    abline(h=h)
    segments(results[[2]],0,results[[2]],dbeta(results[[2]],shape1,shape2))
    segments(results[[1]],0,results[[1]],dbeta(results[[1]],shape1,shape2))
    title(bquote(paste("p(", .(round(results[[1]], 5))," < ", theta, " < ",
                       .(round(results[[2]],5)), " | " , y, ") = ",
                       .(round(results[[3]], 5)), ")")))
  }
  return(c(results[[1]], results[[2]]))}

y=1; n=100; a=1; b=1
# upper <- max(dbeta(seq(0, 1, length=10000),y + a,n-y+b))
p <- seq(0,1,by = 0.001)
solve.HPD.beta(shape1 = y + a, shape2 = n - y + b, plot = TRUE)
# 
vapply(p, blcover, shape1 = y + a, shape2 = n - y + b, plot = FALSE, 0)




