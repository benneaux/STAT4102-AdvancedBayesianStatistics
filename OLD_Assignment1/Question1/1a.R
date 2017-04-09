# Beta-Binomial model
options(digits = 8)
solve.HPD.beta = function(h, y, n, a, b, p, plot=T, ...){
  apost= y + a
  bpost= n - y + b
  mode = (y + a - 1)/(n + a + b - 2)
  lt = uniroot(f=function(x){ dbeta(x,apost, bpost) - h},
               lower=0, upper=mode)$root
  ut = uniroot(f=function(x){ dbeta(x,apost, bpost) - h},
               lower=mode, upper=1)$root
  coverage = pbeta(ut, apost, bpost) - pbeta(lt, apost, bpost)
  
  hpdval = abs(p-coverage)
  
  if (plot) {
    th = seq(0, 1, length=10000)
    plot(th, dbeta(th, apost, bpost),
         t="l", lty=1,xlab=expression(theta),
         ylab="posterior Density", ...)
    abline(h=h)
    segments(ut,0,ut,dbeta(ut,apost,bpost))
    segments(lt,0,lt,dbeta(lt,apost,bpost))
    title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
                       .(round(ut,5)), " | " , y, ") = ",
                       .(round(coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval) }

y=1; n=100; a=1; b=1; p=0.95
upper <- max(dbeta(seq(0, 1, length=10000),y + a,n-y+b))
interval <- seq(0,upper,by = 0.00001)

optimize(solve.HPD.beta,
         interval = interval,
         lower = min(interval),
         tol = .Machine$double.eps,
         y=y,
         n=n,
         a=a,
         b=b,
         p=p)

