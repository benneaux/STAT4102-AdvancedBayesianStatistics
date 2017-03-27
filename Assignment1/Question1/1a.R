# Beta-Binomial model

# solve.HPD.beta = function(h, y, n, a, b, p, plot=T, ...){
#   apost= y + a
#   bpost= n - y + b
#   if (apost < 1 | bpost < 1)
#     warning("code assumes mode is not at 0 or 1")
#   mode = (y + a - 1)/(n + a + b - 2)
#   lt = uniroot(f=function(x){ dbeta(x,apost, bpost) - h},
#                lower=0, upper=mode)$root
#   ut = uniroot(f=function(x){ dbeta(x,apost, bpost) - h},
#                lower=mode, upper=1)$root
#   coverage = pbeta(ut, apost, bpost) - pbeta(lt, apost, bpost)
#   hpdval = abs(p-coverage) 
#   if (plot) {
#     th = seq(0, 1, length=5000)
#     plot(th, dbeta(th, apost, bpost),
#          t="l", lty=1,xlab=expression(theta),
#          ylab="posterior Density", ...)
#     abline(h=h)
#     segments(ut,0,ut,dbeta(ut,apost,bpost))
#     segments(lt,0,lt,dbeta(lt,apost,bpost))
#     title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
#                        .(round(ut,5)), " | " , y, ") = ",
#                        .(round(coverage, 5)), ")")))
#   }
#   results <<- data.frame(lt,ut,coverage,h)
#   return(hpdval) }
# 
y=23; n=80; a=2; b=17; p=0.95
upper <- max(dbeta(seq(0, 1, length=500),y + a,n-y+b))
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

