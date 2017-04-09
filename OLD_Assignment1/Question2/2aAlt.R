posterior=function(theta,data,scale,mu0,tau)
{
  f=function(theta) prod(dcauchy(data,theta,scale))
  likelihood=sapply(theta,f)
  prior=dcauchy(theta, mu0, tau)
  return(10^3*prior*likelihood)
}
p <- seq(0,8,by = 1/1000)
# data=c(11.4, 7.3, 9.8, 13.7, 10.6)
data = c(3,5)
curve(posterior(x,data,1,0,1),from=0,to=8,xlab="THETA",ylab="DENSITY")

data2 <- data.frame(param = p, density = posterior(p,data,1,0,1))
post_mean = max(posterior(p,data,1,0,1))
pmax <- data2[data2$density==post_mean,"param"]

solve.HPD.cauchy <-  function(h, p, data, plot=T, ...){
  
  lt = uniroot(
    f = function(x){
      posterior(x,data,1,0,1) - h
    },
    lower = -5000,
    upper = pmax)$root
  
  ut = uniroot(
    f = function(x){
      posterior(x,data,1,0,1) - h
    },
    lower = pmax,
    upper = 5000)$root
  
  coverage = posterior(ut,data, 1,0, 1) - posterior(lt,data, 1,0, 1)
  hpdval = abs(alpha-coverage) 
  if (plot) {
    th = seq(0,8, length=10000)
    plot(th, posterior(th,data,1,0,1),
         t="l", lty=1,xlab=expression(theta),
         ylab="posterior Density", ...)
    abline(h=h)
    segments(ut,0,ut,posterior(ut,data,1,0,1))
    segments(lt,0,lt,posterior(lt,data,1,0,1))
    title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
                       .(round(ut,5)), " | " , y, ") = ",
                       .(round(coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval)
}

interval <- seq(0,0.015,by = 0.00001)
alpha  <- 0.90
optimize(solve.HPD.cauchy,
         interval = interval,
         lower = min(interval),
         tol = .Machine$double.eps,
         data = data,
         p=p)
