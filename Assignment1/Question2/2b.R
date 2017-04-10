x <- c(3,5)
p <- seq(0,8,1/100)


r <- CauchyHPD(x = x, p = p, alpha = 0.9, tol = 0.0001)

interval <- seq(0,3, by = 1/1000)
g <- vector("numeric",2L)
  g[1] <- optimize(CauchyPercentage,
                   interval = interval,
                   lower = min(interval),
                   tol = .Machine$double.eps,
                   x=x,
                   p=p,
                   alpha = 0.05)$minimum
  
  interval <- seq(5,8, by = 1/1000)
  
  g[2] <- optimize(CauchyPercentage,
                   interval = interval,
                   lower = min(interval),
                   tol = .Machine$double.eps,
                   x=x,
                   p=p,
                   alpha = 0.95)$minimum
  

CauchyPercentage(y = 19, alpha = 0.975, x,p)
