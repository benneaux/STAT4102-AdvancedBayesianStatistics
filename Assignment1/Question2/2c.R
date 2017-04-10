# Box \& Tiao

x <- c(11.4, 7.3, 9.8, 13.7, 10.6)
p <- seq(5,15,by = 1/100)

r <- CauchyHPD(x=x,p=p,alpha = 0.95,tol = 0.00001)

interval <- seq(5,10, by = 1/1000)
h <- vector("numeric",2L)
  h[1] <- optimize(CauchyPercentage,
                interval = interval,
                lower = min(interval),
                tol = .Machine$double.eps,
                x=x,
                p=p,
                alpha = 0.025)$minimum
  
  interval <- seq(10,15, by = 1/1000)
  
  h[2] <- optimize(CauchyPercentage,
                   interval = interval,
                   lower = min(interval),
                   tol = .Machine$double.eps,
                   x=x,
                   p=p,
                   alpha = 0.975)$minimum
  abline(v = h[1], col = "red", lty = 2)
  abline(v = h[2], col = "red", lty = 2)

