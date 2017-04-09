# Box \& Tiao

x <- c(11.4, 7.3, 9.8, 13.7, 10.6)
p <- seq(0,20,by = 1/100)

r <- CauchyHPD(x=x,p=p,alpha = 0.95,tol = 0.00001)

interval <- seq(0,10, by = 1/1000)
h <- vector("numeric",2L)
  h[1] <- optimize(CauchyPercentage,
                interval = interval,
                lower = min(interval),
                tol = .Machine$double.eps,
                x=x,
                p=p,
                alpha = 0.025)$minimum
  
  interval <- seq(10,20, by = 1/1000)
  
  h[2] <- optimize(CauchyPercentage,
                   interval = interval,
                   lower = min(interval),
                   tol = .Machine$double.eps,
                   x=x,
                   p=p,
                   alpha = 0.975)$minimum

#8.9287514,12.270382

# Q2bData <- tbl_df(p) %>%
#   rename(Parameter = value) %>%
#   mutate(Density = c*dens(p))
# 
# Q2bChart <- ggplot(data = Q2bData, aes(x = Parameter, y = Density)) +
#   geom_area() +
#   theme_tufte(base_size = 14)
