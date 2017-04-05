x <- c(3,5)
p <- seq(-10,20,1/100)


r <- CauchyHPD(x = x, p = p, alpha = 0.9, tol = 0.0001)

interval <- seq(0,4, by = 1/1000)
g <- 0
  g[1] <- optimize(CauchyPercentage,
                   interval = interval,
                   lower = min(interval),
                   tol = .Machine$double.eps,
                   x=x,
                   p=p,
                   alpha = 0.025)$minimum
  
  interval <- seq(0,15, by = 1/1000)
  
  g[2] <- optimize(CauchyPercentage,
                   interval = interval,
                   lower = min(interval),
                   tol = .Machine$double.eps,
                   x=x,
                   p=p,
                   alpha = 0.975)$minimum
  
# Q2bData <- tbl_df(p) %>%
#             rename(Parameter = value) %>%
#             mutate(Density = c*dens(p))
# 
# Q2bChart <- ggplot(data = Q2bData, aes(x = Parameter, y = Density)) +
#               geom_area() +
#               theme_tufte(base_size = 14)
