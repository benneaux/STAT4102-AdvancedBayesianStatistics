# Berger


x <- c(4.0, 5.5, 7.5, 4.5, 3.0)
p <- seq(2,8,1/100)

r <- CauchyHPD(x=x,p=p,alpha = 0.95,tol = 0.00001)


interval <- seq(2,5, by = 1/1000)

h <- vector("numeric",2L)
h[1] <- optimize(CauchyPercentage,
                 interval = interval,
                 lower = min(interval),
                 tol = .Machine$double.eps,
                 x=x,
                 p=p,
                 alpha = 0.025)$minimum

interval <- seq(6,8, by = 1/1000)

h[2] <- optimize(CauchyPercentage,
                 interval = interval,
                 lower = min(interval),
                 tol = .Machine$double.eps,
                 x=x,
                 p=p,
                 alpha = 0.975)$minimum

abline(v = h[1], col = "red", lty = 2)
abline(v = h[2], col = "red", lty = 2)
# 
# Q2dData <- tbl_df(p) %>%
#               rename(Parameter = value) %>%
#               mutate(Density = c*dens(p))
# 
# Q2dChart <- ggplot(data = Q2dData, aes(x = Parameter, y = Density)) +
#               geom_area() +
#               theme_tufte(base_size = 14)