x <- c(3,5)
p <- seq(-10,20,1/100)


CauchyHPD(x = x, p = p, alpha = 0.9, tol = 0.0001)


# Q2bData <- tbl_df(p) %>%
#             rename(Parameter = value) %>%
#             mutate(Density = c*dens(p))
# 
# Q2bChart <- ggplot(data = Q2bData, aes(x = Parameter, y = Density)) +
#               geom_area() +
#               theme_tufte(base_size = 14)
