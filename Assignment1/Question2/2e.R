# Berger 2
require(dplyr)
require(ggplot2)
require(ggthemes)
x <- c(-4.0, -5.5, -7.5, -4.5, -3.0)
p <- seq(-10,0,1/100)

CauchyHPD(x,p,alpha = 0.95)
# Q2eData <- tbl_df(p) %>%
#   rename(Parameter = value) %>%
#   mutate(Density = c*dens(p))
# 
# Q2eChart <- ggplot(data = Q2eData, aes(x = Parameter, y = Density)) +
#   geom_area() +
#   theme_tufte(base_size = 14)