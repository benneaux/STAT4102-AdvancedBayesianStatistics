# Berger


x <- c(4.0, 5.5, 7.5, 4.5, 3.0)
p <- seq(0,10,1/100)

CauchyHPD(x,p,alpha = 0.95)
p <- seq(-10,10,1/100)
CauchyHPD(x,p,alpha = 0.95)
# 
# Q2dData <- tbl_df(p) %>%
#               rename(Parameter = value) %>%
#               mutate(Density = c*dens(p))
# 
# Q2dChart <- ggplot(data = Q2dData, aes(x = Parameter, y = Density)) +
#               geom_area() +
#               theme_tufte(base_size = 14)