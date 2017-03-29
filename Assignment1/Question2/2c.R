# Box \& Tiao

x <- c(11.4, 7.3, 9.8, 13.7, 10.6)
p <- seq(0,20,by = 1/100)

CauchyHPD(x=x,p=p,alpha = 0.95)


# Q2bData <- tbl_df(p) %>%
#   rename(Parameter = value) %>%
#   mutate(Density = c*dens(p))
# 
# Q2bChart <- ggplot(data = Q2bData, aes(x = Parameter, y = Density)) +
#   geom_area() +
#   theme_tufte(base_size = 14)
