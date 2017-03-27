library(tidyverse)
x <- c(11.4, 7.3, 9.8, 13.7, 10.6)
p <- seq(0.000001,25,1/100000)

cauchydist <- function(p,x){
  H <- vector(mode = "numeric",length = length(x))
  for(i in 1:length(x)){
    H[i] = (1+(x[i] - p)^2)^(-1)
  }
  return(prod(H))
}
dens <- function(p) vapply(p,cauchydist,0,x=x)
c <- (integrate(dens, -Inf, Inf)$value)^(-1)

cauchyperc <- function(y){
  maxc = max(c*(integrate(dens,-Inf,Inf)$value))
  cauchyql = c*(dens(p))
  difference = y-(cauchyql)
  return(difference)
}

# cauchyperc2 <- function(p) {
#   data = matrix(NA,nrow = length(k))
#   for (k in 1:length(p)) {
#     data[k] = c*dens(p[k])
#   }
#   return(data)
# }

# k = seq(0.00001,0.6, by = 1/10000)
# cauchydata <- cauchyperc2(p)
# cauchydata <- cbind(p,cauchydata)

optalphau <- optimize(cauchyperc,
                     c(0,1),
                     tol = .Machine$double.eps^0.5,
                     maximum = TRUE)$maximum

# 
# 
# data <- tbl_df(p) %>%
#   mutate(density = c*dens(p))
# ggplot(data = data, aes(x = p, y = density)) +
#   geom_line()

