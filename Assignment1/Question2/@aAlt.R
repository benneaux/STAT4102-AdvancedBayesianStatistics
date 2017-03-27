library(tidyverse)
x <- c(11.4, 7.3, 9.8, 13.7, 10.6)
p <- seq(5,15,1/1000)

d <- integrate(dens, -Inf, Inf)
c <- d$value^(-1)

cauchydist <- function(p,x){
  H <- vector(mode = "numeric",length = length(x))
  for(i in 1:length(x)){
    H[i] = (1+(x[i] - p)^2)^(-1)
  }
  return(prod(H))
}

dens <- function(p) vapply(p,cauchydist,0,x=x)
c <- (integrate(dens, -Inf, Inf)$value)^(-1)

data <- data.frame(param = p, density = c*dens(p))
post_mean = c*max(dens(p))
pmax <- data[data$density==post_mean,"param"]
solve.HPD.cauchy <-  function(h, y, p, plot=T, ...){
  
  lt = uniroot(
    f = function(x){
      c*dens(x) - h
    },
    lower = -100,
    upper = pmax)$root
  
  ut = uniroot(
    f = function(x){
      c*dens(x) - h
    },
    lower = pmax,
    upper = 100)$root
  
  coverage = pcauchy(c*dens(ut), pmax, 1) - pcauchy(c*dens(lt), pmax, 1)
  hpdval = abs(alpha-coverage) 
  if (plot) {
    th = seq(5,15, length=10000)
    plot(th, c*dens(th),
         t="l", lty=1,xlab=expression(theta),
         ylab="posterior Density", ...)
    abline(h=h)
    segments(ut,0,ut,c*dens(ut))
    segments(lt,0,lt,c*dens(lt))
    title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
                       .(round(ut,5)), " | " , y, ") = ",
                       .(round(c*coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval)
}

interval <- seq(0,.25,by = 0.00001)
alpha  <- 0.90
optimize(solve.HPD.cauchy,
         interval = interval,
         lower = min(interval),
         tol = .Machine$double.eps^0.25,
         y=x,
         p=p)

Q2bData <- tbl_df(p) %>%
  rename(Parameter = value) %>%
  mutate(Density = c*dens(p))

Q2bChart <- ggplot(data = Q2bData, aes(x = Parameter, y = Density)) +
  geom_area() +
  theme_tufte(base_size = 14)
