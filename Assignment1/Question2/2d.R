# Berger
require(dplyr)
require(ggplot2)
require(ggthemes)
x <- c(4.0, 5.5, 7.5, 4.5, 3.0)
p <- seq(-20,20,1/1000)

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

ncredint(data$param, data$density, level = 0.95, tol = 0.00001, verbose = TRUE)

post_mean = c*max(dens(p))
pmax <- data[data$density==post_mean,"param"]
solve.HPD.cauchy <-  function(h, data, plot=T, ...){
  
  lt = uniroot(
    f = function(x){
      posterior(x) - h
    },
    lower = -10000,
    upper = 11)$root
  
  ut = uniroot(
    f = function(x){
      posterior(x) - h
    },
    lower = 11,
    upper = 10000)$root
  
  coverage = pcauchy(ut,1,1) - pcauchy(lt,1,1)
  hpdval = abs(alpha-coverage) 
  if (plot) {
    th = seq(0,15, length=10000)
    plot(th, c*dens(th),
         t="l", lty=1,xlab=expression(theta),
         ylab="posterior Density", ...)
    abline(h=h)
    segments(ut,0,ut,posterior(ut))
    segments(lt,0,lt,posterior(lt))
    title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
                       .(round(ut,5)), " | " , y, ") = ",
                       .(round(coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval)
}

interval <- seq(0,.25,by = 0.00001)
alpha  <- 0.95
optimize(solve.HPD.cauchy,
         interval = interval,
         lower = min(interval),
         tol = .Machine$double.eps,
         y=x,
         p=p)

Q2dData <- tbl_df(p) %>%
              rename(Parameter = value) %>%
              mutate(Density = c*dens(p))

Q2dChart <- ggplot(data = Q2dData, aes(x = Parameter, y = Density)) +
              geom_area() +
              theme_tufte(base_size = 14)