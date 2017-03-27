# Jaynes.
require(dplyr)
require(ggplot2)
require(ggthemes)
x <- c(3,5)
p <- seq(-5,15,1/100)
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

ncredint(data$param, data$density, level = 0.9, tol = 0.000001, verbose = TRUE)

post_mean = c*max(dens(p))
pmax <- data[data$density==post_mean,"param"]
solve.HPD.cauchy <-  function(h, y, p, plot=T, ...){
  
  lt = uniroot(
    f = function(x){
      c*dens(x) - h
    },
    lower = -50,
    upper = pmax)$root
  
  ut = uniroot(
    f = function(x){
      c*dens(x) - h
    },
    lower = pmax,
    upper = 50)$root
  
  coverage = c*pcauchy(ut, pmax, 1) - c*pcauchy(lt, pmax, 1)
  hpdval = abs(alpha-coverage) 
  if (plot) {
    th = seq(-5,15, length=10000)
    plot(th, c*dens(th),
         t="l", lty=1,xlab=expression(theta),
         ylab="posterior Density", ...)
    abline(h=h)
    segments(ut,0,ut,c*dens(ut))
    segments(lt,0,lt,c*dens(lt))
    title(bquote(paste("p(", .(round(lt, 5))," < ", theta, " < ",
                       .(round(ut,5)), " | " , y, ") = ",
                       .(round(coverage, 5)), ")")))
  }
  results <<- data.frame(lt,ut,coverage,h)
  return(hpdval)
}

interval <- seq(0,.25,by = 0.00001)
alpha  <- 0.90
optimize(solve.HPD.cauchy,
         interval = interval,
         lower = min(interval),
         tol = .Machine$double.eps,
         y=x,
         p=p)

Q2bData <- tbl_df(p) %>%
            rename(Parameter = value) %>%
            mutate(Density = c*dens(p))

Q2bChart <- ggplot(data = Q2bData, aes(x = Parameter, y = Density)) +
              geom_area() +
              theme_tufte(base_size = 14)
