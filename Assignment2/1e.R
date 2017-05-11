## Monte Carlo HPD Interval Approximation 
## # Input posterior parameters x1,x2 

hpd.mc <- function(x1,x2, a){ # Significance level 
  alpha <- 0.05 ### specify significance level here 
# Sample from transformed posterior, a beta density 
n <- 1000   
### sample size
mc.phi <- sort(1/(1-rbeta(n,x1+ a ,x2 + a))-1) 
### sample, transform and sort 
# Compute empirical length & highest lower rank of interval estimates with the probability content of 1-alpha 
hpd.length <- floor((1-alpha)*n) 
hpd.lower.high <- n-hpd.length 
# Form candidates of HPD interval 
can.lower <- mc.phi[1:hpd.lower.high] 
can.upper <- mc.phi[(1+hpd.length):n]
# Determine the order where the shortest interval occurs 
shortest.stack <- 1     
### initial value 
for (i in 2:hpd.lower.high){ 
  if ((can.upper-can.lower)[i] < (can.upper-can.lower)[shortest.stack])
    shortest.stack <- i 
  } 
hpd.rank <- shortest.stack 
# Output the approximated HPD interval 
hpd.appr <- c(can.lower[hpd.rank], can.upper[hpd.rank])
hpd.appr 
}

alpha <- c(0.5,1)
results <- matrix(data = NA, nrow = 200, ncol = 2)
for(i in 1:2){
  for(j in 1:200){
    results[j,i] <- hpd.mc(0,j, alpha[i])[1]
  }}


results[1,] <- hpd.mc(0,j,alpha[1])

ni <- seq(0,100, by = 1)
p <- seq(0.0001,0.9999,1/1000)
a <- 0.05

testdata <- matrix(data = NA, nrow = length(ni),ncol = 3)
testdata[,1] <- sort(ni)
for(i in 1:nrow(testdata)){
  n = testdata[i,1]
  minvalues = vapply(p,blcover0,0, x = 0, one.sided = TRUE, credint = 1-a)
  testdata[i,2] = (which.min(minvalues)-1)/1000
  minvalues = vapply(p,jeffreyscover0, 0, x = 0, one.sided = TRUE, credint = 1-a)
  testdata[i,3] = (which.min(minvalues)-1)/1000
}
blcover0 <- function(p,x = 0,...) {
  dens = dbinom(x,n,p)
  hig = solve.HPD.beta(shape1 = 1, shape2 = n + 1,...)
  intvals = as.numeric(p <= hig)
  sum(intvals * dens)
}

jeffreyscover0 <- function(p,x = 0,...) {
  dens = dbinom(x,n,p)
  hig = solve.HPD.beta(shape1 = 0.5, shape2 = n + 0.5,...)
  intvals = as.numeric(p <= hig)
  sum(intvals * dens)
}
library(ggplot2)
library(tidyverse)
df <- tbl_df(testdata)
ggplot2::ggplot(data = df) +
                  geom_line(aes(x = df$V1,y = df$V2), colour = "blue", group = "BL") +
                  geom_line(aes(x = df$V1,y = df$V3), colour = "red", group = "Jef/Rat")+
  coord_cartesian(ylim = c(0,0.2),xlim =  c(20,100))

solve.HPD.beta = function(shape1, shape2, credint = 0.95, one.sided = FALSE,...){
  if(shape1 <= 1| one.sided == TRUE){
    lt = 0
    ut = qbeta(credint, shape1, shape2)
    coverage = credint
    results = data_frame("Lower"    = lt,
                         "Upper"    = ut,
                         "Coverage" = coverage,
                         "Height"   = ut)
    
    return(c(results[[1]],results[[2]]))
  }
  if(shape1 > n){
    lt = qbeta(1-credint, shape1, shape2)
    ut = 1
    coverage = credint
    results = data_frame("Lower"    = lt,
                         "Upper"    = ut,
                         "Coverage" = coverage,
                         "Height"   = lt)
    
    return(c(results[[1]],results[[2]]))
  } else {
    hpdfunc <- function(h, shape1, shape2){
      mode = (shape1 - 1)/(shape1 + shape2 - 2)
      lt = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
                   lower=0, upper=mode)$root
      ut = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
                   lower=mode, upper=1)$root
      coverage = pbeta(ut, shape1, shape2) - pbeta(lt, shape1, shape2)
      
      hpdval = abs(credint-coverage)
      
      return(hpdval)
    }
    upper = max(dbeta(seq(0,1, by = 0.001), shape1, shape2)) 
    
    h = optimize(hpdfunc,
                 interval = seq(0,upper,by = 0.001),
                 lower = 0,
                 tol = .Machine$double.eps,
                 shape1,
                 shape2)
    
    h = h$minimum
    mode = (shape1 - 1)/(shape1 + shape2 - 2)
    lt = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
                 lower=0, upper=mode)$root
    ut = uniroot(f=function(x){ dbeta(x,shape1, shape2) - h},
                 lower=mode, upper=1)$root
    coverage = pbeta(ut, shape1, shape2) - pbeta(lt, shape1, shape2)
    results = data_frame("Lower"    = lt,
                         "Upper"    = ut,
                         "Coverage" = coverage,
                         "Height"   = h)
    
    return(c(results[[1]],results[[2]]))
  }}
