library(tidyverse)
library(ggthemes)
#source("Assignment1/Question1/1functions.r")


ni <- seq(5,100, by = 1)
nis <- rep(ni, 20)
p <- seq(0.0001,0.5221,1/1000)
a <- rep(c(0.05,0.1,0.25,0.5),5)
           
testdata <- matrix(data = NA, nrow = length(ni)*20,ncol = 7)
testdata[,1] <- sort(nis)
testdata[,7] <- sort(a)
for(i in 1:nrow(testdata)){
  n = testdata[i,1]
  a = testdata[i,7]
  z = abs(qnorm(.5*a,0,1))
  minvalues = vapply(p,adjwaldcover,0, x = 0)
  testdata[i,2] = (which.min(minvalues)-1)/1000
  minvalues = vapply(p,exactcover,0, x = 0)
  testdata[i,3] = (which.min(minvalues)-1)/1000
  minvalues = vapply(p,scorecover,0, x = 0)
  testdata[i,4] = (which.min(minvalues)-1)/1000
  minvalues = vapply(p,blcover0,0, x = 0, one.sided = TRUE, credint = 1-a)
  testdata[i,5] = (which.min(minvalues)-1)/1000
  minvalues = vapply(p,jeffreyscover0, 0, x = 0, one.sided = TRUE, credint = 1-a)
  testdata[i,6] = (which.min(minvalues)-1)/1000
}
  

chartdata <- tbl_df(testdata) %>%
  rename(n = V1,
         AdjWald = V2,
         Exact = V3,
         Score = V4,
         BL = V5,
         Jeffreys = V6,
         a = V7) %>%
  gather(key = covinterval, value, -c(n,a))



q1echart <- ggplot(data = chartdata, aes(x = n, y = value)) +
              geom_line(aes(group = covinterval, colour = covinterval), size  = 1) +
              theme_few(base_size = 14) +
              facet_wrap(~a, ncol = 2) +
              xlim(5,100)
