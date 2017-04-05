
# a <- 0.05
#a <- 0.1
#a <- 0.25
a <- 0.5
ni <- seq(5,100, by = 1)
testdata <- matrix(data = NA, nrow = length(ni),ncol = 6)
testdata[,1] = ni

for(i in 1:length(ni)){
  n = ni[i]
  minvalues = vapply(p,adjwaldcover,0, x = 0)
  testdata[i,2] = which.min(minvalues)/1000
  minvalues = vapply(p,exactcover,0, x = 0)
  testdata[i,3] = which.min(minvalues)/1000
  minvalues = vapply(p,scorecover,0, x = 0)
  testdata[i,4] = which.min(minvalues)/1000
  minvalues = vapply(p,blcover,0, x = 0)
  testdata[i,5] = which.min(minvalues)/1000
  minvalues = vapply(p,jeffreyscover,0, x = 0)
  testdata[i,6] = which.min(minvalues)/1000
}
chartdata <- tbl_df(testdata) %>%
  rename(n = V1,
         AdjWald = V2,
         Exact = V3,
         Score = V4,
         BL = V5,
         Jeffreys = V6) %>%
  gather(key = covinterval, value, -n)

ggplot(data = chartdata, aes(x = n, y = value)) +
  geom_line(aes(group = covinterval, colour = covinterval), size  = 1) + 
  theme_tufte(base_size = 14) + 
  xlim(5,100)
