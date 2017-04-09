# n <- 10
# a <- 0.05
# p <- seq(0.001,0.999,1/1000)
# z <- abs(qnorm(.5*a,0,1))
# # mean coverage the same as Table 1 "Approximate is better than..." for n = 30.
# # Agresti \& Coull 2
# #source("Assignment1/Question1/1functions.R")
# # Coverage Intervals
# 
# Q1ddata10 <- tbl_df(p) %>%
#   rename(p = value) %>%
#   mutate(Wald = vapply(p,waldcover,0)) %>%
#   mutate("Adjusted-Wald" = vapply(p,adjwaldcover, 0)) %>%
#   mutate(Exact = vapply(p,exactcover,0)) %>%
#   mutate(Score = vapply(p,scorecover,0)) %>%
#   mutate("Bayes-Laplace" = vapply(p, blcover,0)) %>%
#   mutate(Jeffreys = vapply(p,jeffreyscover,0)) %>%
#   gather(key = covinterval, value = dens, -p)
# 
# Q1dCoverage10 <- Q1ddata10 %>%
#   group_by(covinterval) %>%
#   summarize(mean(dens), median(dens), min(dens))

n <- 50

# Agresti \& Coull 2
#source("Assignment1/Question1/1functions.R")
# Coverage Intervals

Q1ddata50 <- tbl_df(p) %>%
  rename(p = value) %>%
  mutate(Wald = vapply(p,waldcover,0)) %>%
  mutate("Adjusted-Wald" = vapply(p,adjwaldcover, 0)) %>%
  mutate(Exact = vapply(p,exactcover,0)) %>%
  mutate(Score = vapply(p,scorecover,0)) %>%
  mutate("Bayes-Laplace" = vapply(p, blcover,0)) %>%
  mutate(Jeffreys = vapply(p,jeffreyscover,0)) %>%
  gather(key = covinterval, value = dens, -p)

Q1dCoverage50 <- Q1ddata50 %>%
  group_by(covinterval) %>%
  summarize(mean(dens), median(dens), min(dens))




