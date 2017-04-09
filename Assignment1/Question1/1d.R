

# Agresti \& Coull 2
# Coverage Intervals
n <- 50
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




