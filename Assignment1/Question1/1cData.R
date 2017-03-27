# Agresti \& Coull 2
#source("Assignment1/Question1/1functions.R")
# Coverage Intervals

Q1cdata <- tbl_df(p) %>%
  rename(p = value) %>%
  mutate(Wald = vapply(p,waldcover,0)) %>%
  mutate("Adjusted-Wald" = vapply(p,adjwaldcover, 0)) %>%
  mutate(Exact = vapply(p,exactcover,0)) %>%
  mutate(Score = vapply(p,scorecover,0)) %>%
  mutate("Bayes-Laplace" = vapply(p, blcover,0)) %>%
  mutate(Jeffreys = vapply(p,jeffreyscover,0)) %>%
  gather(key = covinterval, value = dens, -p)