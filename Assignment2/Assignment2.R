#'---
#'title: "Assignment 1"
#'date: "`r Sys.Date()`"
#'output:  pdf_document
#'---

#+ DocSetup, include = FALSE
require(knitr)
opts_chunk$set(echo = FALSE, warning = FALSE)
opts_chunk$set(fig.width = 8, fig.height = 6, fig.align = "center")
#'

#+ Q1Setup, include = FALSE
source("Question1/1setup.R")
#'

#' ## Question 1: Inference for the Poisson parameter $\lambda$:

#' (a) Given a general Gamma(a,b) prior, derive the posterior distribution of $\lambda$, given data $(x_1,\dots,x_n)$, followed by the posterior predictive distribution of $(z_1,\dots,z_m)$.  

#+ Question1a1, echo = TRUE

#'
#' (b) Assuming Gamma(a,b) priors for two Poisson parameters $\lambda_1$ and $\lambda_2$, derive the posterior for $\phi = \lambda_1 / \lambda_2$. (Hint: use nuisance parameter $\mu = \lambda_2$).

#' 
#+ Question1bData, echo = TRUE, cache = TRUE

#'
#' (c) Derive the Jeffreys prior for $(\phi, \mu)$, and the corresponding marginal posterior for $\phi$. 
#'
#+ Question1cData, echo = TRUE, cache = TRUE


#'
#' (d) Derive the reference prior for $(\phi, \mu)$, and the corresponding marginal posterior for $\phi$. 

#+ Question1d 

#'
#' (e) Based on (c), consider the posterior for $\phi$ based on uniform priors for $\lambda_1$ and $\lambda_2$, and comment on when inference based on this posterior could be quite different from that based on the reference posterior from (e). Give a data example, in terms of resulting intervals. (Hint: the above posteriors are related to known pdfs, and by transformation may be simplified even further.).
#'
#+ Question1e, cache = TRUE, echo = FALSE

#'
#+Question1e2, cache = FALSE, fig.height = 10, fig.width  = 8, dependson = "Question1e"

#'

#'\newpage
#'

#+ Q2Setup, include = FALSE
rm(list = ls())
source("Question2/2setup.R")
opts_chunk$set(warning = FALSE, echo = TRUE, cache = TRUE)
opts_chunk$set(fig.width = 8, fig.height = 6, fig.align = "center")
#'
#' ## Question 2: Inference for variance components: 
#' 
#' (a) Use e.g. SAS (PROC VARCOMP) to perform a classical analysis of the data in Table 5.1.4 of Box & Tiao (1973), based on finding point estimates only.
#'
#+ Question2a, eval = FALSE
#'
#' (b)  Use WinBUGS for a Bayesian analysis of (a), and find reasonable point and interval estimates for $\sigma_1^2$ and $\sigma_2^2$. Include graphs, including one of the joint posterior. [6 marks] 
#+ Question2b

#'
#' (c) Box & Tiao also studied a 3-component model (Table 5.3.1). 
#'     i. Derive central credible intervals, for the 3 individual components, based on Table 5.3.3.
#'     ii. Use WinBUGS to do the same, and include graphs.
#'     iii. While not going as far as Box & Tiao’s Figure 5.3.2, produce a graph of the joint posterior of $\sigma_2^2$ and $\sigma_3^2$, and one of $\sigma_2^2 / \sigma_3^2$

#+ Question2c

#'
#' (d) Box, Hunter & Hunter (1976, Chapter 17.3) studied a pigment paste example with three components, focusing on point estimates only. Use WinBUGS again to perform a Bayesian analysis. Include graphs.

#+ Question2d

#'

