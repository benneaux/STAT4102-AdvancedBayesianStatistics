# Generates 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(knitr) # spin
library(rmarkdown) # render
filename <- "Assignment1/Assignment1.R"
spin(filename, 
     knit = FALSE, 
     report = TRUE, 
     format = "Rmd"
)

render("Assignment1/Assignment1.Rmd", clean = FALSE, quiet = FALSE)

rm(list = ls())