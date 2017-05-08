# Generates 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(knitr) # spin
library(rmarkdown) # render
# filename <- "Assignment2/Assignment2.R"
# spin(filename, 
#      knit = FALSE, 
#      report = TRUE, 
#      format = "Rmd"
# )

render("Assignment2/Assignment2.Rmd", clean = FALSE, quiet = FALSE)

rm(list = ls())