library(knitr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rmarkdown)
library(tidyr)
library(ggthemes)
library(grid)

q2Data <- readRDS("Question2/Data/q2Data_514.rds") %>%
  mutate(Batch = as.factor(Batch)) %>%
  mutate(Observation = as.factor(Observation)) %>%
  select(Batch,Observation, Samplevalue)

# obs <- c(1,2,3,4,5)
# B1 <- c(7.298,3.846,2.434,9.566,7.99)
# B2 <- c(5.22,6.556,0.608,11.788,-0.892)
# B3 <- c(0.11,10.386,13.434,5.51,8.166)
# B4 <- c(2.212,4.852,7.092,9.288,4.98)
# B5 <- c(0.282,9.014,4.458,9.446,7.198)
# B6 <- c(1.722,4.782,8.106,0.758,3.758)
# 
# df <- tbl_df(cbind(obs,B1,B2,B3,B4,B5,B6)) %>%
#   gather("Batch","Samplevalue", -obs) %>%
#   select(Observation = obs, everything())
# 
# saveRDS(df, "Assignment2/Question2/Data/q2Data_514.rds")


