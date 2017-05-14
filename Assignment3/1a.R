FMDData <- readRDS("~/R/STAT4102-AdvancedBayesianStatistics/Assignment3/Question1/Data/FMDData.rds")

FMDData2 <- FMDData %>%
  mutate(StCattle = (Cattle-mean(Cattle))/sd(Cattle)) %>%
  select(-c(Province, Cattle))

fm_pois <- glm(FMD1998 ~ .,
               data = FMDData2,
               family = poisson)
fm_poisB <- glm(FMD1998 ~ EasternTurkey*StCattle,
               data = FMDData2,
              family = poisson)
fm_poisC <- glm(FMD1998 ~ StCattle,
                data = FMDData2,
                family = poisson)

(tidy(fm_poisB))
broom::tidy(fm_poisC)
broom::tidy(fm_pois, conf.int = TRUE)
broom::augment(fm_pois)
broom::glance(fm_pois)
broom::glance(fm_poisB)
summary(fm_poisB)
summary(fm_poisC)
lrtest(fm_poisB, fm_poisC)
coeftest(fm_pois, vcov = sandwich)


