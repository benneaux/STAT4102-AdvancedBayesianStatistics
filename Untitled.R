library(lme4)
library(ggplot2)

## EDA for two datasets
dye.plot = qplot(
  Batch, Yield, data=Dyestuff2,
  geom="jitter",
  position=position_jitter(width=.1, height=0)
)


## Standard one-way ANOVA
dmodel.f = lm(Yield~Batch, Dyestuff2)
anova(dmodel.f)

## Random effects model
dmodel.r = lmer(Yield~1|Batch, data=Dyestuff2)
dmodel.r

## Shrinkage
plot(fitted(dmodel.f), fitted(dmodel.r))
abline(0,1)
mean(Dyestuff2$Yield)
points(1527.5, 1527.5, pch=19, col=2)
(fitted(dmodel.r)-mean(Dyestuff2$Yield))/(fitted(dmodel.f)-mean(Dyestuff2$Yield))
dmodel.r
1764/(1764+2451.3/5)

## Diagnostics for dmodel.r and dmodel.f
qqnorm(resid(dmodel.r))
qqnorm(resid(dmodel.f))
qqnorm(ranef(dmodel.r)$Batch[,1])
qplot(Dyestuff2$Batch, resid(dmodel.r))

## AIC for dmodel.r and dmodel.f
extractAIC(dmodel.f)
extractAIC(dmodel.r)
dmodel.rML = lmer(Yield~1|Batch, data=Dyestuff2, REML=FALSE)
extractAIC(dmodel.rML)
-2*logLik(dmodel.f)+2*6
