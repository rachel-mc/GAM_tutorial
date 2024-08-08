## Required packages:

library(DHARMa)
library(FSA)
library(gamlss)
library(gratia)
library(hnp)
library(itsadug)
library(MASS)
library(mgcv)
library(mgcViz)
library(MuMIn)
library(nlme)
library(nlraa)
library(performance)
library(psych)
library(tidyverse)


## Create a function to calculate BICc according to McQuarrie (1999)

BICc <- function (obj) {
  bic <- BIC(obj)
  n <- obj$df.null + 1
  p <- n - obj$df.residual
  bicc <- bic + (log(n) * (p + 1) * p)/(n - p - 1)
  return(bicc)
}


## Loading the three datasets used as case studies in the paper
## Note that these datasets are solely provided to reproduce the GAMM examples
## presented in the related manuscript. If some of these data are of any
## interest for other ecological reasons, please contact the main author.

TROUT <- read_delim("https://raw.githubusercontent.com/rachel-mc/GAM_tutorial/main/TROUT.txt")
WALLEYE <- read_delim("https://raw.githubusercontent.com/rachel-mc/GAM_tutorial/main/WALLEYE.txt")
CHARR <- read_delim("https://raw.githubusercontent.com/rachel-mc/GAM_tutorial/main/CHARR.txt")


##### Example 1: Length-at-age relationship of Aupaluk lake trout females ######

## Non-linear mixed-effects logistic growth model
## Define the logistic growth (LG) equation
  
LG <- function(x, Linf, K, t0) Linf / (1 + exp(-K * (x - t0)))


## Estimate initial parameters (ip) to be used as starting values with FSA 

ip_TROUT <- vbStarts(formula = TL ~ AGE, data = TROUT)
ip_TROUT


## Estimate the growth parameters for each LAKE with nlme

LAKE_PARMS <- nlsList(
  TL ~ LG(AGE, Linf, K, t0) | LAKE,
  data = TROUT,
  start = c(Linf = ip_TROUT$Linf, K = ip_TROUT$K, t0 = ip_TROUT$t0),
  na.action = na.omit)


## Fit the (global) non-linear mixed-effects logistic growth model with nlme

m_TROUT_LG <- nlme(
  TL ~ LG(AGE, Linf, K, t0),
  fixed = Linf + K + t0 ~ 1,
  random = Linf + K + t0 ~ 1 | LAKE,
  start = fixef(LAKE_PARMS),
  data = TROUT)
summary(m_TROUT_LG)


## Calculate the marginal and conditional R2 with nlraa

R2M(m_TROUT_LG)


## Extract the growth parameter estimates 

Linf <- fixef(m_TROUT_LG)[[1]]
K <- fixef(m_TROUT_LG)[[2]]
t0 <- fixef(m_TROUT_LG)[[3]]


## Predict TL (i.e., pred_TL) according to AGE from the logistic growth model

nd_TROUT_LG <- data.frame (AGE = seq(1, 50, by = 0.1))
pred_TROUT_LG <- transform(
  nd_TROUT_LG, pred_TL = Linf / (1 + exp(-K * (AGE - t0))))


## Original and predicted values (Figure 1 of main text, produced with Excel)

plot(TROUT$AGE, TROUT$TL)
plot(pred_TROUT_LG$AGE, pred_TROUT_LG$pred)


## Gaussian GAMM with the default parametrization for s()

TROUT$LAKE <- as.factor(TROUT$LAKE)

m_TROUT_g_k10 <- gam(
  TL ~ s(AGE) + s(LAKE, bs = "re"),
  family = gaussian,
  method = "ML",
  data = TROUT)

summary(m_TROUT_g_k10)


## All the other possible Gaussian GAMMs with k = 3 to 37
## Code has been voluntarily condensed to save space

m_TROUT_g_k3<-gam(TL~s(AGE,k=3)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k4<-gam(TL~s(AGE,k=4)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k5<-gam(TL~s(AGE,k=5)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k6<-gam(TL~s(AGE,k=6)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k7<-gam(TL~s(AGE,k=7)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k8<-gam(TL~s(AGE,k=8)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k9<-gam(TL~s(AGE,k=9)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k11<-gam(TL~s(AGE,k=11)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k12<-gam(TL~s(AGE,k=12)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k13<-gam(TL~s(AGE,k=13)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k14<-gam(TL~s(AGE,k=14)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k15<-gam(TL~s(AGE,k=15)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k16<-gam(TL~s(AGE,k=16)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k17<-gam(TL~s(AGE,k=17)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k18<-gam(TL~s(AGE,k=18)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k19<-gam(TL~s(AGE,k=19)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k20<-gam(TL~s(AGE,k=20)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k21<-gam(TL~s(AGE,k=21)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k22<-gam(TL~s(AGE,k=22)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k23<-gam(TL~s(AGE,k=23)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k24<-gam(TL~s(AGE,k=24)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k25<-gam(TL~s(AGE,k=25)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k26<-gam(TL~s(AGE,k=26)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k27<-gam(TL~s(AGE,k=27)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k28<-gam(TL~s(AGE,k=28)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k29<-gam(TL~s(AGE,k=29)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k30<-gam(TL~s(AGE,k=30)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k31<-gam(TL~s(AGE,k=31)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k32<-gam(TL~s(AGE,k=32)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k33<-gam(TL~s(AGE,k=33)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k34<-gam(TL~s(AGE,k=34)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k35<-gam(TL~s(AGE,k=35)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k36<-gam(TL~s(AGE,k=36)+s(LAKE,bs="re"),method="ML",data=TROUT)

m_TROUT_g_k37<-gam(TL~s(AGE,k=37)+s(LAKE,bs="re"),method="ML",data=TROUT)


## MODEL ADEQUACY
## Diagnostic plots with mgcv: example with m_TROUT_g_k10 (see Figure S2)

gam.check(m_TROUT_g_k10)


## Diagnostic plots with gratia: example with m_TROUT_g_k10 (see Figure S3)

appraise(m_TROUT_g_k10, method = "simulate")


## Worm plot with gratia: example with m_TROUT_g_k10 (see Figure S4) 

worm_plot(m_TROUT_g_k10, method = "simulate")


## Generate scaled quantile residuals from model-derived simulations (DHARMa)
## m_TROUT_g_k10 is again used as an example

simulationOutput <- simulateResiduals(fittedModel = m_TROUT_g_k10)


## Test the uniformity of the quantile residuals relative to theoretical 
## quantiles of the normal distribution (Figure S5)

testUniformity(simulationOutput)


## Test the distribution of the same residuals along the rank-transformed model
## predictions (Figure S6)

testQuantiles(simulationOutput)


## MODEL SELECTION
## Rank candidate models according to AICc (default option with MuMIn)

model.sel(
  m_TROUT_g_k5,
  m_TROUT_g_k6,
  m_TROUT_g_k7,
  m_TROUT_g_k8,
  m_TROUT_g_k9,
  m_TROUT_g_k10,
  m_TROUT_g_k11,
  m_TROUT_g_k12,
  m_TROUT_g_k13,
  m_TROUT_g_k14,
  m_TROUT_g_k15,
  m_TROUT_g_k16,
  m_TROUT_g_k17,
  m_TROUT_g_k18,
  m_TROUT_g_k19,
  m_TROUT_g_k20,
  m_TROUT_g_k21,
  m_TROUT_g_k22,
  m_TROUT_g_k23,
  m_TROUT_g_k24,
  m_TROUT_g_k25,
  m_TROUT_g_k26,
  m_TROUT_g_k27,
  m_TROUT_g_k28,
  m_TROUT_g_k29,
  m_TROUT_g_k30,
  m_TROUT_g_k31,
  m_TROUT_g_k32,
  m_TROUT_g_k33,
  m_TROUT_g_k34,
  m_TROUT_g_k35,
  m_TROUT_g_k36,
  m_TROUT_g_k37)


## Rank candidate models according to BICc

model.sel(
  m_TROUT_g_k5,
  m_TROUT_g_k6,
  m_TROUT_g_k7,
  m_TROUT_g_k8,
  m_TROUT_g_k9,
  m_TROUT_g_k10,
  m_TROUT_g_k11,
  m_TROUT_g_k12,
  m_TROUT_g_k13,
  m_TROUT_g_k14,
  m_TROUT_g_k15,
  m_TROUT_g_k16,
  m_TROUT_g_k17,
  m_TROUT_g_k18,
  m_TROUT_g_k19,
  m_TROUT_g_k20,
  m_TROUT_g_k21,
  m_TROUT_g_k22,
  m_TROUT_g_k23,
  m_TROUT_g_k24,
  m_TROUT_g_k25,
  m_TROUT_g_k26,
  m_TROUT_g_k27,
  m_TROUT_g_k28,
  m_TROUT_g_k29,
  m_TROUT_g_k30,
  m_TROUT_g_k31,
  m_TROUT_g_k32,
  m_TROUT_g_k33,
  m_TROUT_g_k34,
  m_TROUT_g_k35,
  m_TROUT_g_k36,
  m_TROUT_g_k37,
  rank = BICc)


## Compare the best AICc- and BICc-retained models with the proposed
## Chi-squared test of itsadug

compareML(m_TROUT_g_k5, m_TROUT_g_k11, suggest.report = TRUE)


## Generate a single mgcViz diagnostic plot based on 100 model-derived
## simulations (Figure S7)

model <- m_TROUT_g_k10
predictor <- "AGE"

set.seed(2024)
viz <- getViz(model, nsim = 100)
plot <- check1D(viz, predictor) + l_gridCheck1D(n = 100)
diag_plot <- plot$ggObj
est <- diag_plot$layers[[1]]$data
ll <- diag_plot$layers[[3]]$data$ll
ul <- diag_plot$layers[[3]]$data$ul
diag <- cbind(est, ll, ul)
diagnostic <- mutate(diag,
                     ll_diff = y - ll,
                     ul_diff = ul - y,
                     test = ll_diff * ul_diff)
n <- length(diagnostic$test)
i_n <- sum(diagnostic$test >= 0)
mgcViz_perc <- i_n / n * 100


## Print the mgcViz plot

diag_plot


## Obtain the percentage of the binned residuals' means found 
## within the simulated 80% confidence limits

mgcViz_perc


## Calculate the mgcViz score of a given model based on 100 iterations with
## 100 simulations/iteration. The mgcViz score with associated 95% uncertainty
## interval for m_TROUT_g_k10 is shown below

model <- m_TROUT_g_k10
predictor <- "AGE"

set.seed(2024)
viz_fun <- function(model, predictor) {
  viz <- getViz(model, nsim = 100)
  plot <- check1D(viz, predictor) + l_gridCheck1D(n = 100)
  diag_plot <- plot$ggObj
  est <- diag_plot$layers[[1]]$data
  ll <- diag_plot$layers[[3]]$data$ll
  ul <- diag_plot$layers[[3]]$data$ul
  diag <- cbind(est, ll, ul)
  diagnostic <- mutate(diag,
                       ll_diff = y - ll,
                       ul_diff = ul - y,
                       test = ll_diff * ul_diff)
  n <- length(diagnostic$test)
  i_n <- sum(diagnostic$test >= 0)
  mgcViz_perc <- i_n / n * 100
  mgcViz_perc
  }

## This can take several seconds to run
summary_vf <- replicate(100, viz_fun(model, predictor))
mgcViz_score <- mean(summary_vf)
ll_mgcViz_score <- (quantile(summary_vf, probs = 0.025))[[1]]
ul_mgcViz_score <- (quantile(summary_vf, probs = 0.975))[[1]]
cbind(mgcViz_score, ll_mgcViz_score, ul_mgcViz_score)


## Compare the Gaussian GAMM with k = 6 (best mgcViz score) to the AICc-
## and BICc-retained models

compareML(m_TROUT_g_k6, m_TROUT_g_k11, suggest.report = TRUE)
compareML(m_TROUT_g_k6, m_TROUT_g_k5, suggest.report = TRUE)


## Compare the Gaussian GAMM with k = 6 and the non-linear mixed-effects
## logistic growth model based on AICc or BIC

AICc(m_TROUT_LG, m_TROUT_g_k6)
BIC(m_TROUT_LG, m_TROUT_g_k6)


## Check for potential concurvity issues

concurvity(m_TROUT_g_k6)


## Including the additive effect of MASS fitted as a covariate with AGE to
## generate possible concurvity issues since these are linearly correlated
## Pearson's correlation r

cor.test(TROUT$AGE, TROUT$MASS)


## Both covariates are not normally distributed but exhibit a strong
## monotonous, curvilinear statistical association
## Spearman's rank correlation rho 

cor.test(TROUT$AGE, TROUT$MASS, method = "spearman", exact = FALSE)


## Fitting smooth functions to AGE and MASS with default mgcv parametrization

m_TROUT_g_k10_MASS <- gam(
  TL ~ s(AGE) + s(MASS) + s(LAKE, bs = "re"),
  family = gaussian, method = "ML", data = TROUT)

summary(m_TROUT_g_k10_MASS)


## Check for potential concurvity issues

concurvity(m_TROUT_g_k10_MASS)


## Obtain predictions from the best-retained Gaussian GAMM with k = 6
## Figure 1 of main text: same approach for k = 11 (AICc) and k = 5 (BICc)
## by simply specifying the desired model from which predictions should be
## produced on the first line of code

model <- m_TROUT_g_k6

nd_TROUT <- data.frame(
  AGE = seq(1, 50, by = 0.1),
  LAKE = "BRULE")

fitted_TROUT <- predict(
  model,
  nd_TROUT,
  type = "link",	 
  exclude = "s(LAKE)",
  se.fit = TRUE)

pred_TROUT <- fitted_TROUT$fit
LL_TROUT <- fitted_TROUT$fit - 1.96 * fitted_TROUT$se.fit
UL_TROUT <- fitted_TROUT$fit + 1.96 * fitted_TROUT$se.fit
results_TROUT <- cbind(nd_TROUT$AGE, pred_TROUT, LL_TROUT, UL_TROUT)
results_TROUT 


## Predicted values were used in Figure 1 (produced with Excel)
## Here for the central tendency only of k = 6 (best model)

plot(nd_TROUT$AGE, pred_TROUT)


##### Example 2: CPUE temporal trend of small St. Lawrence River walleyes #####

## Observed yearly CPUE data per study AREA (Figure 3a)
## Example for the mean and 95% CI estimated by a profile-likelihood function
## of the CPUE data from AREA = BEBA and YEAR = 2001

WALLEYE_2001 <- WALLEYE[which(WALLEYE$YEAR == 2001),]

WALLEYE_2001_BEBA <- WALLEYE_2001[which(WALLEYE_2001$AREA == "BEBA"),]

m_WALLEYE_2001_BEBA_nb2_NULL <- glm.nb(
  CPUE ~ 1, 
  data = WALLEYE_2001_BEBA)

nd_WALLEYE_2001_BEBA <- data.frame(
  YEAR = 2001,
  AREA = c("BEBA"))

fitted <- predict(
  m_WALLEYE_2001_BEBA_nb2_NULL,
  nd_WALLEYE_2001_BEBA,
  type = "link",
  se.fit = TRUE)

pred<-exp(fitted$fit)

CI <- confint(m_WALLEYE_2001_BEBA_nb2_NULL)
ll <- exp(CI[[1]])
ul <- exp(CI[[2]])
results_WALLEYE_2001_BEBA <- cbind(nd_WALLEYE_2001_BEBA, pred, ll, ul)
results_WALLEYE_2001_BEBA


## Standardization of CPUE data to produce a temporal trend, first with a 
## Poisson GAMM and k = 10

WALLEYE$AREA<-factor(WALLEYE$AREA)

m_WALLEYE_p_k10 <- gam(
  CPUE ~ s(YEAR) + CONDUCT + D1AUG + DEPTH + EFFORT + TEMP + TURBID +
  s(AREA, bs = "re"),
  family = poisson(),
  method = "ML",
  data = WALLEYE)

summary(m_WALLEYE_p_k10)


## MODEL ADEQUACY (Poisson)
## Diagnostic plots based on deviance residuals of m_WALLEYE_p_k10 
## as an example with gratia (Figure S9)

appraise(m_WALLEYE_p_k10, method = "simulate")


## Worm plot based on Pearson residuals of m_WALLEYE_p_k10 
## as an example with gratia (Figure S10)

worm_plot(m_WALLEYE_p_k10, type = "pearson", method = "simulate")


## Rootogram as an alternative diagnostic plot (Figure S11)

rg_WALLEYE_p_k10 <- rootogram(m_WALLEYE_p_k10)
draw(rg_WALLEYE_p_k10)


## Generate scaled randomized quantile residuals from model-derived simulations
## with DHARMa (m_WALLEYE_p_k10)

simulationOutput <- simulateResiduals(fittedModel = m_WALLEYE_p_k10)


## Test the uniformity of the randomized quantile residuals relative to
## the theoretical quantiles of the Poisson distribution (Figure S12)

testUniformity(simulationOutput)


## Test the distribution of the same residuals along the rank-transformed 
## model predictions (Figure S13)

testQuantiles(simulationOutput)


## Diagnostic based on half-normal plots of the same considered model

model <- m_WALLEYE_p_k10

dfun <- function(obj) resid(obj, type = "pearson")

sfun <- function(n, obj){
  y <- rpois(nrow(WALLEYE),
  lambda = predict(model, type = "response"))
  return(y)
  }

ffun <- function(resp) {
  gam(resp ~ s(YEAR) + CONDUCT + D1AUG + DEPTH
  + EFFORT + TEMP + TURBID + s(AREA, bs = "re"), 
  family = poisson, method = "ML", data = WALLEYE)
  }

## This can take several minutes to run
set.seed(2024)
hfun <- list()
  for(i in 1:10) {
  hfun[[i]] <- hnp(
  model, newclass = TRUE, diagfun = dfun, 
  simfun = sfun, fitfun = ffun, how.many.out = TRUE, 
  plot.sim = FALSE)
  }

hnp_summary <- sapply(hfun, function(x) x$out/x$total*100) 
mean(hnp_summary)


## Produce a single hnp diagnostic plot for m_WALLEYE_p_k10 (Figure S14)
## The model and helper functions used were defined in the previous step
## This can take several seconds to run

set.seed(2024)
hnp(model, newclass = TRUE, diagfun = dfun, simfun = sfun,
    fitfun = ffun, how.many.out = TRUE, plot.sim = TRUE, paint = TRUE,
    ylab = "Pearson residuals")


## Test for plausible lack of equidispersion (likely overdispersion)
## The same scaled randomized quantile residuals as defined above are used

testDispersion(simulationOutput, plot = FALSE)


## Compare the overall variance-to-mean ratio in CPUE data

var(WALLEYE$CPUE)
mean(WALLEYE$CPUE)
var(WALLEYE$CPUE) / mean(WALLEYE$CPUE)


## CPUE data are overdispersed, so the negative binomial type-II (NB2)
## extension is used to produce different possible candidate models

WALLEYE$AREA <- factor(WALLEYE$AREA)


## "NULL" NB2 GAMM that does not include YEAR as a predictor
## to mimic temporal stability
m_WALLEYE_nb2_NULL <- gam(
  CPUE ~ CONDUCT + D1AUG + DEPTH + EFFORT + TEMP + TURBID +
  s(AREA, bs = "re"), family = nb, method = "ML", data = WALLEYE)


## "linear" NB2 GAMM which fit YEAR as a parametric component
## No smooth function is applied
m_WALLEYE_nb2 <- gam(
  CPUE ~ YEAR + CONDUCT + D1AUG + DEPTH + EFFORT + TEMP + TURBID +
  s(AREA, bs = "re"), family = nb, method = "ML", data = WALLEYE)


## All possible values for k (minimum of 3 to 20, as 20 years were
## surveyed overall)
m_WALLEYE_nb2_k3<-gam(CPUE~s(YEAR,k=3)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k4<-gam(CPUE~s(YEAR,k=4)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k5<-gam(CPUE~s(YEAR,k=5)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k6<-gam(CPUE~s(YEAR,k=6)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k7<-gam(CPUE~s(YEAR,k=7)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k8<-gam(CPUE~s(YEAR,k=8)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k9<-gam(CPUE~s(YEAR,k=9)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k10<-gam(CPUE~s(YEAR,k=10)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k11<-gam(CPUE~s(YEAR,k=11)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k12<-gam(CPUE~s(YEAR,k=12)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k13<-gam(CPUE~s(YEAR,k=13)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k14<-gam(CPUE~s(YEAR,k=14)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k15<-gam(CPUE~s(YEAR,k=15)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k16<-gam(CPUE~s(YEAR,k=16)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k17<-gam(CPUE~s(YEAR,k=17)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k18<-gam(CPUE~s(YEAR,k=18)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k19<-gam(CPUE~s(YEAR,k=19)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)

m_WALLEYE_nb2_k20<-gam(CPUE~s(YEAR,k=20)+CONDUCT+D1AUG+DEPTH+EFFORT+TEMP+TURBID+
s(AREA,bs="re"),family=nb,method="ML",data=WALLEYE)


## For the model in which YEAR is fitted as a parametric component, one can
## assess for possible multicollinearity issues among the linear predictors

check_collinearity(m_WALLEYE_nb2)


## Check the pairwise correlations among all the covariates when YEAR is still
## fitted as a parametric component with the psych package

KEEP <- c("YEAR", "CONDUCT", "D1AUG", "DEPTH", "EFFORT", "TEMP", "TURBID")
W <- WALLEYE[KEEP]
corr.test(W)


## MODEL ADEQUACY (NB2) for the GAMM with k = 9 ultimately identified as
## the top-ranking one according to a Worm plot based on Pearson residuals
## with gratia (Figure S15)

worm_plot(m_WALLEYE_nb2_k9, type = "pearson", method = "simulate")


## Rootogram of m_WALLEYE_nb2_k9 with gratia (Figure S16) 

rg_WALLEYE_nb2_k9 <- rootogram(m_WALLEYE_nb2_k9)
draw(rg_WALLEYE_nb2_k9)


## Test adequacy with DHARMa using the same model (k = 9) as an example
## (Figures S17-18)
simulationOutput<-simulateResiduals(fittedModel = m_WALLEYE_nb2_k9)
testUniformity(simulationOutput, plot = TRUE)
testQuantiles(simulationOutput, plot = TRUE)


## Diagnostic based on half-normal plots of the same considered model (k = 9)

model <- m_WALLEYE_nb2_k9

dfun <- function(obj) resid(obj, type = "pearson")

sfun <- function(n, obj) {
  y <- rnbinom(nrow(WALLEYE),
  size = model$family$getTheta(TRUE),
  mu = predict(model, type = "response"))
  return(y)
  }

ffun<-function(resp) {
  gam(resp ~ s(YEAR, k = 9) + CONDUCT + D1AUG + DEPTH + EFFORT + TEMP + TURBID +
  s(AREA, bs = "re"), family = nb, method = "ML", data = WALLEYE)
  }

## This can take several minutes to run
set.seed(2024)
hfun <- list()
  for(i in 1:10) {
  hfun[[i]] <- hnp(
  model, newclass = TRUE, diagfun = dfun, simfun = sfun,
  fitfun = ffun, how.many.out = TRUE, plot.sim = FALSE)
  }

hnp_summary <- sapply(hfun, function(x) x$out/x$total*100) 
mean(hnp_summary)


## Produce a single hnp diagnostic plot for m_WALLEYE_nb2_k9 (Figure S19)
## The model and helper functions used were defined in the previous step
## This can take several seconds to run

set.seed(2024)
hnp(model, newclass = TRUE, diagfun = dfun, simfun = sfun,
    fitfun = ffun, how.many.out = TRUE, plot.sim = TRUE, paint = TRUE,
    ylab = "Pearson residuals")


## Zero-inflation test applied to m_WALLEYE_nb2_k9 as an example

simulationOutput <- simulateResiduals(fittedModel = m_WALLEYE_nb2_k9)
testZeroInflation(simulationOutput, plot = FALSE)


## MODEL SELECTION
## Information-theoretic approach based on the default AICc with MuMIn

model.sel(
  m_WALLEYE_nb2_NULL,
  m_WALLEYE_nb2,
  m_WALLEYE_nb2_k3,
  m_WALLEYE_nb2_k4,
  m_WALLEYE_nb2_k5,
  m_WALLEYE_nb2_k6,
  m_WALLEYE_nb2_k7,
  m_WALLEYE_nb2_k8,
  m_WALLEYE_nb2_k9,
  m_WALLEYE_nb2_k10,
  m_WALLEYE_nb2_k11,
  m_WALLEYE_nb2_k12,
  m_WALLEYE_nb2_k13,
  m_WALLEYE_nb2_k14,
  m_WALLEYE_nb2_k15,
  m_WALLEYE_nb2_k16,
  m_WALLEYE_nb2_k17,
  m_WALLEYE_nb2_k18,
  m_WALLEYE_nb2_k19,
  m_WALLEYE_nb2_k20)


## Rank candidate models according to BICc

model.sel(
  m_WALLEYE_nb2_NULL,
  m_WALLEYE_nb2,
  m_WALLEYE_nb2_k3,
  m_WALLEYE_nb2_k4,
  m_WALLEYE_nb2_k5,
  m_WALLEYE_nb2_k6,
  m_WALLEYE_nb2_k7,
  m_WALLEYE_nb2_k8,
  m_WALLEYE_nb2_k9,
  m_WALLEYE_nb2_k10,
  m_WALLEYE_nb2_k11,
  m_WALLEYE_nb2_k12,
  m_WALLEYE_nb2_k13,
  m_WALLEYE_nb2_k14,
  m_WALLEYE_nb2_k15,
  m_WALLEYE_nb2_k16,
  m_WALLEYE_nb2_k17,
  m_WALLEYE_nb2_k18,
  m_WALLEYE_nb2_k19,
  m_WALLEYE_nb2_k20,
  rank = BICc)


## Compare the best AICc- and BICc-retained models with the proposed
## Chi-squared test of itsadug

compareML(m_WALLEYE_nb2_NULL, m_WALLEYE_nb2_k17, suggest.report = TRUE)


## Compare candidate models based on the mgcViz score (Figure S20),
## here computed for m_WALLEYE_nb2_k9 with 95% uncertainty interval

model <- m_WALLEYE_nb2_k9
predictor <- "YEAR"

set.seed(2024)
viz_fun <- function(model, predictor) {
  viz <- getViz(model, nsim = 100)
  plot <- check1D(viz, predictor) + l_gridCheck1D(n = 100)
  diag_plot <- plot$ggObj
  est <- diag_plot$layers[[1]]$data
  ll <- diag_plot$layers[[3]]$data$ll
  ul <- diag_plot$layers[[3]]$data$ul
  diag <- cbind(est, ll, ul)
  diagnostic <- mutate(diag,
                       ll_diff = y - ll,
                       ul_diff = ul - y,
                       test = ll_diff * ul_diff)
  n <- length(diagnostic$test)
  i_n <- sum(diagnostic$test >= 0)
  mgcViz_perc <- i_n / n * 100
  mgcViz_perc
  }

## This can take several seconds to run
summary_vf <- replicate(100, viz_fun(model, predictor))
mgcViz_score <- mean(summary_vf)
ll_mgcViz_score <- (quantile(summary_vf, probs = 0.025))[[1]]
ul_mgcViz_score <- (quantile(summary_vf, probs = 0.975))[[1]]
cbind(mgcViz_score, ll_mgcViz_score, ul_mgcViz_score)


## Assess for plausible temporal autocorrelation with the acf() function
## of the default stats package

acf(resid(m_WALLEYE_nb2_k9), plot = FALSE)$acf[2]


## Refit the same GAMM with the bam() function, method = fREML, and argument
## discrete = TRUE. This model is required to assess for temporal
## autocorrelation with itsadug (next steps)

m_WALLEYE_nb2_k9_fREML <- bam(
  CPUE ~ s(YEAR, k = 9) + CONDUCT + D1AUG + DEPTH + EFFORT + TEMP + TURBID +
  s(AREA, bs = "re"), 
  family = nb, 
  method = "fREML", 
  discrete = TRUE, 
  data = WALLEYE)

summary(m_WALLEYE_nb2_k9_fREML)


## Estimate the first-order autoregressive (AR1) coefficient from the refitted
## model and store it as RHO and produce the related ACF plot (Figure S21)

RHO <- start_value_rho(m_WALLEYE_nb2_k9_fREML, plot = TRUE)
RHO


## Refit the model by now including an autocorrelation matrix, as in a 
## Generalized Estimating Equation (GEE) with the argument rho = RHO

m_WALLEYE_nb2_k9_fREML_AR1 <- bam(
  CPUE ~ s(YEAR, k = 9) + CONDUCT + D1AUG + DEPTH + EFFORT + TEMP + TURBID +
  s(AREA, bs = "re"), 
  family = nb, 
  method = "fREML", 
  discrete = TRUE,
  rho = RHO,
  data = WALLEYE)

summary(m_WALLEYE_nb2_k9_fREML_AR1)


## Check whether adding the AR1 term (RHO) helped to reduce temporal
## autocorrelation (Figure S22) - press "Enter" in RStudio console to see plot

check_resid(m_WALLEYE_nb2_k9_fREML_AR1, select = 3)


## Compare the two models with and without the AR1 term with itsadug

compareML(m_WALLEYE_nb2_k9_fREML, m_WALLEYE_nb2_k9_fREML_AR1)


## Test if the retained model that accounts for temporal autocorrelation
## correctly predicts the number of observed zeros as a further check

simulationOutput <- simulateResiduals(fittedModel = m_WALLEYE_nb2_k9_fREML_AR1)
testZeroInflation(simulationOutput, plot = FALSE)


## Check for possible concurvity issues

concurvity(m_WALLEYE_nb2_k9_fREML_AR1)


## Assess the (adjusted) explanatory power of this retained model

model <- m_WALLEYE_nb2_k9

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
Total_df <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - Total_df) * (100 - D2))
D2_adj


## Predictions

nd_WALLEYE<-data.frame(
  YEAR = seq(2001, 2021, by = 0.1),
  CONDUCT = mean(WALLEYE$CONDUCT),
  D1AUG = mean(WALLEYE$D1AUG),
  DEPTH = mean(WALLEYE$DEPTH),
  EFFORT = mean(WALLEYE$EFFORT),
  TEMP = mean(WALLEYE$TEMP),
  TURBID = mean(WALLEYE$TURBID),
  AREA = "ALSP")

model <- m_WALLEYE_nb2_k9_fREML_AR1

fitted_WALLEYE <- predict(
  model,
  nd_WALLEYE,
  type = "link",
  exclude = "s(AREA)",
  se.fit = TRUE)

pred_WALLEYE <- exp(fitted_WALLEYE$fit)
LL_WALLEYE <- exp(fitted_WALLEYE$fit - 1.96 * fitted_WALLEYE$se.fit)
UL_WALLEYE <- exp(fitted_WALLEYE$fit + 1.96 * fitted_WALLEYE$se.fit)
results_WALLEYE <- cbind(nd_WALLEYE$YEAR, pred_WALLEYE, LL_WALLEYE, UL_WALLEYE)
results_WALLEYE 

## The predicted values were used in Figure 3b (produced with Excel)
## Here for the central tendency only of m_WALLEYE_nb2_k9_fREML_AR1

plot(nd_WALLEYE$YEAR, pred_WALLEYE)

#### Example 3: Prob. to observe developing gonads in Arctic charr females ####

## Compare the three links: logit, probit and complementary log-log (cloglog)
## No smooth function is applied to the fixed effect: GLMM-like GAMM

CHARR$RIVER <- as.factor(CHARR$RIVER)

## logit link
m_CHARR_logit <- gam(
  GONADS ~ FL +s(RIVER, bs = "re"),
  family = binomial(link = logit), 
  method = "ML", 
  data = CHARR)

summary(m_CHARR_logit)


## probit link
m_CHARR_probit <- gam(
  GONADS ~ FL +s(RIVER, bs = "re"),
  family = binomial(link = probit), 
  method = "ML", 
  data = CHARR)

summary(m_CHARR_probit)


## cloglog link
m_CHARR_cloglog <- gam(
  GONADS ~ FL +s(RIVER, bs = "re"),
  family = binomial(link = cloglog), 
  method = "ML", 
  data = CHARR)

summary(m_CHARR_cloglog)


## Assess for a possible lack-of-fit with the test proposed
## by Osius and Rojek (1992) | p-value should be > 0.05

o.r.test <- function(obj) {
  mod <- obj$model
  trials <- rep(1, nrow(mod))
  if(any(colnames(mod) == "(weights)")) 
  trials <- mod[[ncol(mod)]]
  prop <- mod[[1]]
  if(is.factor(prop))
  prop <- as.numeric(prop) == 2
  pi_hat <- obj$fitted.values
  y <- trials * prop
  yhat <- trials * pi_hat
  nu <- yhat * (1 - pi_hat)
  pearson <- sum((y - yhat)^2 / nu)
  c <- (1 - 2 * pi_hat) / nu
  exclude <- c(1, which(colnames(mod) == "(weights)"))
  vars <- data.frame(c, mod[, -exclude]) 
  wlr <- lm(c ~ ., weights = nu, data = vars)
  rss <- sum(nu * residuals(wlr)^2)
  j <- nrow(mod)
  a <- 2 * (j - sum(1 / trials))
  z <- (pearson - (j - ncol(vars) - 1)) / sqrt(a + rss)
  p_value <- 2 * (1 - pnorm(abs(z)))
  cat("z =", z, "with p-value =", p_value, "\n")
}

o.r.test(m_CHARR_logit)
o.r.test(m_CHARR_probit)
o.r.test(m_CHARR_cloglog)


## Assess for possible link misspecification with the proposed test
## of McCullagh and Nelder (1989) | p-value should be > 0.05

eta2 <- predict(m_CHARR_logit)^2
m_CHARR_logit_eta2 <- update(m_CHARR_logit, . ~ . + eta2)
anova(m_CHARR_logit, m_CHARR_logit_eta2, test = "Chisq")

eta2 <- predict(m_CHARR_probit)^2
m_CHARR_probit_eta2 <- update(m_CHARR_probit, . ~ . + eta2)
anova(m_CHARR_probit, m_CHARR_probit_eta2, test = "Chisq")

eta2 <- predict(m_CHARR_cloglog)^2
m_CHARR_cloglog_eta2 <- update(m_CHARR_cloglog, . ~ . + eta2)
anova(m_CHARR_cloglog, m_CHARR_cloglog_eta2, test = "Chisq")


## Model selection based on AICc

model.sel(
  m_CHARR_logit,
  m_CHARR_probit,
  m_CHARR_cloglog)


## Model selection based on BICc

model.sel(
  m_CHARR_logit,
  m_CHARR_probit,
  m_CHARR_cloglog,
  rank = BICc)


## Use the diagnostic plots of gratia for additional
## adequacy assessment (Figure S24)

appraise(m_CHARR_cloglog, type = "pearson", method = "simulate")


## Generate scaled randomized quantile residuals from model-derived
## simulations with DHARMa

simulationOutput <- simulateResiduals(fittedModel = m_CHARR_cloglog)


## Test the uniformity of the randomized quantile residuals relative to
## the theoretical quantiles of the binomial distribution (Figure S25)

testUniformity(simulationOutput)


## Test the distribution of the same residuals along the rank-transformed
## model predictions (Figure S26)

testQuantiles(simulationOutput)


## Can additional non-linearities be incorporated in the cloglog sigmoidal
## curve to to better describe the reproduction ogive?

m_CHARR_cloglog_k3<-gam(GONADS~s(FL,k=3)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k4<-gam(GONADS~s(FL,k=4)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k5<-gam(GONADS~s(FL,k=5)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k6<-gam(GONADS~s(FL,k=6)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k7<-gam(GONADS~s(FL,k=7)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k8<-gam(GONADS~s(FL,k=8)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k9<-gam(GONADS~s(FL,k=9)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k10<-gam(GONADS~s(FL,k=10)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k15<-gam(GONADS~s(FL,k=15)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k20<-gam(GONADS~s(FL,k=20)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k25<-gam(GONADS~s(FL,k=25)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k50<-gam(GONADS~s(FL,k=50)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)

m_CHARR_cloglog_k100<-gam(GONADS~s(FL,k=100)+s(RIVER,bs="re"),
family=binomial(link=cloglog),method="ML",data=CHARR)


## MODEL SELECTION
## Information-theoretic approach based on AICc

model.sel(
  m_CHARR_cloglog,
  m_CHARR_cloglog_k3,
  m_CHARR_cloglog_k4,
  m_CHARR_cloglog_k5,
  m_CHARR_cloglog_k6,
  m_CHARR_cloglog_k7,
  m_CHARR_cloglog_k8,
  m_CHARR_cloglog_k9,
  m_CHARR_cloglog_k10,
  m_CHARR_cloglog_k15,
  m_CHARR_cloglog_k20,
  m_CHARR_cloglog_k25,
  m_CHARR_cloglog_k50,
  m_CHARR_cloglog_k100)


## Information-theoretic approach based on BICc

model.sel(
  m_CHARR_cloglog,
  m_CHARR_cloglog_k3,
  m_CHARR_cloglog_k4,
  m_CHARR_cloglog_k5,
  m_CHARR_cloglog_k6,
  m_CHARR_cloglog_k7,
  m_CHARR_cloglog_k8,
  m_CHARR_cloglog_k9,
  m_CHARR_cloglog_k10,
  m_CHARR_cloglog_k15,
  m_CHARR_cloglog_k20,
  m_CHARR_cloglog_k25,
  m_CHARR_cloglog_k50,
  m_CHARR_cloglog_k100,
  rank = BICc)


## Assess the (adjusted) explanatory power of the retained model

model <- m_CHARR_cloglog

D2 <- 100 * (1 - model$deviance / model$null.deviance)
logLik <- logLik(model)
Total_df <- attributes(logLik)$df
n <- summary(model)$n
D2_adj <- 100 - ((n - 1) / (n - Total_df) * (100 - D2))
D2_adj


## Calculate the mgcViz score of m_CHARR_cloglog and associated 
## 95% uncertainty interval

model <- m_CHARR_cloglog
predictor <- "FL"

set.seed(2024)
viz_fun <- function(model, predictor) {
  viz <- getViz(model, nsim = 100)
  plot <- check1D(viz, predictor) + l_gridCheck1D(n = 100)
  diag_plot <- plot$ggObj
  est <- diag_plot$layers[[1]]$data
  ll <- diag_plot$layers[[3]]$data$ll
  ul <- diag_plot$layers[[3]]$data$ul
  diag <- cbind(est, ll, ul)
  diagnostic <- mutate(diag,
                       ll_diff = y - ll,
                       ul_diff = ul - y,
                       test = ll_diff * ul_diff)
  n <- length(diagnostic$test)
  i_n <- sum(diagnostic$test >= 0)
  mgcViz_perc <- i_n / n * 100
  mgcViz_perc
}

## This can take several seconds to run
summary_vf <- replicate(100, viz_fun(model, predictor))
mgcViz_score <- mean(summary_vf)
ll_mgcViz_score <- (quantile(summary_vf, probs = 0.025))[[1]]
ul_mgcViz_score <- (quantile(summary_vf, probs = 0.975))[[1]]
cbind(mgcViz_score, ll_mgcViz_score, ul_mgcViz_score)


## Predicted values with 95% CI

min(CHARR$FL)
max(CHARR$FL)

nd_CHARR <- data.frame(
  FL = seq(182, 674, by = 1),
  RIVER = "VOLTZ")

model <- m_CHARR_cloglog

fitted_CHARR <- predict(
  model, 
  nd_CHARR,
  type = "link",
  exclude = "s(RIVER)",
  se.fit = TRUE)

pred_CHARR <- 1 - exp(-exp(fitted_CHARR$fit))
LL_pred_CHARR <- 1 - exp(-exp(fitted_CHARR$fit - 1.96 * fitted_CHARR$se.fit))
UL_pred_CHARR <- 1 - exp(-exp(fitted_CHARR$fit + 1.96 * fitted_CHARR$se.fit))
results_CHARR <- cbind(nd_CHARR$FL, pred_CHARR, LL_pred_CHARR, UL_pred_CHARR)
results_CHARR 

## The predicted values were used in Figure 4 (produced with Excel)
## Here for the central tendency only of m_CHARR_cloglog

plot(nd_CHARR$FL, pred_CHARR)


## Estimate the L50 and its SE according to the Delta method (default: p = 0.5)
## The 95% CI of the L50 estimate presented in Figure 4 of the main text
## was produced under normal theory (i.e., ± 1.96 * SE)

dose.p(m_CHARR_cloglog)


## Estimate the L75 and its SE according to the Delta method

dose.p(m_CHARR_cloglog, p = 0.75)