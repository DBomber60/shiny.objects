library(tidyverse)
library(gfoRmula)
library(rstanarm)
# generate some data to validate the gfoRmula function
# TODO: can we compute the true Y^{11} for this example??

set.seed(1)
n = 500
X1 = rnorm(n) # let's assume X is CD4
A1 = rbinom(n, 1, prob = plogis(-1 - X1 )) # higher probability of treatment for low CD4
X2 = rnorm(n, mean = X1 + A1, sd = 1)
A2 = rbinom(n, 1, prob = plogis(-1 - X2 )) # higher probability of treatment for low CD4
Y = rnorm(n, 1, mean = (X2 + A2)) # final outcome

consdat = data.frame(A1 = A1, A2 = A2, X1 = X1, X2 = X2, Y = Y)

dat = consdat %>% pivot_longer(!Y, names_to = c(".value", "set"), names_pattern = "(.)(.)")

dat$id = rep(1:nrow(consdat), each = 2) # to-do.. this automatically (using regex)
dat$set = as.numeric(dat$set) - 1
head(dat)


# use gfoRmula package
id = 'id'
time_name = 'set'
covnames = c('X', 'A')
outcome_name = 'Y'
covtypes = c('normal', 'binary')
histories = c(lagged)
histvars = list(c('X', 'A'))
covparams = list(covmodels=c(X ~ lag1_X + lag1_A,
                             A ~ X + lag1_A))
ymodel = Y ~ A + X
intvars = list('A', 'A')
interventions = list( list(c(static, rep(0,2))),
                      list(c(static, c(1,1))) )  
int_descript = c('Never', 'Always')

b = gformula_continuous_eof(dat, id=id, time_name = time_name, covnames = covnames, covtypes = covtypes,
                        covparams = covparams, histvars = histvars, histories = histories, outcome_name = outcome_name,
                        ymodel = ymodel, intvars = intvars, interventions = interventions, int_descript = int_descript,
                        seed = 1, nsamples = 100, ref_int = 2) # 0.2771665/ 0.2697512

# 0.30787258/ 0.03070607
lower = b$result$`MD lower 95% CI`[2] # lower
upper = b$result$`MD upper 95% CI`[2]
mean = b$result$`Mean difference`[2] # -1.79


#### what is our estimate if we instead consider X2 to be a baseline covariate (not time-varying)


covparams2 = list(covmodels=c(X ~ lag1_X,
                             A ~ X + lag1_A))
ymodel = Y ~ A + X
intvars = list('A', 'A')
interventions = list( list(c(static, rep(0,2))),
                      list(c(static, c(1,1))) )  
int_descript = c('Never', 'Always')

b2 = gformula_continuous_eof(dat, id=id, time_name = time_name, covnames = covnames, covtypes = covtypes,
                            covparams = covparams2, histvars = histvars, histories = histories, outcome_name = outcome_name,
                            ymodel = ymodel, intvars = intvars, interventions = interventions, int_descript = int_descript,
                            seed = 1, nsamples = 100, ref_int = 2) # 0.2771665/ 0.2697512
b2$result

# TODO 
# 1. show sensitivity of causal effect estimate to different modeling assumptions
# 2. prospectus ---> do asthma part!
# 3. prospectus ---> do causal inference part!
# 4. WIHS - more EDA (what centers are these patients from)

# now, let's do Bayesian estimates of this?
# X1 --> A1 --> X2 --> A2 --> Y
# sample from the posterior of relevant parameters, then sample outcome

bd = dat %>% group_by(id) %>% mutate(lag_A = lag(A), lag_X = lag(X)) %>% na.omit
xmod = stan_glm(X ~ lag_A + lag_X, family = gaussian(), data = bd)
ymod = stan_glm(Y ~ A + X, family = gaussian(), data = bd)

# Y ^ 11
nrep = 30
pp.X11 = posterior_predict(xmod, newdata = data.frame(lag_A = 1, lag_X = bd$lag_X), draws = nrep)
pp.X10 = posterior_predict(xmod, newdata = data.frame(lag_A = 0, lag_X = bd$lag_X), draws = nrep)

# these will be 500 * nrep long
pp.Y11 = posterior_predict(ymod, newdata = data.frame(A = 1, X = array(t(pp.X11))), draws = nrep)
pp.Y00 = posterior_predict(ymod, newdata = data.frame(A = 0, X = array(t(pp.X10))), draws = nrep)


mean(pp.Y11 - pp.Y00)

plot.df = data.frame(Y = c( array(t(pp.Y11)), array(t(pp.Y00)) ) ,
                     keep = rep( c(1, rep(0,nrep-1))    , each = n) ) %>% filter(keep == 1)

plot.df$nrep = rep(1:(nrep*2), each = 500)
plot.df$int = c( rep("11", 500 * nrep), rep("00", 500 * nrep))

#Plot.
cbp1 <- c("#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(plot.df, aes(x = Y, group = nrep)) + geom_density(aes(color=int)) + 
  theme_minimal() + theme(legend.position = "none") + scale_colour_manual(values=cbp1)

plot.df %>% group_by(int) %>% summarise(m = mean(Y))

# now show posterior predictive under different assumption
