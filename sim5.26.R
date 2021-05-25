# treatment depends on both covariates/ treatment has an effect on just one of them
# Y is an end of follow up variable

# so, each unit has variabes: (X10, X20, A1, X11, X21, A2, Y)
# compute outcome distribution under: always treat/ never treat
# using various assumed causal graphs
# G0: treatment has no effect on TVC
# G1: treatment effects X1 only
# G2: treatment effects X2
# G3: treatment effects both X1, X2
 
# 1. simulate data for this scenario (where true G is s.t. A --> X1, not X2)
# 2. compute the 'true' effect under this scenario
# 3. confirm that gfoRmula pkg correctly estimates this

n = 500
u = rnorm(n, mean = 0, sd = 2)
Z1 = rnorm(n)
X1 = rnorm(n)

A1 = rbinom(n, 1, prob = plogis(Z1 + X1)) # higher probability of treatment when X's are >

Z2 = rnorm(n, mean = 0.8 * Z1 + u)
X2 = rnorm(n, mean = 0.7 * X1 - A1 + u) # responds to treatment

A2 = rbinom(n, 1, prob = plogis(Z2 + X2)) # higher probability of treatment when X's are >

Y = rnorm(n, mean = 0.7 * X2 - A2 + u)

dat = data.frame(Z1, Z2, X1, X2, A1, A2, Y)

dat = dat %>% pivot_longer(!Y, names_to = c(".value", "set"), names_pattern = "(.)(.)")

dat$id = rep(1:n, each = 2) # to-do.. this automatically (using regex)
dat$set = as.numeric(dat$set) - 1
head(dat)

# what is the true effect? - .7, - 1 (-1.7)
# use gfoRmula package
id = 'id'
time_name = 'set'
covnames = c('X', 'Z', 'A')
outcome_name = 'Y'
covtypes = c('normal', 'normal', 'binary')
histories = c(lagged)
histvars = list(c('X', 'Z', 'A'))
covparams = list(covmodels=c(X ~ lag1_X + lag1_A,
                             Z ~ lag1_Z,
                             A ~ X + Z))
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


############# Bayesian version: true G ##################


bd = dat %>% group_by(id) %>% mutate(lag_A = lag(A), lag_X = lag(X), lag_Z = lag(Z)) %>% na.omit
xmod = stan_glm(X ~ lag_A + lag_X, family = gaussian(), data = bd)
ymod = stan_glm(Y ~ A + X, family = gaussian(), data = bd)

# Y ^ 11
nrep = 30
pp.X21 = posterior_predict(xmod, newdata = data.frame(lag_A = 1, lag_X = bd$lag_X), draws = nrep)
pp.X20 = posterior_predict(xmod, newdata = data.frame(lag_A = 0, lag_X = bd$lag_X), draws = nrep)

# these will be 500 * nrep long
pp.Y11 = posterior_predict(ymod, newdata = data.frame(A = 1, X = array(t(pp.X21))), draws = nrep)
pp.Y00 = posterior_predict(ymod, newdata = data.frame(A = 0, X = array(t(pp.X20))), draws = nrep)


mean(pp.Y11 - pp.Y00)

############# Bayesian version: wrong G ##################


bd = dat %>% group_by(id) %>% mutate(lag_A = lag(A), lag_X = lag(X), lag_Z = lag(Z)) %>% na.omit
xmod = stan_glm(X ~ lag_X, family = gaussian(), data = bd)
ymod = stan_glm(Y ~ A + X, family = gaussian(), data = bd)

# Y ^ 11
nrep = 30
pp.X21 = posterior_predict(xmod, newdata = data.frame(lag_X = bd$lag_X), draws = nrep)
pp.X20 = posterior_predict(xmod, newdata = data.frame(lag_X = bd$lag_X), draws = nrep)

# these will be 500 * nrep long
pp.Y11 = posterior_predict(ymod, newdata = data.frame(A = 1, X = array(t(pp.X21))), draws = nrep)
pp.Y00 = posterior_predict(ymod, newdata = data.frame(A = 0, X = array(t(pp.X20))), draws = nrep)


mean(pp.Y11 - pp.Y00)


# draw graphs
# plot posterior predictive in each


