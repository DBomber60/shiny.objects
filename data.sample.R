library(tidyverse)
# sample data
#set.seed(1)
#source('~/Documents/shiny')


# start with a single iteration of simulating data, generate bayesian g formula results
# msmsim = function(outfile, from=1, to=1000, nobs=500, ncov=5, coeffa, coeffb, coeffc, effect, a0, nboot=100) {
# 2. generate parametric g-formula

ncov = 1 # number of (time-varying) covariates
a0=-1
coeffa=1; coeffb=-2; coeffc=1; coeffd = 0; effect=-0.25; su = 1.2
nobs=2000
ntot = 10000 # total number of 'true' observations generated (observed is a subset of this)
m = 2 # time periods

data.sample = function(coeffa, coeffb, coeffc, coeffd, su, effect, ntot, nobs) {

x <- array(NA, c(2, ntot, m * ncov)) # theoretical matrix/ 3 time periods
xobs <- matrix(NA, ntot, m * ncov) # observed matrix
y <- array(NA, c(2, 2, ntot)) # counter-factual response
z <- matrix(NA, ntot, 2) # treatment
xdiffsd <- 0.25
a <- rep(coeffa, ncov)
b <- rep(coeffb, ncov)
c <- rep(coeffc, ncov)
d <- rep(coeffd, ncov)
#s <- matrix(0.25, ncov, ncov)
s <- 1.0
#su <- matrix(0.05, ncov, ncov)
# su <- 0.1


# 1. intervention:

u <- rnorm(ntot, 0, su)
xobs[,1:ncov] <- x[1, , 1:ncov] <- x[2, , 1:ncov] <- rnorm(ntot, 0, s)
pz = plogis( a0 + a * xobs[,1])
#pz <- 1.0/(1.0+exp(-(a0 + xobs[,1:ncov] %*% a)))
z[,1] <- (runif(ntot) < pz)

# 2. intervention:

# untreated at time 1
x[1, ,(ncov+1):(2*ncov)] = rnorm(ntot, 0, (xdiffsd^2) * s)
x[1, ,(ncov+1):(2*ncov)] = x[1, , (ncov+1):(2*ncov)] + x[1, ,1:ncov] + u

# treated at time 1
x[2, ,(ncov+1):(2*ncov)] = rnorm(ntot, b, (xdiffsd^2) * s)
x[2, ,(ncov+1):(2*ncov)] = x[2, , (ncov+1):(2*ncov)] + x[1, ,1:ncov] + u


for (j in (ncov+1):(2*ncov)) {
  xobs[,j] <- x[cbind(z[,1]+1, 1:ntot, j)]
}

pz = plogis(a0 + 2 * z[,1] + c * xobs[,2] )
z[,2] <- (runif(ntot) < pz)
zobs = z #rowSums(z)


# Outcome:

ylong <- NULL
zlong <- NULL
counter <- 1
for (j in 0:1) {  # treatment time 1  
  for (k in 0:1) { # treatment time 2
    py = plogis(-1 + effect * (j + k) + x[(j+1), , (ncov+1):(2*ncov)] * d + u)
    y[j+1,k+1, ] <- (runif(ntot) < py)
    ylong <- c(ylong, y[cbind(j+1,k+1,1:ntot)])
    zlong <- c(zlong, rep(j+k, ntot))
    }
}
yobs <- y[cbind(z[,1]+1,z[,2]+1,1:ntot)]
table(yobs)

true.model <- glm(ylong ~ zlong + rep(u,4), family=binomial(link=logit))

#coef(truemodel)[2]
#vcov(truemodel)[2,2]

# true 'causal' parameter

#predict(truemodel, newdata = data.frame(zlong = 3), type = "response") - 
#  predict(truemodel, newdata = data.frame(zlong = 0), type = "response")

y <- y[,,1:nobs]
yobs <- yobs[1:nobs]
#z <- z[1:nobs,]
zobs <- zobs[1:nobs,]
x <- x[, 1:nobs, ]
xobs <- xobs[1:nobs,]

observed.data = data.frame(cbind(y=yobs, xobs, zobs))
names(observed.data) = c("Y", str_c("X",1:2), str_c("Z",1:2))
observed.data$Z = with(observed.data, Z1 + Z2)
observed.data$u = u[1:nobs]

# conditional on a 'high' X2, 
#summary(glm(Z1 ~ u, data = observed.data)) # unconditional
#cond = filter(observed.data, X2 > mean(X2))
#summary(glm(Z1 ~ u, data = cond)) # conditional

naive.model = glm(Y ~ Z + X1 + X2, data = observed.data, family=binomial)

# iptw
tmodel1 <- glm(Z1 ~ X1, data = observed.data, family=binomial(link=logit))
tmodel2 <- glm(Z2 ~ X2 + Z1, data = observed.data, family=binomial(link=logit))
p1 <- dbinom(observed.data$Z1, 1, prob = (predict(tmodel1, type = "response")))
p2 <- dbinom(observed.data$Z2, 1, prob = (predict(tmodel2, type = "response")))
pt <- p1 * p2

smodel1 <- glm(Z1 ~ 1, data = observed.data, family=binomial(link=logit))
smodel2 <- glm(Z2 ~ Z1, data = observed.data, family=binomial(link=logit))
sp1 <- dbinom(observed.data$Z1, 1, prob = (predict(smodel1, type = "response")))
sp2 <- dbinom(observed.data$Z2, 1, prob = (predict(smodel2, type = "response")))
sc <- sp1 * sp2

# never treat

# always treat
iptw <- 1.0/pt
iptws <- sc * iptw
#observed.data$s = with(observed.data, Z1 + Z2)

fit.iptw = (glm(Y ~ Z, data = observed.data, weights = iptws, family = quasibinomial))



return(list(obs.dat = observed.data,
            true.model = true.model,
            naive.model = naive.model,
            fit.iptw = fit.iptw,
            weights = iptws))

}


# 
data.fit = function(nsim, coeffa, coeffb, coeffc, coeffd, su, effect, ntot, nobs) {
  
  # initialize a data frame to hold results
  # bring this up if looping through
  sim = str_c("sim",1:nsim)
  sim = factor(sim, levels = sim)
  
  
  df.plot.naive = data.frame(sim = sim, mod = "naive", theta.hat = NA, se = NA, theta = NA)
  df.plot.iptw = data.frame(sim = sim, mod = "iptw", theta.hat = NA, se = NA, theta = NA)

  for(sim in 1:nsim) {
    # simulate data / fit model
    dat = data.sample(coeffa, coeffb, coeffc, coeffd, su, effect, ntot, nobs)
    
    df.plot.naive[sim,3:4] = summary(dat$naive.model)$coefficients[2,1:2]
    df.plot.iptw[sim,3:4] = summary(dat$fit.iptw)$coefficients[2,1:2]
    df.plot.naive[sim, 5] = df.plot.iptw[sim, 5] = coef(dat$true.model)[2]
    
  }
  return(df.plot=rbind(df.plot.naive, df.plot.iptw))
}

coeffa=1; coeffb=-2; coeffc=1; coeffd = 0; effect=-0.25; su = 1.2;

a = data.fit(nsim = 10, coeffa = coeffa,
             coeffb = coeffb,
             coeffc = coeffc,
             coeffd = coeffd,
             su = su,
             effect = effect,
             ntot = ntot,
             nobs = nobs)

par(mfrow=c(2,1))

ggplot(a, aes(x = sim, y = theta.hat, group = mod, color = mod)) +
  geom_errorbar(aes(ymin = theta.hat - 2 * se, ymax = theta.hat + 2 * se),
                position = position_dodge()) +
  geom_hline(yintercept=mean(a[,5])) + theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00'))

# can we see the non-causal association?






