## Simulating from an HMM (slide 25)

## Choose parameters
mu <- c(2, 4) # state-dependent means
sigma <- c(0.5, 0.8) # state-dependent standard deviations
delta <- c(0.5, 0.5) # initial state distribution
Gamma <- matrix(c(0.95, 0.05, 0.05, 0.95), nrow = 2, byrow = TRUE) # t.p.m.

## Simulate states and observations
T <- 1000 # length of the time series
s <- y <- rep(NA, times = T) # empty state and observation vectors
set.seed(123)
s[1] <- sample(x = 1:2, size = 1, prob = delta) # state at time 1
y[1] <- rnorm(n = 1, mu[s[1]], sigma[s[1]]) # observation at time 1
for(t in 2:T){
  s[t] <- sample(x = 1:2, size = 1, prob = Gamma[s[t - 1], ]) # state at time t
  y[t] <- rnorm(n = 1, mu[s[t]], sigma[s[t]]) # observation at time t
}

## Plot states and observations
pal <- c("orange","deepskyblue") # colour palette
par(mfrow = c(2, 1))
plot(y, ylab = "y", xlab = "t", col = pal[s], bty = "n", pch = 20) # observations
plot(s, ylab = "state", xlab = "t", col = pal[s], bty = "n", pch = 20) # states



## Fitting an HMM to data (slide 44)

## Install and load the LaMa package
install.packages("LaMa") # install package
library(LaMa) # load package

## Load data
data(elephant)
elephant <- elephant[-c(1, nrow(elephant)), -4]

## Get an overview of the data
dim(elephant)
head(elephant)
summary(elephant)
hist(elephant$step, xlab = "step length", ylab = "density", main = "", breaks = 20, bty = "n", prob = TRUE)
plot(elephant$step[1:1000], xlab = "hour", ylab = "step length", type = "l", bty = "n")

## Function that returns the negative log-likelihood of a 2-state gamma HMM
negative_log_likelihood <- function(theta_star, y) {
  # parameter transformations for unconstrained optimisation
  Gamma = tpm(theta_star[1:2]) # applies the inverse logit link
  delta = stationary(Gamma) # stationary HMM
  mu = exp(theta_star[3:4]) # mu has to be positive
  sigma = exp(theta_star[5:6]) # sigma has to be positive
  allprobs <- matrix(1, nrow = length(y), ncol = 2)
  for(i in 1:2){
    allprobs[, i] = dgamma2(x = y, mu = mu[i], sigma = sigma[i]) # state-dependent densities
  }
  llh = forward(delta, Gamma, allprobs) # applies the forward algorithm
  return(-llh) # return negative for minimisation
}

## Numerical maximisation of the likelihood
theta_star0 <- c(qlogis(0.2), qlogis(0.2), log(0.5), log(2), log(0.3), log(1.5)) # choose starting values
negative_log_likelihood(theta_star = theta_star0, y = elephant$step) # check the likelihood function
mod <- nlm(negative_log_likelihood, theta_star0, y = elephant$step, print.level = 2) # model fitting

## Extract estimates
mod
Gamma <- tpm(mod$estimate[1:2])
delta <- stationary(Gamma = Gamma)
mu <- exp(mod$estimate[3:4])
sigma <- exp(mod$estimate[5:6])

## Plot estimated state-dependent distributions
hist(elephant$step, xlab = "step length", ylab = "density", main = "", breaks = 20, bty = "n", prob = TRUE)
for(i in 1:2) {
  curve(delta[i] * dgamma2(x, mu = mu[i], sigma = sigma[i]), lwd = 2, col = pal[i], n = 500, add = TRUE)
}
curve(delta[1] * dgamma2(x, mu = mu[1], sigma = sigma[1]) + delta[2] * dgamma2(x, mu = mu[2], sigma = sigma[2]), lwd = 2, lty = 2, n = 500, add = TRUE)
legend("topright", col = c(pal[1], pal[2], "black"), lwd = 2, bty = "n", lty = c(1, 1, 2), legend = c("state 1", "state 2", "marginal"))



## Model selection (slide 61)

llh <- -mod$minimum # log-likelihood
n_parameters <- length(mod$estimate) # number of parameters
T <- length(elephant$step) # number of observations
-2 * llh + 2 * n_parameters # AIC
-2 * llh + log(T) * n_parameters # BIC



## Model checking (slide 61)

## Function that computes the pseudo-residuals for a 2-state gamma HMM
pseudo_residuals <- function(y, delta, Gamma, mu, sigma) {
  T <- length(y)
  log_allprobs <- lC <- la <- matrix(NA, nrow = T, ncol = 2)
  for(i in 1:2) {
    log_allprobs[, i] <- dgamma2(x = y, mu = mu[i], sigma = sigma[i], log = TRUE)
    lC[, i] <- pgamma2(y, mu = mu[i], sigma = sigma[i])
  }
  c <- max(log_allprobs[1, ])
  foo <- t(delta) * exp(log_allprobs[1, ] - c)
  sumfoo <- sum(foo)
  lscale <- log(sumfoo) + c
  la[1, ] <- lscale
  foo <- foo / sumfoo
  for (t in 2:T) {
    c <- max(log_allprobs[t, ])
    foo <- (foo %*% Gamma) * exp(log_allprobs[t, ] - c)
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo) + c
    foo <- foo / sumfoo
    la[t, ] <- foo
  }
  residuals <- rep(NA, times = length(y))
  for (t in 2:T) {
    c <- max(la[t - 1,])
    a <- exp(la[t - 1,] - c)
    res[t] <- qnorm(t(a) %*% (Gamma / sum(a)) %*% lC[t,])
  }
  return(residuals)
}

## Plot pseudo-residuals
residuals <- pseudo_residuals(y = elephant$step, delta = delta, Gamma = Gamma, mu = mu, sigma = sigma)
qqnorm(residuals, pch = 20, bty = "n", xlab = "theoretical quantiles", ylab = "sample quantiles", main = "")
abline(a = 0, b = 1)
acf(residuals[-1], xlab = "lag", ylab = "autocorrelation", main = "", bty = "n", lwd = 2, ci.col = "black")



## State decoding (slide 61)

## Compute state-dependent densities
allprobs <- matrix(1, nrow = length(elephant$step), ncol = 2)
for(j in 1:2) { 
  allprobs[, j] = dgamma2(x = elephant$step, mu = mu[j], sigma = sigma[j]) 
}

## Global decoding
states <- viterbi(delta = delta, Gamma = Gamma, allprobs = allprobs) # Viterbi algorithm

## Local decoding
state_probs <- stateprobs(delta = delta, Gamma = Gamma, allprobs = allprobs) # compute state probabilities
states <- apply(state_probs, 1, which.max) # get the states with the highest probabilities

## Plot decoded time series
plot(elephant$step[1:1000], xlab = "hour", ylab = "step length", col = pal[states[1:1000]], pch = 20, bty = "n")
legend("topright", col = c(pal[1], pal[2]), lwd = 2, bty = "n", lty = c(1, 1), legend = c("state 1", "state 2"))
