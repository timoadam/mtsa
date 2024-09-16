## Multivariate HMMs (slide 67)

## Plot the data
head(elephant)
summary(elephant)
par(mfrow = c(2, 2))
hist(elephant$step, xlab = "step length", ylab = "density", main = "", breaks = 30, bty = "n", prob = TRUE)
hist(elephant$angle, xlab = "turning angle", ylab = "density", main = "", breaks = 30, bty = "n", prob = TRUE)
plot(elephant$step[1:1000], xlab = "hour", ylab = "step length", type = "l", bty = "n")
plot(elephant$angle[1:1000], xlab = "hour", ylab = "turning angle", type = "l", bty = "n")

## Function that returns the negative log-likelihood of a bivariate gamma-HMM
negative_log_likelihood <- function(theta_star, y1, y2) {
  # parameter transformations for unconstrained optimisation
  Gamma = tpm(theta_star[1:2]) # applies the inverse logit link
  delta = stationary(Gamma) # stationary HMM
  mu_step = exp(theta_star[3:4]) # mu has to be positive
  sigma = exp(theta_star[5:6]) # sigma has to be positive
  mu_turn = theta_star[7:8]
  kappa = exp(theta_star[9:10]) # kappa has to be positive
  # calculate all state-dependent probabilities assuming contemporaneous conditional independence
  allprobs = matrix(1, nrow = length(y1), ncol = 2)
  for(i in 1:2){
    allprobs[, i] = dgamma2(x = y1, mu = mu_step[i], sigma = sigma[i]) * dvm(x = y2, mu = mu_turn[i], kappa = kappa[i]) # state-dependent densities
  }
  llh = forward(delta, Gamma, allprobs) # applies the forward algorithm
  return(-llh) # return negative for minimisation
}

## Numerical maximisation of the likelihood
theta_star0 <- c(qlogis(0.1), qlogis(0.1), log(2), log(4), log(1.5), log(3), 0, 0, log(1), log(2)) 
negative_log_likelihood(theta_star = theta_star0, y1 = elephant$step, y2 = elephant$angle)
mod <- nlm(negative_log_likelihood, theta_star0,  y1 = elephant$step, y2 = elephant$angle, print.level = 2)

## Extract estimates
mod
Gamma <- tpm(mod$estimate[1:2])
delta <- stationary(Gamma = Gamma)
mu_step <- exp(mod$estimate[3:4])
sigma <- exp(mod$estimate[5:6])
mu_angle <- mod$estimate[7:8]
kappa <- exp(mod$estimate[9:10])

## Plot estimated state-dependent distributions
par(mfrow = c(2, 2))
hist(elephant$step, xlab = "step length", ylab = "density", main = "", breaks = 20, bty = "n", prob = TRUE)
for(i in 1:2) {
  curve(delta[i] * dgamma2(x, mu = mu_step[i], sigma = sigma[i]), lwd = 2, col = pal[i], n = 500, add = TRUE)
}
curve(delta[1] * dgamma2(x, mu = mu_step[1], sigma = sigma[1]) + delta[2] * dgamma2(x, mu = mu_step[2], sigma = sigma[2]), lwd = 2, lty = 2, n = 500, add = TRUE)
legend("topright", col = c(pal[1], pal[2], "black"), lwd = 2, bty = "n", lty = c(1, 1, 2), legend = c("state 1", "state 2", "marginal"))

hist(elephant$angle, xlab = "turning angle", ylab = "density", main = "", breaks = 20, bty = "n", prob = TRUE)
for(i in 1:2) {
  curve(delta[i] * dvm(x, mu = mu_angle[i], kappa = kappa[i]), lwd = 2, col = pal[i], n = 500, add = TRUE)
}
curve(delta[1] * dvm(x, mu = mu_angle[1], kappa = kappa[1]) + delta[2] * dvm(x, mu = mu_angle[2], kappa = kappa[2]), lwd = 2, lty = 2, n = 500, add = TRUE)
legend("topright", col = c(pal[1], pal[2], "black"), lwd = 2, bty = "n", lty = c(1, 1, 2), legend = c("state 1", "state 2", "marginal"))

## Compute state-dependent densities
allprobs <- matrix(1, nrow = length(elephant$step), ncol = 2)
for(i in 1:2) { 
  allprobs[, i] = dgamma2(x = elephant$step, mu = mu_step[i], sigma = sigma[i]) * dvm(elephant$angle, mu = mu_angle[i], kappa = kappa[i])
}

## Global decoding
states <- viterbi(delta = delta, Gamma = Gamma, allprobs = allprobs) # Viterbi algorithm

## Plot decoded time series
plot(elephant$step[1:1000], xlab = "hour", ylab = "step length", col = pal[states[1:1000]], pch = 20, bty = "n")
legend("topright", col = c(pal[1], pal[2]), lwd = 2, bty = "n", lty = c(1, 1), legend = c("state 1", "state 2"))
plot(elephant$angle[1:1000], xlab = "hour", ylab = "turning angle", col = pal[states[1:1000]], pch = 20, bty = "n")
legend("topright", col = c(pal[1], pal[2]), lwd = 2, bty = "n", lty = c(1, 1), legend = c("state 1", "state 2"))



## Covariate-dependent state processes (slide 78)

beta <- matrix(c(-2, 0, 0, -2, 0, 0), nrow = 2, byrow = TRUE)
Gamma <- tpm_p(tod = 1:24, L = 24, beta)

## Function that returns the negative log-likelihood of an inhomogeneous, bivariate gamma-HMM
negative_log_likelihood <- function(theta_star, x, y1, y2){
  beta = matrix(theta_star[1:6], nrow = 2) # matrix of coefficients
  Gamma = tpm_p(tod = 1:24, L = 24, beta = beta, degree = 1) # calculating all L tpms
  delta = stationary_p(Gamma, t = x[1]) # periodically stationary start
  mu_step = exp(theta_star[7:8]) # mu has to be positive
  sigma = exp(theta_star[9:10]) # sigma has to be positive
  mu_turn = theta_star[11:12]
  kappa = exp(theta_star[13:14]) # kappa has to be positive
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(y1), 2)
  for(i in 1:2){ 
    allprobs[, i] = dgamma2(x = y1, mu = mu_step[i], sigma = sigma[i]) * dvm(x = y2, mu = mu_turn[i], kappa = kappa[i]) # state-dependent densities
  }
  llh <- forward_p(delta, Gamma, allprobs, x) # applies the forward algorithm
  return(-llh) # return negative for minimisation
}

## Numerical maximisation of the likelihood
theta_star0 <- c(-2, 0, 0, -2, 0, 0, # starting values state process
                 log(2), log(4), log(1.5), log(3), 0, 0, log(1), log(2)) # starting values state-dependent process

negative_log_likelihood(theta_star = theta_star0, x = elephant$tod, y1 = elephant$step, y2 = elephant$angle)
mod <- nlm(negative_log_likelihood, theta_star0, x = elephant$tod, y1 = elephant$step, y2 = elephant$angle, print.level = 2)

## Extract estimates
mod
beta <- matrix(mod$estimate[1:6], nrow = 2)
Gamma <- tpm_p(tod = 1:24, L = 24, beta = beta, degree = 1)
delta <- stationary_p(Gamma)
mu_step <- exp(mod$estimate[7:8])
sigma <- exp(mod$estimate[9:10])
mu_angle <- mod$estimate[11:12]
kappa <- exp(mod$estimate[13:14])
delta <- apply(delta, 2, mean)

## Plot estimated state-dependent distributions
par(mfrow = c(2, 2))
hist(elephant$step, xlab = "step length", ylab = "density", main = "", breaks = 20, bty = "n", prob = TRUE)
for(i in 1:2) {
  curve(delta[i] * dgamma2(x, mu = mu_step[i], sigma = sigma[i]), lwd = 2, col = pal[i], n = 500, add = TRUE)
}
curve(delta[1] * dgamma2(x, mu = mu_step[1], sigma = sigma[1]) + delta[2] * dgamma2(x, mu = mu_step[2], sigma = sigma[2]), lwd = 2, lty = 2, n = 500, add = TRUE)
legend("topright", col = c(pal[1], pal[2], "black"), lwd = 2, bty = "n", lty = c(1, 1, 2), legend = c("state 1", "state 2", "marginal"))

hist(elephant$angle, xlab = "turning angle", ylab = "density", main = "", breaks = 20, bty = "n", prob = TRUE)
for(i in 1:2) {
  curve(delta[i] * LaMa::dvm(x, mu = mu_angle[i], kappa = kappa[i]), lwd = 2, col = pal[i], n = 500, add = TRUE)
}
curve(delta[1] * LaMa::dvm(x, mu = mu_angle[1], kappa = kappa[1]) + delta[2] * LaMa::dvm(x, mu = mu_angle[2], kappa = kappa[2]), lwd = 2, lty = 2, n = 500, add = TRUE)
legend("topright", col = c(pal[1], pal[2], "black"), lwd = 2, bty = "n", lty = c(1, 1, 2), legend = c("state 1", "state 2", "marginal"))

## Compute state-dependent densities
allprobs <- matrix(1, nrow = length(elephant$step), ncol = 2)
for(i in 1:2) { 
  allprobs[, i] = dgamma2(x = elephant$step, mu = mu_step[i], sigma = sigma[i]) * LaMa::dvm(elephant$angle, mu = mu_angle[i], kappa = kappa[i])
}

## Global decoding
states <- viterbi(delta = delta, Gamma = Gamma, allprobs = allprobs) # Viterbi algorithm

## Plot decoded time series
plot(elephant$step[1:1000], xlab = "hour", ylab = "step length", col = pal[states[1:1000]], pch = 20, bty = "n")
legend("topright", col = c(pal[1], pal[2]), lwd = 2, bty = "n", lty = c(1, 1), legend = c("state 1", "state 2"))
plot(elephant$angle[1:1000], xlab = "hour", ylab = "turning angle", col = pal[states[1:1000]], pch = 20, bty = "n")
legend("topright", col = c(pal[1], pal[2]), lwd = 2, bty = "n", lty = c(1, 1), legend = c("state 1", "state 2"))

## Plot t.p.m.
par(mfrow = c(2, 2))
plot(Gamma[1, 1, ], type = "l", lwd = 2, bty = "n", xlab = "time of day", ylab = "Pr(state 1 -> 1)", ylim = c(0, 1))
plot(Gamma[1, 2, ], type = "l", lwd = 2, bty = "n", xlab = "time of day", ylab = "Pr(state 1 -> 2)", ylim = c(0, 1))
plot(Gamma[2, 1, ], type = "l", lwd = 2, bty = "n", xlab = "time of day", ylab = "Pr(state 2 -> 1)", ylim = c(0, 1))
plot(Gamma[2, 2, ], type = "l", lwd = 2, bty = "n", xlab = "time of day", ylab = "Pr(state 2 -> 2)", ylim = c(0, 1))

## Plot stationary distributions
par(mfrow = c(1, 1))
plot(delta[, 1], type = "l", lwd = 2, col = pal[1], bty = "n", xlab = "time of day", ylab = "Pr(state 1)", ylim = c(0, 1))
lines(delta[, 2], lwd = 2, col = pal[2])
