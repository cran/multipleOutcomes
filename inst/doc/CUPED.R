## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
## generating data
library(multipleOutcomes)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)

genData <- function(seed = NULL){
  set.seed(seed)
  n <- 500
  sigma <- matrix(.7, 5, 5)
  diag(sigma) <- 1.
  
  x <- rmvnorm(n, sigma = sigma)
  x[, 1] # pre-experiment data of the targeted metric
  x[, 2] # continuous metric
  x[, 3] <- log(abs(x[, 3])) # continuous and skewed metric
  x[, 4] <- ifelse(x[, 4] > .2, 1, 0) # binary metric
  a <- rbinom(n, 1, .5) # 1:1 randomization
  
  # y is post-experiment data of the targeted metric
  y <- x[, 5]
  y[a == 0] <- y[a == 0] + rnorm(sum(a == 0), sd = .1)
  y[a == 1] <- y[a == 1] + rnorm(sum(a == 1), mean = .1, sd = .1) # Delta = .1
  
  ## x5 and x6 are metrics (random noise) but uncorrelated with the targeted metric
  data.frame(
    y, a = a, 
    x1 = x[, 1], x2 = x[, 2], x3 = x[, 3], x4 = x[, 4], 
    x5 = rnorm(n), x6 = rbinom(n, 1, .5))
  
}

dat <- genData(1)
head(dat) %>% print()
xtabs(~x4, dat) %>% print()

plot(
  dat %>% 
    select(-x4, -x6, -a) %>% 
    pivot_longer(everything(), names_to = "indicator", values_to = "value") %>% 
    ggplot(aes(value)) + 
    geom_histogram() + facet_wrap(~indicator)
)

## function to compute CUPED estimate
## ...: formulas for each of the metrics
## family: a vector of family for each of the formulas
## data: a data frame of experiment data
cuped <- function(..., family, data){
  fit <- multipleOutcomes(..., family = family, data = data, data_index = NULL)
  id <- sapply(fit$id_map, function(x){x['a']})
  id1 <- id[1]
  id2 <- id[-1]
  Delta <- coef(fit)[id1]
  delta <- coef(fit)[id2]
  mcov <- vcov(fit)
  opt_c <- as.vector(- solve(mcov[id2, id2, drop = FALSE]) %*% mcov[id2, id1, drop = FALSE])
  
  estimate <- Delta + sum(delta * opt_c)
  stderr <- (
    mcov[id1, id1] + 
    2 * as.vector(t(opt_c) %*% mcov[id2, id1, drop = FALSE]) + 
    as.vector(t(opt_c) %*% mcov[id2, id2, drop = FALSE] %*% opt_c)
  ) %>% sqrt()
  
  list(estimate = estimate, stderr = stderr)
}

## -----------------------------------------------------------------------------
## An example
## y ~ a: treatment effect of post-experiment data for the targeted metric
## x1 ~ a: delta (= 0) using pre-experiment data for the targeted metric
## x2 ~ a: delta (= 0) using a continuous pre-experiment metric
## x3 ~ a: delta (= 0) using a skewed pre-experiment metric
## abs(x3) ~ a: delta (= 0) using a transformation of a pre-experiment metric
## x4 ~ a: delta (= 0) using a binary pre-experiment metric
## x5 ~ a: delta (= 0) using a metric that is a random noise
## family: linear regression model is fitted for all working models (even if x4 is binary)
cuped_fit <- 
  cuped(
    y ~ a, 
    x1 ~ a, 
    x2 ~ a, 
    x3 ~ a, 
    abs(x3) ~ a, 
    x4 ~ a, 
    x5 ~ a, 
    family = 'gaussian', 
    data = dat)

## CUPED estimate and its standard error
print(cuped_fit)

## treatment effect on targeted metric (post-experiment data only)
summary(lm(y ~ a, dat))$coef


## -----------------------------------------------------------------------------
n_simu <- 1000
all_estimate <- NULL
all_stderr <- NULL
for(i in 1:n_simu){
  dat <- genData(i)
  
  # use post- data only
  post_only <- summary(lm(y ~ a, dat))$coef
  
  # add pre- data of targeted metric
  pre_target <- cuped(y ~ a, x1 ~ a, family = 'gaussian', data = dat)
  
  # add more pre- data
  pre_more <- cuped(y ~ a, x1 ~ a, x3 ~ a, x4 ~ a, family = 'gaussian', data = dat)
  
  # multiple delta (transformation) for the same metric
  # transform_delta <- cuped(y ~ a, x1 ~ a, x3 ~ a, exp(x3) ~ a, x4 ~ a, x4 ~ a, family = c('gaussian', 'gaussian', 'gaussian', 'gaussian', 'gaussian', 'binomial'), data = dat)
  
  # add all correlated pre- data
  all_corr <- cuped(y ~ a, x1 ~ a, x2 ~ a, x3 ~ a, x4 ~ a, family = 'gaussian', data = dat)
  
  # add noise data
  add_noise <- cuped(y ~ a, x1 ~ a, x2 ~ a, x3 ~ a, x4 ~ a, x5 ~ a, x6 ~ a, family = 'gaussian', data = dat)
  
  # noise only
  noise_only <- cuped(y ~ a, x5 ~ a, x6 ~ a, family = c('gaussian', 'gaussian', 'binomial'), data = dat)
  
  extractEstimate <- function(mod){
    mod$estimate
  }
  
  extractStderr <- function(mod){
    mod$stderr
  }
  
  estimate <- 
    data.frame(
      post_only = post_only['a', 'Estimate'], 
      pre_target = extractEstimate(pre_target), 
      pre_more = extractEstimate(pre_more), 
      # transform_delta = extractEstimate(transform_delta), 
      all_corr = extractEstimate(all_corr), 
      add_noise = extractEstimate(add_noise), 
      noise_only = extractEstimate(noise_only)
    )
  
  stderr <- 
    data.frame(
      post_only = post_only['a', 'Std. Error'], 
      pre_target = extractStderr(pre_target), 
      pre_more = extractStderr(pre_more), 
      # transform_delta = extractStderr(transform_delta), 
      all_corr = extractStderr(all_corr), 
      add_noise = extractStderr(add_noise), 
      noise_only = extractStderr(noise_only)
    )
  
  all_estimate <- rbind(all_estimate, estimate)
  all_stderr <- rbind(all_stderr, stderr)

}

## mean of estimates of simulation (true Delta = .1)
apply(all_estimate, 2, mean) %>% print(digits = 3)

## empirical standard deviation
apply(all_estimate, 2, sd) %>% print(digits = 3)

## theoretical standard error
apply(all_stderr, 2, mean) %>% print(digits = 3)

