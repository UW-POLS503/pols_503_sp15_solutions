## ----init,results='hide', echo = FALSE, message = FALSE------------------
set.seed(135125)

## ----setup, echo = TRUE, eval = FALSE------------------------------------
## knitr::opts_chunk$set(cache = TRUE)

## ----load,message = FALSE------------------------------------------------
library("ggplot2")
library("dplyr")
library("broom")
library("assertthat")
library("tidyr")
library("pols503")

## ----r-pols503,eval = FALSE----------------------------------------------
## devtools::install_github("POLS503/r-pols503")

## ----sdcor2cov-----------------------------------------------------------
sdcor2cov <- function(s, r = diag(length(s))) {
  s <- diag(s, nrow = length(s), ncol = length(s))
  s %*% r %*% s
}

## ------------------------------------------------------------------------
n <- 128
mu_X <- rep(0, 3)
s_X <- rep(1, 3)
R_X <- diag(3)
Sigma_X <- sdcor2cov(s_X, R_X)

## ------------------------------------------------------------------------
X <- MASS::mvrnorm(n, mu_X, Sigma_X, empirical = TRUE)

## ------------------------------------------------------------------------
round(cor(X), 1)
round(apply(X, 2, mean), 1)

## ------------------------------------------------------------------------
beta <- c(0, 1, 1, 1)

## ------------------------------------------------------------------------
sigma <- 1.7

## ------------------------------------------------------------------------
mu_y <- cbind(1, X) %*% beta

## ------------------------------------------------------------------------
epsilon_y <- rnorm(n, mean = 0, sd = sigma)

## ------------------------------------------------------------------------
y <- mu_y + epsilon_y

## ------------------------------------------------------------------------
mod <- lm(y ~ X)

## ------------------------------------------------------------------------
mod_df <- tidy(mod)

## ------------------------------------------------------------------------
mod_df[["std.error.robust"]] <- sqrt(diag(car::hccm(mod)))

## ------------------------------------------------------------------------
mod_df

## ----sim_lin_norm--------------------------------------------------------
sim_lin_norm <- function(iter, n, mu_X, s_X, R_X, beta, sigma) {
  # Error checking so that bugs are caught quicker :-)
  assert_that(length(s_X) == length(mu_X),
              ncol(R_X) == nrow(R_X),
              ncol(R_X) == length(mu_X),
              length(beta) == (length(mu_X) + 1))
  # Generate an X
  X <- MASS::mvrnorm(n, mu = mu_X, Sigma = sdcor2cov(s_X, R_X),
                     empirical = TRUE)
  # Create a list to stor the results
  simulations <- list()
  # Create a progress bar because we're impatient
  p <- progress_estimated(iter, min_time = 2)
  # Loop over the simulation runs
  for (j in 1:iter) {
    # Draw y
    mu <- cbind(1, X) %*% beta
    epsilon <- rnorm(n, mean = 0, sd = sigma)
    y <- mu + epsilon
    # Run a regression
    mod <- lm(y ~ X)
    # Save the coefficients in a data frame
    mod_df <- tidy(mod) %>%
      # Add a column indicating the simulation number
      mutate(.iter = j)
    # Add hetroskedasticity consistent se to the data
    mod_df[["std.error.robust"]] <- sqrt(diag(car::hccm(mod)))
    # Save these results as the next element in the storage list
    simulations[[j]] <- mod_df
    # Update the progress bar
    p$tick()$print()
  }
  # Combine the list of data frames into a single data frame
  bind_rows(simulations)
}


## ----sim0----------------------------------------------------------------
iter <- 4096
sim0 <- sim_lin_norm(iter, n, mu_X, s_X, R_X, beta, sigma)
head(sim0)

## ----summarize_sim-------------------------------------------------------
summarize_sim <- function(.data, beta) {
  ret <- .data %>%
    group_by(term) %>%
    summarize(estimate_mean = mean(estimate),
              estimate_sd = sd(estimate),
              se_mean = sqrt(mean(std.error) ^ 2),
              se_robust_mean = sqrt(mean(std.error.robust) ^ 2),
              samples = length(estimate))
  ret[["beta_true"]] <- beta
  ret
}

## ----sim0_summary--------------------------------------------------------
sim0_summary <- summarize_sim(sim0, beta)
sim0_summary

## ------------------------------------------------------------------------
select(sim0_summary, estimate_mean, beta_true)

## ------------------------------------------------------------------------
select(sim0_summary, estimate_sd)

## ------------------------------------------------------------------------
select(sim0_summary, estimate_sd, se_mean, se_robust_mean)

## ----sim_lin_norm_omitted------------------------------------------------
sim_lin_norm_omitted <- function(iter, n, mu_X, s_X, R_X, beta, sigma,
                                 omit = integer(0)) {
  assert_that(length(s_X) == length(mu_X),
              ncol(R_X) == nrow(R_X),
              ncol(R_X) == length(mu_X),
              length(beta) == (length(mu_X) + 1))
  # Generate an X
  k <- length(mu_X)
  X <- MASS::mvrnorm(n, mu_X, sdcor2cov(s_X, R_X), empirical = TRUE)
  # ------
  # NEW: ensure colnames of X are consistent despite omitting some in lm
  colnames(X) <- paste("X", 1:k, sep = "")
  # ------
  simulations <- list()
  p <- progress_estimated(iter, min_time = 2)
  for (j in 1:iter) {
    mu <- cbind(1, X) %*% beta
    epsilon <- rnorm(n, mean = 0, sd = sigma)
    y <- mu + epsilon
    # ---------
    # NEW: omit columns of X
    # Look up paste and setdiff function to see what they does
    Xomit <- as.data.frame(X)[ , setdiff(1:k, omit)]
    # ~ . means use all variables from `data` on the RHS of the formula
    mod <- lm(y ~ . , data = Xomit)
    # ---------
    mod_df <- tidy(mod) %>%
      mutate(.iter = j)
    mod_df[["std.error.robust"]] <- sqrt(diag(car::hccm(mod)))
    simulations[[j]] <- mod_df
    p$tick()$print()
  }
  # Combine the list of data frames into a single data frame
  bind_rows(simulations)
}

## ------------------------------------------------------------------------
sim_lin_norm_heterosked <- function(iter, n, mu_X, s_X, R_X, beta, sigma, gamma) {
  assert_that(length(s_X) == length(mu_X),
              ncol(R_X) == nrow(R_X),
              ncol(R_X) == length(mu_X),
              length(beta) == (length(mu_X) + 1),
              length(gamma) == (length(mu_X) + 1))
  X <- MASS::mvrnorm(n, mu_X, sdcor2cov(s_X, R_X), empirical = TRUE)
  simulations <- list()
  p <- progress_estimated(iter, min_time = 2)
  for (j in 1:iter) {
    mu <- cbind(1, X) %*% beta
    # ------------
    # NEW: variance varies by each observation
    sigma <- sqrt(exp(cbind(1, X) %*% gamma))
    # ------------
    epsilon <- rnorm(n, mean = 0, sd = sigma)
    y <- mu + epsilon
    # Run a regression
    mod <- lm(y ~ X)
    # Save the coefficients in a data frame
    # and Add a column indicating the simulation number
    mod_df <- tidy(mod) %>%
      mutate(.iter = j)
    mod_df[["std.error.robust"]] <- sqrt(diag(car::hccm(mod)))
    simulations[[j]] <- mod_df
    p$tick()$print()
  }
  bind_rows(simulations)
}

## ----sim_lin_norm_truncated----------------------------------------------
sim_lin_norm_truncated <- function(iter, n, mu_X, s_X, R_X, beta, sigma, truncation = 0.5) {
  X <- MASS::mvrnorm(n, mu_X, sdcor2cov(s_X, R_X), empirical = TRUE)
  simulations <- list()
  p <- progress_estimated(iter, min_time = 2)
  for (j in 1:iter) {
    mu <- cbind(1, X) %*% beta
    epsilon <- rnorm(n, mean = 0, sd = sigma)
    y <- mu + epsilon
    # -------
    # NEW: drop cases in which y > mean(y)
    is_obs <- (y < quantile(y, prob = truncation))
    yobs <- y[is_obs, ]
    Xobs <- X[is_obs, ]
    # -------
    # Run a regression
    mod <- lm(yobs ~ Xobs)
    # Save the coefficients in a data frame
    # and Add a column indicating the simulation number
    mod_df <- tidy(mod) %>%
      mutate(.iter = j)
    mod_df[["std.error.robust"]] <- sqrt(diag(car::hccm(mod)))
    simulations[[j]] <- mod_df
    p$tick()$print()
  }
  bind_rows(simulations)
}

## ----msim_lin_norm_ma1---------------------------------------------------
sim_lin_norm_ma1 <- function(iter, n, mu_X, s_X, R_X, beta, sigma,
                             rho = 0.5, rho_X = 0) {
  assert_that(length(s_X) == length(mu_X),
              ncol(R_X) == nrow(R_X),
              ncol(R_X) == length(mu_X),
              length(beta) == (length(mu_X) + 1),
              length(rho_X) == length(mu_X),
              all(rho >= 0 & rho <= 1),
              all(rho_X >= 0 & rho_X <= 1))
  k <- length(mu_X)
  X <- mvrnorm_ma(n, mu_X, sdcor2cov(s_X, R_X), rho = rho_X, empirical = TRUE)
  simulations <- list()
  p <- progress_estimated(iter, min_time = 2)
  for (j in 1:iter) {
    mu <- cbind(1, X) %*% beta
    # ---------
    # NEW: generate y with serially correlated errors
    lag_epsilson <- rep(0, k)
    for (i in 1:n) {
      epsilon <- rnorm(1, mean = 0, sd = sigma)
      y <- mu[i, ] + rho * lag_epsilon + epsilon
      lag_epsilon <- epsilon
    }
    # ---------
    mod <- lm(y ~ X)
    mod_df <- tidy(mod) %>%
      mutate(.iter = j)
    mod_df[["std.error.robust"]] <- sqrt(diag(car::hccm(mod)))
    simulations[[j]] <- mod_df
    p$tick()$print()
  }
  # Combine the list of data frames into a single data frame
  bind_rows(simulations)
}

