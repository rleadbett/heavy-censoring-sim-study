library(dplyr)

# likelihood function for the Weibull lifetime model
# (used for MLE in GetHyperparams())
ll_weibull <- function(
  params,
  lifetimes
) {
  beta <- params[1]
  eta <- params[2]

  ll <- sum(log(dweibull(lifetimes, shape = beta, scale = eta)))

  return(-ll)
}

# calculates beta parameters given a mean and standard deviation
# (uses the same parameterisation of the Beta distribution as R)
BetaParameters <- function(mu, sigma) {
  a <- (mu^2 * (1 - mu)) / sigma
  b <- (a / mu) - a
  return(c(a, b))
}

# aux function used to calculate Weibull params in KaminskyJointPrior()
g <- function(x) log(-log(1 - x))

KaminskyJointPrior <- function(t1, t2, t1_mean, t1_sd, t2_mean, t2_sd) {
  # calculate the parameters of the true beta priors
  beta_params_t1 <- BetaParameters(t1_mean, t1_sd)
  beta_params_t2 <- BetaParameters(t2_mean, t2_sd)

  # sample from the beta deistributions
  cdf_t1_samples <- rbeta(
    10000, beta_params_t1[1], beta_params_t1[2]
  )
  cdf_t2_samples <- rbeta(
    10000, beta_params_t2[1], beta_params_t2[2]
  )
  # apply constraint that CDF must be monotonic
  valid_samples <- cdf_t2_samples > cdf_t1_samples
  # calculate the corresponding draws of the Weibull parameters
  draws_df <- data.frame(
    cdf_t1_draw = cdf_t1_samples[valid_samples],
    cdf_t2_draw = cdf_t2_samples[valid_samples]
  ) %>%
    mutate(
      beta_draw = (g(cdf_t2_draw) - g(cdf_t1_draw)) / log(t2 / t1),
      eta_draw = exp(log(t1) - (g(cdf_t1_draw) / beta_draw))
    )
  # return a data frame of the draws for the CDf at t1 and t2 and the parameters
  return(draws_df)
}

Misspecify <- function(t1_mean_true, t1_sd, t2_mean_true, t2_sd, phi) {
  # calculate the beta distribution parameters
  beta_params_t1 <- BetaParameters(t1_mean_true, t1_sd)
  beta_params_t2 <- BetaParameters(t2_mean_true, t2_sd)
  # mispesification is defined in terms of the quantile of phi minus the median
  t1_misspec <- qbeta(phi, beta_params_t1[1], beta_params_t1[2]) -
    qbeta(0.5, beta_params_t1[1], beta_params_t1[2])
  t2_misspec <- qbeta(phi, beta_params_t2[1], beta_params_t2[2]) -
    qbeta(0.5, beta_params_t2[1], beta_params_t2[2])

  return(
    c(t1_mean_true + t1_misspec, t2_mean_true + t2_misspec)
  )
}

GetHyperparams <- function(
  beta,    # shape parameter of Weibull dist
  eta,     # scale parameter of Weibull dist
  N,       # sample size that uncertainty is calculated from
  t1,      # first exposure time that information is elicited at
  t2,      # second exposure time that information is elicited at
  phi      # misspecification(0.5 is none, <0.5 is underspecified, >0.5 is over)
) {
  # simulate 10000 Weibull lifetime datasets each with N observations,
  # estimate the Weibull parameters using MLE, and calculate the value of
  # the CDFs at t1 and t2
  N_sample_estimates <- lapply(
    1:10000,
    function(i) {
      samples <- rweibull(N, shape = beta, scale = eta)
      mel_fit <- optim(
        c(1, 1),
        fn = function(params) {
          ll_weibull(
            params,
            lifetimes = samples
          )
        }
      )
      df <- data.frame(
        beta_est = mel_fit$par[1],
        eta_est = mel_fit$par[2]
      )
      return(df)
    }
  ) %>%
    bind_rows() %>%
    mutate(
      cdf_t1 = pweibull(t1, beta_est, eta_est),
      cdf_t2 = pweibull(t2, beta_est, eta_est)
    )

  # get mean and sd of CDFs at t1 and t2
  joint_hyperparams <- data.frame(
    cdf_t1_mean = mean(N_sample_estimates$cdf_t1),
    cdf_t1_sd = sd(N_sample_estimates$cdf_t1),
    cdf_t2_mean = mean(N_sample_estimates$cdf_t2),
    cdf_t2_sd = sd(N_sample_estimates$cdf_t2)
  )
  # apply misspecification parameter
  misspecified_means <- Misspecify(
    t1_mean_true = joint_hyperparams$cdf_t1_mean,
    t1_sd = joint_hyperparams$cdf_t1_sd,
    t2_mean_true = joint_hyperparams$cdf_t2_mean,
    t2_sd = joint_hyperparams$cdf_t2_sd,
    phi = phi
  )
  joint_hyperparams$cdf_t1_mean_misspec <- misspecified_means[1]
  joint_hyperparams$cdf_t2_mean_misspec <- misspecified_means[2]
  # simulate from Kaminsky joint prior
  joint_prior_draws <- KaminskyJointPrior(
    t1 = t1,
    t2 = t2,
    t1_mean = joint_hyperparams$cdf_t1_mean_misspec,
    t1_sd = joint_hyperparams$cdf_t1_sd,
    t2_mean = joint_hyperparams$cdf_t2_mean_misspec,
    t2_sd = joint_hyperparams$cdf_t2_sd
  ) %>%
    filter(
      (eta_draw < 100) & (beta_draw < 20)
    )

  # determine the hyperparameter values of the equivalent Normal
  # marginal prior for beta and eta.
  joint_hyperparams <- joint_hyperparams %>%
    mutate(
      beta_eq_mean = mean(joint_prior_draws$beta_draw),
      beta_eq_sd = sd(joint_prior_draws$beta_draw),
      eta_eq_mean = mean(joint_prior_draws$eta_draw),
      eta_eq_sd = sd(joint_prior_draws$eta_draw),
    )

  return(joint_hyperparams)
}

GetHyperparams(
  beta = 1.1,
  eta = 1,
  N = 500,
  t1 = 1,
  t2 = 2,
  phi = 0.7
)
