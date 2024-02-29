library(dplyr)
library(testthat)

SimData <- function(
  beta,
  eta,
  n_units,
  t_start,
  t_end
){
  weibull_05_quant <- qweibull(0.05, shape = beta, scale = eta)

  n_lifetimes_per_unit <- ceiling(t_end / weibull_05_quant) * 2

  sim_df <- data.frame(
      unit = factor(
        rep(1:n_units, each = n_lifetimes_per_unit),
        1:n_units
      ),
      lifetime = rweibull(
        n_units * n_lifetimes_per_unit,
        shape = beta,
        scale = eta
      )
    ) %>% 
    group_by(unit) %>%
    mutate(
      failure_time = cumsum(lifetime),
      install_time = lag(failure_time),
      lifetime_id = 1:n()
    ) %>%
    ungroup() %>%
    replace(is.na(.), 0) %>%
    filter(
      between(install_time, t_start, t_end) |
      between(failure_time, t_start, t_end)
    ) %>%
    mutate(
      int_censored = !between(install_time, t_start, t_end),
      right_censored = !between(failure_time, t_start, t_end)
    )

  sim_df$install_time[sim_df$int_censored] <- t_start
  sim_df$failure_time[sim_df$right_censored] <- t_end
  sim_df$observed_lifetime <- sim_df$failure_time - sim_df$install_time

  return(sim_df)
}
