# A function to generate a fictitious simulated data set.
library(dplyr)

genCensoredDataset <- function(
t_end,
t_start = 0,
replacement = TRUE,
n_units = 100,
n_samples_per_unit = 99,
beta,
eta
) {

  if (!replacement) {
    if ((t_start != 0) | !(n_samples_per_unit %in% c(1, 99))) {
      message("
        !For replacement = FALSE the start time must be zero and number 
        of samples per unit must be one. Changing t_start to zero and 
        n_samples_per_unit to one.
      ")
    }
    t_start  <- 0
    n_samples_per_unit <- 1
  }

  # Sample lifetimes from Weibull distribution and arrange in a matrix.
  sample_size <- n_units * n_samples_per_unit

  samples <- rweibull(sample_size, shape = beta, scale = eta)

  unit_lifetimes <- matrix(
    samples,
    nrow = n_units
  )

  # Calculate the replacement times.
  unit_replacement_times <- apply(
    unit_lifetimes,
    1,
    cumsum
  ) %>% t()

  # Censor the observations according to the specified parameters.
  if (!replacement) {
    censored_logical <- unit_replacement_times > t_end
    unit_replacement_times[censored_logical] <- t_end

    replacements_df <- data.frame(
      install_time = 0,
      replacement_time = unit_replacement_times,
      censoring = ifelse(censored_logical, "right", NA),
      unit = 1:n_units
    )
  } else {

    if (
      sum(unit_replacement_times[, ncol(unit_replacement_times)] < t_end) >= 1
    ) {
      message("
        !For the inputs you have selected the repeated replacements of some of 
        the units stops before t_end.
      ")
    }

    replacements_df <- lapply(
      1:n_units,
      function(i) {
        censorUnitReplacements(
          unit_replacement_times[i, ],
          t_start,
          t_end
        ) %>% mutate(unit = i)
      }
    ) %>% bind_rows()
  }

  # Return a data frame with columns for; the install time of the unit, the
  # replacement times, if the lifetime is right or interval censored, and which
  # unit.
  return(replacements_df)

}

censorUnitReplacements <- function(
  unit_replacements_vec,
  t_start,
  t_end
) {
  in_window_logical <- between(
    unit_replacements_vec,
    left = t_start,
    right = t_end
  )

  install_times  <- lag(unit_replacements_vec[in_window_logical])
  replacement_times  <- unit_replacements_vec[in_window_logical]
  install_times[1] <- t_start

  last_life_included <- tail(replacement_times, 1) ==
    tail(unit_replacements_vec, 1)

  if (!last_life_included) {
    install_times <- append(
      install_times,
      replacement_times[length(replacement_times)]
    )
    replacement_times <- append(
      replacement_times,
      t_end
    )
  }

  df <- data.frame(
    install_time = install_times,
    replacement_time = replacement_times,
    censoring = c(
      ifelse(t_start == 0, NA, "int"),
      rep(NA, (length(install_times) - 2)),
      ifelse(last_life_included, NA, "right")
    )
  )

  return(df)
}