library(dplyr)
library(knitr)
library(kableExtra)
library(ggplot2)

set.seed(461)

n_units <- 3
n_lifetimes_per_unit <- 50
beta <- 1.1
eta <- 1

small_sim <- data.frame(
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
  mutate(observed_failure = between(failure_time, 5, 12.5))

png(
  file.path(".", "figures_and_tables", "simulation_method_example.png"),
  width = 500, height = 200, units = "px"
)
small_sim %>%
  ggplot() +
  geom_point(
    aes(x = failure_time, y = unit, colour = observed_failure),
    size = 2
  ) +
  scale_colour_manual(values = c("grey", "black")) +
  geom_segment(
    data = data.frame(
      x = c(5, 5, 5), y = c(3, 2, 1),
      xend = c(12.5, 12.5, 12.5), yend = c(3, 2, 1)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "black",
    linetype = 1
  ) +
  geom_segment(
    data = data.frame(
      x = c(0, 0, 0, 12.5, 12.5, 12.5), y = c(3, 2, 1, 3, 2, 1),
      xend = c(5, 5, 5, 15, 15, 15), yend = c(3, 2, 1, 3, 2, 1)
    ),
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey",
    linetype = 2
  ) +
  geom_point(
    data = data.frame(
      t = c(5, 5, 5, 12.5, 12.5, 12.5),
      unit = c(3, 2, 1, 3, 2, 1)
    ),
    aes(x = t, y = unit),
    colour = "black",
    shape = 1,
    size = 2
  ) +
  geom_vline(xintercept = c(5, 12.5), colour = "red", linetype = 2) +
  xlim(0, 15) +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate(
    "text",
    x = 5.1, y = 3.2,
    label = "start of observation",
    hjust = -0.01,
    size = 2.5,
    colour = "red"
  ) +
  annotate(
    "text",
    x = 12.6, y = 3.2,
    label = "end of observation",
    hjust = -0.01,
    size = 2.5,
    colour = "red"
  ) +
  xlab("operational time")
dev.off()

small_sim_obs <- small_sim %>%
  filter(
    between(install_time, 5, 12.5) |
    between(failure_time, 5, 12.5)
  ) %>%
  mutate(
    int_censored = !between(install_time, 5, 12.5),
    right_censored = !between(failure_time, 5, 12.5)
  )

small_sim_obs$install_time[small_sim_obs$int_censored] <- 5
small_sim_obs$failure_time[small_sim_obs$right_censored] <- 12.5
small_sim_obs$observed_lifetime <- small_sim_obs$failure_time -
  small_sim_obs$install_time

sim_table <- small_sim_obs %>%
  select(
    unit,
    true_lifetime = lifetime,
    install_time,
    failure_time,
    int_censored,
    right_censored,
    observed_lifetime
  ) %>%
  kbl() %>%
  kable_styling(bootstrap_options = "striped") %>%
  scroll_box(width = "500px", height = "200px")

readr::write_file(
  sim_table,
  file.path(".", "figures_and_tables", "simulation_example_table.html")
)
