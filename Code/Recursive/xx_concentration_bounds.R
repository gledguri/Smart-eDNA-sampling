# Concentration bounds from prediction error ------------------------------------
#
# GOAL
# Given an estimated eDNA concentration and the survey's prediction error,
# return the lower and upper concentration you can defend.
#
# This is the operational version of Supplementary Fig. 3 (`sf_3` in
# 03_Plots.R), which drew, for a true concentration of 100 copies/L:
#     lower = exp(log(100) - epsilon)      upper = exp(log(100) + epsilon)
#
# WHY IT WORKS THIS WAY
# The Gaussian process is fitted to LOG concentration (01_GP_sim.R:96), so the
# prediction error epsilon (an RMSE) lives on the natural-log scale. On the log
# scale a symmetric band is  log(C) +/- epsilon. Back on the concentration
# scale that becomes a multiplicative "fold" band:
#     C / exp(epsilon)   ...   C * exp(epsilon)
# epsilon behaves like a standard deviation on the log scale, so a +/-1*epsilon
# band covers ~68% of cases and +/-1.96*epsilon ~95% (the `n_sd` argument).
# This is a multiplicative error scale, NOT a guaranteed per-site interval.

library(dplyr)
library(tidyr)
library(ggplot2)


# 1. CORE FUNCTION -------------------------------------------------------------
# concentration    : estimated eDNA concentration (copies/L, ordinary scale)
# prediction_error : epsilon, the log-scale RMSE from the survey
# n_sd             : how many multiples of epsilon to span (1 = the sf_3 band)
#
# Give it single numbers or vectors; it returns every combination as a tibble.

concentration_bounds <- function(concentration, prediction_error, n_sd = 1) {
  expand_grid(
    concentration    = concentration,
    prediction_error = prediction_error,
    n_sd             = n_sd
  ) %>%
    mutate(
      fold  = exp(n_sd * prediction_error),   # multiplicative half-width
      lower = concentration / fold,
      upper = concentration * fold,
      coverage = 2 * pnorm(n_sd) - 1          # ~0.68 at n_sd = 1, ~0.95 at 1.96
    )
}


# 2. GET epsilon FROM A SURVEY DESIGN (manuscript Eq. 6) -----------------------
# If you don't already have epsilon, work it out from the sampling density and
# the species' spatial range using the manuscript regression:
#     log(epsilon) = omega + beta*log(E) + theta*log(rho) + gamma*log(E)*log(rho)
# density : samples per 10,000 km^2 per stratum (E)
# rho_km  : spatial range in km

epsilon_from_density <- function(density, rho_km) {
  omega <- 1.65; beta <- -0.176; theta <- -0.228; gamma <- -0.0684
  a <- omega + theta * log(rho_km)
  b <- beta  + gamma * log(rho_km)
  exp(a + b * log(density))
}

# Go straight from a survey design to concentration bounds.
bounds_from_design <- function(concentration, density, rho_km, n_sd = 1) {
  expand_grid(density = density, rho_km = rho_km) %>%
    mutate(prediction_error = epsilon_from_density(density, rho_km)) %>%
    rowwise() %>%
    mutate(band = list(concentration_bounds(concentration, prediction_error, n_sd))) %>%
    ungroup() %>%
    unnest(band, names_sep = "_") %>%
    transmute(density, rho_km,
              concentration = band_concentration,
              prediction_error = band_prediction_error,
              n_sd = band_n_sd,
              lower = band_lower, upper = band_upper)
}


# 3. SIMULATION WITH FAKE DATA -------------------------------------------------
# Does a +/- epsilon band really contain the truth as often as it claims?
# We fake a survey of many sites that all share one true concentration. The GP
# predicts each site's LOG concentration with a random error whose typical size
# (standard deviation) is epsilon. Then we check how often the true value falls
# inside the band drawn around each prediction.

simulate_coverage <- function(true_concentration = 100,
                              prediction_error = 0.8,
                              n_sd = 1,
                              n_sites = 5000,
                              seed = 2026) {
  set.seed(seed)

  fake_survey <- tibble(site = 1:n_sites) %>%
    mutate(
      true_log      = log(true_concentration),
      # prediction = truth + a random log-scale error of size epsilon
      predicted_log = true_log + rnorm(n_sites, mean = 0, sd = prediction_error),
      predicted     = exp(predicted_log),
      # band around each prediction
      lower = predicted / exp(n_sd * prediction_error),
      upper = predicted * exp(n_sd * prediction_error),
      covered = true_concentration >= lower & true_concentration <= upper
    )

  fake_survey %>%
    summarise(
      prediction_error   = prediction_error,
      n_sd               = n_sd,
      empirical_coverage = mean(covered),      # how often truth was inside
      nominal_coverage   = 2 * pnorm(n_sd) - 1 # what we expected
    )
}


# 4. PLOT (the sf_3 figure) ----------------------------------------------------
# Concentration bounds against sampling density, one line per spatial range.

plot_concentration_bounds <- function(concentration, density, rho_km, n_sd = 1) {
  bounds_from_design(concentration, density, rho_km, n_sd) %>%
    mutate(rho_km = factor(rho_km)) %>%
    ggplot(aes(x = density, colour = rho_km, fill = rho_km)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(aes(y = concentration), linewidth = 0.8) +
    geom_point(aes(y = concentration), size = 1.5) +
    scale_y_log10() +
    labs(
      x = expression('Sampling density (N per 10,000 km'^2*' per stratum)'),
      y = paste0("Concentration bounds around ", concentration, " copies/L"),
      colour = expression(rho*' (km)'), fill = expression(rho*' (km)'),
      title = "Defensible concentration range from prediction error"
    ) +
    theme_bw(base_size = 12)
}


# EXAMPLES (edit and rerun) ----------------------------------------------------

# A. The headline use case: concentration + error -> bounds
concentration_bounds(concentration = 100, prediction_error = 0.8)
# -> lower ~ 45, upper ~ 223 copies/L (the sf_3 band)

# Ask for a wider, ~95% band as well:
concentration_bounds(100, 0.8, n_sd = c(1, 1.96))

# B. Work epsilon out from the survey design instead of supplying it
bounds_from_design(
  concentration = 100,
  density = c(2, 5, 10, 20),   # samples per 10,000 km^2 per stratum
  rho_km  = c(50, 100, 500)
)

# C. Redraw the sf_3 figure
bounds_plot <- plot_concentration_bounds(
  concentration = 100,
  density = seq(2, 20, by = 2),
  rho_km  = c(25, 50, 100, 200, 500, 1000)
)
bounds_plot

# D. Check the band with fake data: does +/-1*epsilon cover ~68%, +/-1.96 ~95%?
bind_rows(
  simulate_coverage(prediction_error = 0.8, n_sd = 1),
  simulate_coverage(prediction_error = 0.8, n_sd = 1.96)
)
