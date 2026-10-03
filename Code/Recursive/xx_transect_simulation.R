# Illustrative transect simulation ("Hey, do you want to see what is happening?")
#
# GOAL
# Recreate Figure 1 from the eDNA survey planner, with epsilon as a single,
# honest error.
#
# THE IDEA
#   1. Build a hidden TRUE concentration field along a 1,000 km transect.
#   2. SAMPLE from that truth at random locations, each observed with NOISE_SD.
#   3. PROJECT (krige) a concentration surface back from those samples.
#   4. epsilon = RMSE between the truth and the projection, over the n_grid
#      points, on the log scale. It is one number.
#   5. Choose the NUMBER OF SAMPLES so that this RMSE equals your target epsilon.
#      The exact samples that achieve it are the dots we plot.
#
# small epsilon -> many samples -> tight projection
# large epsilon -> few samples  -> loose projection
#
# WHAT YOU CAN CHANGE (see the PARAMETERS block at the end)
#   n_grid  : number of rows building the TRUE transect (the RMSE points). 400.
#   epsilon : target error (a single RMSE). The sample count follows from it.
#   rho_km  : spatial range in km (small = patchy, large = smooth).

library(dplyr)
library(tidyr)
library(ggplot2)

# Fixed illustrative constants (as in the planner HTML).
MU_LOG      <- log(100)   # mean log-concentration (~100 copies/L)
NOISE_SD    <- 0.15       # observation noise on the log scale
TRANSECT_KM <- 1000       # length of the illustrative transect


# 1. THE HIDDEN TRUE FIELD -----------------------------------------------------
# A squared-exponential GP approximated by random cosine/sine "features".
# Returns a function giving log-concentration at any distance x (km).

make_true_field <- function(rho_km, n_features = 128, field_seed = 1) {
  set.seed(field_seed)
  w <- rnorm(n_features) / rho_km   # feature frequencies
  a <- rnorm(n_features)            # cosine weights
  b <- rnorm(n_features)            # sine weights
  function(x) {
    MU_LOG + colSums(a * cos(outer(w, x)) + b * sin(outer(w, x))) / sqrt(n_features)
  }
}

# Squared-exponential covariance between two sets of locations.
sq_exp_cov <- function(x1, x2, rho_km) {
  exp(-0.5 * (outer(x1, x2, "-") / rho_km)^2)
}


# 2. PROJECT A SURFACE FROM SOME SAMPLES ---------------------------------------
# Krige log-concentration over grid_x from the sampled locations and values.
# Returns, at every grid point:
#   est_log : the projected (posterior mean) log-concentration
#   post_sd : the projection's own uncertainty (posterior SD, log scale). This
#             is what the blue ribbon is made of -- it comes straight from the
#             samples: ~0 right at a sample, widening in the gaps between them.

project_from_samples <- function(sample_x, sample_log, grid_x, rho_km) {
  K <- sq_exp_cov(sample_x, sample_x, rho_km)
  diag(K) <- diag(K) + NOISE_SD^2 + 1e-9          # observation noise + jitter
  Kinv  <- solve(K)
  cross <- sq_exp_cov(grid_x, sample_x, rho_km)    # grid-to-sample covariance

  est_log  <- MU_LOG + as.vector(cross %*% (Kinv %*% (sample_log - MU_LOG)))
  # Posterior variance = prior variance (1) minus the part explained by samples.
  post_var <- pmax(1 - rowSums((cross %*% Kinv) * cross), 0)

  list(est_log = est_log, post_sd = sqrt(post_var))
}


# 3. PLOTS ---------------------------------------------------------------------

# Panel 1: the hidden true field with the samples drawn from it.
plot_truth <- function(grid, samples, epsilon, n_samples) {
  ggplot() +
    geom_line(data = grid, aes(x, truth), linewidth = 1, colour = "#102a2e") +
    geom_point(data = samples, aes(x, observed),
               colour = "#365b26", fill = "#8bc34a", shape = 21, size = 2.4) +
    labs(x = "Distance along illustrative transect (km)",
         y = "Concentration (copies/L)",
         title = "1. Hidden true concentration and the samples drawn from it",
         subtitle = sprintf("epsilon = %.2f  ->  N = %d samples (observed with noise SD %.2f)",
                            epsilon, n_samples, NOISE_SD)) +
    theme_bw(base_size = 12)
}

# Panel 2: the projection and the +/- epsilon band (RMSE = epsilon).
plot_reconstruction <- function(grid, samples, epsilon, n_samples, rmse) {
  ggplot(grid, aes(x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#006d77", alpha = 0.17) +
    geom_line(aes(y = truth), colour = "#102a2e", linewidth = 0.6, alpha = 0.4) +
    geom_line(aes(y = estimate), colour = "#d45b24", linewidth = 1,
              linetype = "dashed") +
    geom_point(data = samples, aes(x, observed),
               colour = "#365b26", fill = "#8bc34a", shape = 21, size = 2) +
    labs(x = "Distance along illustrative transect (km)",
         y = "Concentration (copies/L)",
         title = "2. Projection from the samples, with its uncertainty ribbon",
         subtitle = sprintf(
           "N = %d samples, RMSE(truth vs projection) = %.2f. Solid = truth, dashed = projection, ribbon = +/- 1 posterior SD from the samples",
           n_samples, rmse)) +
    theme_bw(base_size = 12)
}

# Optional: how the projection RMSE falls as the sample count grows.
plot_rmse_curve <- function(rmse_curve, epsilon, n_samples) {
  pick <- rmse_curve %>% filter(n_samples == !!n_samples)
  ggplot(rmse_curve, aes(n_samples, rmse)) +
    geom_line(colour = "#006d77", linewidth = 1) +
    geom_hline(yintercept = epsilon, linetype = "dashed", colour = "#d45b24") +
    geom_point(data = pick, aes(n_samples, rmse), colour = "#d45b24", size = 3) +
    labs(x = "Number of samples",
         y = "RMSE(truth vs projection), log scale",
         title = sprintf("epsilon = %.2f is reached at about N = %d samples",
                         epsilon, n_samples)) +
    theme_bw(base_size = 12)
}


# =============================================================================
# PARAMETERS -- change these and rerun
# =============================================================================
n_grid  <- 400             # rows building the true transect / RMSE points
epsilon <- 0.3             # target error (a single RMSE); sample count follows
rho_km  <- 100             # spatial range (km): small = patchy, large = smooth

max_samples <- 200         # largest sample count to search over
field_seed  <- 1           # seed for the hidden true field
sample_seed <- 1           # seed for where the samples land

# --- Build the true field and the grid it is evaluated on. --------------------
field     <- make_true_field(rho_km, field_seed = field_seed)
grid_x    <- seq(0, TRANSECT_KM, length.out = n_grid)
truth_log <- field(grid_x)

# --- Sample from the truth (with NOISE_SD). ----------------------------------
# One pool of samples we reuse: taking the first n of them adds samples one by
# one, so the projection RMSE falls smoothly as n grows.
set.seed(sample_seed)
sample_pool <- tibble(x = runif(max_samples, 0, TRANSECT_KM)) %>%
  mutate(log_obs  = field(x) + rnorm(n(), 0, NOISE_SD),
         observed = exp(log_obs))

# --- RMSE = truth vs projection, for each possible sample count. --------------
rmse_curve <- tibble(n_samples = 1:max_samples) %>%
  rowwise() %>%
  mutate(rmse = {
    s   <- sample_pool[1:n_samples, ]
    est <- project_from_samples(s$x, s$log_obs, grid_x, rho_km)$est_log
    sqrt(mean((truth_log - est)^2))
  }) %>%
  ungroup()

# --- Choose the sample count whose projection RMSE equals epsilon. -----------
n_samples <- rmse_curve$n_samples[which.min(abs(rmse_curve$rmse - epsilon))]

# --- The samples that achieve it, and the projection they give. --------------
samples    <- sample_pool[1:n_samples, ]
projection <- project_from_samples(samples$x, samples$log_obs, grid_x, rho_km)
est_log    <- projection$est_log
rmse       <- sqrt(mean((truth_log - est_log)^2))

# The blue ribbon is the projection's own uncertainty (+/- 1 posterior SD),
# so it comes from the samples: it pinches in at each sample and bulges in gaps.
grid <- tibble(x = grid_x, truth = exp(truth_log),
               estimate = exp(est_log),
               lower = exp(est_log - projection$post_sd),
               upper = exp(est_log + projection$post_sd))

# --- Draw the figures. -------------------------------------------------------
plot_truth(grid, samples, epsilon, n_samples)
plot_reconstruction(grid, samples, epsilon, n_samples, rmse)
plot_rmse_curve(rmse_curve, epsilon, n_samples)

cat(sprintf("epsilon = %.2f -> N = %d samples (RMSE truth vs projection = %.3f)\n",
            epsilon, n_samples, rmse))

# Tips:
#   - Patchier species need more samples: set rho_km <- 40 and rerun.
#   - Different sample locations: change sample_seed. Different truth: field_seed.
