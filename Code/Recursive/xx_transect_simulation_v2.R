# Transect simulation v2 -------------------------------------------------------
#
# GOAL
#   1. Simulate a smooth 1-D "true" LATENT concentration field along a transect,
#      as a realization from a Gaussian process (inputs: n_field, rho, spatial_sd,
#      mean concentration).
#   2. OBSERVE a randomly selected subset of locations. Each observation carries
#      non-spatial noise (the nugget, sigma_sd), so observing the same location
#      twice gives two different values.
#   3. Reconstruct the field using ONLY the selected (noisy) observations.
#   4. Evaluate the reconstruction at ALL of the original locations, so every
#      location has both a TRUE value and a RECONSTRUCTED value.
#   5. Summarise the mismatch as a single RMSE across all original locations.
#
# TWO VARIANCE COMPONENTS (this is the fix from the previous version)
#   spatial_sd (alpha) : amplitude of the SPATIAL field. Controls how high/low
#                        the true peaks and troughs are. Two points at the SAME
#                        location are perfectly correlated, so alpha alone cannot
#                        make a repeated observation differ.
#   sigma_sd   (sigma) : the NUGGET / non-spatial observation noise added to each
#                        observation. THIS is what makes the same location wiggle
#                        when you observe it more than once.

library(dplyr)
library(tidyr)
library(ggplot2)

# Fixed constants.
MU_LOG      <- log(100)   # mean log-concentration (~100 copies/L)
TRANSECT_KM <- 1000       # length of the transect


# 1. SPATIAL COVARIANCE --------------------------------------------------------
# Squared-exponential. spatial_sd is the amplitude; rho_km the range.
sq_exp_cov <- function(x1, x2, rho_km, spatial_sd) {
  spatial_sd^2 * exp(-0.5 * (outer(x1, x2, "-") / rho_km)^2)
}


# 2. SIMULATE THE TRUE (LATENT) FIELD ------------------------------------------
# One draw from a Gaussian process at the locations x. This is the smooth truth,
# WITHOUT observation noise. Returns log-concentration.
simulate_true_field <- function(x, rho_km, spatial_sd, seed = 1) {
  set.seed(seed)
  K <- sq_exp_cov(x, x, rho_km, spatial_sd)
  L <- chol(K + diag(1e-8, length(x)))          # small jitter for stability
  MU_LOG + as.vector(t(L) %*% rnorm(length(x)))
}


# 3. OBSERVE WITH NON-SPATIAL NOISE (the nugget) -------------------------------
# Add independent noise of size sigma_sd to each latent value. Call it twice on
# the same location and you get two different observations.
observe <- function(latent_log, sigma_sd) {
  latent_log + rnorm(length(latent_log), mean = 0, sd = sigma_sd)
}


# 4. RECONSTRUCT FROM THE SELECTED OBSERVATIONS --------------------------------
# GP regression: the nugget sigma_sd^2 goes on the covariance diagonal, so the
# fit smooths through the noisy observations instead of interpolating them
# exactly. Evaluated at grid_x (all original locations). Returns log-concentration.
reconstruct_field <- function(sample_x, sample_log, grid_x, rho_km, spatial_sd,
                              sigma_sd) {
  K <- sq_exp_cov(sample_x, sample_x, rho_km, spatial_sd)
  diag(K) <- diag(K) + sigma_sd^2 + 1e-8        # nugget + jitter
  weights <- solve(K, sample_log - MU_LOG)
  MU_LOG + as.vector(sq_exp_cov(grid_x, sample_x, rho_km, spatial_sd) %*% weights)
}


# 5. PLOTS ---------------------------------------------------------------------

# The true field: every one of the n_field locations plotted as an orange point
# (behind the black line), with the black line drawn on top.
plot_true_field <- function(field, n_field) {
  ggplot(field, aes(x)) +
    geom_point(aes(y = observed), colour = "orange", size = 1.6) +   # behind (drawn first)
    geom_line(aes(y = true), colour = "#102a2e", linewidth = 1) + # on top
    scale_y_log10()+
    labs(x = "Distance along transect (km)", y = "Concentration (copies/L)",
         title = "True field",
         subtitle = sprintf("%d locations drawn uniformly along the transect", n_field)) +
    theme_bw(base_size = 12)
}

# True vs reconstructed at every original location, with the RMSE.
plot_reconstruction <- function(field, rmse) {
  ggplot(field, aes(x)) +
    geom_line(aes(y = true), colour = "#102a2e", linewidth = 1) +
    geom_line(aes(y = reconstructed), colour = "#d45b24", linewidth = 1,
              linetype = "dashed") +
    geom_point(data = filter(field, selected), aes(y = observed),
               colour = "#365b26", fill = "#8bc34a", shape = 21, size = 2.2) +
    scale_y_log10()+
    labs(x = "Distance along transect (km)", y = "Concentration (copies/L)",
         title = "True vs reconstructed field",
         subtitle = sprintf("Solid = true, dashed = reconstructed, dots = noisy observations. RMSE across all %d locations = %.3f (log scale)",
                            nrow(field), rmse)) +
    theme_bw(base_size = 12)
}


# =============================================================================
# PARAMETERS -- change these and rerun
# =============================================================================
n_field    <- 500          # number of locations building the TRUE latent field
rho_km     <- 250          # spatial range (km): small = patchy, large = smooth
spatial_sd <- 3            # alpha: amplitude of the spatial field
sigma_sd   <- 1          # sigma: non-spatial observation noise (the nugget)

field_seed  <- round(runif(n = 1,min=1,max=100000))           # seed for the true field
select_seed <- round(runif(n = 1,min=1,max=100000))            # seed for which locations are selected
obs_seed    <- round(runif(n = 1,min=1,max=100000))            # seed for the observation noise
# field_seed  <- 1           # seed for the true field
# select_seed <- 1            # seed for which locations are selected
# obs_seed    <- 1            # seed for the observation noise

# --- Step 1: the true latent field at n_field locations. ----------------------
# Locations are drawn uniformly at random along the transect (not a grid).
# Sorted only so the line and reconstruction plot in x-order.
set.seed(field_seed)
x_all    <- sort(runif(n_field, 0, TRANSECT_KM))
true_log <- simulate_true_field(x_all, rho_km, spatial_sd, seed = field_seed)

# --- Step 2: randomly select locations and OBSERVE them (with nugget noise). ---

set.seed(obs_seed)
observed_log <- observe(true_log, sigma_sd)

# =============================================================================
# 6. HOW MANY SAMPLES DO I NEED FOR A TARGET RMSE?  (added on top; nothing above
#    is changed) Reuses the SAME generated field and the SAME reconstruct_field().
# =============================================================================
# Idea: for a candidate number of samples, randomly thin the observations to that
# many, reconstruct, and measure the RMSE against the true field. Repeat the
# thinning a few times and average, so the answer is not tied to one lucky draw.
# Then report the smallest number of samples whose mean RMSE reaches the target.

target_rmse <- 0.8          # <-- the RMSE you want to achieve (change me)
search_reps <- 10           # random thinnings averaged per candidate sample count
search_seed <- 1            # reproducible search on this field

# Mean RMSE when we deploy `n` samples (averaged over `reps` random thinnings).
mean_rmse_for_n <- function(n, reps = search_reps) {
  rmses <- sapply(seq_len(reps), function(r) {
    idx <- sample(n_field, n)
    rec <- reconstruct_field(x_all[idx], observed_log[idx], x_all,
                             rho_km, spatial_sd, sigma_sd)
    sqrt(mean((true_log - rec)^2))
  })
  mean(rmses)
}

# Sweep a range of sample counts and build the RMSE-vs-effort curve.
set.seed(search_seed)
candidate_n <- unique(round(seq(2, n_field, length.out = 100)))
rmse_vs_effort <- tibble(n_samples = candidate_n) %>%
  mutate(mean_rmse = sapply(n_samples, mean_rmse_for_n))

# Smallest number of samples whose mean RMSE reaches (is <=) the target.
reached   <- rmse_vs_effort %>% filter(mean_rmse <= target_rmse)
n_needed  <- if (nrow(reached) > 0) min(reached$n_samples) else NA_integer_


# reconstructing ------------------------------------------------------------------------------


n_selected <- n_needed
set.seed(select_seed)
selected_idx <- sort(sample(n_field, n_selected))
sample_x <- x_all[selected_idx]

sample_log <- observed_log[selected_idx]

# --- Step 3 & 4: reconstruct from the observations, evaluate at ALL locations. -
recon_log <- reconstruct_field(sample_x, sample_log, x_all, rho_km, spatial_sd,
                               sigma_sd)

field <- tibble(
  x             = x_all,
  true          = exp(true_log),
  reconstructed = exp(recon_log),
  observed      = exp(observed_log),
  true_log      = true_log,
  recon_log     = recon_log,
  selected      = seq_len(n_field) %in% selected_idx
)


# --- Step 5: RMSE between true and reconstructed across all locations. --------
rmse <- sqrt(mean((field$true_log - field$recon_log)^2))

# --- Output. -----------------------------------------------------------------
plot_true_field(field, n_field)
plot_reconstruction(field, rmse)
# 
cat(sprintf("Selected %d of %d locations -> RMSE (log scale) = %.3f\n",
            n_selected, n_field, rmse))
