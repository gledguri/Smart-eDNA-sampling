# eDNA survey planner -----------------------------------------------------------
# Standalone prototype: source this file in RStudio. Editable scenarios and
# executable checks are at the VERY END. No downloads or files are written.
# To load functions only: options(edna.run_examples = FALSE); source("...")
# Dependencies: install.packages(c("dplyr", "ggplot2", "tidyr", "stringr"))
#
# Main functions:
#   calibrate_survey() : manuscript illustration, supplied GP parameters, or pilot
#   plan_samples()    : epsilon + area -> samples PER STRATUM
#   plan_area()       : samples PER STRATUM + epsilon -> conditional area
# Methods: predict(calibration, ...), plot(calibration), plot(plan), print(...).
#
# Units: coordinates and rho in km; area in km^2; alpha/sigma are STANDARD
# DEVIATIONS on the supplied response scale. Nothing is automatically logged.
# One sample = one observation at a distinct site within one stratum.
# Plan different strata separately. This does not optimise sampling locations.
#
# Two deliberately distinct error targets:
# * manuscript: RMSE versus a fitted reference map; rounded published regression
#   coefficients. Exploratory ONLY: archived sigma/count provenance unresolved.
# * GP simulation: RMSE versus known latent truth over a supplied grid; refits
#   by maximum likelihood (not the manuscript's Bayesian Stan fitting).
#   Results condition on the chosen/estimated generating parameters. Simulation
#   variability is NOT a confidence interval or a target-attainment guarantee.
# All area conversions assume the density relationship transfers to the new
# extent/geometry. They do not identify a geographic boundary of good inference.
#
# Local design brief: eDNA_sampling_package_brainstorm.md
# GP model: https://mc-stan.org/docs/stan-users-guide/gaussian-processes.html
# Paper relationship: log(epsilon) = omega + beta*log(E) + theta*log(rho)
#                                      + gamma*log(E)*log(rho)

.edna_packages <- c("dplyr", "ggplot2", "tidyr", "stringr")
.edna_missing <- .edna_packages[!vapply(.edna_packages, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(.edna_missing)) {
  stop("Install required packages first: install.packages(c(",
       paste(sprintf('"%s"', .edna_missing), collapse = ", "), "))",
       call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

# Internal helpers -------------------------------------------------------------
.edna_positive <- function(z, name, scalar = FALSE, integer = FALSE,
                           allow_zero = FALSE) {
  if (!is.numeric(z) || !length(z) || any(!is.finite(z)) ||
      any(if (allow_zero) z < 0 else z <= 0) ||
      (scalar && length(z) != 1L) ||
      (integer && any(abs(z - round(z)) > 1e-8))) {
    stop(name, " must contain ", if (scalar) "one " else "",
         if (allow_zero) "non-negative" else "positive", " finite ",
         if (integer) "integer(s)." else "number(s).", call. = FALSE)
  }
}

.edna_seed <- function(seed, expr) {
  .edna_positive(seed, "seed", scalar = TRUE, integer = TRUE, allow_zero = TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(expr)
}

.edna_xy <- function(data, x, y, name) {
  if (!is.data.frame(data) || !all(c(x, y) %in% names(data)))
    stop(name, " must be a data frame with columns ", x, " and ", y, ".",
         call. = FALSE)
  if (x == y || !all(vapply(data[c(x, y)], is.numeric, logical(1))))
    stop(name, " needs two different numeric coordinate columns in km.",
         call. = FALSE)
  xy <- as.matrix(data[c(x, y)])
  if (nrow(xy) < 5L || any(!is.finite(xy)) || anyDuplicated(as.data.frame(xy)))
    stop(name, " needs at least 5 finite, distinct coordinate pairs. Analyse ",
         "strata separately; do not count replicates as independent sites.",
         call. = FALSE)
  xy
}

.edna_d2 <- function(x, y = x) {
  pmax(outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y), 0)
}
.edna_cov <- function(d2, rho_km, spatial_sd) {
  spatial_sd^2 * exp(-d2 / (2 * rho_km^2))
}
.edna_solve <- function(R, b) backsolve(R, forwardsolve(t(R), b))

# Profile the constant mean using GLS; optimise free log covariance parameters.
# Bounded multistart MLE keeps the standalone script independent of Stan.
.edna_fit_gp <- function(xy, response, fixed = list()) {
  d2 <- .edna_d2(xy)
  distances <- sqrt(d2[upper.tri(d2)])
  sy <- stats::sd(response)
  if (!is.finite(sy) || sy < 1e-10)
    stop("The response has no usable variation for fitting a spatial GP.",
         call. = FALSE)
  lower <- c(rho_km = min(distances[distances > 0]) / 20,
             spatial_sd = sy * 1e-3, noise_sd = sy * 1e-3)
  upper <- c(rho_km = max(distances) * 10,
             spatial_sd = sy * 10, noise_sd = sy * 10)
  fixed <- unlist(fixed)
  free <- setdiff(names(lower), names(fixed))
  assemble <- function(par) {
    out <- setNames(numeric(3), names(lower))
    out[names(fixed)] <- fixed
    out[free] <- exp(par)
    out
  }
  evaluate <- function(par, details = FALSE) {
    p <- assemble(par)
    K <- .edna_cov(d2, p[["rho_km"]], p[["spatial_sd"]])
    diag(K) <- diag(K) + p[["noise_sd"]]^2 +
      1e-9 * (p[["spatial_sd"]]^2 + p[["noise_sd"]]^2)
    R <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(R)) return(1e100)
    inv_one <- .edna_solve(R, rep(1, nrow(xy)))
    mu <- sum(inv_one * response) / sum(inv_one)
    residual <- response - mu
    value <- sum(log(diag(R))) +
      sum(residual * .edna_solve(R, residual)) / 2 +
      length(response) * log(2 * pi) / 2
    if (details) list(parameters = p, mean = mu, R = R, value = value)
    else value
  }
  if (!length(free)) {
    result <- evaluate(numeric(), details = TRUE)
    result$at_boundary <- FALSE
    return(result)
  }
  starts <- lapply(c(0.1, 0.35, 0.8), function(fraction) {
    initial <- c(rho_km = max(distances) * fraction,
                 spatial_sd = sy * 0.8, noise_sd = sy * 0.5)
    initial <- pmin(pmax(initial, lower * 1.01), upper / 1.01)
    tryCatch(stats::optim(log(initial[free]), evaluate, method = "L-BFGS-B",
                         lower = log(lower[free]), upper = log(upper[free]),
                         control = list(maxit = 150)),
             error = function(e) NULL)
  })
  good <- vapply(starts, function(s) !is.null(s) && s$convergence == 0 &&
                   is.finite(s$value) && s$value < 1e99, logical(1))
  if (!any(good)) stop("GP optimisation did not converge. Check the data and ",
                       "parameter assumptions.", call. = FALSE)
  fits <- starts[good]
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  result <- evaluate(best$par, details = TRUE)
  result$at_boundary <- any(abs(best$par - log(lower[free])) < 0.01 |
                             abs(best$par - log(upper[free])) < 0.01)
  result
}

.edna_cal_check <- function(calibration) {
  if (!inherits(calibration, "edna_calibration"))
    stop("calibration must be returned by calibrate_survey().", call. = FALSE)
}

# Calibration ------------------------------------------------------------------
# Supply exactly one of:
# 1. preset = "manuscript", rho_km = ... (explicit exploratory preset)
# 2. rho_km, spatial_sd, noise_sd, domain, area_km2 (parameter simulation)
# 3. data, response, domain, area_km2 (pilot fit, then simulation).
#    Optional rho_km/spatial_sd/noise_sd FIX these parameters during pilot fit.
#
# domain: data frame of prediction AND candidate sampling locations, x/y in km.
# Optional domain$weight gives positive area weights; equal-area cells otherwise.
# Samples are selected uniformly among candidate locations without replacement.
# Pilot observations estimate parameters; they are NOT treated as legacy samples
# in the planned survey. Simulation refits all GP covariance parameters when
# refit = TRUE. refit = FALSE assumes parameters known, often optimistically.
#
# n_sim controls Monte Carlo accuracy, not statistical confidence. Use >= 100
# for substantive exploration, then assess convergence and independent accuracy.
# Dense GP calculations are capped at 500 grid points for this prototype.
calibrate_survey <- function(data = NULL, x = "x_km", y = "y_km",
                             response = "log_edna", domain = NULL,
                             area_km2 = NULL, rho_km = NULL,
                             spatial_sd = NULL, noise_sd = NULL,
                             preset = NULL, n_samples = c(10, 20, 40, 80),
                             n_sim = 30, seed = 42,
                             response_scale = "natural_log", refit = TRUE) {
  if (!is.character(response_scale) || length(response_scale) != 1L ||
      is.na(response_scale) || !nzchar(str_trim(response_scale)))
    stop("response_scale must name the scale of your response.", call. = FALSE)
  supplied <- list(rho_km = rho_km, spatial_sd = spatial_sd, noise_sd = noise_sd)
  for (nm in names(supplied)) if (!is.null(supplied[[nm]]))
    .edna_positive(supplied[[nm]], nm, scalar = TRUE, allow_zero = nm == "noise_sd")
  if (!is.null(preset)) {
    if (!identical(preset, "manuscript"))
      stop('The only preset is "manuscript".', call. = FALSE)
    if (!is.null(data) || !is.null(domain) || !is.null(spatial_sd) ||
        !is.null(noise_sd) || !is.null(area_km2))
      stop("The manuscript preset cannot accept data, domain, area, or altered ",
           "variance parameters. Use simulation calibration for those.", call. = FALSE)
    .edna_positive(rho_km, "rho_km", scalar = TRUE)
    if (rho_km < 25 || rho_km > 1000)
      stop("The manuscript illustration supports rho_km from 25 to 1000.",
           call. = FALSE)
    coefficients <- c(omega = 1.65, beta = -0.176, theta = -0.228, gamma = -0.0684)
    a <- unname(coefficients["omega"] + coefficients["theta"] * log(rho_km))
    b <- unname(coefficients["beta"] + coefficients["gamma"] * log(rho_km))
    if (b >= -1e-8) stop("The error curve cannot be safely inverted.", call. = FALSE)
    density <- exp(seq(log(3.58), log(20), length.out = 200))
    result <- list(engine = "manuscript", rho_km = rho_km,
                   spatial_sd = NA_real_, noise_sd = NA_real_,
                   coefficients = coefficients, a = a, b = b,
                   area_km2 = 200000, response_scale = "natural_log",
                   error_target = "RMSE against a fitted reference map",
                   source = "Exploratory manuscript equation (rounded coefficients)",
                   uncertainty = "Unavailable: no coefficient covariance supplied",
                   caveat = "Sigma/count provenance unresolved; not validated survey advice.",
                   curve = tibble(density = density, central = exp(a + b * log(density)),
                                  lower = NA_real_, upper = NA_real_))
    message("Using the exploratory manuscript equation; variance assumptions ",
            "and calibration provenance are unresolved.")
  } else {
    xy_grid <- .edna_xy(domain, x, y, "domain")
    if (nrow(xy_grid) > 500)
      stop("Use <= 500 domain cells for this dense-GP prototype. Check grid ",
           "resolution sensitivity before substantive use.", call. = FALSE)
    .edna_positive(area_km2, "area_km2", scalar = TRUE)
    .edna_positive(n_samples, "n_samples", integer = TRUE)
    .edna_positive(n_sim, "n_sim", scalar = TRUE, integer = TRUE)
    if (!is.logical(refit) || length(refit) != 1L || is.na(refit))
      stop("refit must be TRUE or FALSE.", call. = FALSE)
    sizes <- sort(unique(n_samples))
    if (length(sizes) < 3 || min(sizes) < 5 || max(sizes) > nrow(xy_grid))
      stop("Use at least 3 distinct sample sizes, all >= 5 and <= domain rows.",
           call. = FALSE)
    if (n_sim < 10) stop("Use n_sim >= 10 (prefer >= 100 for substantive work).",
                        call. = FALSE)
    weights <- if ("weight" %in% names(domain)) domain$weight else rep(1, nrow(domain))
    .edna_positive(weights, "domain$weight")
    weights <- weights / sum(weights)
    pilot_fit <- NULL
    fixed <- supplied[!vapply(supplied, is.null, logical(1))]
    if (!is.null(data)) {
      xy <- .edna_xy(data, x, y, "data")
      if (nrow(xy) > 500) stop("Use <= 500 pilot sites for this prototype.", call. = FALSE)
      if (!response %in% names(data) || !is.numeric(data[[response]]) ||
          any(!is.finite(data[[response]])))
        stop("response must name a finite numeric column. Transform and handle ",
             "zeros/detection limits explicitly before calibration.", call. = FALSE)
      pilot_fit <- .edna_fit_gp(xy, data[[response]], fixed)
      pars <- pilot_fit$parameters
      mu <- pilot_fit$mean
      if (pilot_fit$at_boundary)
        warning("Pilot covariance fit reached an optimisation bound; inspect ",
                "identifiability before using this calibration.", call. = FALSE)
    } else {
      if (length(fixed) != 3L)
        stop("Without pilot data, supply rho_km, spatial_sd, AND noise_sd.", call. = FALSE)
      pars <- unlist(fixed)
      mu <- 0
    }
    if (n_sim < 100) message("Small n_sim: quick exploration only; inspect Monte Carlo stability.")
    d2 <- .edna_d2(xy_grid)
    K <- .edna_cov(d2, pars[["rho_km"]], pars[["spatial_sd"]])
    R_latent <- chol(K + diag(1e-9 * pars[["spatial_sd"]]^2, nrow(K)))
    runs <- .edna_seed(seed, {
      rows <- vector("list", n_sim * length(sizes))
      k <- 0L
      for (iteration in seq_len(n_sim)) {
        truth <- mu + as.vector(t(R_latent) %*% rnorm(nrow(K)))
        observed <- truth + rnorm(nrow(K), sd = pars[["noise_sd"]])
        permutation <- sample.int(nrow(K))
        for (n in sizes) {
          idx <- permutation[seq_len(n)]
          fitted <- tryCatch(.edna_fit_gp(xy_grid[idx, , drop = FALSE], observed[idx],
                                         if (refit) list() else as.list(pars)),
                             error = function(e) stop("Simulation ", iteration,
                               ", n = ", n, ": ", conditionMessage(e), call. = FALSE))
          p <- fitted$parameters
          cross <- .edna_cov(d2[, idx, drop = FALSE], p[["rho_km"]], p[["spatial_sd"]])
          prediction <- fitted$mean + as.vector(cross %*%
            .edna_solve(fitted$R, observed[idx] - fitted$mean))
          k <- k + 1L
          rows[[k]] <- tibble(iteration = iteration, n_per_stratum = n,
                              density = 10000 * n / area_km2,
                              epsilon = sqrt(sum(weights * (prediction - truth)^2)),
                              boundary_fit = fitted$at_boundary)
        }
      }
      bind_rows(rows)
    })
    if (any(!is.finite(runs$epsilon) | runs$epsilon <= 0))
      stop("Simulation produced invalid RMSE values.", call. = FALSE)
    curve <- runs %>%
      group_by(density) %>%
      summarise(central_raw = exp(mean(log(epsilon))),
                lower_raw = as.numeric(quantile(epsilon, 0.1)),
                upper_raw = as.numeric(quantile(epsilon, 0.9)), .groups = "drop") %>%
      arrange(density)
    # Log-scale isotonic smoothing: an explicitly non-increasing planning curve.
    # Retain unsmoothed errors/runs so users can inspect Monte Carlo artefacts.
    decreasing <- function(z) exp(-stats::isoreg(log(curve$density), -log(z))$yf)
    curve <- curve %>% mutate(central = decreasing(central_raw),
                              lower = decreasing(lower_raw), upper = decreasing(upper_raw))
    boundary_rate <- mean(runs$boundary_fit)
    if (boundary_rate > 0.1)
      warning(round(100 * boundary_rate), "% of simulation fits reached an optimisation ",
              "bound. Small samples may not identify covariance parameters.", call. = FALSE)
    result <- list(engine = "simulation", rho_km = unname(pars["rho_km"]),
                   spatial_sd = unname(pars["spatial_sd"]),
                   noise_sd = unname(pars["noise_sd"]), area_km2 = area_km2,
                   response_scale = response_scale,
                   error_target = "Area-weighted RMSE against simulated latent truth",
                   source = if (is.null(data)) "GP simulation from supplied parameters"
                            else "GP simulation from pilot MLE parameters",
                   uncertainty = "10th-90th percentiles across simulated surveys, conditional on generating parameters",
                   caveat = "No pilot-parameter uncertainty or independent validation; area transfer is conditional.",
                   curve = curve, runs = runs, pilot_fit = pilot_fit,
                   diagnostics = list(boundary_fit_fraction = boundary_rate,
                                      pilot_boundary = !is.null(pilot_fit) && pilot_fit$at_boundary),
                   n_sim = n_sim, seed = seed, refit = refit,
                   domain = domain, grid_weights = weights,
                   parameters_fixed_in_pilot = names(fixed))
  }
  class(result) <- "edna_calibration"
  result
}

print.edna_calibration <- function(x, ...) {
  cat(x$source, "\n", x$error_target, " [", x$response_scale, "]\n", sep = "")
  print(tibble(rho_km = x$rho_km, spatial_sd = x$spatial_sd, noise_sd = x$noise_sd,
               reference_area_km2 = x$area_km2,
               min_density = min(x$curve$density), max_density = max(x$curve$density)))
  cat(x$uncertainty, "\n", x$caveat, "\n", sep = "")
  invisible(x)
}

.edna_curve_at <- function(calibration, density, column = "central") {
  if (calibration$engine == "manuscript") {
    if (column != "central") return(rep(NA_real_, length(density)))
    return(exp(calibration$a + calibration$b * log(density)))
  }
  exp(stats::approx(log(calibration$curve$density), log(calibration$curve[[column]]),
                    xout = log(density), rule = 1)$y)
}

.edna_required <- function(calibration, epsilon, criterion) {
  .edna_cal_check(calibration)
  .edna_positive(epsilon, "epsilon")
  if (criterion == "upper" && calibration$engine == "manuscript")
    stop("The manuscript preset has no uncertainty estimates. Use central, or ",
         "a simulation calibration for an upper-quantile scenario.", call. = FALSE)
  lo <- min(calibration$curve$density)
  hi <- max(calibration$curve$density)
  column <- if (criterion == "upper") "upper" else "central"
  bind_rows(lapply(epsilon, function(target) {
    minimum_error <- .edna_curve_at(calibration, hi, column)
    low_density_error <- .edna_curve_at(calibration, lo, column)
    if (target < minimum_error * (1 - 1e-10)) {
      density <- NA_real_
      status <- "target_beyond_calibrated_density"
    } else if (target >= low_density_error) {
      density <- lo
      status <- "lowest_calibrated_density_suffices"
    } else {
      # Leftmost density meeting the target, including isotonic flat segments.
      # Bisection is valid for both the analytic curve and monotone interpolation.
      left <- log(lo); right <- log(hi)
      for (i in seq_len(60)) {
        middle <- (left + right) / 2
        if (.edna_curve_at(calibration, exp(middle), column) <= target)
          right <- middle else left <- middle
      }
      density <- exp(right)
      status <- "within_density_calibration"
    }
    tibble(epsilon_target = target, required_density = density, status = status)
  }))
}

.edna_plan <- function(table, calibration, kind, criterion) {
  table <- table %>% mutate(criterion = criterion, response_scale = calibration$response_scale,
                            error_target = calibration$error_target)
  class(table) <- c("edna_plan", class(table))
  attr(table, "calibration") <- calibration
  attr(table, "kind") <- kind
  table
}

# Target error -> sampling effort -----------------------------------------------
# Vectors of epsilon and area produce ALL combinations (not pairwise recycling).
# criterion = "central": geometric mean RMSE curve.
# criterion = "upper": 90th percentile simulation curve, NOT validated assurance.
# Counts refer to ONE stratum. A rounded count outside density support is NA.
plan_samples <- function(calibration, epsilon, area_km2 = 10000,
                          criterion = c("central", "upper")) {
  criterion <- match.arg(criterion)
  .edna_positive(area_km2, "area_km2")
  requirements <- .edna_required(calibration, epsilon, criterion)
  table <- tidyr::expand_grid(requirements, area_km2 = unique(area_km2)) %>%
    mutate(n_per_stratum = pmax(1, ceiling(area_km2 * required_density / 10000 - 1e-10)),
           realised_density = 10000 * n_per_stratum / area_km2,
           rounded_outside = !is.na(realised_density) &
             realised_density > max(calibration$curve$density) * (1 + 1e-9),
           status = if_else(rounded_outside, "integer_count_outside_density_support", status),
           n_per_stratum = if_else(rounded_outside, NA_real_, n_per_stratum),
           realised_density = if_else(rounded_outside, NA_real_, realised_density),
           epsilon_central = .edna_curve_at(calibration, realised_density),
           epsilon_q10 = .edna_curve_at(calibration, realised_density, "lower"),
           epsilon_q90 = .edna_curve_at(calibration, realised_density, "upper"),
           area_scaling = if_else(abs(area_km2 - calibration$area_km2) < 1e-8,
                                 "reference_area; geometry_must_match",
                                 "conditional_transfer_to_different_area")) %>%
    select(-rounded_outside)
  .edna_plan(table, calibration, "samples", criterion)
}

# Available effort -> density-equivalent area ----------------------------------
# This reports a conditional size, NOT the map of a trustworthy region.
plan_area <- function(calibration, n_per_stratum, epsilon,
                       criterion = c("central", "upper")) {
  criterion <- match.arg(criterion)
  .edna_positive(n_per_stratum, "n_per_stratum", integer = TRUE)
  requirements <- .edna_required(calibration, epsilon, criterion)
  table <- tidyr::expand_grid(requirements, n_per_stratum = unique(n_per_stratum)) %>%
    mutate(area_km2 = 10000 * n_per_stratum / required_density,
           area_interpretation = "conditional density-equivalent area; geometry not validated")
  .edna_plan(table, calibration, "area", criterion)
}

# Existing effort and area -> predicted error -----------------------------------
predict.edna_calibration <- function(object, n_per_stratum, area_km2 = 10000, ...) {
  .edna_cal_check(object)
  .edna_positive(n_per_stratum, "n_per_stratum", integer = TRUE)
  .edna_positive(area_km2, "area_km2")
  table <- tidyr::expand_grid(n_per_stratum = unique(n_per_stratum),
                             area_km2 = unique(area_km2)) %>%
    mutate(density = 10000 * n_per_stratum / area_km2,
           supported = density >= min(object$curve$density) * (1 - 1e-9) &
             density <= max(object$curve$density) * (1 + 1e-9),
           epsilon_central = if_else(supported, .edna_curve_at(object, density), NA_real_),
           epsilon_q10 = if_else(supported, .edna_curve_at(object, density, "lower"), NA_real_),
           epsilon_q90 = if_else(supported, .edna_curve_at(object, density, "upper"), NA_real_),
           status = if_else(supported, "within_density_calibration", "outside_density_calibration"),
           area_scaling = if_else(abs(area_km2 - object$area_km2) < 1e-8,
                                 "reference_area; geometry_must_match",
                                 "conditional_transfer_to_different_area")) %>%
    select(-supported)
  .edna_plan(table, object, "prediction", "central")
}

print.edna_plan <- function(x, ...) {
  cal <- attr(x, "calibration")
  cat(cal$source, "\n", cal$error_target, "\n", sep = "")
  cat("All counts are per stratum; area conversions assume transferable geometry.\n")
  if (cal$engine == "simulation")
    cat("q10-q90 describes simulated survey variability, not a confidence interval.\n")
  else cat(cal$caveat, "\n")
  print(tibble::as_tibble(x), ...)
  invisible(x)
}

# Plotting ---------------------------------------------------------------------
plot.edna_calibration <- function(x, ...) {
  .edna_cal_check(x)
  # Draw the SAME log-interpolated curve used by the planner, rather than
  # straight segments between sparse simulation sizes on an arithmetic axis.
  draw_density <- exp(seq(log(min(x$curve$density)), log(max(x$curve$density)),
                          length.out = 300))
  draw_density[c(1, length(draw_density))] <- range(x$curve$density)
  draw_curve <- tibble(density = draw_density,
                       central = .edna_curve_at(x, draw_density),
                       lower = .edna_curve_at(x, draw_density, "lower"),
                       upper = .edna_curve_at(x, draw_density, "upper"))
  p <- ggplot(draw_curve, aes(x = density, y = central))
  if (x$engine == "simulation") {
    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2B8CBE", alpha = 0.18) +
      geom_point(data = x$runs, aes(x = density, y = epsilon),
                 inherit.aes = FALSE, colour = "#466B7A", alpha = 0.15, size = 1)
  }
  p + geom_line(colour = "#0868AC", linewidth = 1) +
    labs(title = "Sampling effort and spatial prediction error",
         subtitle = stringr::str_wrap(paste(x$source, "| rho =", round(x$rho_km, 1), "km"), 90),
         x = "Samples per 10,000 km2 per stratum",
         y = paste0("RMSE (", x$response_scale, ")"),
         caption = stringr::str_wrap(paste(x$error_target, x$uncertainty,
                                           x$caveat, sep = "; "), 105)) +
    theme_bw(base_size = 12) + theme(plot.caption = element_text(hjust = 0))
}

plot.edna_plan <- function(x, ...) {
  cal <- attr(x, "calibration")
  p <- plot(cal)
  if (cal$engine == "simulation" && any(x$criterion == "upper"))
    p <- p + geom_line(aes(y = upper), colour = "#0868AC", linetype = "longdash")
  if (attr(x, "kind") == "prediction") {
    points <- tibble::as_tibble(x) %>% filter(is.finite(epsilon_central))
    p <- p + geom_point(data = points, aes(x = density, y = epsilon_central),
                         inherit.aes = FALSE, colour = "#D95F0E", size = 3)
  } else {
    targets <- tibble::as_tibble(x) %>% distinct(epsilon_target)
    positions <- tibble::as_tibble(x) %>% filter(is.finite(required_density)) %>%
      distinct(required_density)
    p <- p + geom_hline(data = targets, aes(yintercept = epsilon_target),
                         colour = "#D95F0E", linetype = "dashed") +
      geom_vline(data = positions, aes(xintercept = required_density),
                 colour = "#D95F0E", linetype = "dotted")
  }
  p
}

# ==============================================================================
# PLAYGROUND AND RUNNABLE CHECKS -- edit this section and run it again
# ==============================================================================
# Source the whole file once. Then change variables below and rerun scenarios.
# Set options(edna.run_examples = FALSE) before sourcing to load functions only.
# No plots/files are automatically saved. In RStudio, plots appear in Plots.
# For a console-only run, plot objects remain available after source().

RUN_EXAMPLES <- getOption("edna.run_examples", TRUE)
RUN_GP_EXAMPLES <- TRUE  # set FALSE for only the instant manuscript scenarios
SHOW_PLOTS <- interactive()

if (RUN_EXAMPLES) {
  # A. Change acceptable error, patchiness, and survey area -----------------------
  rho_to_try <- 100                   # km: manuscript illustration range 25-1000
  epsilon_to_try <- c(0.6, 0.8, 1.0)  # natural-log RMSE, NOT percentage error
  area_to_try <- 50000                # km^2 for ONE stratum
  available_samples <- c(25, 50, 100) # per stratum, not PCR replicates

  paper_cal <- calibrate_survey(preset = "manuscript", rho_km = rho_to_try)
  sample_scenarios <- plan_samples(paper_cal, epsilon_to_try, area_to_try)
  area_scenarios <- plan_area(paper_cal, available_samples, epsilon = 0.8)
  current_survey <- predict(paper_cal, n_per_stratum = available_samples,
                            area_km2 = area_to_try)
  print(sample_scenarios)
  print(area_scenarios)
  print(current_survey)
  sample_plot <- plot(sample_scenarios)
  if (SHOW_PLOTS) print(sample_plot)

  # B. Compare patchy and smooth species -----------------------------------------
  rho_scenarios <- c(50, 100, 500)
  patchiness_comparison <- bind_rows(lapply(rho_scenarios, function(rho) {
    cal <- calibrate_survey(preset = "manuscript", rho_km = rho)
    plan_samples(cal, epsilon = 0.8, area_km2 = area_to_try) %>%
      tibble::as_tibble() %>% mutate(rho_km = rho)
  })) %>% mutate(scenario = str_glue("rho = {rho_km} km"))
  print(patchiness_comparison %>%
          select(scenario, epsilon_target, required_density, n_per_stratum, status))

  if (RUN_GP_EXAMPLES) {
    # C. Supplied parameters: change sigma and see its actual effect -------------
    # A 100 x 100 km domain divided into equal-area cells; coords are centroids.
    # Change grid resolution as a sensitivity check (<= 500 total grid points).
    grid_side <- 8L
    domain_width_km <- 100
    domain_height_km <- 100
    demo_domain <- tidyr::expand_grid(
      x_km = (seq_len(grid_side) - 0.5) * domain_width_km / grid_side,
      y_km = (seq_len(grid_side) - 0.5) * domain_height_km / grid_side)
    demo_area <- domain_width_km * domain_height_km
    demo_rho <- 30
    demo_spatial_sd <- 1.0             # alpha, not alpha^2
    noise_levels <- c(0.2, 0.8)        # sigma, not sigma^2; fixed per calibration
    sample_sizes <- c(10, 20, 35, 50)
    simulation_repeats <- 12L          # quick demo; increase >= 100 for exploration
    demo_seed <- 2026

    noise_calibrations <- lapply(noise_levels, function(sigma) {
      calibrate_survey(domain = demo_domain, area_km2 = demo_area,
                       rho_km = demo_rho, spatial_sd = demo_spatial_sd,
                       noise_sd = sigma, n_samples = sample_sizes,
                       n_sim = simulation_repeats, seed = demo_seed, refit = TRUE)
    })
    noise_comparison <- bind_rows(lapply(seq_along(noise_levels), function(i) {
      noise_calibrations[[i]]$curve %>% mutate(noise_sd = noise_levels[i])
    })) %>% mutate(scenario = str_glue("sigma = {noise_sd}"))
    noise_plot <- ggplot(noise_comparison, aes(density, central, colour = scenario)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      labs(x = "Samples per 10,000 km2 per stratum", y = "Latent-field RMSE",
           colour = "Noise SD", title = "Explore changing observation noise",
           caption = "Simulation calibration; generating parameters fixed. Small demonstration run.") +
      theme_bw(base_size = 12)
    print(noise_comparison %>% select(scenario, density, central, lower, upper))
    if (SHOW_PLOTS) print(noise_plot)

    # D. Raw pilot data: synthetic example, replace pilot with YOUR data ----------
    # Your table needs x_km, y_km, log_edna (or specify other column names).
    # pilot <- read.csv("my_pilot.csv")
    # Supply the real domain and its area; do not treat lon/lat degrees as km.
    pilot <- .edna_seed(81, {
      xy <- as.matrix(demo_domain[c("x_km", "y_km")])
      K <- .edna_cov(.edna_d2(xy), demo_rho, demo_spatial_sd)
      z <- as.vector(t(chol(K + diag(1e-9, nrow(K)))) %*% rnorm(nrow(K)))
      demo_domain %>% mutate(log_edna = 2 + z + rnorm(n(), sd = 0.4)) %>%
        slice_sample(n = 50)
    })
    pilot_cal <- calibrate_survey(data = pilot, domain = demo_domain,
                                  area_km2 = demo_area, n_samples = sample_sizes,
                                  n_sim = simulation_repeats, seed = demo_seed)
    print(pilot_cal)
    pilot_epsilon <- 0.6               # edit this acceptable-error target
    pilot_plan <- plan_samples(pilot_cal, epsilon = pilot_epsilon, area_km2 = demo_area)
    upper_plan <- plan_samples(pilot_cal, epsilon = pilot_epsilon,
                               area_km2 = demo_area, criterion = "upper")
    print(pilot_plan)
    print(upper_plan)
    pilot_plot <- plot(pilot_plan)
    if (SHOW_PLOTS) print(pilot_plot)
  }

  # E. Deterministic self-checks: separate fixed fixtures from your scenarios ----
  # These check code behaviour, NOT ecological validity of the manuscript preset.
  check_cal <- calibrate_survey(preset = "manuscript", rho_km = 100)
  check_n <- plan_samples(check_cal, epsilon = 0.8, area_km2 = 50000)
  check_area <- plan_area(check_cal, n_per_stratum = 100, epsilon = 0.8)
  check_roundtrip <- predict(check_cal, n_per_stratum = 100,
                             area_km2 = check_area$area_km2)
  check_strict <- plan_samples(check_cal, epsilon = c(0.6, 0.8), area_km2 = 50000)
  check_unreachable <- plan_samples(check_cal, epsilon = 0.001)
  check_loose <- plan_samples(check_cal, epsilon = 100, area_km2 = 50000)
  stopifnot(
    abs(check_n$required_density - 5.34696083244) < 1e-8,
    check_n$n_per_stratum == 27,
    check_n$epsilon_central <= 0.8,
    abs(check_roundtrip$epsilon_central - 0.8) < 1e-8,
    check_strict$n_per_stratum[1] >= check_strict$n_per_stratum[2],
    is.na(check_unreachable$n_per_stratum),
    check_unreachable$status == "target_beyond_calibrated_density",
    check_loose$status == "lowest_calibrated_density_suffices",
    nrow(plan_samples(check_cal, c(0.7, 0.8), c(10000, 50000))) == 4,
    inherits(try(plan_samples(check_cal, epsilon = -1), silent = TRUE), "try-error"),
    inherits(try(plan_area(check_cal, n_per_stratum = 1.5, epsilon = 0.8),
                 silent = TRUE), "try-error")
  )
  if (RUN_GP_EXAMPLES) {
    stopifnot(
      inherits(pilot_cal, "edna_calibration"),
      all(is.finite(pilot_cal$runs$epsilon)),
      all(vapply(seq_along(noise_levels), function(i)
        noise_calibrations[[i]]$noise_sd == noise_levels[i], logical(1)))
    )
  }
  message("Planner self-checks passed. Edit the scenario variables above and rerun.")
}
