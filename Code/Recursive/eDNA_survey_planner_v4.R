# =============================================================================
# eDNA survey planner (R version)
# -----------------------------------------------------------------------------
# Use this script in five steps. Copy the complete example at the bottom and
# replace its coordinates and parameter values with your own.
#   I.   create_survey_area(coordinates)     GPS corners -> study area
#   II.  plan_samples(area, rho, epsilon)    effort per 10,000 km2 and whole area
#   III. generate_gp_field(area, rho, sigma, alpha, mu)  full reference field
#   IV.  select_survey_samples(field, plan)  random recommended survey
#   V.   Four plot functions, one picture each, shown in RStudio.
#
# Plotting uses ggplot2. The maths matches the web tool exactly. The one difference
# is the random number generator: the web page has its own; here we use R's, so
# results are reproducible with set.seed() but a given seed will not reproduce the
# exact random field shown in the browser.
#
# A note on wording. Three ideas appear throughout:
#   * rho (km) - how quickly spatial similarity fades.
#   * epsilon - target root mean squared prediction error on the log scale.
#               Smaller epsilon = sharper maps = more samples.
#   * alpha - spatial variance; sigma - observation noise SD; mu - log mean.
#   * "sampling density" (E) - samples per 10,000 km^2 per depth layer.
# =============================================================================


# ---------------------------------------------------------------------------
# Fixed settings (these match the web tool)
# ---------------------------------------------------------------------------

# Coefficients of the fitted equation that links precision, effort and patchiness:
#   log(precision) = omega + beta*log(E) + theta*log(rho) + gamma*log(E)*log(rho)
DEFAULT_COEFFICIENTS <- list(omega = 1.65, beta = -0.176,
                             theta = -0.228, gamma = -0.0684)

# The equation was fitted within these ranges. Results outside them still
# calculate, but are flagged as extrapolations.
CALIBRATED_PATCHINESS_MIN <- 25     # km
CALIBRATED_PATCHINESS_MAX <- 1000   # km
CALIBRATED_DENSITY_MIN    <- 1.5    # samples per 10,000 km^2
CALIBRATED_DENSITY_MAX    <- 20     # samples per 10,000 km^2

EARTH_RADIUS_KM        <- 6371.0088
DEGREES_TO_RADIANS     <- pi / 180
FIELD_POINTS_PER_10000 <- 40        # simulated field points per 10,000 km^2
MAX_FIELD_POINTS       <- 4000      # cap -> largest area is 1,000,000 km^2
MIN_AREA_KM2           <- 1000      # smallest area allowed
GRID_SPACING_KM        <- 5         # evaluation grid: 5 x 5 km cells


# ---------------------------------------------------------------------------
# Part I. Create a survey area from GPS coordinates
# ---------------------------------------------------------------------------

# Convert latitude/longitude to flat x/y coordinates in kilometres, using an
# equal-area projection centred on the study region. Works on a whole matrix of
# points at once (one row = one point, columns = latitude, longitude).
project_to_km <- function(lat_lon, center_lat_lon) {
  lat    <- lat_lon[, 1] * DEGREES_TO_RADIANS
  lon    <- lat_lon[, 2]
  c_lat  <- center_lat_lon[1] * DEGREES_TO_RADIANS
  d_lon  <- (lon - center_lat_lon[2]) * DEGREES_TO_RADIANS

  denominator <- 1 + sin(c_lat) * sin(lat) + cos(c_lat) * cos(lat) * cos(d_lon)
  if (any(denominator < 1e-8)) stop("This area is too wide for a regional projection.")
  scale <- sqrt(2 / denominator)

  x <- EARTH_RADIUS_KM * scale * cos(lat) * sin(d_lon)
  y <- EARTH_RADIUS_KM * scale * (cos(c_lat) * sin(lat) -
                                  sin(c_lat) * cos(lat) * cos(d_lon))
  cbind(x = x, y = y)
}

# Convert a single flat x/y point (km) back to latitude/longitude.
km_to_lat_lon <- function(x, y, center_lat_lon) {
  distance <- sqrt(x * x + y * y)
  if (distance < 1e-12) return(center_lat_lon)
  angle  <- 2 * asin(min(1, distance / (2 * EARTH_RADIUS_KM)))
  c_lat  <- center_lat_lon[1] * DEGREES_TO_RADIANS
  lat <- asin(cos(angle) * sin(c_lat) +
              y * sin(angle) * cos(c_lat) / distance) / DEGREES_TO_RADIANS
  lon <- center_lat_lon[2] +
         atan2(x * sin(angle),
               distance * cos(c_lat) * cos(angle) -
               y * sin(c_lat) * sin(angle)) / DEGREES_TO_RADIANS
  c(lat, lon)
}

# Area of a polygon (km^2) from its x/y corners, using the shoelace formula.
polygon_area_km2 <- function(polygon_xy) {
  x <- polygon_xy[, 1]
  y <- polygon_xy[, 2]
  next_corner <- c(2:nrow(polygon_xy), 1)   # wrap back to the first corner
  abs(sum(x * y[next_corner] - x[next_corner] * y)) / 2
}

# Is a set of points inside the polygon? Classic "ray casting" test: a point is
# inside if a ray going right crosses the boundary an odd number of times.
# points_xy: matrix of points (columns x, y). polygon_xy: the boundary corners.
points_inside_polygon <- function(points_xy, polygon_xy) {
  if (is.null(dim(points_xy))) points_xy <- matrix(points_xy, ncol = 2)
  px <- points_xy[, 1]
  py <- points_xy[, 2]
  inside <- logical(length(px))

  n_corners <- nrow(polygon_xy)
  previous <- n_corners
  for (current in seq_len(n_corners)) {
    ax <- polygon_xy[current, 1]; ay <- polygon_xy[current, 2]
    bx <- polygon_xy[previous, 1]; by <- polygon_xy[previous, 2]
    crosses <- ((ay > py) != (by > py)) &
               (px < (bx - ax) * (py - ay) / (by - ay) + ax)
    crosses[is.na(crosses)] <- FALSE
    inside <- xor(inside, crosses)
    previous <- current
  }
  inside
}

# --- helpers for checking a polygon does not cross itself ---

# Orientation of the turn from a->b->c (positive = left turn, negative = right).
turn_direction <- function(a, b, c) {
  (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
}

# Does point c lie on the line segment a-b?
point_on_segment <- function(a, b, c) {
  abs(turn_direction(a, b, c)) < 1e-8 &&
    c[1] >= min(a[1], b[1]) - 1e-8 && c[1] <= max(a[1], b[1]) + 1e-8 &&
    c[2] >= min(a[2], b[2]) - 1e-8 && c[2] <= max(a[2], b[2]) + 1e-8
}

# Do the two line segments a-b and c-d cross (or touch)?
segments_cross <- function(a, b, c, d) {
  proper_crossing <-
    turn_direction(a, b, c) * turn_direction(a, b, d) < 0 &&
    turn_direction(c, d, a) * turn_direction(c, d, b) < 0
  proper_crossing ||
    point_on_segment(a, b, c) || point_on_segment(a, b, d) ||
    point_on_segment(c, d, a) || point_on_segment(c, d, b)
}

# Does the polygon boundary cross itself anywhere? (Checks every pair of edges
# that do not already share a corner.)
polygon_crosses_itself <- function(polygon_xy) {
  n <- nrow(polygon_xy)
  for (i in seq_len(n)) {
    edge1_start <- polygon_xy[i, ]
    edge1_end   <- polygon_xy[(i %% n) + 1, ]
    for (j in seq_len(n)) {
      if (j <= i) next
      shares_a_corner <- (j == i + 1) || (i == 1 && j == n)
      if (shares_a_corner) next
      edge2_start <- polygon_xy[j, ]
      edge2_end   <- polygon_xy[(j %% n) + 1, ]
      if (segments_cross(edge1_start, edge1_end, edge2_start, edge2_end)) return(TRUE)
    }
  }
  FALSE
}

# Build a study area from boundary coordinates.
#   coordinates: a matrix or data frame with two columns - latitude, longitude -
#                one corner per row, in order around the boundary, and WITHOUT
#                repeating the first corner at the end.
# Returns a list describing the area (corners, centre, projected polygon, area in
# km^2, number of field points, and the bounding box).
create_survey_area <- function(coordinates) {
  coordinates <- as.matrix(coordinates)
  storage.mode(coordinates) <- "double"
  if (ncol(coordinates) != 2) stop("Give coordinates as two columns: latitude, longitude.")

  n_corners <- nrow(coordinates)
  if (n_corners < 3 || n_corners > 200) stop("Use between 3 and 200 boundary corners.")
  if (any(!is.finite(coordinates)) || any(abs(coordinates[, 1]) > 80))
    stop("Latitudes must be between 80 S and 80 N and all values must be finite.")

  # Keep longitudes consistent even if the area crosses the +/-180 line.
  first_longitude <- coordinates[1, 2]
  wrapped_longitude <- first_longitude +
    ((coordinates[, 2] - first_longitude + 180) %% 360 + 360) %% 360 - 180
  corners_lat_lon <- cbind(coordinates[, 1], wrapped_longitude)

  center    <- c(mean(corners_lat_lon[, 1]), mean(corners_lat_lon[, 2]))
  polygon   <- project_to_km(corners_lat_lon, center)

  # Basic checks that the shape is usable.
  next_corner <- c(2:n_corners, 1)
  edge_lengths <- sqrt(rowSums((polygon - polygon[next_corner, ])^2))
  if (any(edge_lengths < 1e-6)) stop("Remove repeated or identical neighbouring corners.")
  if (max(dist(polygon)) > 2000) stop("Keep the area within a 2,000 km span.")
  if (polygon_crosses_itself(polygon)) stop("The boundary crosses itself; use a simple shape.")

  area <- polygon_area_km2(polygon)
  if (area < MIN_AREA_KM2)
    stop(sprintf("This area is only %s km2. Use an area of at least 1,000 km2.",
                 format(round(area), big.mark = ",")))

  field_points <- max(1, ceiling(area * FIELD_POINTS_PER_10000 / 10000))
  if (field_points > MAX_FIELD_POINTS)
    stop(sprintf("This area needs %s field points; the limit is %d. Use an area of at most 1,000,000 km2.",
                 format(field_points, big.mark = ","), MAX_FIELD_POINTS))

  bounding_box <- c(min(polygon[, 1]), max(polygon[, 1]),
                    min(polygon[, 2]), max(polygon[, 2]))  # x_min, x_max, y_min, y_max

  list(corners_lat_lon = corners_lat_lon, center = center, polygon = polygon,
       area_km2 = area, field_points = field_points, bounding_box = bounding_box)
}

# Older name retained for existing scripts.
build_study_area <- create_survey_area


# ---------------------------------------------------------------------------
# Part II. Plan the number of samples from rho and epsilon
# ---------------------------------------------------------------------------

# For a given patchiness, the equation becomes a straight line (on log scales):
#   log(precision) = intercept + slope * log(sampling_density)
# This returns that intercept and slope.
equation_line <- function(patchiness, coefficients = DEFAULT_COEFFICIENTS) {
  list(intercept = coefficients$omega + coefficients$theta * log(patchiness),
       slope     = coefficients$beta  + coefficients$gamma * log(patchiness))
}

# Predicted precision at a chosen sampling density.
precision_at_density <- function(sampling_density, patchiness,
                                 coefficients = DEFAULT_COEFFICIENTS) {
  line <- equation_line(patchiness, coefficients)
  exp(line$intercept + line$slope * log(sampling_density))
}

# The reverse: the sampling density needed to reach a target precision.
# Returns the density plus a note on whether it falls inside the tested range.
density_for_precision <- function(precision, patchiness,
                                  coefficients = DEFAULT_COEFFICIENTS) {
  line <- equation_line(patchiness, coefficients)
  # The line must slope downward (more samples -> better precision) to invert it.
  if (!is.finite(line$intercept) || !is.finite(line$slope) || line$slope >= -1e-8)
    return(list(density = NaN, range = "equation does not decrease"))
  density <- exp((log(precision) - line$intercept) / line$slope)
  range <- if (density > CALIBRATED_DENSITY_MAX) "above tested range"
           else if (density < CALIBRATED_DENSITY_MIN) "below tested range"
           else "within tested range"
  list(density = density, range = range)
}

# Main Part II function. The density is samples per 10,000 km2 for one depth.
# The whole-area count is rounded up to a whole sample.
plan_samples <- function(study_area, rho, epsilon,
                         coefficients = DEFAULT_COEFFICIENTS) {
  if (!is.finite(rho) || rho <= 0 || !is.finite(epsilon) || epsilon <= 0)
    stop("rho and epsilon must be positive, finite numbers.")
  needed <- density_for_precision(epsilon, rho, coefficients)
  if (!is.finite(needed$density) || needed$density <= 0)
    stop("The planning equation cannot produce a finite sample density.")
  count <- max(1, ceiling(needed$density * study_area$area_km2 / 10000 - 1e-10))
  list(rho = rho, epsilon = epsilon, area_km2 = study_area$area_km2,
       samples_per_10000_km2 = needed$density, samples_in_area = count,
       density_range = needed$range,
       is_extrapolation = needed$range != "within tested range" ||
         rho < CALIBRATED_PATCHINESS_MIN || rho > CALIBRATED_PATCHINESS_MAX)
}

# =====================  QUESTION 1: how many samples?  =====================
# Given the boundary coordinates of a study area, a species patchiness and a
# target precision, return (and print) the number of samples needed.
#   coordinates  : matrix/data frame of latitude, longitude rows
#   patchiness   : rho, in km
#   precision    : epsilon (smaller = sharper, more samples)
#   depth_layers : multiplies only the total number of samples
how_many_samples <- function(coordinates, patchiness = 100, precision = 0.8,
                             depth_layers = 1,
                             coefficients = DEFAULT_COEFFICIENTS,
                             show_summary = TRUE) {
  stopifnot(patchiness > 0, precision > 0,
            depth_layers >= 1, depth_layers == round(depth_layers))
  study_area <- build_study_area(coordinates)

  needed  <- density_for_precision(precision, patchiness, coefficients)
  density <- needed$density

  samples_per_depth <- if (is.finite(density))
    max(1, ceiling(study_area$area_km2 * density / 10000 - 1e-10)) else NA_real_
  total_samples <- samples_per_depth * depth_layers

  patchiness_ok <- patchiness >= CALIBRATED_PATCHINESS_MIN &&
                   patchiness <= CALIBRATED_PATCHINESS_MAX
  is_extrapolation <- needed$range != "within tested range" || !patchiness_ok

  answer <- list(study_area = study_area, patchiness = patchiness,
                 precision = precision, depth_layers = depth_layers,
                 sampling_density = density, samples_per_depth = samples_per_depth,
                 total_samples = total_samples, density_range = needed$range,
                 is_extrapolation = is_extrapolation)

  if (show_summary) {
    meaning <- precision_in_copies_per_litre(precision)
    cat("===================  How many samples?  ===================\n")
    cat(sprintf("Study area        : %s km2 (%s field points)\n",
                format(round(study_area$area_km2), big.mark = ","),
                format(study_area$field_points, big.mark = ",")))
    cat(sprintf("Species patchiness: %s km\n", format(patchiness, big.mark = ",")))
    cat(sprintf("Map precision     : %.3f  (around 100 copies/L -> %.1f to %.1f)\n",
                precision, meaning$lower, meaning$upper))
    cat(sprintf("Sampling density  : %.2f samples / 10,000 km2 / depth\n", density))
    cat(sprintf("Samples per depth : %s\n", format(samples_per_depth, big.mark = ",")))
    cat(sprintf("Total samples     : %s  (across %d depth layers)\n",
                format(total_samples, big.mark = ","), depth_layers))
    cat(sprintf("Tested range      : %s%s\n", needed$range,
                if (is_extrapolation) "  [EXTRAPOLATION - validate with a pilot]" else ""))
    cat("===========================================================\n")
  }
  invisible(answer)
}

# =====================  QUESTION 2: how much area?  ========================
# (Optional - this inverse calculation is not part of the manuscript.)
# Given a fixed sample budget, a patchiness and a target precision, return (and
# print) how much area that budget can cover.
how_much_area <- function(available_samples, patchiness = 100, precision = 0.8,
                          coefficients = DEFAULT_COEFFICIENTS,
                          show_summary = TRUE) {
  stopifnot(available_samples > 0, available_samples == round(available_samples),
            patchiness > 0, precision > 0)

  needed  <- density_for_precision(precision, patchiness, coefficients)
  density <- needed$density
  area_covered <- if (is.finite(density)) 10000 * available_samples / density else NA_real_

  patchiness_ok <- patchiness >= CALIBRATED_PATCHINESS_MIN &&
                   patchiness <= CALIBRATED_PATCHINESS_MAX
  is_extrapolation <- needed$range != "within tested range" || !patchiness_ok

  answer <- list(available_samples = available_samples, patchiness = patchiness,
                 precision = precision, sampling_density = density,
                 area_covered_km2 = area_covered, density_range = needed$range,
                 is_extrapolation = is_extrapolation)

  if (show_summary) {
    cat("=============  How much area can my budget cover?  ============\n")
    cat(sprintf("Available samples : %s per depth\n",
                format(available_samples, big.mark = ",")))
    cat(sprintf("Species patchiness: %s km\n", format(patchiness, big.mark = ",")))
    cat(sprintf("Map precision     : %.3f\n", precision))
    cat(sprintf("Sampling density  : %.2f samples / 10,000 km2 / depth\n", density))
    cat(sprintf("Area covered      : %s km2 per depth\n",
                format(round(area_covered), big.mark = ",")))
    cat(sprintf("Tested range      : %s%s\n", needed$range,
                if (is_extrapolation) "  [EXTRAPOLATION - validate with a pilot]" else ""))
    cat("==============================================================\n")
  }
  invisible(answer)
}


# ---------------------------------------------------------------------------
# Additional Part II calculation: what does epsilon mean in concentrations?
# ---------------------------------------------------------------------------

# Translate a precision value into a concentration range around a reference value
# (default 100 copies/L). The range is reference x/ exp(precision), the reference
# is the median, and the distribution of predictions is right-skewed (lognormal).
precision_in_copies_per_litre <- function(precision, reference = 100) {
  stopifnot(precision > 0, reference > 0)
  list(reference        = reference,
       lower            = exp(log(reference) - precision),
       median           = reference,
       upper            = exp(log(reference) + precision),
       average          = reference * exp(precision^2 / 2),
       spread_factor    = exp(precision))
}


# ---------------------------------------------------------------------------
# Part III. Generate the reference GP field from rho, sigma, alpha and mu
# ---------------------------------------------------------------------------
# The idea: build a realistic "true" concentration map, pretend we only sampled
# some of it (the recommended survey), rebuild the map from those samples, and
# measure how far the rebuilt map is from the true one.

# Distance-based correlation between locations. Nearby points are more alike;
# "amplitude" is the overall variance and "patchiness" sets how fast similarity
# fades with distance. (This is a squared-exponential / Gaussian correlation.)
correlation_matrix <- function(points_xy, amplitude, patchiness,
                               noise_sd = 0, add_noise = TRUE) {
  squared_distance <- outer(points_xy[, 1], points_xy[, 1], "-")^2 +
                      outer(points_xy[, 2], points_xy[, 2], "-")^2
  covariance <- amplitude * exp(-squared_distance / (2 * patchiness^2))
  if (add_noise) {
    tiny_stabiliser <- 1e-10 * max(amplitude + noise_sd^2, 1)
    diag(covariance) <- diag(covariance) + noise_sd^2 + tiny_stabiliser
  }
  covariance
}

# Predict concentrations at new locations from a set of sampled locations.
# Each prediction is the mean plus a weighted blend of the samples, where the
# weights come from how correlated each new point is with each sample.
predict_concentrations <- function(target_points, sample_points, weights,
                                   settings) {
  prediction <- rep(settings$mean, nrow(target_points))
  two_patch_squared <- 2 * settings$patchiness^2
  for (s in seq_len(nrow(sample_points))) {
    squared_distance <- (target_points[, 1] - sample_points[s, 1])^2 +
                        (target_points[, 2] - sample_points[s, 2])^2
    correlation <- settings$amplitude * exp(-squared_distance / two_patch_squared)
    prediction <- prediction + weights[s] * correlation
  }
  prediction
}

# Scatter a given number of points uniformly at random inside the area.
random_points_in_area <- function(study_area, how_many, max_tries = 2e6) {
  box <- study_area$bounding_box
  chosen <- matrix(NA_real_, nrow = how_many, ncol = 2)
  found <- 0L
  tries <- 0L
  while (found < how_many && tries < max_tries) {
    batch_size <- max(1L, how_many - found) * 4L
    candidates <- cbind(runif(batch_size, box[1], box[2]),
                        runif(batch_size, box[3], box[4]))
    keep <- which(points_inside_polygon(candidates, study_area$polygon))
    if (length(keep) > 0) {
      take <- min(length(keep), how_many - found)
      chosen[(found + 1):(found + take), ] <- candidates[keep[seq_len(take)], , drop = FALSE]
      found <- found + take
    }
    tries <- tries + batch_size
  }
  if (found < how_many) stop("This area is too thin to place the field points.")
  chosen
}

# Step 4a. Create the "true" reference field: scatter field points, draw a
# correlated random surface over them, and work out the weights needed to
# redraw that surface anywhere.
# settings = list(mean, amplitude, noise_sd, patchiness).
make_true_field <- function(study_area, settings, seed = 317) {
  set.seed(seed)
  points <- random_points_in_area(study_area, study_area$field_points)

  covariance <- correlation_matrix(points, settings$amplitude, settings$patchiness,
                                   settings$noise_sd, add_noise = TRUE)
  lower_triangle <- tryCatch(t(chol(covariance)),
    error = function(e) stop("The surface is numerically unstable; try gentler settings."))

  true_values <- settings$mean + as.vector(lower_triangle %*% rnorm(nrow(points)))
  if (any(!is.finite(true_values))) stop("The simulated surface overflowed; try gentler settings.")

  weights <- solve(covariance, true_values - settings$mean)
  list(points = points, true_values = true_values, weights = weights)
}

# Main Part III function. The reference GP uses 40 randomly placed observations
# per 10,000 km2, rounded up. Values and surfaces are ln(copies/L).
# K(distance) = alpha * exp(-distance^2 / (2 * rho^2)); sigma is noise SD.
generate_gp_field <- function(study_area, rho, sigma = 3, alpha = 4,
                              mu = log(100), seed = 317) {
  if (any(!is.finite(c(rho, sigma, alpha, mu))) ||
      rho <= 0 || sigma < 0 || alpha <= 0)
    stop("Use positive rho and alpha, nonnegative sigma, and finite mu.")
  settings <- list(patchiness = rho, noise_sd = sigma,
                   amplitude = alpha, mean = mu)
  field <- make_true_field(study_area, settings, seed)
  grid <- make_comparison_grid(study_area)
  grid$true_surface <- predict_concentrations(grid$points, field$points,
                                              field$weights, settings)
  list(study_area = study_area, settings = settings, field = field,
       grid = grid, field_seed = seed)
}

# Part IV. Randomly select the planned survey observations and reconstruct

# The survey size the equation recommends for this area.
recommended_survey_size <- function(settings, area_km2,
                                    coefficients = DEFAULT_COEFFICIENTS) {
  needed <- density_for_precision(settings$precision, settings$patchiness, coefficients)
  if (!is.finite(needed$density)) stop("The equation cannot be inverted for these settings.")
  count <- max(1, ceiling(needed$density * area_km2 / 10000 - 1e-10))
  list(density = needed$density, count = count)
}

# Step 4b. Reconstruct the field from only the recommended survey: pick that many
# field points at random, and redraw the map using just those.
reconstruct_from_survey <- function(field, settings, survey_size, seed = 911) {
  n_points <- nrow(field$points)
  if (survey_size < 1 || survey_size > n_points)
    stop("The recommended survey needs more points than the field has. Use a larger precision value.")

  set.seed(seed)
  chosen <- sample.int(n_points)[seq_len(survey_size)]
  sample_points <- field$points[chosen, , drop = FALSE]

  covariance <- correlation_matrix(sample_points, settings$amplitude,
                                   settings$patchiness, settings$noise_sd, add_noise = TRUE)
  weights <- solve(covariance, field$true_values[chosen] - settings$mean)

  list(chosen = chosen, sample_points = sample_points, weights = weights)
}

# A regular 5 x 5 km grid covering the area, keeping only cells inside it.
# This is where the true and reconstructed maps are compared.
make_comparison_grid <- function(study_area, spacing = GRID_SPACING_KM) {
  box <- study_area$bounding_box
  width  <- box[2] - box[1]
  height <- box[4] - box[3]
  n_cols <- max(1, ceiling(width  / spacing))
  n_rows <- max(1, ceiling(height / spacing))
  if (n_cols * n_rows > 3e6) stop("The comparison grid is too large; use a smaller area.")

  col_centres <- if (width  < spacing) rep((box[1] + box[2]) / 2, n_cols)
                 else box[1] + spacing / 2 + (0:(n_cols - 1)) * spacing
  row_centres <- if (height < spacing) rep((box[3] + box[4]) / 2, n_rows)
                 else box[3] + spacing / 2 + (0:(n_rows - 1)) * spacing

  all_cells <- expand.grid(col = seq_len(n_cols), row = seq_len(n_rows))
  cell_x <- col_centres[all_cells$col]
  cell_y <- row_centres[all_cells$row]
  inside <- points_inside_polygon(cbind(cell_x, cell_y), study_area$polygon)
  if (!any(inside)) stop("The area is too small or narrow to hold a 5 km grid point.")

  list(points = cbind(cell_x[inside], cell_y[inside]),
       col = all_cells$col[inside], row = all_cells$row[inside],
       col_centres = col_centres, row_centres = row_centres,
       n_cols = n_cols, n_rows = n_rows, spacing = spacing)
}

# Main Part IV function. The survey is a random subset of the reference points.
# The error at each 5 km grid center is |reconstructed - reference|.
select_survey_samples <- function(gp_field, sample_plan, seed = 911) {
  if (!isTRUE(all.equal(sample_plan$area_km2, gp_field$study_area$area_km2)))
    stop("The sample plan and GP field must use the same study area.")
  if (!isTRUE(all.equal(sample_plan$rho, gp_field$settings$patchiness)))
    stop("Use the same rho for the plan and GP field.")
  count <- sample_plan$samples_in_area
  if (count > nrow(gp_field$field$points))
    stop("The plan needs more samples than the 40-per-10,000-km2 reference field contains. Increase epsilon.")
  reconstruction <- reconstruct_from_survey(gp_field$field, gp_field$settings,
                                            count, seed)
  grid <- gp_field$grid
  reconstruction$surface <- predict_concentrations(
    grid$points, reconstruction$sample_points,
    reconstruction$weights, gp_field$settings)
  difference <- reconstruction$surface - grid$true_surface
  reconstruction$error_at_each_point <- abs(difference)
  measured_epsilon <- sqrt(mean(difference^2))
  list(study_area = gp_field$study_area, settings = gp_field$settings,
       survey = list(density = sample_plan$samples_per_10000_km2, count = count),
       sample_plan = sample_plan, field = gp_field$field, grid = grid,
       reconstruction = reconstruction, target_precision = sample_plan$epsilon,
       measured_precision = measured_epsilon,
       field_seed = gp_field$field_seed, survey_seed = seed,
       message = sprintf("%d survey samples; measured epsilon = %.3f (target %.3f).",
                         count, measured_epsilon, sample_plan$epsilon))
}

# Put the whole simulation together for one study area.
# Returns the true map, the reconstructed map, the error at every grid point, and
# the overall prediction error.
simulate_and_reconstruct <- function(study_area, patchiness, precision,
                                     mean = log(100), amplitude = 4, noise_sd = 3,
                                     field_seed = 317, survey_seed = 911,
                                     coefficients = DEFAULT_COEFFICIENTS) {
  settings <- list(mean = mean, amplitude = amplitude, noise_sd = noise_sd,
                   patchiness = patchiness, precision = precision)

  survey <- recommended_survey_size(settings, study_area$area_km2, coefficients)
  field  <- make_true_field(study_area, settings, field_seed)

  grid <- make_comparison_grid(study_area)
  grid$true_surface <- predict_concentrations(grid$points, field$points,
                                              field$weights, settings)

  result <- list(study_area = study_area, settings = settings, survey = survey,
                 field = field, grid = grid, target_precision = precision,
                 field_seed = field_seed, survey_seed = survey_seed)

  # If the recommended survey is bigger than the field, we cannot reconstruct.
  if (survey$count > nrow(field$points)) {
    result$reconstruction <- NULL
    result$measured_precision <- NA_real_
    result$message <- sprintf(
      "The equation asks for %s samples, more than the %s field points. Use a larger precision value.",
      format(survey$count, big.mark = ","), format(nrow(field$points), big.mark = ","))
    return(result)
  }

  reconstruction <- reconstruct_from_survey(field, settings, survey$count, survey_seed)
  reconstructed_surface <- predict_concentrations(grid$points,
                                                  reconstruction$sample_points,
                                                  reconstruction$weights, settings)
  difference <- reconstructed_surface - grid$true_surface

  reconstruction$surface <- reconstructed_surface
  reconstruction$error_at_each_point <- abs(difference)      # |reconstructed - true|
  measured_precision <- sqrt(mean(difference^2))             # overall error

  result$reconstruction <- reconstruction
  result$measured_precision <- measured_precision
  result$message <- sprintf(
    "%s survey samples / %s field points; measured precision = %.3f (%s target %.3f) over %s grid points.",
    format(survey$count, big.mark = ","), format(nrow(field$points), big.mark = ","),
    measured_precision, ifelse(measured_precision <= precision, "meets", "exceeds"),
    precision, format(length(grid$true_surface), big.mark = ","))
  result
}


# ---------------------------------------------------------------------------
# Part V. Four plots matching "Simulate & Reconstruct" in the HTML tool
# ---------------------------------------------------------------------------

concentration_palette <- grDevices::colorRampPalette(c("#193466", "#197589", "#5DB59C", "#EFDD6C"))
error_palette         <- grDevices::colorRampPalette(c("#FFF4D2", "#F9A45B", "#CA4237", "#6B163B"))

# Draw just the study area (the polygon from the coordinates), in km. No basemap.
plot_study_area <- function(study_area) {
  polygon <- as.data.frame(study_area$polygon)
  names(polygon) <- c("east", "north")
  polygon <- rbind(polygon, polygon[1, ])
  ggplot2::ggplot(polygon, ggplot2::aes(east, north)) +
    ggplot2::geom_polygon(fill = "#006d77", alpha = 0.08) +
    ggplot2::geom_path(color = "#006d77", linewidth = 1) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "East (km)", y = "North (km)",
      title = sprintf("Study area from coordinates: %s km2",
                      format(round(study_area$area_km2), big.mark = ","))) +
    ggplot2::theme_minimal()
}

# The "precision vs. sampling density" curve.
plot_precision_curve <- function(patchiness, precision,
                                 coefficients = DEFAULT_COEFFICIENTS) {
  needed <- density_for_precision(precision, patchiness, coefficients)
  densities <- seq(max(0.1, CALIBRATED_DENSITY_MIN * 0.5),
                   CALIBRATED_DENSITY_MAX * 1.5, length.out = 300)
  precisions <- precision_at_density(densities, patchiness, coefficients)

  curve_data <- data.frame(density = densities, precision = precisions)
  p <- ggplot2::ggplot(curve_data, ggplot2::aes(density, precision)) +
    ggplot2::annotate("rect", xmin = CALIBRATED_DENSITY_MIN,
      xmax = CALIBRATED_DENSITY_MAX, ymin = -Inf, ymax = Inf,
      fill = "#006d77", alpha = 0.07) +
    ggplot2::geom_line(color = "#006d77", linewidth = 1) +
    ggplot2::geom_hline(yintercept = precision, color = "#d45b24", linetype = 2) +
    ggplot2::labs(x = "Sampling density (samples / 10,000 km2)",
      y = "Map precision (lower = sharper)",
      title = sprintf("Precision vs. sampling density (patchiness = %g km)", patchiness)) +
    ggplot2::theme_minimal()
  if (is.finite(needed$density)) {
    p <- p + ggplot2::geom_vline(xintercept = needed$density,
      color = "#d45b24", linetype = 3) +
      ggplot2::geom_point(data = data.frame(density = needed$density,
        precision = precision), color = "#d45b24", size = 3)
  }
  print(p)
  invisible(needed)
}

# Turn a list of grid values into a matrix the image() function can draw.
grid_values_to_matrix <- function(grid, values) {
  z <- matrix(NA_real_, nrow = grid$n_cols, ncol = grid$n_rows)
  z[cbind(grid$col, grid$row)] <- values
  z
}

# A common 5 km grid for both GP surfaces and the error map.
gp_map_data <- function(grid, values) {
  data.frame(east = grid$points[, 1], north = grid$points[, 2],
             value = values)
}

# 1. GP surface inferred from all 40 observations per 10,000 km2.
plot_reference_field <- function(gp_field) {
  grid <- gp_field$grid
  points <- data.frame(east = gp_field$field$points[, 1],
                       north = gp_field$field$points[, 2])
  ggplot2::ggplot(gp_map_data(grid, grid$true_surface),
                  ggplot2::aes(east, north, fill = value)) +
    ggplot2::geom_tile(width = grid$spacing, height = grid$spacing) +
    ggplot2::geom_point(data = points, ggplot2::aes(east, north),
                        inherit.aes = FALSE, shape = 1, color = "white", size = 1) +
    ggplot2::scale_fill_gradientn(colors = concentration_palette(64),
                                  name = "ln(copies/L)") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "1. Reference GP field", x = "East (km)", y = "North (km)") +
    ggplot2::theme_minimal()
}

# 2. GP surface inferred from only the recommended survey observations.
plot_reconstructed_field <- function(survey) {
  grid <- survey$grid
  points <- data.frame(east = survey$reconstruction$sample_points[, 1],
                       north = survey$reconstruction$sample_points[, 2])
  ggplot2::ggplot(gp_map_data(grid, survey$reconstruction$surface),
                  ggplot2::aes(east, north, fill = value)) +
    ggplot2::geom_tile(width = grid$spacing, height = grid$spacing) +
    ggplot2::geom_point(data = points, ggplot2::aes(east, north),
                        inherit.aes = FALSE, shape = 1, color = "white", size = 1) +
    ggplot2::scale_fill_gradientn(colors = concentration_palette(64),
                                  name = "ln(copies/L)") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "2. GP field from recommended samples",
                  x = "East (km)", y = "North (km)") +
    ggplot2::theme_minimal()
}

# 3. Absolute pointwise difference between the two GP surfaces.
plot_prediction_error_map <- function(survey) {
  grid <- survey$grid
  ggplot2::ggplot(gp_map_data(grid, survey$reconstruction$error_at_each_point),
                  ggplot2::aes(east, north, fill = value)) +
    ggplot2::geom_tile(width = grid$spacing, height = grid$spacing) +
    ggplot2::scale_fill_gradientn(colors = error_palette(64),
      limits = c(0, max(1e-9, survey$reconstruction$error_at_each_point)),
      name = "Absolute log error") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "3. Prediction error map",
                  x = "East (km)", y = "North (km)") +
    ggplot2::theme_minimal()
}

# 4. Distribution of all pointwise errors on the same 5 km grid.
plot_prediction_error_distribution <- function(survey) {
  errors <- data.frame(error = survey$reconstruction$error_at_each_point)
  ggplot2::ggplot(errors, ggplot2::aes(error)) +
    ggplot2::geom_histogram(bins = 30, fill = "#d98050", color = "white") +
    ggplot2::geom_vline(xintercept = survey$target_precision,
                        color = "#006d77", linetype = 2) +
    ggplot2::geom_vline(xintercept = survey$measured_precision,
                        color = "#8e2437") +
    ggplot2::labs(title = "4. Distribution of pointwise prediction error",
                  subtitle = sprintf("Dashed: target epsilon %.3f; solid: measured epsilon %.3f",
                    survey$target_precision, survey$measured_precision),
                  x = "Absolute log error", y = "Number of grid cells") +
    ggplot2::theme_minimal()
}

# The four simulation pictures: true map, reconstructed map, error map, and the
# distribution of errors.
plot_simulation <- function(simulation) {
  reference <- list(grid = simulation$grid, field = simulation$field)
  print(plot_reference_field(reference))
  if (is.null(simulation$reconstruction)) {
    message(simulation$message)
    return(invisible(NULL))
  }
  print(plot_reconstructed_field(simulation))
  print(plot_prediction_error_map(simulation))
  print(plot_prediction_error_distribution(simulation))
  invisible(NULL)
}


# The precision-to-concentration picture (a right-skewed lognormal curve).
plot_precision_meaning <- function(precision, reference = 100) {
  values <- precision_in_copies_per_litre(precision, reference)
  x <- seq(reference * exp(-3.5 * precision), reference * exp(2.8 * precision),
           length.out = 400)
  density <- dlnorm(x, meanlog = log(reference), sdlog = precision)

  curve_data <- data.frame(concentration = x, likelihood = density)
  p <- ggplot2::ggplot(curve_data, ggplot2::aes(concentration, likelihood)) +
    ggplot2::geom_line(color = "#006d77", linewidth = 1) +
    ggplot2::geom_vline(xintercept = values$lower, color = "#587074", linetype = 3) +
    ggplot2::geom_vline(xintercept = values$median, color = "#d45b24") +
    ggplot2::geom_vline(xintercept = values$upper, color = "#587074", linetype = 3) +
    ggplot2::labs(x = "Concentration (copies/L)", y = "Relative likelihood",
      title = sprintf("What precision = %.2f means (median %g copies/L)", precision, reference)) +
    ggplot2::theme_minimal()
  print(p)
  invisible(values)
}


# ---------------------------------------------------------------------------
# Example: replace these GPS corners and parameter values with your own
# ---------------------------------------------------------------------------
# A simple box off the U.S. West Coast, used by all examples below.
# Each row is latitude, longitude.
EXAMPLE_COORDINATES <- rbind(
  c(42.0, -126.0),
  c(42.0, -124.5),
  c(45.0, -124.5),
  c(45.0, -126.0)
)

# Source the file, then select and run these lines in RStudio one step at a time.
# The block also runs if this file is entered directly in an interactive session.
if (interactive() && sys.nframe() == 0L) {
  # I. GPS corners are latitude, longitude in boundary order.
  area <- create_survey_area(coordinates = EXAMPLE_COORDINATES)
  print(area$area_km2)
  print(plot_study_area(area))

  # II. The result includes density and the rounded whole-area count.
  plan <- plan_samples(study_area = area, rho = 100, epsilon = 0.8)
  print(plan$samples_per_10000_km2)
  print(plan$samples_in_area)

  # III. These parameters shape the visual example, not the sample count.
  reference <- generate_gp_field(
    study_area = area, rho = 100, sigma = 3, alpha = 4,
    mu = log(100), seed = 317
  )

  # IV. Randomly select the planned count from the reference observations.
  survey <- select_survey_samples(gp_field = reference,
                                  sample_plan = plan, seed = 911)

  # V. Run each line to display one plot in RStudio's Plots pane.
  print(plot_reference_field(reference))
  print(plot_reconstructed_field(survey))
  print(plot_prediction_error_map(survey))
  print(plot_prediction_error_distribution(survey))
}
