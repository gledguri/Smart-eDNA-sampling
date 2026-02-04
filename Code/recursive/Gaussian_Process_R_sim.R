# Simulate coordinates
n_points <- 900
x <- rep(seq(0,100,length.out=sqrt(n_points)),sqrt(n_points))
y <- rep(seq(0,100,length.out=sqrt(n_points)),each=sqrt(n_points))

coords <- cbind(x, y) %>% as.data.frame()

# Create spatial covariance matrix using exponential or Gaussian kernel
# Choose kernel - Gaussian: exp(- (d^2) / (2 * l^2))
alpha <- 1        # Variance
length_scale <- 40  # Controls smoothness

# Compute distance matrix
dist_mat <- rdist(coords)
cov_mat <- alpha * exp(- (dist_mat^2) / (2 * length_scale^2))

center_point <- c(50, 50)
center_effect <- function(point, center, spread = 15) {
	dist_to_center <- sqrt(sum((point - center)^2))
	return(3 * exp(-(dist_to_center^2) / (2 * spread^2)))
}

# Generate mean values for each point based on distance to center
mu <- apply(coords, 1, function(point) center_effect(point, center_point))

# Simulate GP with the specified mean and covariance
z <- mvrnorm(n = 1,mu = mu, Sigma = cov_mat)

data.frame(x = coords$x, y = coords$y, abundance = z) %>%
		ggplot(aes(x = x, y = y, fill = abundance)) +
		geom_raster() +  # or use geom_tile()
		scale_fill_viridis_c() +
		coord_equal() +
		theme_minimal() +
		labs(title = paste0('Gaussian Process, μ=3; ',"α=",alpha,'; d=',length_scale),
				 fill = "Abundance")
