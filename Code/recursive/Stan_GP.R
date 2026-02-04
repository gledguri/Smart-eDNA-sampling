# R code to fit and visualize a Gaussian Process model for species abundance data
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(here)
library(fields)
library(MASS)

# Set Stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

cond_raw <- expand.grid(a = c(6), d = c(1,1.1,1.5,2,3,6), n = c(10,20))

# cond_raw <- data.frame(a=1,d=c(5,15,30,40),n=1600)
cond <- bind_rows(replicate(20, cond_raw, simplify = FALSE))

# 1. Example/simulated data (replace with your own data)
# Example function: simulate species abundance with spatial pattern
# n_obs <- 100  # number of observations
n_obs <- 10^2  # number of observations

# Create x-y coordinates in a 2D space
coords <- data.frame(x = runif(n_obs, 0, 10),	y = runif(n_obs, 0, 10))
# coords <- expand.grid(
# x = seq(0,10,length.out=sqrt(n_obs)),
# y = seq(0,10,length.out=sqrt(n_obs)))

# Create spatial covariance matrix using exponential or Gaussian kernel
# Choose kernel - Gaussian: exp(- (d^2) / (2 * l^2))
alpha <- 6
length_scale <- 1.1
sigma_sim <- 0.5
mu_g <- 0

# Compute distance matrix
dist_mat <- rdist(coords)
cov_mat <- (alpha^2) * exp(- (dist_mat^2) / (2 * length_scale^2))

# Simulate a species abundance pattern with some spatial structure
# Here we use a simple function with two Gaussian "hotspots"
# true_abundance <- function(x, y) {
# 	exp(-(x - 3)^2/8 - (y - 3)^2/8) + 
# 		exp(-(x - 7)^2/5 - (y - 8)^2/5) + 
# 		0.2 * sin(x/2) * cos(y/2)
# }


# Calculate true abundance and add noise
# coords$abundance <- true_abundance(coords$x, coords$y) + rnorm(n_obs, 0, 0.1)
coords$Conc <- mvrnorm(mu = rep(mu_g,n_obs), Sigma = cov_mat)+ rnorm(n_obs, 0, sigma_sim)

mu_sim <- coords$Conc %>% mean()

coords_raw <- coords
p0 <- coords_raw %>% 
	ggplot(aes(x = x, y = y, color = Conc)) +
	geom_point()+
	scale_color_viridis(option = "D") +
	theme_bw() +
	coord_fixed()

thinning_50 <- sample(1:n_obs,n_obs*0.5)
coords_50 <- coords %>% slice(-thinning_50)

thinning_80 <- sample(1:n_obs,n_obs*0.8)
coords_80 <- coords %>% slice(-thinning_80)

p1 <- coords_50 %>% 
	ggplot(aes(x = x, y = y, color = Conc)) +
	geom_point()+
	scale_color_viridis(option = "D") +
	theme_bw() +
	coord_fixed()

cowplot::plot_grid(p0,p1)
# p0 <- coords %>% 
# ggplot(aes(x = x, y = y, fill = Conc)) +
# 	geom_raster() +  # or use geom_tile()
# 	scale_fill_viridis_c() +
# 	coord_equal() +
# 	theme_minimal() +
# 	labs(title = paste0('Gaussian Process, μ=3; ',"α=",alpha,'; d=',length_scale),
# 			 fill = "Conc")


# 2. Create prediction grid
grid_size <- 30  # 30x30 grid for predictions
pred_grid <- expand.grid(
	x = seq(0, 10, length.out = grid_size),
	y = seq(0, 10, length.out = grid_size)
)


# pred_grid %>% 
# 	ggplot()+
# 	geom_point(aes(x=x,y=y))

# 3. Prepare data for Stan
stan_data_thinned_50 <- list(
	N = nrow(coords_50),
	X = cbind(coords_50$x, coords_50$y),
	y = coords_50$Conc,
	N_pred = nrow(pred_grid),
	X_pred = cbind(pred_grid$x, pred_grid$y),
	alpha_prior_scale = 1.0,
	rho_prior_scale = 1.0,
	sigma_prior_scale = 0.5
)

stan_data_thinned_80 <- list(
	N = nrow(coords_80),
	X = cbind(coords_80$x, coords_80$y),
	y = coords_80$Conc,
	N_pred = nrow(pred_grid),
	X_pred = cbind(pred_grid$x, pred_grid$y),
	alpha_prior_scale = 1.0,
	rho_prior_scale = 1.0,
	sigma_prior_scale = 0.5
)

stan_data <- list(
	N = nrow(coords_raw),
	X = cbind(coords_raw$x, coords_raw$y),
	y = coords_raw$Conc,
	N_pred = nrow(pred_grid),
	X_pred = cbind(pred_grid$x, pred_grid$y),
	alpha_prior_scale = 1.0,
	rho_prior_scale = 1.0,
	sigma_prior_scale = 0.5
)

# 4. Fit the model
# Load the Stan model (assuming the model is saved as "species_gp.stan")
stan_model <- stan_model(here('code',"gp.stan"))

# Run MCMC
fit_thin_80 <- sampling(stan_model,data = stan_data_thinned_80,chains = 4,iter = 2000,warmup = 1000)
fit_thin_50 <- sampling(stan_model,data = stan_data_thinned_50,chains = 4,iter = 2000,warmup = 1000)
fit <- sampling(stan_model,data = stan_data,chains = 4,iter = 2000,warmup = 1000)

# # 5. Extract and process results
# # Get parameter summaries
# print(fit_thin, pars = c("alpha", "rho", "sigma", "mu"))
# data.frame(alpha=alpha,rho=length_scale,sigma=sigma_sim,mu=mu_sim) %>% 
# 	t() %>% 
# 	as.data.frame() %>% setNames('mean')

g_sigma_list <- readRDS(here('Data','g_sigma_list.RDS'))
pred_grid_list <- readRDS(here('Data','pred_grid_list.RDS'))
stan_param_list <- readRDS(here('Data','stan_param_list.RDS'))

g_sigma <- bind_rows(g_sigma_list,.id='id') %>% 
	left_join(.,cond %>% rownames_to_column('id'),by='id')
pred_grid <- bind_rows(pred_grid_list,.id='id') %>% 
	left_join(.,cond %>% rownames_to_column('id'),by='id')
stan_param <- bind_rows(stan_param_list,.id='id') %>% 
	left_join(.,cond %>% rownames_to_column('id'),by='id')

g_sigma %>% 
	mutate(delta_0=round(1/((n^2)/10^2),2)) %>% 
	mutate(delta_50=round(1/(((n^2)*0.5)/10^2),2)) %>% 
	mutate(delta_80=round(1/(((n^2)*0.2)/10^2),2)) %>% 
	# filter(n==30) %>% 
	ggplot()+
	geom_point(aes(x=d,y=g_sigma_80,color=as.factor(delta_80)))+
	geom_smooth(aes(x=d,y=g_sigma_80,color=as.factor(delta_80)),lty=2,se=F,span=0.5)+
	geom_point(aes(x=d,y=g_sigma_50,color=as.factor(delta_50)))+
	geom_smooth(aes(x=d,y=g_sigma_50,color=as.factor(delta_50)),se=F,span=0.5)+
	# geom_point(aes(x=d,y=g_sigma_80,color=as.factor(delta_50)))+
	# geom_smooth(aes(x=d,y=g_sigma_80,color=as.factor(n)),se=F,span=0.5)+
	# geom_point(aes(x=d,y=g_sigma_50,color=as.factor(n)))+
	facet_wrap(~delta_0)+
	scale_y_log10()+
	theme_bw()

# Extract posterior predictions
pred_samples <- rstan::extract(fit, "f_pred")$f_pred
pred_mean <- apply(pred_samples, 2, mean)
pred_samples_thin_50 <- rstan::extract(fit_thin_50, "f_pred")$f_pred
pred_mean_thin_50 <- apply(pred_samples_thin_50, 2, mean)
pred_samples_thin_80 <- rstan::extract(fit_thin_80, "f_pred")$f_pred
pred_mean_thin_80 <- apply(pred_samples_thin_80, 2, mean)
# pred_sd <- apply(pred_samples, 2, sd)
# pred_lower <- apply(pred_samples, 2, quantile, 0.025)
# pred_upper <- apply(pred_samples, 2, quantile, 0.975)

# Add predictions to the grid
pred_grid$Conc <- pred_mean
pred_grid$Conc_thin_50 <- pred_mean_thin_50
pred_grid$Conc_thin_80 <- pred_mean_thin_80
# pred_grid$abundance_sd <- pred_sd
# pred_grid$abundance_lower <- pred_lower
# pred_grid$abundance_upper <- pred_upper

# 6. Visualize results
# Plot original data points

# Plot predicted surface
p2 <- ggplot(pred_grid, aes(x = x, y = y, fill = Conc)) +
	geom_tile() +
	scale_fill_viridis(option = "D") +
	geom_point(data = coords_raw, aes(x = x, y = y), size = 3, shape = 4, 
						 fill = "white", color = "tomato3") +
	ggtitle("Predicted Concentration (all samples)") +
	theme_minimal() +
	coord_fixed()

p3 <- ggplot(pred_grid, aes(x = x, y = y, fill = Conc_thin_80)) +
	geom_tile() +
	scale_fill_viridis(option = "D") +
	geom_point(data = coords, aes(x = x, y = y), size = 3, shape = 4, 
						 fill = "white", color = "tomato3") +
	ggtitle("Predicted Concentration (thinned samples)") +
	theme_minimal() +
	coord_fixed()

cowplot::plot_grid(p2,p3)

sd(pred_grid$Conc-pred_grid$Conc_thin_50)*3
sd(pred_grid$Conc-pred_grid$Conc_thin_80)*3


# p4 <- 
ggplot(pred_grid, aes(x = x, y = y, fill = Conc-Conc_thin_50)) +
	geom_tile() +
	scale_fill_viridis(option = "D") +
	geom_point(data = anti_join(coords_raw, coords), aes(x = x, y = y), size = 2, shape = 21, 
						 fill = 'white',color = "black", alpha = 0.5) +
	ggtitle("ΔConcentration (raw-thinned)") +
	theme_minimal() +
	coord_fixed()
