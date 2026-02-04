# .libPaths(c('/mnt/nfs/home/KellyCEG/ggu/Rpack','/mnt/nfs/home/KellyCEG/R/x86_64-pc-linux-gnu-library/4.3',.libPaths()))

# # R code to fit and visualize a Gaussian Process model for species abundance data
# library(rstan)
# library(dplyr)
# library(tidyr)
# library(fields)
# library(MASS)
# library(tibble)
# 
# # Set Stan options
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
# 
# extract_param <- function (model, par) {
# 	fit <- (methods::selectMethod("summary", signature = "stanfit"))(object = model, par = par)
# 	fit <- fit$summary
# 	return(fit %>% unlist() %>% as.data.frame %>% round(., 9))}
# 
# # cond_raw <- expand.grid(a = c(6), d = c(3), n = c(10))
# # cond <- cond_raw
# cond_raw <- expand.grid(a = c(6), d = c(1,1.5,3,4,6), n = c(10,20,30))
# cond <- bind_rows(replicate(10, cond_raw, simplify = FALSE))
# 
# stan_param_list <- list()
# pred_grid_list <- list()
# g_sigma_list <- list()
# mu_sim_list <- list()
# coords_raw_list <- list()
# coords_10_list <- list()
# coords_30_list <- list()
# coords_50_list <- list()
# coords_80_list <- list()

for (i in st_iter:end_iter) {

n_obs <- cond$n[i]^2  # number of observations
coords <- data.frame(x = runif(n_obs, 0, 10),	y = runif(n_obs, 0, 10))

# Choose kernel - Gaussian: exp(- (d^2) / (2 * l^2))
# alpha <- 6
# length_scale <- 1.1
# sigma_sim <- 0.5
# mu_g <- 0
alpha <- cond$a[i]
length_scale <- cond$d[i]
sigma_sim <- 0.5
mu_g <- 0

# Compute distance matrix
dist_mat <- rdist(coords)
cov_mat <- (alpha^2) * exp(- (dist_mat^2) / (2 * length_scale^2))

# Simulate a species abundance pattern with some spatial structure
coords$Conc <- mvrnorm(mu = rep(mu_g,n_obs), Sigma = cov_mat)+ rnorm(n_obs, 0, sigma_sim)

mu_sim_list[[i]] <- coords$Conc %>% mean()

coords_raw_list[[i]] <- coords

# Spatial thinning
thinning_10 <- sample(1:n_obs,n_obs*0.1)
coords_10 <- coords %>% slice(-thinning_10)

thinning_30 <- sample(1:n_obs,n_obs*0.3)
coords_30 <- coords %>% slice(-thinning_30)

thinning_50 <- sample(1:n_obs,n_obs*0.5)
coords_50 <- coords %>% slice(-thinning_50)

thinning_80 <- sample(1:n_obs,n_obs*0.8)
coords_80 <- coords %>% slice(-thinning_80)

coords_10_list[[i]] <- coords_10
coords_30_list[[i]] <- coords_30
coords_50_list[[i]] <- coords_50
coords_80_list[[i]] <- coords_80

# 2. Create prediction grid
grid_size <- 30  # 30x30 grid for predictions
pred_grid <- expand.grid(
	x = seq(0, 10, length.out = grid_size),
	y = seq(0, 10, length.out = grid_size)
)


# 3. Prepare data for Stan
stan_data_thinned_10 <- list(
	N = nrow(coords_10),
	X = cbind(coords_10$x, coords_10$y),
	y = coords_10$Conc,
	N_pred = nrow(pred_grid),
	X_pred = cbind(pred_grid$x, pred_grid$y),
	alpha_prior_scale = 1.0,
	rho_prior_scale = 1.0,
	sigma_prior_scale = 0.5
)

stan_data_thinned_30 <- list(
	N = nrow(coords_30),
	X = cbind(coords_30$x, coords_30$y),
	y = coords_30$Conc,
	N_pred = nrow(pred_grid),
	X_pred = cbind(pred_grid$x, pred_grid$y),
	alpha_prior_scale = 1.0,
	rho_prior_scale = 1.0,
	sigma_prior_scale = 0.5
)

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
	N = nrow(coords),
	X = cbind(coords$x, coords$y),
	y = coords$Conc,
	N_pred = nrow(pred_grid),
	X_pred = cbind(pred_grid$x, pred_grid$y),
	alpha_prior_scale = 1.0,
	rho_prior_scale = 1.0,
	sigma_prior_scale = 0.5
)

# 4. Fit the model
# Load the Stan model (assuming the model is saved as "species_gp.stan")
# stan_model <- stan_model(here('code',"gp.stan"))
stan_model <- stan_model("GP.stan")

# Run MCMC
cat('\n');cat('N=');cat(stan_data_thinned_80$N);cat('; iter',i,'/',end_iter);cat('\n')
fit_thin_80 <- sampling(stan_model,data = stan_data_thinned_80,chains = 4,iter = 2000,warmup = 1000)
cat('\n');cat('N=');cat(stan_data_thinned_50$N);cat('; iter',i,'/',end_iter);cat('\n')
fit_thin_50 <- sampling(stan_model,data = stan_data_thinned_50,chains = 4,iter = 2000,warmup = 1000)
cat('\n');cat('N=');cat(stan_data_thinned_30$N);cat('; iter',i,'/',end_iter);cat('\n')
fit_thin_30 <- sampling(stan_model,data = stan_data_thinned_30,chains = 4,iter = 2000,warmup = 1000)
cat('\n');cat('N=');cat(stan_data_thinned_10$N);cat('; iter',i,'/',end_iter);cat('\n')
fit_thin_10 <- sampling(stan_model,data = stan_data_thinned_10,chains = 4,iter = 2000,warmup = 1000)
cat('\n');cat('N=');cat(stan_data$N);cat('; iter',i,'/',end_iter);cat('\n')
fit <- sampling(stan_model,data = stan_data,chains = 4,iter = 2000,warmup = 1000)

# # 5. Extract and process results
# Get parameter summaries
stan_param_list[[i]] <- 
	bind_rows(
extract_param(fit, c("alpha", "rho", "sigma", "mu")) %>% 
	rownames_to_column('param') %>% 
	mutate(N=cond$n[i]^2),
extract_param(fit_thin_50, c("alpha", "rho", "sigma", "mu")) %>% 
	rownames_to_column('param') %>% 
	mutate(N=(cond$n[i]^2)*0.5)) %>% 
	bind_rows(.,
extract_param(fit_thin_80, c("alpha", "rho", "sigma", "mu")) %>% 
rownames_to_column('param') %>% 
	mutate(N=(cond$n[i]^2)*0.2)) %>% 
	bind_rows(.,
extract_param(fit_thin_10, c("alpha", "rho", "sigma", "mu")) %>% 
rownames_to_column('param') %>% 
	mutate(N=(cond$n[i]^2)*0.9)) %>% 
	bind_rows(.,
extract_param(fit_thin_30, c("alpha", "rho", "sigma", "mu")) %>% 
rownames_to_column('param') %>% 
	mutate(N=(cond$n[i]^2)*0.7))
	

# data.frame(alpha=alpha,rho=length_scale,sigma=sigma_sim,mu=mu_sim) %>% 
# 	t() %>% 
# 	as.data.frame() %>% setNames('mean')

# Extract posterior predictions
pred_samples <- rstan::extract(fit, "f_pred")$f_pred
pred_mean <- apply(pred_samples, 2, mean)
pred_samples_thin_10 <- rstan::extract(fit_thin_10, "f_pred")$f_pred
pred_mean_thin_10 <- apply(pred_samples_thin_10, 2, mean)
pred_samples_thin_30 <- rstan::extract(fit_thin_30, "f_pred")$f_pred
pred_mean_thin_30 <- apply(pred_samples_thin_30, 2, mean)
pred_samples_thin_50 <- rstan::extract(fit_thin_50, "f_pred")$f_pred
pred_mean_thin_50 <- apply(pred_samples_thin_50, 2, mean)
pred_samples_thin_80 <- rstan::extract(fit_thin_80, "f_pred")$f_pred
pred_mean_thin_80 <- apply(pred_samples_thin_80, 2, mean)

# Add predictions to the grid
pred_grid$Conc <- pred_mean
pred_grid$Conc_thin_10 <- pred_mean_thin_10
pred_grid$Conc_thin_30 <- pred_mean_thin_30
pred_grid$Conc_thin_50 <- pred_mean_thin_50
pred_grid$Conc_thin_80 <- pred_mean_thin_80

pred_grid_list[[i]] <- pred_grid

g_sigma_list[[i]] <-
sd(pred_grid$Conc-pred_grid$Conc_thin_50)*3 %>% 
	as.data.frame() %>% 
	setNames('g_sigma_50') %>% 
	cbind(.,
sd(pred_grid$Conc-pred_grid$Conc_thin_80)*3 %>% 
	as.data.frame() %>% 
	setNames('g_sigma_80') ) %>% 
	cbind(.,
sd(pred_grid$Conc-pred_grid$Conc_thin_10)*3 %>% 
	as.data.frame() %>% 
	setNames('g_sigma_10') ) %>% 
	cbind(.,
sd(pred_grid$Conc-pred_grid$Conc_thin_30)*3 %>% 
	as.data.frame() %>% 
	setNames('g_sigma_30') )

cat('#######################################')
cat('\n');cat(((i-st_iter+1)/end_iter)*100,'% complete');cat('\n')
cat('#######################################');cat('\n')
}