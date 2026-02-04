library(dplyr)
library(ggplot2)
library(here)
library(fields)
library(MASS)
library(PNWColors)
library(purrr)
library(tidyr)
library(tibble)
library(rstan);options(mc.cores = parallel::detectCores()); rstan_options(auto_write = TRUE)

select <- dplyr::select


# Functions -----------------------------------------------------------------------------------

extract_param <- function(model, par) {
	fit <- (methods::selectMethod("summary", signature = "stanfit"))(object = model, pars = par)
	fit <- fit$summary
	fit_df <- as.data.frame(fit)
	fit_df <- round(fit_df, 9)
	return(fit_df)
}

make_index <- function(data, index_variable, index_name) {
	data %>%
		group_by(across(all_of(index_variable))) %>%
		mutate(!!index_name := cur_group_id()) %>%
		ungroup()
}

simulate_spatial_gp <- function(coords, alpha = 4, length_scale = 1, sigma_sim = 3, mu_sim = 0) {
	# Compute distance matrix
	dist_mat <- rdist(coords)
	
	# Create covariance matrix using Gaussian kernel
	cov_mat <- alpha * exp(-(dist_mat^2) / (2 * length_scale^2))
	
	# Add noise to diagonal (matching Stan's approach)
	diag(cov_mat) <- diag(cov_mat) + sigma_sim^2
	
	# Create mean vector
	mu_vec <- rep(mu_sim, nrow(coords))
	
	# Simulate GP (single stage, matching Stan model)
	z <- mvrnorm(n = 1, mu = mu_vec, Sigma = cov_mat)
	
	return(as.numeric(z))
}


# Data ----------------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))


coords <- pred_data %>% 
	mutate(x=X_utm/100,
				 y=Y_utm/100) %>% 
	slice_sample(n = prop_samp_thinn) %>% 
	arrange(depth_cat)
	# select(depth_cat,x,y) %>% 

# dim(coords)

#Bring this outside of the script so that I have it easier to control the simulations
# sim_param_1 <- data.frame(alpha=4,
# 												rho=1,
# 												sigma_sim=1)

for (y in 1:10) {
simulated_data <- coords %>%
	as.data.frame() %>% 
	group_by(depth_cat) %>%
	nest() %>%
	mutate(
		coords = map(data, ~ select(.x, x, y)),
		mu_sim = map_dbl(coords, ~ rnorm(1, 0, 3)),
		z = map2(coords, mu_sim, ~ simulate_spatial_gp(.x, 
																									 alpha = sim_param_1$alpha, 
																									 length_scale = sim_param_1$rho, 
																									 sigma_sim = sim_param_1$sigma_sim, 
																									 mu_sim = .y))
	) %>%
	unnest(c(data, z, mu_sim))

mu_depth <- simulated_data %>% 
	distinct(depth_cat,mu_sim) %>% 
	rename(parameter='depth_cat',
				 value='mu_sim')


sim_param <- bind_rows(sim_param_1 %>% 
											 	pivot_longer(cols = everything(), names_to = "parameter", values_to = "value"),
											 mu_depth)

# simulated_data %>%
# ggplot() +
# 	geom_point(aes(x = x, y = y, color = z),size=3)+
# 	coord_equal() +
# 	scale_color_gradientn(colors = pnw_palette("Bay", 18, type = "continuous")) +
# 	facet_grid(~depth_cat)+
# 	theme_minimal()

pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

sample_100 <- c(1:prop_samp_thinn)
sample_90 <- sample(1:prop_samp_thinn, prop_samp_thinn*0.9,replace = F)
sample_80 <- sample(sample_90, prop_samp_thinn*0.8,replace = F)
sample_70 <- sample(sample_80, prop_samp_thinn*0.7,replace = F)
sample_60 <- sample(sample_70, prop_samp_thinn*0.6,replace = F)
sample_50 <- sample(sample_60, prop_samp_thinn*0.5,replace = F)

thinning_it <- seq(1,0.5,by=-0.1)

samp_lsit <- list(sample_100,sample_90,sample_80,sample_70,sample_60,sample_50)

for (oo in 1:6) {
	samp_size=samp_lsit[[oo]]
	th=thinning_it[oo]

sim_data <- simulated_data %>% 
	ungroup() %>% 
	rename(mean='z') %>% 
	select(-coords) %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	slice(samp_size) %>% 
	# slice_sample(prop = th) %>% 
	arrange(depth,X_utm,Y_utm)

saveRDS(sim_data,here('Output',paste0(ww,'_',y,'_th',th,'_simulated_data.rds')))


# Priors --------------------------------------------------------------------------------------
priors <- list(
	alpha_prior = c(4,1),
	mu_prior = c(0,5),
	rho_prior = c(0,1),
	# mag_rho_prior = c(0,1),
	sigma_prior = c(0,1))


# Model stan -------------------------------------------------------------------------------------

stan_data_7 <- list(
	N_total = nrow(sim_data),
	N_depths = nrow(sim_data %>% distinct(depth)),
	X = cbind(sim_data$x, sim_data$y),
	y = sim_data$mean,
	N_by_depth = sim_data %>% count(depth) %>% pull(n),
	start_idx = c(1,which(c(FALSE, diff(sim_data$depth) != 0))),
	# Pred
	N_pred = nrow(pred_data_by_depth),
	X_pred = cbind(pred_data_by_depth$X_utm/100, pred_data_by_depth$Y_utm/100),
	N_depths_pred = nrow(pred_data_by_depth %>% distinct(depth)),
	N_by_depth_pred = pred_data_by_depth %>% count(depth) %>% pull(n),
	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0)))
)

stan_data_7 <- c(stan_data_7,priors);str(stan_data_7)

fit_7 <- sampling(stan_model_7,
										data = stan_data_7,
										chains = 4,
										iter = 2000,
										warmup = 1000)

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output',paste0(ww,'_',y,'_th',th,'_est_GP_param.csv')))
write.csv(sim_param,file = here('Output',paste0(ww,'_',y,'_th',th,'_sim_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output',paste0(ww,'_',y,'_th',th,'_pred_GP.rds')))
}
}

# pred_data_stan %>% 
# 	ggplot() +
# 	geom_point(aes(x = X_utm, y = Y_utm, color = conc),size=2)+
# 	coord_equal() +
# 	scale_color_gradientn(colors = pnw_palette("Bay", 18, type = "continuous")) +
# 	facet_grid(~depth)+
# 	theme_minimal()
