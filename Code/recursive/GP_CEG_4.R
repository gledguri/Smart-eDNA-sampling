# utils::savehistory(file = NULL)

# Full simulation  ----------------------------------------------------------------------------
library(here)
library(dplyr)
library(ggplot2)
library(tibble)
library(MoMAColors)
library(PNWColors)
library(tidyr)
# library(plotly)
# library(sf)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



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


# Import multifish data -----------------------------------------------------------------------

# pred_data <- readRDS(here('Data','edna_data.rds'))
# readRDS(here('Data','pred_data.rds'))

sp <- 'Engraulis mordax'
# sp <- 'Sardinops sagax'
# sp <- 'Sebastes entomelas'
# sp <- 'Stenobrachius leucopsarus'
# sp <- 'Tarletonbeania crenularis'

edna_data <-
	readRDS(here('Test_rho_sale','Data','edna_data.rds')) %>% filter(species==sp) 

pred_data <- readRDS(here('Test_rho_sale','Data','pred_data.rds'))

obs_data <- edna_data %>%
	# filter(depth%in%c(0,50)) %>%
	arrange(depth)

pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

# Folder path
folder_path <- here("Test_rho_sale", "Plots", paste0("3_", sp))

# Create only if it doesn't exist
if (!dir.exists(folder_path)) {
	dir.create(folder_path, recursive = TRUE)
}

# Priors --------------------------------------------------------------------------------------
priors <- list(
	alpha_prior = c(0,1),
	mu_prior = c(0,1),
	rho_prior = c(0,1),
	mag_rho_prior = c(0,1),
	sigma_prior = c(0,1))

stan_model_4 <- stan_model(here('Test_rho_sale','Code','GP_3_CEG_3.stan'))


# Model 4.1 -------------------------------------------------------------------------------------

stan_data_4_1 <- list(
	N_total = nrow(obs_data),
	N_depths = nrow(obs_data %>% distinct(depth)),
	X = cbind(obs_data$X_utm/100, obs_data$Y_utm/100),
	y = obs_data$mean,
	N_by_depth = obs_data %>% count(depth) %>% pull(n),
	start_idx = c(1,which(c(FALSE, diff(obs_data$depth) != 0))),
	# Pred
	N_pred = nrow(pred_data_by_depth),
	X_pred = cbind(pred_data_by_depth$X_utm/100, pred_data_by_depth$Y_utm/100),
	N_depths_pred = nrow(pred_data_by_depth %>% distinct(depth)),
	N_by_depth_pred = pred_data_by_depth %>% count(depth) %>% pull(n),
	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0)))
)

stan_data_4_1 <- c(stan_data_4_1,priors)

fit_4_1 <- sampling(stan_model_4,
											 data = stan_data_4_1,
											 chains = 4,
											 iter = 2000,
											 warmup = 1000)

saveRDS(stan_data_4_1,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_data_1.rds')))
saveRDS(fit_4_1,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_fit_1.rds')))

param_list <- fit_4_1@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

x1_4 <- extract_param(fit_4_1,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`)

pred_data_stan_4_1 <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_4_1,c('y_pred'))) %>% 
	rename(conc=mean)
saveRDS(fit_4_1,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_pred_1.rds')))
saveRDS(x1_4,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_param_1.rds')))

# Thinning 10% --------------------------------------------------------------------------------

obs_data_filtered <- obs_data %>%
	group_by(depth) %>%
	slice_sample(prop = 0.9) %>%
	ungroup()

stan_data_4_2 <- list(
	N_total = nrow(obs_data_filtered),
	N_depths = nrow(obs_data_filtered %>% distinct(depth)),
	X = cbind(obs_data_filtered$X_utm/100, obs_data_filtered$Y_utm/100),
	y = obs_data_filtered$mean,
	N_by_depth = obs_data_filtered %>% count(depth) %>% pull(n),
	start_idx = c(1,which(c(FALSE, diff(obs_data_filtered$depth) != 0))),
	# Pred
	N_pred = nrow(pred_data_by_depth),
	X_pred = cbind(pred_data_by_depth$X_utm/100, pred_data_by_depth$Y_utm/100),
	N_depths_pred = nrow(pred_data_by_depth %>% distinct(depth)),
	N_by_depth_pred = pred_data_by_depth %>% count(depth) %>% pull(n),
	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0)))
)

stan_data_4_2 <- c(stan_data_4_2,priors)

fit_4_2 <- sampling(stan_model_4,
											 data = stan_data_4_2,
											 chains = 4,
											 iter = 2000,
											 warmup = 1000)

saveRDS(stan_data_4_2,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_data_2.rds')))
saveRDS(fit_4_2,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_fit_2.rds')))

param_list <- fit_4_2@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

x2_4 <- extract_param(fit_4_2,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`)

pred_data_stan_4_2 <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_4_2,c('y_pred'))) %>% 
	rename(conc=mean)
saveRDS(fit_4_2,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_pred_2.rds')))
saveRDS(x2_4,here('Test_rho_sale','Plots',paste0('4_',sp),paste0(sp,'4_stan_param_2.rds')))
