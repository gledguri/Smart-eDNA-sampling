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

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))
pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

h=20
for (j in 2:7) { #Loop around different rho values
# Read simulated GP data
name_file <- paste0(j,'_',h)
   
sim_data_raw_up <- readRDS(here('Output','Raw_GP_fileds_simulated',paste0(name_file,'.rds')))

sim_data_raw_down <- readRDS(here('Output','GP_3000',paste0(name_file,'_simulated_data.rds'))) %>% 
 select(-8) %>% as_tibble()

sim_data_raw <- anti_join(sim_data_raw_up, sim_data_raw_down, 
 by = c("X_utm", "Y_utm", "depth_cat", "x", "y", "mu_sim")) %>% 
 slice_sample(n=1000) %>% 
 rbind(sim_data_raw_down %>% rename(z='mean'))
 
# Final formatting of the simulated data before stan model fitting
sim_data <- sim_data_raw %>% 
 rename(mean='z') %>% 
 mutate(depth=as.numeric(as.character(depth_cat))) %>% 
 arrange(depth,x,y)

saveRDS(sim_data,here('Output','GP_4000',paste0(name_file,'_simulated_data.rds')))

# Priors 
priors <- list(
	alpha_prior = c(4,1),
	mu_prior = c(0,5),
	rho_prior = c(0,1),
	# mag_rho_prior = c(0,1),
	sigma_prior = c(0,1))

# Model stan 

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
									warmup = 1000,
									refresh=10)

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output','GP_4000',paste0(name_file,'_est_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_4000',paste0(name_file,'_pred_GP.rds')))

}