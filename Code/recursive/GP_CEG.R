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

# sp <- 'Engraulis mordax'
# sp <- 'Sardinops sagax'
sp <- 'Stenobrachius leucopsarus'
edna_data <-
	readRDS(here('Data','edna_data.rds')) %>% filter(species==sp) 

pred_data <- readRDS(here('Data','pred_data.rds'))

obs_data <- edna_data %>%
	# filter(depth%in%c(0,50)) %>%
	arrange(depth)

pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)


# Priors --------------------------------------------------------------------------------------
priors <- list(
	alpha_prior = c(0,5),
	mu_prior = c(0,1),
	rho_prior = c(0,00),
	rho_sd_prior = c(1,0.5),
	sigma_prior = c(3,1))

# stan_model_3 <- stan_model(here('Code','GP_3_CEG.stan'))

# Model 3.1 -------------------------------------------------------------------------------------

stan_data_3_1 <- list(
	N_total = nrow(obs_data),
	N_depths = nrow(obs_data %>% distinct(depth)),
	X = cbind(obs_data$X_utm*100, obs_data$Y_utm*100),
	y = obs_data$mean,
	N_by_depth = obs_data %>% count(depth) %>% pull(n),
	start_idx = c(1,which(c(FALSE, diff(obs_data$depth) != 0))),
	# Pred
	N_pred = nrow(pred_data_by_depth),
	X_pred = cbind(pred_data_by_depth$X_utm*100, pred_data_by_depth$Y_utm*100),
	N_depths_pred = nrow(pred_data_by_depth %>% distinct(depth)),
	N_by_depth_pred = pred_data_by_depth %>% count(depth) %>% pull(n),
	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0)))
)

stan_data_3_1 <- c(stan_data_3_1,priors)

fit_3_1 <- sampling(stan_model_3,
											 data = stan_data_3_1,
											 chains = 4,
											 iter = 2000,
											 warmup = 1000)

param_list <- fit_3_1@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]
# 
extract_param(fit_3_1,c('rho'))
extract_param(fit_3_1,param_list)
# 
# p1_1 <- extract_param(fit_3_1,param_list) %>%
# 	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
# 	# filter(param!='lp__') %>%
# 	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
# 	theme_bw()
# 
# ggsave(here('Plots','param_1.jpg'),p1_1)


pred_data_stan_3_1 <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_3_1,c('y_pred'))) %>% 
	rename(conc=mean)

p2_1 <- obs_data %>% #filter(depth%in%c(0,50)) %>% 
	ggplot() +
	geom_point(data = pred_data_stan_3_1, aes(x = X_utm, y = Y_utm, colour = conc)) +
	# geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
	scale_color_gradientn(colors = pnw_palette("Bay", 18, type = "continuous")) +
	facet_grid(~depth)+
	coord_equal() +
	theme_bw()
p2_1
# ggsave(here('Plots','maps_1.jpg'),p2_1,width = 8,height = 6)



# Priors --------------------------------------------------------------------------------------
priors <- list(
	alpha_prior = c(0,5),
	mu_prior = c(0,1),
	rho_prior = c(0,00),
	rho_sd_prior = c(5,00),
	mag_rho_prior = c(0,2),
	sigma_prior = c(3,1))

# stan_model_4 <- stan_model(here('Code','GP_4_CEG.stan'))

# Model 4.1 -------------------------------------------------------------------------------------

stan_data_4_1 <- list(
	N_total = nrow(obs_data),
	N_depths = nrow(obs_data %>% distinct(depth)),
	X = cbind(obs_data$X_utm/10, obs_data$Y_utm/10),
	y = obs_data$mean,
	N_by_depth = obs_data %>% count(depth) %>% pull(n),
	start_idx = c(1,which(c(FALSE, diff(obs_data$depth) != 0))),
	# Pred
	N_pred = nrow(pred_data_by_depth),
	X_pred = cbind(pred_data_by_depth$X_utm/10, pred_data_by_depth$Y_utm/10),
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

param_list <- fit_4_1@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]
# 
extract_param(fit_4_1,c('rho','mag_rho'))
# 
p1_2 <- extract_param(fit_4_1,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
	# filter(param!='lp__') %>%
	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
	theme_bw()
# 
# ggsave(here('Plots','param_1.jpg'),p1_1)


pred_data_stan_4_1 <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_4_1,c('y_pred'))) %>% 
	rename(conc=mean)

p2_2 <- obs_data %>% #filter(depth%in%c(0,50)) %>% 
	ggplot() +
	geom_point(data = pred_data_stan_4_1, aes(x = X_utm, y = Y_utm, colour = conc)) +
	# geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
	scale_color_gradientn(colors = pnw_palette("Bay", 18, type = "continuous")) +
	facet_grid(~depth)+
	coord_equal() +
	theme_bw()
p2_2
# ggsave(here('Plots','maps_1.jpg'),p2_1,width = 8,height = 6)




# 
# # Model 3.2 -------------------------------------------------------------------------------------
# 
# stan_data_3_2 <- list(
# 	N_total = nrow(obs_data),
# 	N_depths = nrow(obs_data %>% distinct(depth)),
# 	X = cbind(obs_data$X_utm/10, obs_data$Y_utm/10),
# 	y = obs_data$mean,
# 	N_by_depth = obs_data %>% count(depth) %>% pull(n),
# 	start_idx = c(1,which(c(FALSE, diff(obs_data$depth) != 0))),
# 	# Pred
# 	N_pred = nrow(pred_data_by_depth),
# 	X_pred = cbind(pred_data_by_depth$X_utm/10, pred_data_by_depth$Y_utm/10),
# 	N_depths_pred = nrow(pred_data_by_depth %>% distinct(depth)),
# 	N_by_depth_pred = pred_data_by_depth %>% count(depth) %>% pull(n),
# 	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0)))
# )
# 
# stan_data_3_2 <- c(stan_data_3_2,priors)
# 
# fit_3_2 <- sampling(stan_model_3,
# 											 data = stan_data_3_2,
# 											 chains = 4,
# 											 iter = 2000,
# 											 warmup = 1000)
# 
# 
# param_list <- fit_3_2@model_pars
# param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]
# 
# likelihood <- extract_param(fit_3_2,c('lp__'))
# 
# p1_2 <- extract_param(fit_3_2,param_list) %>%
# 	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
# 	# filter(param!='lp__') %>%
# 	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
# 	theme_bw()
# 
# ggsave(here('Plots','param_2.jpg'),p1_2)
# 
# 
# pred_data_stan_3_2 <-
# 	pred_data_by_depth %>% bind_cols(.,extract_param(fit_3_2,c('y_pred'))) %>% 
# 	rename(conc=mean)
# 
# p2_2 <- obs_data %>% #filter(depth%in%c(0,50)) %>% 
# 	ggplot() +
# 	geom_point(data = pred_data_stan_3_2, aes(x = X_utm, y = Y_utm, colour = conc)) +
# 	# geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
# 	scale_color_gradientn(colors = pnw_palette("Bay", 18, type = "continuous")) +
# 	facet_grid(~depth)+
# 	coord_equal() +
# 	theme_bw()
# 
# ggsave(here('Plots','maps_2.jpg'),p2_2,width = 8,height = 6)
# 
# # Model 3.3 -------------------------------------------------------------------------------------
# 
# stan_data_3_3 <- list(
# 	N_total = nrow(obs_data),
# 	N_depths = nrow(obs_data %>% distinct(depth)),
# 	X = cbind(obs_data$X_utm/100, obs_data$Y_utm/100),
# 	y = obs_data$mean,
# 	N_by_depth = obs_data %>% count(depth) %>% pull(n),
# 	start_idx = c(1,which(c(FALSE, diff(obs_data$depth) != 0))),
# 	# Pred
# 	N_pred = nrow(pred_data_by_depth),
# 	X_pred = cbind(pred_data_by_depth$X_utm/100, pred_data_by_depth$Y_utm/100),
# 	N_depths_pred = nrow(pred_data_by_depth %>% distinct(depth)),
# 	N_by_depth_pred = pred_data_by_depth %>% count(depth) %>% pull(n),
# 	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0)))
# )
# 
# stan_data_3_3 <- c(stan_data_3_3,priors)
# 
# fit_3_3 <- sampling(stan_model_3,
# 											 data = stan_data_3_3,
# 											 chains = 4,
# 											 iter = 2000,
# 											 warmup = 1000)
# 
# 
# param_list <- fit_3_3@model_pars
# param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]
# 
# likelihood <- extract_param(fit_3_3,c('lp__'))
# 
# p1_3 <- extract_param(fit_3_3,param_list) %>%
# 	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
# 	# filter(param!='lp__') %>%
# 	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
# 	theme_bw()
# 
# ggsave(here('Plots','param_3.jpg'),p1_3)
# 
# 
# pred_data_stan_3_3 <-
# 	pred_data_by_depth %>% bind_cols(.,extract_param(fit_3_3,c('y_pred'))) %>% 
# 	rename(conc=mean)
# 
# p2_3 <- obs_data %>% #filter(depth%in%c(0,50)) %>% 
# 	ggplot() +
# 	geom_point(data = pred_data_stan_3_3, aes(x = X_utm, y = Y_utm, colour = conc)) +
# 	# geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
# 	scale_color_gradientn(colors = pnw_palette("Bay", 18, type = "continuous")) +
# 	facet_grid(~depth)+
# 	coord_equal() +
# 	theme_bw()
# 
# ggsave(here('Plots','maps_3.jpg'),p2_3,width = 8,height = 6)



