# Libraries --------------------------------------------------------------
# srun --cpus-per-task=4 --mem=40G --time=14-00:00:00 --pty bash
# .libPaths(c('/mnt/nfs/home/KellyCEG/ggu/Rpack','/mnt/nfs/home/KellyCEG/R/x86_64-pc-linux-gnu-library/4.3',.libPaths()))

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








# Functions --------------------------------------------------------------

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

	# Add noise to diagonal
	diag(cov_mat) <- diag(cov_mat) + sigma_sim^2

	# Create mean vector
	mu_vec <- rep(mu_sim, nrow(coords))

	# Simulate GP (single stage)
	z <- mvrnorm(n = 1, mu = mu_vec, Sigma = cov_mat)

	return(as.numeric(z))
}









# Data import and declaration --------------------------------------------
# stan_model <- stan_model(here('Code','GP_7.stan')) # Load stan model
stan_model <- stan_model(here('Code','GP.stan')) # Load stan model
pred_data <- readRDS(here('Data','pred_data.rds')) # Load the coordinates where the prediction will be made
edna_data <- readRDS(here('Data','edna_data.rds')) # Load the real eDNA data

rho_sim <- c(0.05,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20) # Length scales parameter to be simulated
alpha_sim <- 4 # Amplitude of the GP (this will be fixed throughtout the simulations)
sigma_sim <- 0.1 # Standard deviation of the GP (this will be fixed throughtout the simulations)
depth <- c(0,50,150,300,500) # Depth categories

# Dividing the coordinates by 100 to have center the GP parameters around 0.
coords <- pred_data %>% mutate(x=X_utm/100,y=Y_utm/100)

# Prepare prediction data
pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

# Declare prior distributions for stan model 
priors <- list(
	alpha_prior = c(4,1),
	mu_prior = c(0,5),
	rho_prior = c(0,1),
	sigma_prior = c(0,1))











# # Simulate the raw GP fields ---------------------------------------------
# Initialize lists for simulation data
c <- vector("list", 5) # Simulated points at each depth
z <- vector("list", 5) # Simulated eDNA concentration (log) through GP at each depth
x <- vector("list", 5) # Combined coordinates and eDNA concentration data

# for (h in 1:20) { # Repeat for 20 iterations
h=6
	for (j in 3:12) { # Vary length-scale parameter (rho)
		for (i in 1:5) {  # Loop over each depth as a different layer (with the same rho)
			mu_selected <- rnorm(1, 0, 3)

   c[[i]] <- coords %>% filter(depth_cat==depth[i])

   z[[i]] <- simulate_spatial_gp(c[[i]] %>% select(x,y),
																		alpha = alpha_sim,
																		length_scale = rho_sim[j],
																		sigma_sim = sigma_sim,
																		mu_sim = mu_selected) %>% as.data.frame() %>% setNames('z')

			x[[i]] <- bind_cols(c[[i]],z[[i]],as.data.frame(mu_selected) %>% setNames('mu_sim'))

   cat('depth=',i,' - rho=',j, ' - iteration=',h);cat('\n')
		}

		saveRDS(bind_rows(x),here('Output',
															'Raw_GP_fileds_simulated_alternative',
															paste0(j,'_',h,'.rds')))
	}
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# GP 1-st thinning estimations loops -------------------------------------

scenarios_1 <- c(4000,2000,1000,600,350,200) # Declare number of sampling points scenarios

#--- --- --- loops --- --- ---
# for (k in 1:length(scenarios_1)){ # Loop across different sampling size scenarios
k=5
# j=3
 # for (h in 1:20){ # Loop across 20 different iterations
h=1
  for (j in 3:12){ # Loop across different length-scale parameters
   
  # Read simulated GP data
  name_file <- paste0(j,'_',h)
   
  if (k==1) {sim_data_raw <- readRDS(here('Output','Raw_GP_fileds_simulated',paste0(name_file,'.rds'))) %>% rename(mean='z')}
  if (k!=1) {sim_data_raw <- readRDS(here('Output',paste0('GP_',scenarios_1[k-1]),paste0(name_file,'_simulated_data.rds')))}

  # Summarize simulated mu by depth category
  mu_depth <- sim_data_raw %>% 
   group_by(depth_cat) %>% 
    summarise(est_mu=mean(mean)) %>% 
   rename(parameter='depth_cat', value='est_mu')
  
  sim_param <- bind_rows(data.frame(alpha=alpha_sim,rho=rho_sim[j],sigma_sim=sigma_sim) %>% 
               pivot_longer(cols = everything(), names_to = "parameter", values_to = "value"),mu_depth)

  # Final formatting of the simulated data before stan model fitting
  sim_data <- sim_data_raw %>% 
   slice_sample(n=scenarios_1[k]) %>% 
  #  rename(mean='z') %>% 
   mutate(depth=as.numeric(as.character(depth_cat))) %>% 
   arrange(depth,x,y)

  saveRDS(sim_data,here('Output',paste0('GP_',scenarios_1[k],'x'),paste0(name_file,'_simulated_data.rds')))

  # Model stan 
  stan_data_raw <- list(
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

  stan_data <- c(stan_data_raw,priors);str(stan_data)

  fit <- sampling(stan_model,
           data = stan_data,
           chains = 4,
           iter = 2000,
           warmup = 1000,
          open_progress = FALSE)

  param_list <- fit@model_pars
  param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

  GP_param <- extract_param(fit,param_list) %>%
   as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

  write.csv(GP_param,file = here('Output',paste0('GP_',scenarios_1[k],'x'),paste0(name_file,'_est_GP_param.csv')))
  if (k==1) {write.csv(sim_param,file = here('Output',paste0('GP_',scenarios_1[k]),paste0(name_file,'_sim_GP_param.csv')))}

  pred_data_stan <- pred_data_by_depth %>% bind_cols(.,extract_param(fit,c('y_pred'))) %>% rename(conc=mean)

  saveRDS(pred_data_stan, here('Output',paste0('GP_',scenarios_1[k],'x'),paste0(name_file,'_pred_GP.rds')))
  # } # end loop across different sampling size scenarios
 # } # end loop across 20 different iterations
} # end loop across different length-scale parameters













# GP 2-nd thinning estimations loops -------------------------------------

scenarios_2 <- c(330,300,260,220,150,100) # Declare number of sampling points scenarios

#--- --- --- loops --- --- ---
# for (k in 1:length(scenarios_2)){ # Loop across different sampling size scenarios
k=6
 for (h in 2:20){ # Loop across 20 different iterations
  for (j in 2:15){ # Loop across different length-scale parameters
   
  # Read simulated GP data
  name_file <- paste0(j,'_',h)
   
  if (k==1) {sim_data_raw <- readRDS(here('Output','GP_350',paste0(name_file,'_simulated_data.rds')))}
  if (k!=1) {sim_data_raw <- readRDS(here('Output','GP_sth',paste0('GP_',scenarios_2[k-1]),paste0(name_file,'_simulated_data.rds')))}

  # Summarize simulated mu by depth category
  mu_depth <- sim_data_raw %>% 
   group_by(depth_cat) %>% 
   summarise(est_mu=mean(mean)) %>% 
   rename(parameter='depth_cat', value='est_mu')

  sim_param <- bind_rows(data.frame(alpha=4,rho=rho_sim[j],sigma_sim=3) %>% 
               pivot_longer(cols = everything(), names_to = "parameter", values_to = "value"),mu_depth)

  # Final formatting of the simulated data before stan model fitting
  sim_data <- sim_data_raw %>% 
   slice_sample(n=scenarios_2[k]) %>% 
   rename(mean='mean') %>% 
   mutate(depth=as.numeric(as.character(depth_cat))) %>% 
   arrange(depth,x,y)

  saveRDS(sim_data,here('Output','GP_sth',paste0('GP_',scenarios_2[k]),paste0(name_file,'_simulated_data.rds')))

  # Model stan 
  stan_data_raw <- list(
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

  stan_data <- c(stan_data_raw,priors);str(stan_data_raw)

  fit <- sampling(stan_model,
           data = stan_data,
           chains = 4,
           iter = 2000,
           warmup = 1000,
           open_progress = FALSE)

  param_list <- fit@model_pars
  param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

  GP_param <- extract_param(fit,param_list) %>%
   as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

  write.csv(GP_param,file = here('Output','GP_sth',paste0('GP_',scenarios_2[k]),paste0(name_file,'_est_GP_param.csv')))

  pred_data_stan <- pred_data_by_depth %>% bind_cols(.,extract_param(fit,c('y_pred'))) %>% rename(conc=mean)

  saveRDS(pred_data_stan, here('Output','GP_sth',paste0('GP_',scenarios_2[k]),paste0(name_file,'_pred_GP.rds')))
  } # end loop across different sampling size scenarios
 } # end loop across 20 different iterations
# } # end loop across different length-scale parameters


# Real world data --------------------------------------------------------
spp <- edna_data %>% distinct(species) %>% slice(-n()) %>% pull(species)
scenarios_2 <- c(330,300,260,220,150,100) # Declare number of sampling points scenarios

#--- --- --- loops --- --- ---
for (k in 1:length(scenarios_2)) { # Loop across different sampling size scenarios
for (h in 1:20) { # Loop across 20 different iterations
for (i in 1:length(spp)) { # Loop across different species
if(k==1){
 obs_data_raw <- edna_data %>%
  filter(species==spp[i]) %>% 
  arrange(depth) %>% 
  mutate(x=X_utm/100,
      y=Y_utm/100) %>% 
  arrange(depth,x,y)
 } else {
  obs_data_raw <- readRDS(here('Output','GP_multifish',paste0('Multifish_',scenarios_2[k-1]),
  paste0(spp[i],'_',h,'_thinned_data.rds')))
 }

if (sum(colnames(obs_data_raw)=='depth_cat')){
 obs_data_raw <- obs_data_raw %>% select(-depth_cat)
}

# Final formatting of the observed data before stan model fittin
obs_data <- obs_data_raw %>% rename(depth_cat='depth') %>% 
	slice_sample(n=scenarios_2[k]) %>%
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(obs_data,here('Output','GP_multifish',paste0('Multifish_',scenarios_2[k]),
paste0(spp[i],'_',h,'_thinned_data.rds')))


# Model stan 
  stan_data_raw <- list(
   N_total = nrow(obs_data),
   N_depths = nrow(obs_data %>% distinct(depth)),
   X = cbind(obs_data$x, obs_data$y),
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

  stan_data <- c(stan_data_raw,priors);str(stan_data)

  fit <- sampling(
   stan_model,
   data = stan_data,
   chains = 4,
   iter = 2000,
   warmup = 1000,
   refresh = 10,
   open_progress = FALSE
  )

	param_list <- fit@model_pars
 param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

 GP_param <- extract_param(fit,param_list) %>% 
   as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output','GP_multifish',paste0('Multifish_',scenarios_2[k]),
paste0(spp[i],'_',h,'_est_GP_param.csv')))

pred_data_stan <-
		pred_data_by_depth %>% bind_cols(.,extract_param(fit,c('y_pred'))) %>% 
		rename(conc=mean)
	
	saveRDS(pred_data_stan, here('Output','GP_multifish',paste0('Multifish_',scenarios_2[k]),
 paste0(spp[i],'_',h,'_pred_GP.rds')))

 } # end loop across different species
} # end loop across 20 different iterations
} # end loop across different sampling size scenarios