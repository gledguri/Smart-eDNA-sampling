# Server code ---------------------------------------------------------------------------------

# srun --cpus-per-task=4 --mem=40G --time=30-00:00:00 --pty bash
# 
.libPaths(c('/mnt/nfs/home/KellyCEG/ggu/Rpack','/mnt/nfs/home/KellyCEG/R/x86_64-pc-linux-gnu-library/4.3',.libPaths()))

# Libraries -----------------------------------------------------------------------------------
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



# Simulate the raw GP fields ------------------------------------------------------------------
rho_sim <- c(0.01,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20)
alpha_sim <- 4
sigma_sim <- 3
depth <- c(0,50,150,300,500)

coords <- pred_data %>% 
	mutate(x=X_utm/100,
				 y=Y_utm/100)


i=1
j=1

z <- vector("list", 5)
c <- vector("list", 5)
x <- vector("list", 5)

for (h in 1:20) {
	for (j in 1:15) {
		for (i in 1:5) {
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
		# bind_rows(x) %>% group_by(depth_cat) %>% summarise(z_mean=mean(z),mu_sim=mean(mu_sim))
		
		saveRDS(bind_rows(x),here('Output',
															'Raw_GP_fileds_simulated',
															paste0(j,'_',h,'.rds')))
	}
}

# Pred 3000 -----------------------------------------------------------------------------------

# Data 
pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))
rho_sim <- c(0.01,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20)

coords <- pred_data %>% mutate(x=X_utm/100,y=Y_utm/100)

j=1
h=1
for (j in 1:15) {
name_file <- paste0(j,'_',h)
sim_data_raw <- readRDS(here('Output','Raw_GP_fileds_simulated',
														 paste0(name_file,'.rds'))) %>% 
	cbind(.,coords %>% select(depth_cat)) %>% select(-8) 

sim_param_1 <- data.frame(alpha=4,rho=rho_sim[j],sigma_sim=3)

mu_depth <-
sim_data_raw %>% as_tibble() %>% 
	group_by(depth_cat) %>% 
	summarise(est_mu=mean(z)) %>% 
	rename(parameter='depth_cat',
				 value='est_mu')

sim_param <- bind_rows(sim_param_1 %>% 
											 	pivot_longer(cols = everything(), names_to = "parameter", values_to = "value"),
											 mu_depth)

pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

	
sim_data <- sim_data_raw %>% 
	slice_sample(n=3000) %>% 
	rename(mean='z') %>% 
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)


saveRDS(sim_data,here('Output','GP_3000',paste0(name_file,'_simulated_data.rds')))
	
	
# Priors 
priors <- list(
	alpha_prior = c(4,1),
	mu_prior = c(0,5),
	rho_prior = c(0,1),
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

write.csv(GP_param,file = here('Output','GP_3000',paste0(name_file,'_est_GP_param.csv')))
write.csv(sim_param,file = here('Output','GP_3000',paste0(name_file,'_sim_GP_param.csv')))

pred_data_stan <- pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_3000',paste0(name_file,'_pred_GP.rds')))
}


# Pred 2000 ------------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

j=1
h=5
for (j in 1:15) {
	name_file <- paste0(j,'_',h)
	sim_data_raw <- readRDS(here('Output','GP_3000',
															 paste0(name_file,'_simulated_data.rds')))
	
	
	pred_data_by_depth <- pred_data %>% 
		rename(depth='depth_cat') %>% 
		mutate(depth=as.numeric(as.character(depth))) %>% 
		arrange(depth)
	
	
	sim_data <- sim_data_raw %>% 
		slice_sample(n=2000) %>% 
		mutate(depth=as.numeric(as.character(depth_cat))) %>% 
		arrange(depth,x,y)
	
	
	saveRDS(sim_data,here('Output','GP_2000',paste0(name_file,'_simulated_data.rds')))
	
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
	
	write.csv(GP_param,file = here('Output','GP_2000',paste0(name_file,'_est_GP_param.csv')))
	# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))
	
	pred_data_stan <-
		pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
		rename(conc=mean)
	
	saveRDS(pred_data_stan, here('Output','GP_2000',paste0(name_file,'_pred_GP.rds')))
}

# Pred 1000 ------------------------------------------------------------------------------------
	
	pred_data <- readRDS(here('Data','pred_data.rds'))
	stan_model_7 <- stan_model(here('Code','GP_7.stan'))
	
	j=2
	h=17
	# for (j in 14:15) {
		name_file <- paste0(j,'_',h)
		sim_data_raw <- readRDS(here('Output','GP_2000',
																 paste0(name_file,'_simulated_data.rds')))
		
		
		pred_data_by_depth <- pred_data %>% 
			rename(depth='depth_cat') %>% 
			mutate(depth=as.numeric(as.character(depth))) %>% 
			arrange(depth)
		
		
		sim_data <- sim_data_raw %>% 
			slice_sample(n=1000) %>% 
			mutate(depth=as.numeric(as.character(depth_cat))) %>% 
			arrange(depth,x,y)
		
		
		saveRDS(sim_data,here('Output','GP_1000',paste0(name_file,'_simulated_data.rds')))
		
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
		
		write.csv(GP_param,file = here('Output','GP_1000',paste0(name_file,'_est_GP_param.csv')))
		# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))
		
		pred_data_stan <-
			pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
			rename(conc=mean)
		
		saveRDS(pred_data_stan, here('Output','GP_1000',paste0(name_file,'_pred_GP.rds')))
	# }

# Pred 600 ------------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

j=11
h=11

# for (j in 14:15) {
	name_file <- paste0(j,'_',h)
	sim_data_raw <- readRDS(here('Output','GP_1000',
															 paste0(name_file,'_simulated_data.rds')))
	
	
	pred_data_by_depth <- pred_data %>% 
		rename(depth='depth_cat') %>% 
		mutate(depth=as.numeric(as.character(depth))) %>% 
		arrange(depth)
	
	
	sim_data <- sim_data_raw %>% 
		slice_sample(n=600) %>% 
		mutate(depth=as.numeric(as.character(depth_cat))) %>% 
		arrange(depth,x,y)
	
	
	saveRDS(sim_data,here('Output','GP_600',paste0(name_file,'_simulated_data.rds')))
	
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
	
	write.csv(GP_param,file = here('Output','GP_600',paste0(name_file,'_est_GP_param.csv')))
	# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))
	
	pred_data_stan <-
		pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
		rename(conc=mean)
	
	saveRDS(pred_data_stan, here('Output','GP_600',paste0(name_file,'_pred_GP.rds')))
# }


# Pred 350 ------------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

j=11
h=11
# for (j in 1:15) {
	name_file <- paste0(j,'_',h)
	sim_data_raw <- readRDS(here('Output','GP_600',
															 paste0(name_file,'_simulated_data.rds')))
	
	
	pred_data_by_depth <- pred_data %>% 
		rename(depth='depth_cat') %>% 
		mutate(depth=as.numeric(as.character(depth))) %>% 
		arrange(depth)
	
	
	sim_data <- sim_data_raw %>% 
		slice_sample(n=350) %>% 
		mutate(depth=as.numeric(as.character(depth_cat))) %>% 
		arrange(depth,x,y)
	
	
	saveRDS(sim_data,here('Output','GP_350',paste0(name_file,'_simulated_data.rds')))
	
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
	
	write.csv(GP_param,file = here('Output','GP_350',paste0(name_file,'_est_GP_param.csv')))
	# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))
	
	pred_data_stan <-
		pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
		rename(conc=mean)
	
	saveRDS(pred_data_stan, here('Output','GP_350',paste0(name_file,'_pred_GP.rds')))
# }
# Pred 200 ------------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

j=11
h=11
# for (j in 1:15) {
	name_file <- paste0(j,'_',h)
	sim_data_raw <- readRDS(here('Output','GP_350',
															 paste0(name_file,'_simulated_data.rds')))
	
	
	pred_data_by_depth <- pred_data %>% 
		rename(depth='depth_cat') %>% 
		mutate(depth=as.numeric(as.character(depth))) %>% 
		arrange(depth)
	
	
	sim_data <- sim_data_raw %>% 
		slice_sample(n=200) %>% 
		mutate(depth=as.numeric(as.character(depth_cat))) %>% 
		arrange(depth,x,y)
	
	
	saveRDS(sim_data,here('Output','GP_200',paste0(name_file,'_simulated_data.rds')))
	
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
	
	write.csv(GP_param,file = here('Output','GP_200',paste0(name_file,'_est_GP_param.csv')))
	# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))
	
	pred_data_stan <-
		pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
		rename(conc=mean)
	
	saveRDS(pred_data_stan, here('Output','GP_200',paste0(name_file,'_pred_GP.rds')))
# }





# Import data ---------------------------------------------------------------------------------
spp <- c(
	'Sebastes entomelas',
	'Stenobrachius leucopsarus',
	'Engraulis mordax',
	'Sardinops sagax',
	'Tarletonbeania crenularis',
	'Clupea pallasii',
	'Leuroglossus stilbius',
	'Microstomus pacificus',
	'Scomber japonicus',
	'Tactostoma macropus',
	'Thaleichthys pacificus',
	'Trachurus symmetricus'
)

est_gp_param_sp_list <- vector("list", length(spp))

for (j in 1:length(spp)) {
	
	est_gp_param_sp <- list.files(here('Output','multifish_spp'), pattern = paste0("^", spp[j], ".*\\.csv$"), full.names = TRUE)
	
	df_list <- lapply(est_gp_param_sp, read.csv)
	
	# library(dplyr)
	combined_df <- bind_rows(
		lapply(seq_along(df_list), function(i) {
			df_list[[i]] %>% mutate(source_file = basename(est_gp_param_sp[i]))
		})
	)
	
	est_gp_param_sp_list[[j]] <- combined_df %>% 
		mutate(iteration=gsub(".*_(\\d+)_th.*", "\\1", source_file)) %>% 
		mutate(iteration=as.numeric(iteration)) %>% 
		mutate(thinning=gsub(".*_th([0-9.]+)_est.*", "\\1", source_file)) %>% 
		mutate(thinning=as.numeric(thinning)) %>% 
		arrange(iteration,thinning) %>% 
		mutate(sp=spp[j]) %>% 
		select(-source_file,-X)
}	


est_gp_param_sp <- bind_rows(est_gp_param_sp_list)

est_rho_sp <- est_gp_param_sp %>% filter(thinning==1) %>% 
	filter(param=='rho') %>% 
	group_by(sp) %>% 
	summarise(est_rho=exp(mean(mean))*100)

sim_param_files <- list.files(here('Output','GP_3000'),pattern = 'sim_GP_param.csv')
my_colors <- c("#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d")
my_sp_colors <- c('#6a3d9a','#cab2d6','#ff7f01','#fdbf6f','#e31b1d','#fb9a99','#33a02c','#b2df8a','#1f78b4','#009999','#999900', '#a6cee3')

sim_raw <- here("Output", "GP_3000", sim_param_files) %>% 
	setNames(basename(.)) %>% 
	lapply(read.csv) %>% 
	bind_rows(.id = "source") %>% 
	select(-X) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% 
	arrange(rho,it) %>% as_tibble()
rho_sim_df <- sim_raw %>% filter(parameter=='rho') %>% distinct(rho,value) %>% 
	mutate(value=value*100)

est_param_files_3000 <- list.files(here('Output','GP_3000'),pattern = 'est_GP_param.csv')
est_param_files_2000 <- list.files(here('Output','GP_2000'),pattern = 'est_GP_param.csv')
est_param_files_1000 <- list.files(here('Output','GP_1000'),pattern = 'est_GP_param.csv')
est_param_files_600 <- list.files(here('Output','GP_600'),pattern = 'est_GP_param.csv')
est_param_files_350 <- list.files(here('Output','GP_350'),pattern = 'est_GP_param.csv')
est_param_files_200 <- list.files(here('Output','GP_200'),pattern = 'est_GP_param.csv')

est_3000 <- here("Output", "GP_3000", est_param_files_3000) %>% 
	setNames(basename(.)) %>% 
	lapply(read.csv) %>% 
	bind_rows(.id = "source") %>% 
	select(-X) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% 
	arrange(rho,it) %>% 
	mutate(N=3000) %>% 
	as_tibble()

est_2000 <- here("Output", "GP_2000", est_param_files_2000) %>% 
	setNames(basename(.)) %>% 
	lapply(read.csv) %>% 
	bind_rows(.id = "source") %>% 
	select(-X) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% 
	arrange(rho,it) %>% 
	mutate(N=2000) %>% 
	as_tibble()

est_1000 <- here("Output", "GP_1000", est_param_files_1000) %>% 
	setNames(basename(.)) %>% 
	lapply(read.csv) %>% 
	bind_rows(.id = "source") %>% 
	select(-X) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% 
	arrange(rho,it) %>% 
	mutate(N=1000) %>% 
	as_tibble()

est_600 <- here("Output", "GP_600", est_param_files_600) %>% 
	setNames(basename(.)) %>% 
	lapply(read.csv) %>% 
	bind_rows(.id = "source") %>% 
	select(-X) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% 
	arrange(rho,it) %>% 
	mutate(N=600) %>% 
	as_tibble()

est_350 <- here("Output", "GP_350", est_param_files_350) %>% 
	setNames(basename(.)) %>% 
	lapply(read.csv) %>% 
	bind_rows(.id = "source") %>% 
	select(-X) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% 
	arrange(rho,it) %>% 
	mutate(N=350) %>% 
	as_tibble()

est_200 <- here("Output", "GP_350", est_param_files_200) %>% 
	setNames(basename(.)) %>% 
	lapply(read.csv) %>% 
	bind_rows(.id = "source") %>% 
	select(-X) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% 
	arrange(rho,it) %>% 
	mutate(N=200) %>% 
	as_tibble()

est <- bind_rows(est_3000,est_2000,est_1000,est_600,est_350,est_200) %>% 
	rename(parameter='param') %>% 
	rename(value='mean') %>% 
	mutate(parameter = if_else(parameter == "mu[1]", "mu_0", parameter)) %>% 
	mutate(parameter = if_else(parameter == "mu[2]", "mu_50", parameter)) %>% 
	mutate(parameter = if_else(parameter == "mu[3]", "mu_150", parameter)) %>% 
	mutate(parameter = if_else(parameter == "mu[4]", "mu_300", parameter)) %>% 
	mutate(parameter = if_else(parameter == "mu[5]", "mu_500", parameter)) %>% 
	mutate(value=if_else(parameter=='rho',exp(value),value)) %>% 
	mutate(lo_q=if_else(parameter=='rho',exp(X2.5.),X2.5.)) %>% 
	mutate(up_q=if_else(parameter=='rho',exp(X97.5.),X97.5.))

sim <- sim_raw %>% 
	mutate(N=9999) %>% 
	mutate(parameter = if_else(parameter == "0", "mu_0", parameter)) %>% 
	mutate(parameter = if_else(parameter == "50", "mu_50", parameter)) %>% 
	mutate(parameter = if_else(parameter == "150", "mu_150", parameter)) %>% 
	mutate(parameter = if_else(parameter == "300", "mu_300", parameter)) %>% 
	mutate(parameter = if_else(parameter == "500", "mu_500", parameter)) %>% 
	mutate(parameter = if_else(parameter == "sigma_sim", "sigma", parameter))

# Analysis ------------------------------------------------------------------------------------

my_theme <- function() {
		theme(
			axis.title = element_text(size=16),
			axis.text = element_text(size=15),
			strip.text = element_text(size = 15),
			legend.title = element_text(size = 15),
			legend.text  = element_text(size = 14))
}

p1 <- est %>% 
	left_join(.,sim %>% select(-N,-source) %>%rename(sim_val='value'),
						by=c('parameter','rho','it')) %>% 
	filter(parameter=='rho') %>% 
	# filter(!rho==1) %>% 
	mutate(value=value*100) %>% 
	mutate(sim_val=sim_val*100) %>% 
	mutate(lo_q=lo_q*100) %>% 
	mutate(up_q=up_q*100) %>%
	mutate(
		N_dens = case_when(N == 3000 ~ "N==3000*' (r=18 km; s=350 km'^2*')'",N == 2000 ~ "N==2000*' (r=22 km; s=500 km'^2*')'",N == 1000 ~ "N==1000*' (r=32 km; s=1000 km'^2*')'",N == 600  ~ "N==600*' (r=40 km; s=1600 km'^2*')'",N == 350  ~ "N==350*' (r=53 km; s=3000 km'^2*')'",N == 200  ~ "N==200*' (r=70 km; s=5000 km'^2*')'"),
		N_dens = factor(N_dens,levels = c("N==200*' (r=70 km; s=5000 km'^2*')'","N==350*' (r=53 km; s=3000 km'^2*')'","N==600*' (r=40 km; s=1600 km'^2*')'","N==1000*' (r=32 km; s=1000 km'^2*')'","N==2000*' (r=22 km; s=500 km'^2*')'","N==3000*' (r=18 km; s=350 km'^2*')'"))) %>%
	ggplot() +
	# geom_errorbar(aes(x = sim_val, y = value, ymin = lo_q, ymax = up_q),
	# 							color = 'grey', width = 0.02) +
	geom_point(aes(x = sim_val, y = value), alpha = 0.4, size = 2) +
	geom_smooth(aes(x = sim_val, y = value), alpha = 0.4, linewidth = 0) +
	scale_color_manual(values = my_sp_colors)+
	scale_x_log10() +
	scale_y_log10() +
	geom_abline(intercept = 0, slope = 1, color = 'red', lty = 2, size = 0.5) +
	facet_wrap(~ N_dens, labeller = label_parsed) +
	labs(
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		y = expression("Estimated spatial autocorrelation - " * rho * " (km)")) +
	theme_bw() +
	my_theme()
# ggsave(here('Plots','Estimation of rho.jpg'),p1,width=16,height = 10)
p1

p1_1 <- est %>% 
	left_join(.,sim %>% select(-N,-source) %>%rename(sim_val='value'),
						by=c('parameter','rho','it')) %>% 
	filter(parameter=='rho') %>% 
	filter(!rho==1) %>% 
	mutate(value=value*100) %>% 
	mutate(sim_val=sim_val*100) %>% 
	mutate(lo_q=lo_q*100) %>% 
	mutate(up_q=up_q*100) %>%
	filter(N==350) %>% 
	mutate(
		N_dens = case_when(N == 350  ~ "N==350*' (r=53 km; s=3000 km'^2*')'"),
		N_dens = factor(N_dens,levels = c("N==350*' (r=53 km; s=3000 km'^2*')'"))) %>%
	ggplot() +
	geom_point(aes(x = sim_val, y = value), alpha = 0.2, size = 2) +
	geom_smooth(aes(x = sim_val, y = value), alpha = 0.4, size = 0) +
	geom_point(data=est_rho_sp %>% mutate(est_rho_y=est_rho*c(1,1.1,1.05,0.95,0.85,1.0,1,1.2,1,1,1,0.8)) %>% 
						 	mutate(N_dens="N==350*' (r=53 km; s=3000 km'^2*')'"),
						 aes(x=est_rho,y=est_rho_y,colour = sp),size=5)+
	scale_color_manual(values = my_sp_colors)+
	scale_x_log10() +
	scale_y_log10() +
	geom_abline(intercept = 0, slope = 1, color = 'red', lty = 2, size = 0.5,alpha=0.4) +
	facet_wrap(~ N_dens, labeller = label_parsed) +
	labs(
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		y = expression("Estimated spatial autocorrelation - " * rho * " (km)")) +
	theme_bw() +
	my_theme()+
	guides(shape = "none")
# ggsave(here('Plots','Estimation of rho with species.jpg'),p1_1,width=16,height = 10)


p2 <-
est %>% 
	left_join(.,sim %>% select(-N,-source) %>%rename(sim_val='value'),
						by=c('parameter','rho','it')) %>% 
	filter(!parameter%in%c('alpha','rho','sigma')) %>% 
	mutate(parameter = case_when(
		parameter == "mu_0"   ~ "0 m",
		parameter == "mu_50"  ~ "50 m",
		parameter == "mu_150" ~ "150 m",
		parameter == "mu_300" ~ "300 m",
		parameter == "mu_500" ~ "500 m"
	)) %>% 
	ggplot()+
	geom_abline(intercept = 0, slope = 1, color = 'red', lty = 2, size = 0.5) +
	geom_point(aes(x=sim_val,y=value),alpha=0.24)+
	labs(
		x = expression("Simulated global mean - " * mu * " (log copies/L)"),
		y = expression("Estimated global mean - " * mu * " (log copies/L)")
	) +
	facet_wrap(~parameter)+
	theme_bw()+
	my_theme()

# ggsave(here('Plots','Estimation of mu.jpg'),p2,width=16,height = 10)


sim_rds_files <- list.files(here('Output','GP_3000'),pattern = 'simulated_data.rds')

sim_3000 <- sim_rds_files %>%
	file.path(here("Output", "GP_3000"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble()

sim_2000 <- sim_rds_files %>%
	file.path(here("Output", "GP_2000"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble()

sim_1000 <- sim_rds_files %>%
	file.path(here("Output", "GP_1000"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble()

sim_600 <- sim_rds_files %>%
	file.path(here("Output", "GP_600"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble()

sim_350 <- sim_rds_files %>%
	file.path(here("Output", "GP_350"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble()

sim_200 <- sim_rds_files %>%
	file.path(here("Output", "GP_200"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble()

# Check that all points belong to the data frame previous to that (insuring thinning)
# sim_list <- list(sim_3000, sim_2000, sim_1000, sim_600, sim_350, sim_200)
# v <- as.data.frame(matrix(NA,15,20))
# l <- vector()
# for (k in 1:5) {
# 	for (i in 1:15) {
# 		for (j in 1:20) {
# 			v[i, j] <- sim_list[[k+1]] %>%
# 				filter(rho == i, it == j) %>%
# 				select(X_utm, Y_utm) %>%
# 				anti_join(
# 					sim_list[[k]] %>% filter(rho == i, it == j) %>% select(X_utm, Y_utm),
# 					by = c("X_utm", "Y_utm")) %>%
# 				nrow() == 0
			
# 		}
# 	}
# 	l[k] <- sum(colSums(v)==15)==20
# }
# l

est_rds_files <- list.files(here('Output','GP_3000'),pattern = 'pred_GP.rds')

pred_3000 <- est_rds_files %>%
	file.path(here("Output", "GP_3000"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=3000)

pred_2000 <- est_rds_files %>%
	file.path(here("Output", "GP_2000"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=2000)

pred_1000 <- est_rds_files %>%
	file.path(here("Output", "GP_1000"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=1000)

pred_600 <- est_rds_files %>%
	file.path(here("Output", "GP_600"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=600)

pred_350 <- est_rds_files %>%
	file.path(here("Output", "GP_350"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=350)

pred_200 <- est_rds_files %>%
	file.path(here("Output", "GP_200"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=200)

pred <- bind_rows(pred_2000 %>% cbind(pred_3000 %>% select(conc) %>% setNames('Conc_0')),
									pred_1000 %>% cbind(pred_3000 %>% select(conc) %>% setNames('Conc_0')),
									pred_600 %>% cbind(pred_3000 %>% select(conc) %>% setNames('Conc_0')),
									pred_350 %>% cbind(pred_3000 %>% select(conc) %>% setNames('Conc_0')),
									pred_200 %>% cbind(pred_3000 %>% select(conc) %>% setNames('Conc_0'))) %>% 
	mutate(diff=abs(Conc_0-conc))

pred_summ <- pred %>% 
	group_by(it,rho,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim %>% rename(rho_sim='value'),by=c('rho'))


p3 <- pred_summ %>% 
	filter(!rho%in%c(1,2)) %>% 
	mutate(mean_diff=exp(mean_diff)) %>%
	ggplot()+
	geom_point(aes(x=rho_sim,y=mean_diff,color=factor(N)))+
	geom_smooth(aes(x=rho_sim,y=mean_diff,color=factor(N)),se=F)+
	labs(
		y = expression("Absolute mean ratio of error - " * epsilon * " (" * frac(predicted, simulated) * ")"),
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		color='Number of\nsamples\ndeployed')+
	scale_x_log10()+
	scale_color_manual(values=my_colors)+
	theme_bw()+
	my_theme()

ggsave(here('Plots','Estimation of epsilon.jpg'),p3,width=16,height = 10)

mod <- pred_summ %>%
	filter(N == 350, !rho %in% c(1, 2)) %>%
	# lm(mean_diff ~ log(rho_sim), data = .) %>% 
	lm(mean_diff ~ I((log(rho_sim))^(-1.5)), data = .)

pred_summ %>% 
	filter(N==350) %>% 
	filter(!rho%in%c(1,2)) %>% 
	# mutate(pred=predict(lm(mean_diff ~ sqrt(sqrt(sqrt(log(rho_sim)))), data = .),newdata = tibble(rho_sim=rho_sim))) %>% 
	mutate(pred=predict(lm(mean_diff ~ I((log(rho_sim))^(-1.5)), data = .),newdata = tibble(rho_sim=rho_sim))) %>% 
	ggplot()+
	geom_point(aes(x=rho_sim,y=mean_diff))+
	geom_line(aes(x=rho_sim,y=pred))+
	theme_bw()+
	my_theme()
	

est_rho_sp_pred <- est_rho_sp %>% mutate(pred=predict(mod,newdata = tibble(rho_sim=est_rho)))

p4 <- pred_summ %>% 
	filter(!rho%in%c(1,2)) %>% 
	mutate(mean_diff=exp(mean_diff)) %>%
	ggplot()+
	geom_point(aes(x=N,y=mean_diff,color=log(rho_sim)),alpha=0.8)+
	geom_smooth(aes(x=N,y=mean_diff,group=rho_sim,color=log(rho_sim)),se=F,size=0.4)+
	labs(
		y = expression("Absolute mean ratio of error - " * epsilon * " (" * frac(predicted, simulated) * ")"),
		x = 'Number of samples deployed',
		color=expression("" * rho * " (log km)"))+
	scale_x_log10()+
	theme_bw()+
	my_theme()
ggsave(here('Plots','Estimation of epsilon_2.jpg'),p4,width=16,height = 10)

pred_summ %>% 
	filter(!rho%in%c(1,2)) %>% 
	mutate(mean_diff=exp(mean_diff)) %>%
	ggplot()+
	geom_point(data=est_rho_sp_pred %>% mutate(pred=exp(pred)),
						 aes(x=est_rho,y=pred))+
	geom_smooth(aes(x=rho_sim,y=mean_diff,color=factor(N)),se=F)+
	scale_x_log10()+
	scale_color_manual(values=my_colors)+
	theme_bw()+
	my_theme()

library(ggnewscale)

p5 <- pred_summ %>%
	filter(!rho %in% c(1, 2)) %>%
	mutate(mean_diff = exp(mean_diff)) %>%
	ggplot() +
	geom_point(aes(x = rho_sim, y = mean_diff, color = factor(N)),alpha = 0.3) +
	geom_smooth(aes(x = rho_sim, y = mean_diff, color = factor(N)),se = FALSE) +
	scale_color_manual(values = my_colors) +
	labs(
		y = expression("Absolute mean ratio of error - " * epsilon * " (" * frac(predicted, simulated) * ")"),
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		color='Number of samples\ndeployed')+
	ggnewscale::new_scale_color() +  # reset color scale
	
	geom_point(data = est_rho_sp_pred %>% mutate(pred = exp(pred)),aes(x = est_rho, y = pred, color = sp),size = 5) +
	scale_color_manual(values = my_sp_colors) +
	labs(color='Species')+
	scale_x_log10() +
	
	theme_bw() +
	my_theme()
ggsave(here('Plots','Estimation of epsilon with species.jpg'),p5,width=16,height = 10)
