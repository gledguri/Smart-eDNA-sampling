# Server code ---------------------------------------------------------------------------------

# srun --cpus-per-task=4 --mem=40G --time=14-00:00:00 --pty bash
# 
.libPaths(c('/mnt/nfs/home/KellyCEG/ggu/Rpack','/mnt/nfs/home/KellyCEG/R/x86_64-pc-linux-gnu-library/4.3',.libPaths()))

# Libraries -----------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(fields)
library(MASS)
library(PNWColors)
library(purrr)
library(tidyr)
library(tibble)
library(rstan);options(mc.cores = parallel::detectCores()); rstan_options(auto_write = TRUE)

my_theme <- function() {
	theme(
		axis.title = element_text(size=16),
		axis.text = element_text(size=15),
		strip.text = element_text(size = 15),
		legend.title = element_text(size = 15),
		legend.text  = element_text(size = 14))
}

select <- dplyr::select

rho_sim <- c(0.01,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20)
my_colors <- c('#6a3d9a',"#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d")
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

# GP pred 330 ---------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

for (h in 8:20) {
for (j in 2:15) {
name_file <- paste0(j,'_',h)
sim_data_raw <- readRDS(here('Output','GP_350',
														 paste0(name_file,'_simulated_data.rds')))


pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

sim_data <- sim_data_raw %>% 
	slice_sample(n=330) %>% 
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(sim_data,here('Output','GP_sth','GP_330',paste0(name_file,'_simulated_data.rds')))

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

write.csv(GP_param,file = here('Output','GP_sth','GP_330',paste0(name_file,'_est_GP_param.csv')))
# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_sth','GP_330',paste0(name_file,'_pred_GP.rds')))
}
}


# GP pred 300 ---------------------------------------------------------------------------------
pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

j=1
for (h in 11:20) {
for (j in 2:15) {
name_file <- paste0(j,'_',h)
sim_data_raw <- readRDS(here('Output','GP_sth','GP_330',
														 paste0(name_file,'_simulated_data.rds')))


pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

sim_data <- sim_data_raw %>% 
	slice_sample(n=300) %>% 
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(sim_data,here('Output','GP_sth','GP_300',paste0(name_file,'_simulated_data.rds')))

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
									refresh=500)

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output','GP_sth','GP_300',paste0(name_file,'_est_GP_param.csv')))
# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_sth','GP_300',paste0(name_file,'_pred_GP.rds')))
}
}

# GP pred 260 ---------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

j=1
for (h in 16:20) {
for (j in 2:15) {
name_file <- paste0(j,'_',h)
sim_data_raw <- readRDS(here('Output','GP_sth','GP_300',
														 paste0(name_file,'_simulated_data.rds')))


pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

sim_data <- sim_data_raw %>% 
	slice_sample(n=260) %>% 
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(sim_data,here('Output','GP_sth','GP_260',paste0(name_file,'_simulated_data.rds')))

# Priors
priors <- list(
	alpha_prior = c(4, 1),
	mu_prior = c(0, 5),
	rho_prior = c(0, 1),
	# mag_rho_prior = c(0,1),
	sigma_prior = c(0, 1)
)


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
									refresh=2000,
									open_progress = FALSE)

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output','GP_sth','GP_260',paste0(name_file,'_est_GP_param.csv')))
# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_sth','GP_260',paste0(name_file,'_pred_GP.rds')))
}
}

# GP pred 220 ---------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

for (h in 8:15) {
for (j in 2:15) {
name_file <- paste0(j,'_',h)
sim_data_raw <- readRDS(here('Output','GP_sth','GP_260',
														 paste0(name_file,'_simulated_data.rds')))


pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

sim_data <- sim_data_raw %>% 
	slice_sample(n=220) %>% 
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(sim_data,here('Output','GP_sth','GP_220',paste0(name_file,'_simulated_data.rds')))

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
									refresh=2000,
									open_progress = FALSE)

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output','GP_sth','GP_220',paste0(name_file,'_est_GP_param.csv')))
# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_sth','GP_220',paste0(name_file,'_pred_GP.rds')))
}
}

# GP pred 150 ---------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

j=2
for (h in 3:7) {
	for (j in 2:15) {
name_file <- paste0(j,'_',h)
sim_data_raw <- readRDS(here('Output','GP_sth','GP_220',
														 paste0(name_file,'_simulated_data.rds')))


pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

sim_data <- sim_data_raw %>% 
	slice_sample(n=150) %>% 
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(sim_data,here('Output','GP_sth','GP_150',paste0(name_file,'_simulated_data.rds')))

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
									refresh=500)

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output','GP_sth','GP_150',paste0(name_file,'_est_GP_param.csv')))
# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_sth','GP_150',paste0(name_file,'_pred_GP.rds')))
}
}

# GP pred 100 ---------------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))

for (h in 3:7) {
for (j in 2:15) {
name_file <- paste0(j,'_',h)
sim_data_raw <- readRDS(here('Output','GP_sth','GP_150',
														 paste0(name_file,'_simulated_data.rds')))


pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

sim_data <- sim_data_raw %>% 
	slice_sample(n=100) %>% 
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(sim_data,here('Output','GP_sth','GP_100',paste0(name_file,'_simulated_data.rds')))

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
									refresh=500)

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output','GP_sth','GP_100',paste0(name_file,'_est_GP_param.csv')))
# write.csv(sim_param,file = here('Output',paste0(name_file,'_sim_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output','GP_sth','GP_100',paste0(name_file,'_pred_GP.rds')))
}
}




# GP multifish species ------------------------------------------------------------------------

pred_data <- readRDS(here('Data','pred_data.rds'))
stan_model_7 <- stan_model(here('Code','GP_7.stan'))
edna_data <- readRDS(here('Data','edna_data.rds'))

spp <- edna_data %>% distinct(species) %>% slice(-n()) %>% pull(species)
thinning <- c(330,300,260,220,150,100)

for (i in 1:12) {
for (t in 1:6) {

obs_data_raw <- edna_data %>%
	filter(species==spp[i]) %>% 
	arrange(depth) %>% 
	mutate(x=X_utm/100,
				 y=Y_utm/100) %>% 
	arrange(depth,x,y)

obs_data <- obs_data_raw %>% rename(depth_cat='depth') %>% 
	slice_sample(n=thinning[t]) %>%
	mutate(depth=as.numeric(as.character(depth_cat))) %>% 
	arrange(depth,x,y)

saveRDS(obs_data,here('Output','GP_multifish',paste0(spp[i],'_',thinning[t],'_thinned_data.rds')))

	pred_data_by_depth <- pred_data %>% 
		rename(depth='depth_cat') %>% 
		mutate(depth=as.numeric(as.character(depth))) %>% 
		arrange(depth)
	
	# Priors 
	priors <- list(
		alpha_prior = c(4,1),
		mu_prior = c(0,5),
		rho_prior = c(0,1),
		# mag_rho_prior = c(0,1),
		sigma_prior = c(0,1))
	
	
	# Model stan 
	
	stan_data_7 <- list(
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
	
	stan_data_7 <- c(stan_data_7,priors);str(stan_data_7)
	
	fit_7 <- sampling(stan_model_7,
										data = stan_data_7,
										chains = 4,
										iter = 2000,
										warmup = 1000,
										refresh=500)
	
	param_list <- fit_7@model_pars
	param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]
	
	GP_param <- extract_param(fit_7,param_list) %>%
		as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)
	
	write.csv(GP_param,file = here('Output','GP_multifish',paste0(spp[i],'_',thinning[t],'_est_GP_param.csv')))
	# write.csv(GP_param,file = here('Output','GP_multifish',paste0(spp[i],'_est_GP_param.csv')))
	
	pred_data_stan <-
		pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
		rename(conc=mean)
	
	saveRDS(pred_data_stan, here('Output','GP_multifish',paste0(spp[i],'_',thinning[t],'_pred_GP.rds')))
	# saveRDS(pred_data_stan, here('Output','GP_multifish',paste0(spp[i],'_pred_GP.rds')))
}
}



# Analysis ------------------------------------------------------------------------------------

pred_350_sim <- list.files(here('Output','GP_350'),pattern = 'pred_GP.rds') %>%
	file.path(here("Output", "GP_350"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=350)

pred_330_sim <- list.files(here('Output','GP_sth','GP_330'),pattern = 'pred_GP.rds') %>%
	file.path(here("Output", "GP_sth",'GP_330'), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=330)

pred_300_sim <- list.files(here('Output','GP_sth','GP_300'),pattern = 'pred_GP.rds') %>%
	file.path(here("Output", "GP_sth",'GP_300'), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=300)

pred_260_sim <- list.files(here('Output','GP_sth','GP_260'),pattern = 'pred_GP.rds') %>%
	file.path(here("Output", "GP_sth",'GP_260'), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=260)

pred_220_sim <- list.files(here('Output','GP_sth','GP_220'),pattern = 'pred_GP.rds') %>%
	file.path(here("Output", "GP_sth",'GP_220'), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=220)

pred_150_sim <- list.files(here('Output','GP_sth','GP_150'),pattern = 'pred_GP.rds') %>%
	file.path(here("Output", "GP_sth",'GP_150'), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=150)

pred_100_sim <- list.files(here('Output','GP_sth','GP_100'),pattern = 'pred_GP.rds') %>%
	file.path(here("Output", "GP_sth",'GP_100'), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(rho = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho=as.numeric(rho),
				 it=as.numeric(it)) %>% select(-source) %>% 
	as_tibble() %>% 
	mutate(N=100)

pred_0_sim <- pred_350_sim %>% select(rho,it,X_utm,Y_utm,depth,conc) %>% rename(Conc_0='conc')

pred_sim <- bind_rows(
pred_330_sim %>% left_join(pred_0_sim,by=c('X_utm','Y_utm','depth','it','rho')),
pred_300_sim %>% left_join(pred_0_sim,by=c('X_utm','Y_utm','depth','it','rho')),
pred_260_sim %>% left_join(pred_0_sim,by=c('X_utm','Y_utm','depth','it','rho')),
pred_220_sim %>% left_join(pred_0_sim,by=c('X_utm','Y_utm','depth','it','rho')),
pred_150_sim %>% left_join(pred_0_sim,by=c('X_utm','Y_utm','depth','it','rho')),
pred_100_sim %>% left_join(pred_0_sim,by=c('X_utm','Y_utm','depth','it','rho'))
) %>% 
	mutate(diff=abs(Conc_0-conc)) %>% 
	filter(!it==1)

pred_summ_sim <- pred_sim %>% 
	group_by(it,rho,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim %>% 
							as.data.frame() %>% 
							setNames('rho_sim') %>% 
							mutate(rho=row_number(),
										 rho_sim=rho_sim*100),by=c('rho'))

spp

est_rho_sp <- list.files(here('Output','GP_multifish'),pattern = paste0('_350_est_GP_param.csv$')) %>%
  map_df(~ read_csv(here('Output','GP_multifish',.x)) %>%
           mutate(file = .x)) %>%
  separate(file,into = c('species','n_fish','type','est','ext1','ext2'),sep = '_',remove = FALSE) %>%
  select(-ext1,-ext2,-est,-type,-`...1`) |> 
  mutate(N = as.numeric(n_fish)) |> 
  filter(param=='rho') |> 
  group_by(species,param) %>%
  summarise(
    est_rho = mean(mean,na.rm=TRUE)) |> 
  select(-param)


pred_350_obs <- list.files(here('Output','GP_multifish'),pattern = paste0('350_pred_GP.rds')) %>% 
	file.path(here("Output", "GP_multifish"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(N = as.numeric(sub(".*_(\\d+)_pred_GP.*", "\\1", source)),
				 species = sub("_.*", "", source)) %>% 
	as_tibble()

pred_0_obs <- pred_350_obs %>% select(X_utm,Y_utm,depth,conc,species) %>% rename(Conc_0='conc') %>% 
	left_join(., est_rho_sp,by='species')

pred_330_obs <- list.files(here('Output','GP_multifish'),pattern = paste0('330_pred_GP.rds')) %>% 
	file.path(here("Output", "GP_multifish"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(N = as.numeric(sub(".*_(\\d+)_pred_GP.*", "\\1", source)),
				 species = sub("_.*", "", source)) %>% 
	as_tibble()

pred_300_obs <- list.files(here('Output','GP_multifish'),pattern = paste0('300_pred_GP.rds')) %>% 
	file.path(here("Output", "GP_multifish"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(N = as.numeric(sub(".*_(\\d+)_pred_GP.*", "\\1", source)),
				 species = sub("_.*", "", source)) %>% 
	as_tibble()

pred_260_obs <- list.files(here('Output','GP_multifish'),pattern = paste0('260_pred_GP.rds')) %>% 
	file.path(here("Output", "GP_multifish"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(N = as.numeric(sub(".*_(\\d+)_pred_GP.*", "\\1", source)),
				 species = sub("_.*", "", source)) %>% 
	as_tibble()

pred_220_obs <- list.files(here('Output','GP_multifish'),pattern = paste0('220_pred_GP.rds')) %>% 
	file.path(here("Output", "GP_multifish"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(N = as.numeric(sub(".*_(\\d+)_pred_GP.*", "\\1", source)),
				 species = sub("_.*", "", source)) %>% 
	as_tibble()

pred_150_obs <- list.files(here('Output','GP_multifish'),pattern = paste0('150_pred_GP.rds')) %>% 
	file.path(here("Output", "GP_multifish"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(N = as.numeric(sub(".*_(\\d+)_pred_GP.*", "\\1", source)),
				 species = sub("_.*", "", source)) %>% 
	as_tibble()

pred_100_obs <- list.files(here('Output','GP_multifish'),pattern = paste0('100_pred_GP.rds')) %>% 
	file.path(here("Output", "GP_multifish"), .) %>%
	map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>% 
	mutate(N = as.numeric(sub(".*_(\\d+)_pred_GP.*", "\\1", source)),
				 species = sub("_.*", "", source)) %>% 
	as_tibble()


pred_obs <- bind_rows(
pred_330_obs %>% left_join(pred_0_obs,by=c('X_utm','Y_utm','depth','species')),
pred_300_obs %>% left_join(pred_0_obs,by=c('X_utm','Y_utm','depth','species')),
pred_260_obs %>% left_join(pred_0_obs,by=c('X_utm','Y_utm','depth','species')),
pred_220_obs %>% left_join(pred_0_obs,by=c('X_utm','Y_utm','depth','species')),
pred_150_obs %>% left_join(pred_0_obs,by=c('X_utm','Y_utm','depth','species')),
pred_100_obs %>% left_join(pred_0_obs,by=c('X_utm','Y_utm','depth','species'))
) %>% 
	mutate(diff=abs(Conc_0-conc)) %>%

# Plotting species thinning from 350
pred_summ_obs <- pred_obs %>% 
  group_by(species, N, est_rho) %>%
  summarise(mean_diff = exp(mean(diff)), .groups = "drop") %>%
  mutate(
    rho = exp(est_rho) * 100,
    source = "Observed"
  )


pred_summ_sim <- pred_sim %>% rename(rho_it='rho') %>%
	group_by(it,rho_it,N) %>% 
	summarise(mean_diff=exp(mean(diff))) %>% 
	left_join(.,rho_sim %>% 
							as.data.frame() %>% 
							setNames('rho_sim') %>% 
							mutate(rho_it=row_number(),
										 rho=rho_sim*100),by=c('rho_it')) %>% 
  mutate(source = "Simulated")

pred_summ <- bind_rows(pred_summ_obs, pred_summ_sim)

# Plot both together
p1 <- pred_summ %>%
  ggplot(., 
	aes(x = rho, y = mean_diff, color = factor(N))) +
  geom_point() +
  geom_smooth(se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = my_colors) +
  labs(
    y = expression("Absolute mean ratio of error - " * epsilon * " (" * frac(predicted, simulated) * ")"),
    x = expression("Spatial autocorrelation - " * rho * " (km)"),
    color = "Number of\nsamples\ndeployed",
    shape = "Data source"
  ) +
  facet_wrap(~source) +
  theme_bw() +
  my_theme()

ggsave(here('Plots','thinning_analysis_plot.jpg'),p1,
	width = 14, height = 8)
