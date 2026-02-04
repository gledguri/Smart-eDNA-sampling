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

# sp <- 'Engraulis mordax'
# sp <- 'Sardinops sagax'
# sp <- 'Sebastes entomelas'
# sp <- 'Stenobrachius leucopsarus'
# sp <- 'Tarletonbeania crenularis'

edna_data <- readRDS(here('Data','edna_data.rds'))

pred_data <- readRDS(here('Data','pred_data.rds'))

stan_model_7 <- stan_model(here('Code','GP_7.stan'))

obs_data_raw <- edna_data %>%
	filter(species==sp) %>% 
	# filter(depth%in%c(0,50)) %>%
	arrange(depth) %>% 
	mutate(x=X_utm/100,
				 y=Y_utm/100)

pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth) %>% 
	mutate(x=X_utm/100,
				 y=Y_utm/100)

for (y in 1:2) {

sample_100 <- c(1:nrow(obs_data_raw))
sample_90 <- sample(1:nrow(obs_data_raw), nrow(obs_data_raw)*0.9,replace = F)
sample_80 <- sample(sample_90, nrow(obs_data_raw)*0.8,replace = F)
sample_70 <- sample(sample_80, nrow(obs_data_raw)*0.7,replace = F)
sample_60 <- sample(sample_70, nrow(obs_data_raw)*0.6,replace = F)
sample_50 <- sample(sample_60, nrow(obs_data_raw)*0.5,replace = F)

thinning_it <- seq(1,0.5,by=-0.1)
samp_lsit <- list(sample_100,sample_90,sample_80,sample_70,sample_60,sample_50)

for (oo in 1:6) {
	samp_size=samp_lsit[[oo]]
	th=thinning_it[oo]
	
	obs_data <- obs_data_raw %>% 
		slice(samp_size) %>% 
		arrange(depth,X_utm,Y_utm)

saveRDS(obs_data,here('Output','multifish_spp',paste0(sp,'_',y,'_th',th,'_real_data.rds')))
	
# Priors --------------------------------------------------------------------------------------
priors <- list(
	alpha_prior = c(4,1),
	mu_prior = c(0,5),
	rho_prior = c(0,1),
	sigma_prior = c(0,1))

# Model stan -------------------------------------------------------------------------------------

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

# Capture warnings specifically
warning_file <- here('Output','multifish_spp',paste0(sp,'_',y,'_th',th,'_stan_warnings.txt'))

# Create a custom warning handler
warnings_captured <- character(0)

withCallingHandlers({
	fit_7 <- sampling(stan_model_7,
										data = stan_data_7,
										chains = 4,
										iter = 2000,
										warmup = 1000)
}, warning = function(w) {
	warnings_captured <<- c(warnings_captured, conditionMessage(w))
	writeLines(conditionMessage(w), warning_file, sep = "\n")
	invokeRestart("muffleWarning")
})

param_list <- fit_7@model_pars
param_list <- param_list[-which(param_list%in%c('y_pred','lp__'))]

GP_param <- extract_param(fit_7,param_list) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`25%`,`50%`,`75%`,`97.5%`)

write.csv(GP_param,file = here('Output', 'multifish_spp', paste0(sp,'_',y,'_th',th,'_est_GP_param.csv')))

pred_data_stan <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_7,c('y_pred'))) %>% 
	rename(conc=mean)

saveRDS(pred_data_stan, here('Output', 'multifish_spp', paste0(sp,'_',y,'_th',th,'_pred_GP.rds')))
}
}