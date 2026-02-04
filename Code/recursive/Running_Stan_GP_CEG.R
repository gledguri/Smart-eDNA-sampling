.libPaths(c('/mnt/nfs/home/KellyCEG/ggu/Rpack','/mnt/nfs/home/KellyCEG/R/x86_64-pc-linux-gnu-library/4.3',.libPaths()))

library(rstan)
library(dplyr)
library(tidyr)
library(fields)
library(MASS)
library(tibble)
library(here)
library(purrr)
library(ggplot2)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

extract_param <- function (model, par) {
	fit <- (methods::selectMethod("summary", signature = "stanfit"))(object = model, par = par)
	fit <- fit$summary
	return(fit %>% unlist() %>% as.data.frame %>% round(., 9))}

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
# 
# st_iter <- 136
# end_iter <- 150
# 
# source('Stan_GP_CEG.R')
# # source(here('Code','Stan_GP_CEG.R'))

# Run this individually in each screen
saveRDS(stan_param_list,'stan_param_list_10.RDS')
saveRDS(pred_grid_list,  'pred_grid_list_10.RDS')
saveRDS(g_sigma_list,      'g_sigma_list_10.RDS')
saveRDS(mu_sim_list,        'mu_sim_list_10.RDS')
saveRDS(coords_raw_list,'coords_raw_list_10.RDS')
saveRDS(coords_10_list,  'coords_10_list_10.RDS')
saveRDS(coords_30_list,  'coords_30_list_10.RDS')
saveRDS(coords_50_list,  'coords_50_list_10.RDS')
saveRDS(coords_80_list,  'coords_80_list_10.RDS')

# After saving each rds pull them together
files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'stan_param_list'))
stan_param_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()


files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'pred_grid_list'))
files_1 <- files[c(1:13,23,25:30)]
files_2 <- files[c(14:22,24)]
pred_grid_list_1 <- files_1 %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
pred_grid_list_2 <- files_2 %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()

pred_grid_list_1 <- pred_grid_list_1 %>% bind_rows(,.id = 'id') %>% 
	group_by(id) %>%
	summarize(deviation_10=mean(abs(Conc_thin_10-Conc)),
						deviation_30=mean(abs(Conc_thin_30-Conc)),
						deviation_50=mean(abs(Conc_thin_50-Conc)),
						deviation_80=mean(abs(Conc_thin_80-Conc))) %>% 
	mutate(id=as.numeric(id)) %>% arrange(id) %>% 
	pivot_longer(cols= -id,
							 names_to = 'thinning',
							 values_to = 'dev') %>% 
	mutate(th=gsub('deviation_','',thinning)) %>% 
	mutate(th=as.numeric(th)/100)

pred_grid_list_2 <- pred_grid_list_2 %>% bind_rows(,.id = 'id') %>% 
	rename(Conc_thin_20='Conc_thin_10') %>% 
	rename(Conc_thin_40='Conc_thin_30') %>% 
	rename(Conc_thin_60='Conc_thin_50') %>% 
	rename(Conc_thin_70='Conc_thin_80') %>% 
	group_by(id) %>%
	summarize(deviation_20=mean(abs(Conc_thin_20-Conc)),
						deviation_40=mean(abs(Conc_thin_40-Conc)),
						deviation_60=mean(abs(Conc_thin_60-Conc)),
						deviation_70=mean(abs(Conc_thin_70-Conc))) %>% 
	mutate(id=as.numeric(id)) %>% arrange(id) %>% 
	pivot_longer(cols= -id,
							 names_to = 'thinning',
							 values_to = 'dev') %>% 
	mutate(th=gsub('deviation_','',thinning)) %>% 
	mutate(th=as.numeric(th)/100) %>% 
	mutate(id=id+max(pred_grid_list_1$id))


pred_grid_list <- pred_grid_list_1 %>% 
	rbind(.,pred_grid_list_2)

files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'g_sigma_list'))
files <- files[c(14:22,24)]
g_sigma_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
# for (i in c(1,12,14:20,2:11,13)) {
# 	v <- c(1,12,14:20,2:11,13)
# 	j <- which(v==i)
# 	l <- readRDS(files[[i]])
# 	g_sigma_list[(15*(j-1)+1):(15*(j))] <- l[(15*(j-1)+1):(15*(j))]
# }

files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'mu_sim_list'))
mu_sim_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
# for (i in c(1,12,14:20,2:11,13)) {
# 	v <- c(1,12,14:20,2:11,13)
# 	j <- which(v==i)
# 	l <- readRDS(files[[i]])
# 	mu_sim_list[(15*(j-1)+1):(15*(j))] <- l[(15*(j-1)+1):(15*(j))]
# }

files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'coords_raw_list'))
coords_raw_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
# for (i in c(1,12,14:20,2:11,13)) {
# 	v <- c(1,12,14:20,2:11,13)
# 	j <- which(v==i)
# 	l <- readRDS(files[[i]])
# 	coords_raw_list[(15*(j-1)+1):(15*(j))] <- l[(15*(j-1)+1):(15*(j))]
# }

files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'coords_10_list'))
coords_10_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
# for (i in c(1,12,14:20,2:11,13)) {
# 	v <- c(1,12,14:20,2:11,13)
# 	j <- which(v==i)
# 	l <- readRDS(files[[i]])
# 	coords_10_list[(15*(j-1)+1):(15*(j))] <- l[(15*(j-1)+1):(15*(j))]
# }

files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'coords_30_list'))
coords_30_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
# for (i in c(1,12,14:20,2:11,13)) {
# 	v <- c(1,12,14:20,2:11,13)
# 	j <- which(v==i)
# 	l <- readRDS(files[[i]])
# 	coords_30_list[(15*(j-1)+1):(15*(j))] <- l[(15*(j-1)+1):(15*(j))]
# }


files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'coords_50_list'))
coords_50_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
# for (i in c(1,12,14:20,2:11,13)) {
# 	v <- c(1,12,14:20,2:11,13)
# 	j <- which(v==i)
# 	l <- readRDS(files[[i]])
# 	coords_50_list[(15*(j-1)+1):(15*(j))] <- l[(15*(j-1)+1):(15*(j))]
# }

files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'coords_80_list'))
coords_80_list <- files %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c()
# for (i in c(1,12,14:20,2:11,13)) {
# 	v <- c(1,12,14:20,2:11,13)
# 	j <- which(v==i)
# 	l <- readRDS(files[[i]])
# 	coords_80_list[(15*(j-1)+1):(15*(j))] <- l[(15*(j-1)+1):(15*(j))]
# }

# Save the pooled list
saveRDS(stan_param_list,here('Data','full_sim','stan_param_list.RDS'))
saveRDS(pred_grid_list,  here('Data','full_sim','pred_grid_list.RDS'))
saveRDS(g_sigma_list,      here('Data','full_sim','g_sigma_list.RDS'))
saveRDS(mu_sim_list,        here('Data','full_sim','mu_sim_list.RDS'))
saveRDS(coords_raw_list,here('Data','full_sim','coords_raw_list.RDS'))
saveRDS(coords_10_list,  here('Data','full_sim','coords_10_list.RDS'))
saveRDS(coords_30_list,  here('Data','full_sim','coords_30_list.RDS'))
saveRDS(coords_50_list,  here('Data','full_sim','coords_50_list.RDS'))
saveRDS(coords_80_list,  here('Data','full_sim','coords_80_list.RDS'))
