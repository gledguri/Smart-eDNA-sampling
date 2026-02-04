# Full simulation  ----------------------------------------------------------------------------
library(here)
library(dplyr)
library(ggplot2)
library(tibble)
library(MoMAColors)
library(PNWColors)
library(tidyr)
library(plotly)
library(sf)
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

# Import files --------------------------------------------------------------------------------
cond_raw <- expand.grid(a = c(6), d = c(1,1.5,3,4,6), n = c(10,20,30))
cond_raw_2 <- 
	cond_raw %>% 
	mutate(d=d/sqrt(100),
				 n=n/sqrt(100))

files <- here('Data','sim_out_partitioned',here('Data','sim_out_partitioned') %>% list.files(pattern = 'pred_grid_list'))
# files_1 <- files[c(1:13,23,25:30)]
files_1 <- files[c(1:13,23,34,36:40)]

pred_grid_list_1 <- files_1 %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c() %>% bind_rows(,.id = 'id') %>% 
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

# files_2 <- files[c(14:22,24)]
files_2 <- files[c(14:22,24:33,35)]
pred_grid_list_2 <- files_2 %>%
	map(readRDS) %>%
	map(compact) %>%
	list_c() %>% bind_rows(,.id = 'id') %>% 
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

deviation <-
	pred_grid_list %>%
	left_join(.,bind_rows(replicate(40, cond_raw, simplify = FALSE)) %>% rowid_to_column('id'),
						by='id') %>%
	mutate(N_ini=n^2,
				 N=(n^2)-(th*(n^2))) %>%
	mutate(point_density_sim_gp=N/100) %>%
	group_by(thinning,a,d,n) %>%
	mutate(mean_deviation=mean(dev)) %>%
	as.data.frame()

# # write.csv(deviation,here('Data','error_df.csv'),row.names = F)

p1 <- deviation %>%
	# filter(N_ini==100) %>% 
	ggplot()+
	# geom_point(aes(x=th,y=dev,color=as.factor(d)))+
	# geom_line(aes(x=th,y=mean_deviation,color=as.factor(d)))+
	geom_point(aes(x=th,y=dev,color=(d)))+
	geom_line(aes(x=th,y=mean_deviation,group=d,color=(d)))+
	scale_color_gradientn(colors = rev(pnw_palette("Bay", 11, type = "continuous")))+
	scale_x_continuous(breaks=seq(0.1,0.8,by=0.1))+
	# geom_point(aes(x=d,y=dev,color=as.factor(th)))+
	# geom_line(aes(x=d,y=mean_deviation,color=as.factor(th)))+
	# scale_x_continuous(breaks=c(1,2,3,4,5,6))+
	facet_wrap(~N_ini)+
	# geom_point(aes(x=N_ini,y=dev,color=as.factor(d)))+
	# geom_line(aes(x=N_ini,y=mean_deviation,color=as.factor(d)))+
	# scale_x_continuous(breaks=c(1,2,3,4,5,6))+
	# facet_wrap(~th)+
	ylab('Precision term (ε)')+
	scale_y_log10()+
	theme_bw()
# 
# 
# # Plotly 3D plot ------------------------------------------------------------------------------
# 
# 
# dd <- deviation %>%
# 	group_by(th,d,N_ini) %>%
# 	mutate(mean_dev=mean(dev)) %>%
# 	mutate(N_ini = as.factor(N_ini))
# 
# plot_ly(
# 	data = dd,
# 	x = ~th,
# 	y = ~d,
# 	# z = ~log(dev),
# 	z = ~log(mean_dev),
# 	color = ~N_ini,
# 	colors = "Set1",  # You can change the palette
# 	type = "scatter3d",
# 	mode = "markers"
# ) %>%
# 	layout(
# 		scene = list(
# 			xaxis = list(title = "T"),
# 			yaxis = list(title = "ρ"),
# 			zaxis = list(title = "ε")
# 		)
# 	)
# 
# 
# # Stan predicting the relationship ------------------------------------------------------------
# 
# stan_data <- list(
# 	N = nrow(deviation),
# 	dev = deviation$dev,
# 	th = deviation$th,
# 	d = deviation$d,
# 	N_ini = deviation$N_ini
# )
# 
# # Compile and run the Stan model
# # Save the Stan model code above as "regression_model.stan"
# fit <- stan(file = here('Code','Error_relationship_with_d_4.stan'), 
# 						data = stan_data,
# 						chains = 4,
# 						iter = 2000,
# 						warmup = 1000,
# 						cores = 4)
# 
# extract_param(fit,'y_pred') %>% 
# 	dplyr::select(mean) %>% rename(pred_dev='mean') %>% 
# 	cbind(.,deviation) %>% 
# 	ggplot()+
# 	# geom_line(aes(x=th,y=pred_dev,color=as.factor(d)))+
# 	# geom_point(aes(x=th,y=dev,color=as.factor(d)))+
# 	# # geom_line(aes(x=th,y=mean_deviation,color=as.factor(d)))+
# 	# scale_x_continuous(breaks=seq(0.1,0.8,by=0.1))+
# 	geom_point(aes(x=d,y=dev,color=as.factor(th)))+
# 	geom_line(aes(x=d,y=pred_dev,color=as.factor(th)))+
# 	# geom_line(aes(x=d,y=mean_deviation,color=as.factor(th)))+
# 	scale_x_continuous(breaks=c(1,2,3,4,5,6))+
# 	facet_wrap(~N_ini)+
# 	# geom_point(aes(x=d,y=pred_dev,color=as.factor(N_ini)))+
# 	# # geom_point(aes(x=d,y=dev,color=as.factor(N_ini)))+
# 	# geom_line(aes(x=d,y=mean_deviation,color=as.factor(N_ini)))+
# 	# scale_x_continuous(breaks=c(1,2,3,4,5,6))+
# 	# facet_wrap(~th)+
# 	ylab('Precision term (ε)')+
# 	scale_y_log10()+
# 	# scale_x_log10()+
# 	theme_bw()



# Model 2 -------------------------------------------------------------------------------------
deviation <-
	pred_grid_list %>% 
	left_join(.,bind_rows(replicate(30, cond_raw_2, simplify = FALSE)) %>% rowid_to_column('id'),
						by='id') %>% 
	# mutate(N_ini=n^2,
	# 			 N=(n^2)-(th*(n^2))) %>% 
	# mutate(point_density_sim_gp=N/100) %>% 
	group_by(thinning,a,d,n) %>%
	mutate(mean_deviation=mean(dev)) %>% 
	as.data.frame()

# write.csv(deviation,here('Data','error_df.csv'),row.names = F)

deviation %>% 
	ggplot()+
	geom_point(aes(x=th,y=dev,color=as.factor(d)))+
	# geom_line(aes(x=th,y=mean_deviation,color=as.factor(d)))+
	scale_x_continuous(breaks=seq(0.1,0.8,by=0.1))+
	# geom_point(aes(x=d,y=dev,color=as.factor(th)))+
	# geom_line(aes(x=d,y=mean_deviation,color=as.factor(th)))+
	# scale_x_continuous(breaks=c(1,2,3,4,5,6))+
	facet_wrap(~n)+
	# geom_point(aes(x=N_ini,y=dev,color=as.factor(d)))+
	# geom_line(aes(x=N_ini,y=mean_deviation,color=as.factor(d)))+
	# scale_x_continuous(breaks=c(1,2,3,4,5,6))+
	# facet_wrap(~th)+
	ylab('Precision term (ε)')+
	scale_y_log10()+
	theme_bw()

stan_data <- list(
	N = nrow(deviation),
	dev = deviation$dev,
	th = deviation$th,
	d = deviation$d,
	N_ini = deviation$n
)

# Compile and run the Stan model
# Save the Stan model code above as "regression_model.stan"
fit <- stan(file = here('Code','Error_relationship_with_d_5.stan'), 
						data = stan_data,
						chains = 4,
						iter = 2000,
						warmup = 1000,
						cores = 4)

extract_param(fit,'y_pred') %>% 
	dplyr::select(mean) %>% rename(pred_dev='mean') %>% 
	cbind(.,deviation) %>% 
	ggplot()+
	# geom_line(aes(x=th,y=pred_dev,color=as.factor(d)))+
	# geom_point(aes(x=th,y=dev,color=as.factor(d)))+
	# # geom_line(aes(x=th,y=mean_deviation,color=as.factor(d)))+
	# scale_x_continuous(breaks=seq(0.1,0.8,by=0.1))+
	geom_point(aes(x=d,y=dev,color=as.factor(th)))+
	geom_line(aes(x=d,y=pred_dev,color=as.factor(th)))+
	# geom_line(aes(x=d,y=mean_deviation,color=as.factor(th)))+
	scale_x_continuous(breaks=c(1,2,3,4,5,6))+
	facet_wrap(~n)+
	# geom_point(aes(x=d,y=pred_dev,color=as.factor(N_ini)))+
	# # geom_point(aes(x=d,y=dev,color=as.factor(N_ini)))+
	# geom_line(aes(x=d,y=mean_deviation,color=as.factor(N_ini)))+
	# scale_x_continuous(breaks=c(1,2,3,4,5,6))+
	# facet_wrap(~th)+
	ylab('Precision term (ε)')+
	scale_y_log10()+
	# scale_x_log10()+
	theme_bw()

extract_param(fit,c('alpha','beta_th','beta_d','beta_N_ini',
										'beta_th_d','beta_th_N','beta_d_N','beta_th_d_N')) %>% 
	rownames_to_column('param') %>%
	ggplot()+
	geom_point(aes(y=param,x=mean))+
	geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
	theme_bw()
	

# Import multifish data -----------------------------------------------------------------------

multifish <- readRDS(here('Data','Multifish_data','Log_D_est.rds'))
pred_data_all <- readRDS(here('Data','Multifish_data','Log_D_est_smoothed.rds'))

multifish <- multifish %>% 
	bind_cols(.,
	st_as_sf(multifish %>% select(lon,lat), coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = 32610) %>% 
	st_coordinates(multifish_utm) %>% 
	as.data.frame() %>% 
	mutate(X_utm=X/1000) %>% 
	mutate(Y_utm=Y/1000))

# sp <- 'Engraulis mordax'
sp <- 'Sardinops sagax'
edna_data <-
	multifish %>% 
		filter(species==sp) %>% 
		# filter(depth%in%c(0,50)) %>%
		select(mean,Y_utm,X_utm,depth) %>% 
	mutate(depth=if_else(depth<50,0,depth)) %>% 
	mutate(depth=if_else(depth<150&depth>49,50,depth)) %>% 
	mutate(depth=if_else(depth<300&depth>149,150,depth)) %>% 
	mutate(depth=if_else(depth<500&depth>299,300,depth)) %>% 
	mutate(depth=if_else(depth<500&depth>300,300,depth)) 

# pred_data <- pred_data_all[[1]] %>% select(utm.lon.km,utm.lat.km) %>% 
# 	as.data.frame() %>% select(-x) %>% rename(Y_utm='utm.lat.km',
# 																						X_utm='utm.lon.km') %>% 
# 	slice_sample(n = 4000)
pred_data <- pred_data_all[[1]] %>% select(utm.lon.km,utm.lat.km,depth_cat) %>%
	as.data.frame() %>% select(-x) %>% rename(Y_utm='utm.lat.km',
																						X_utm='utm.lon.km') %>%
	group_by(depth_cat) %>%
  slice_sample(n = 4000) %>%
  ungroup()


# Model 1 -------------------------------------------------------------------------------------

	
clup_0 <- edna_data %>% filter(depth%in%c(0)) 



# clup_0_stan_data_1 <- list(
# 	N = nrow(clup_0),
# 	X = cbind(clup_0$X_utm, clup_0$Y_utm),
# 	y = clup_0$mean,
# 	N_pred = nrow(pred_data),
# 	X_pred = cbind(pred_data$X_utm, pred_data$Y_utm),
# 	alpha_prior = c(0,1),
# 	rho_prior = c(0,2),
# 	sigma_prior = c(4,0.1),
# 	mu_prior = c(0,2)
# )
# 	
# pred_data_stan <- pred_data %>% 
# 	bind_cols(.,extract_param(fit_clup_0,'f_pred') %>% 
# 							select(mean) %>% 
# 							rename(conc_pred='mean'))
# 
# stan_model <- stan_model(here('Code','GP.stan'))

# fit_clup_0 <- sampling(stan_model,
# 											 data = clup_0_stan_data,
# 											 chains = 4,
# 											 iter = 2000,
# 											 warmup = 1000)
#  
# # library(PNWColors)
# 
# clup_0 %>% 
# 	ggplot() +
# 	geom_point(data = pred_data_stan, aes(x = X_utm, y = Y_utm, colour = conc_pred)) +
# 	# geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
# 	scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous")) +
# 	coord_equal() +
# 	theme_bw()


# Model 2 -------------------------------------------------------------------------------------

clup_50 <- clup %>% filter(depth%in%c(50)) 

clup_stan_data_2 <- list(
	N_0 = nrow(clup_0),
	N_50 = nrow(clup_50),
	X_0 = cbind(clup_0$X_utm, clup_0$Y_utm),
	X_50 = cbind(clup_50$X_utm, clup_50$Y_utm),
	y_0 = clup_0$mean,
	y_50 = clup_50$mean,
	N_pred = nrow(pred_data),
	X_pred = cbind(pred_data$X_utm, pred_data$Y_utm),
	# alpha_prior = c(5,1),
	# rho_prior = c(0,1),
	# sigma_prior = c(4,1),
	# mu_prior = c(0,2)
	alpha_prior = c(2,1),
	rho_prior = c(0,1),
	sigma_prior = c(4,1),
	mu_prior = c(-2,2)
)

# stan_model_2 <- stan_model(here('Code','GP_2.stan'))

fit_clup_2 <- sampling(stan_model_2,
											 data = clup_stan_data_2,
											 chains = 4,
											 iter = 2000,
											 warmup = 1000)

extract_param(fit_clup_2,c('alpha','rho','sigma','mu','lp__'))

pred_data_stan_2 <-
pred_data %>% bind_cols(.,
	extract_param(fit_clup_2,'f_pred_0') %>% select(mean) %>% rename(conc_pred_0='mean'),
	extract_param(fit_clup_2,'f_pred_50') %>% select(mean) %>% rename(conc_pred_50='mean')) %>% 
	as.data.frame() %>% 
	pivot_longer(cols = conc_pred_0:conc_pred_50,
							 names_to = 'depth',
							 values_to = 'conc') %>% 
	mutate(depth=gsub('conc_pred_','',depth)) %>% 
	mutate(depth=as.numeric(depth))


# library(PNWColors)

clup %>% filter(depth%in%c(0,50)) %>% 
	ggplot() +
	geom_point(data = pred_data_stan_2, aes(x = X_utm, y = Y_utm, colour = conc)) +
	# geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
	scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous")) +
	facet_grid(~depth)+
	coord_equal() +
	theme_bw()


# Model 3 -------------------------------------------------------------------------------------

obs_data <- edna_data %>%
	# filter(depth%in%c(0,50)) %>%
	arrange(depth)

pred_data_by_depth <- pred_data %>% 
	rename(depth='depth_cat') %>% 
	mutate(depth=as.numeric(as.character(depth))) %>% 
	arrange(depth)

stan_data_3 <- list(
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
	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0))),
	alpha_prior = c(5,1),
	mu_prior = c(0,5),
	rho_prior = c(0,4.5),
	mag_rho_prior = c(0,4.5),
	sigma_prior = c(3,1)
);str(stan_data_3)

# stan_model_3 <- stan_model(here('Code','GP_3_CEG.stan'))
# stan_model_3 <- stan_model(here('Code','GP_3.stan'))

fit_clup_3 <- sampling(stan_model_3,
										 data = stan_data_3,
										 chains = 4,
										 iter = 2000,
										 warmup = 1000)

# extract_param(fit_clup_2,c('alpha','rho','sigma','mu','lp__')) %>% 
# 	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
# 	filter(param!='lp__') %>%
# 	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
# 	xlim(-5,7)+theme_bw()
# 
extract_param(fit_clup_3,c('lp__'))
extract_param(fit_clup_3,c('rho_sd'))

p7 <- extract_param(fit_clup_3,c('alpha','rho','sigma','mu','lp__')) %>%
	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
	filter(param!='lp__') %>%
	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
	theme_bw()

pred_data_stan_3 <-
	pred_data_by_depth %>% bind_cols(.,extract_param(fit_clup_3,c('y_pred'))) %>% 
		rename(conc=mean)

obs_data %>% #filter(depth%in%c(0,50)) %>% 
	ggplot() +
	geom_point(data = pred_data_stan_3, aes(x = X_utm, y = Y_utm, colour = conc)) +
	# geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
	scale_color_gradientn(colors = pnw_palette("Bay", 18, type = "continuous")) +
	facet_grid(~depth)+
	coord_equal() +
	theme_bw()

# Model 4 -------------------------------------------------------------------------------------
# 
# clup_050 <- clup %>%
# 	filter(depth%in%c(0,50)) %>%
# 	arrange(depth) 
# 
# pred_data_by_depth <- pred_data %>% #slice_sample(n=2000) %>% 
# 	expand_grid(depth = unique(clup_050$depth)) %>%
# 	arrange(depth)
# 
# clup_stan_data_4 <- list(
# 	N_total = nrow(clup_050),
# 	N_depths = nrow(clup_050 %>% distinct(depth)),
# 	X = cbind(clup_050$X_utm, clup_050$Y_utm),
# 	y = clup_050$mean,
# 	des_mat = model.matrix(~ depth - 1, data = clup_050 %>% mutate(depth=as.factor(depth))),
# 	d_idx = make_index(clup_050,'depth','d_idx') %>% pull(d_idx),
# 	N_by_depth = clup_050 %>% count(depth) %>% pull(n),
# 	start_idx = c(1,which(c(FALSE, diff(clup_050$depth) != 0))),
# 	# Pred
# 	N_pred = nrow(pred_data_by_depth),
# 	X_pred = cbind(pred_data_by_depth$X_utm, pred_data_by_depth$Y_utm),
# 	N_depths_pred = nrow(pred_data_by_depth %>% distinct(depth)),
# 	N_by_depth_pred = pred_data_by_depth %>% count(depth) %>% pull(n),
# 	start_idx_pred = c(1,which(c(FALSE, diff(pred_data_by_depth$depth) != 0))),
# 	alpha_prior = c(5,1),
# 	mu_prior = c(0,1),
# 	rho_prior = c(0,1),
# 	sigma_prior = c(0,1)
# )
# 
# # stan_model_3 <- stan_model(here('Code','GP_4.stan'))
# 
# fit_clup_4 <- sampling(stan_model_4,
# 										 data = clup_stan_data_4,
# 										 chains = 4,
# 										 iter = 2000,
# 										 warmup = 1000)
# 
# # extract_param(fit_clup_2,c('alpha','rho','sigma','mu','lp__')) %>% 
# # 	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
# # 	filter(param!='lp__') %>%
# # 	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
# # 	xlim(-5,7)+theme_bw()
# # 
# extract_param(fit_clup_4,c('alpha','rho','sigma','mu','lp__')) %>%
# 	as.data.frame() %>% rownames_to_column('param') %>% select(param,mean,`2.5%`,`97.5%`) %>%
# 	filter(param!='lp__') %>%
# 	ggplot()+geom_point(aes(x=mean,y=param))+geom_errorbar(aes(y=param,xmin=`2.5%`,xmax=`97.5%`),width=0.2)+
# 	theme_bw()
# 
# pred_data_stan_4 <-
# 	pred_data_by_depth %>% bind_cols(.,extract_param(fit_clup_4,c('y_pred'))) %>% 
# 		rename(conc=mean)
# 
# clup %>% #filter(depth%in%c(0,50)) %>% 
# 	ggplot() +
# 	geom_point(data = pred_data_stan_4, aes(x = X_utm, y = Y_utm, colour = conc)) +
# 	geom_point(aes(x = X_utm, y = Y_utm, colour = mean),size=4) +
# 	scale_color_gradientn(colors = pnw_palette("Bay", 8, type = "continuous")) +
# 	facet_grid(~depth)+
# 	coord_equal() +
# 	theme_bw()
