
# Libraries --------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(here)
library(fields)
library(MASS)
library(PNWColors)
library(purrr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(spatstat.geom)
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

my_theme <- function() {
	theme(
		axis.title = element_text(size=16),
		axis.text = element_text(size=15),
		strip.text = element_text(size = 15),
		legend.title = element_text(size = 15),
		legend.text  = element_text(size = 14))
}



compute_NNI <- function(data, x, y, window) {
	x_vec <- dplyr::pull(data, {{ x }})
	y_vec <- dplyr::pull(data, {{ y }})
	
	if (inherits(window, "owin")) {win <- window} else if (is.numeric(window) && length(window) == 4) {
		win <- owin(xrange = c(window[1], window[2]),
								yrange = c(window[3], window[4]))} else {stop("`window` must be either an 'owin' object or numeric c(xmin, xmax, ymin, ymax).")}
	
	# Create point pattern
	pp <- ppp(x_vec, y_vec, window = win)
	
	# Observed mean nearest-neighbor distance
	nn <- nndist(pp)
	observed <- mean(nn)
	
	# Intensity (points per unit area)
	lambda <- pp$n / area.owin(pp$window)
	
	# Expected mean distance under CSR
	expected <- 1 / (2 * sqrt(lambda))
	
	# Nearest Neighbor Index
	NNI <- observed / expected
	
	NNI
}




# Import data ---------------------------------------------------------------------------------

## Declare variables ---------------------------------------------------------------------------------
my_colors <- c("#046dedff","#007fa2ff","#33a02c","#ff7f01", "#e31b1d")
my_sp_colors <- c('#6a3d9a','#cab2d6','#ff7f01','#fdbf6f','#e31b1d','#fb9a99','#33a02c','#b2df8a','#1f78b4','#009999','#999900', '#a6cee3')
# my_colors <- c("#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d")
my_colors_2 <- c("#0660f1ff","#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d")
N_values <- c(3000, 2000, 1000, 600, 350, 200)
N_values_2 <- c(350, 330, 300, 260, 220, 150, 100)
rho_sim <- c(0.01,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20) # Length scales parameter to be simulated
rho_sim_df <- rho_sim %>% as.data.frame() %>% setNames('rho') %>% mutate(rho_it=row_number(.)) 
edna_data <- readRDS(here('Data','edna_data.rds'))
pred_data <- readRDS(here('Data','pred_data.rds'))
spp <- edna_data %>% distinct(species) %>% slice(-n()) %>% pull(species)
trophic <- c('Surface','Surface','Midwater','Deepwater','Surface','Surface','Midwater','Midwater','Deepwater','Midwater','Surface','Surface')
sp_trophic <- data.frame(species=spp,trophic=trophic)

sim_param <- read.csv(here('Output','Complete_db','sim_GP_param_sim.csv'))

est_sim_param_1 <- read.csv(here('Output','Complete_db','est_GP_param_sim_1.csv'))
est_sim_param_2 <- read.csv(here('Output','Complete_db','est_GP_param_sim_2.csv'))
est_obs_param <- read.csv(here('Output','Complete_db','est_obs_param.csv'))

sim_conc_1 <- read_rds(here('Output','Complete_db','xy_conc_sim_1.rds'))
sim_conc_2 <- read_rds(here('Output','Complete_db','xy_conc_sim_2.rds'))
obs_conc <- read_rds(here('Output','Complete_db','xy_conc_obs.rds'))

sim_conc_pred_1 <- read_rds(here('Output','Complete_db','xy_conc_pred_sim_1.rds'))
sim_conc_pred_2 <- read_rds(here('Output','Complete_db','xy_conc_pred_sim_2.rds'))
obs_conc_pred <- read_rds(here('Output','Complete_db','xy_conc_pred_obs.rds'))

NNI <- read_rds(here('Output','Complete_db','NNI_obs.rds'))


# Plots ---------------------------------------------------------------------------------------


est_obs_param %>%
	# mutate(it=if_else(is.na(it),0,it)) %>%
	filter(param=='rho') %>%
	ggplot()+
	geom_point(aes(x=N,y=exp(mean)*100,colour = species),alpha=0.6)+
	geom_smooth(aes(x=N,y=exp(mean)*100,colour = species),se=F)+
	facet_wrap(~species,scales='free_y')+
	ylab('rho')+
	scale_y_log10()+
	scale_color_manual(values=my_sp_colors)+
	theme_bw()

ggsave(here('Plots','expl_1.jpg'))

est_sim_param_2 %>%
	filter(rho_it>2&rho_it<13) %>%
	# mutate(it=if_else(is.na(it),0,it)) %>%
	filter(parameter=='rho') %>%
	ggplot()+
	geom_point(aes(x=N,y=exp(value)*100,colour = factor(rho_it)),alpha=0.5)+
	geom_smooth(aes(x=N,y=exp(value)*100,colour = factor(rho_it)),se=F)+
	facet_wrap(~rho_it)+
	ylab('rho')+
	scale_y_log10()+
	theme_bw()
ggsave(here('Plots','expl_2.jpg'))

est_sim_param_1 %>%
	filter(rho_it>2&rho_it<13) %>% 
	# mutate(it=if_else(is.na(it),0,it)) %>%
	filter(parameter=='rho') %>%
bind_rows(.,
					est_sim_param_2 %>%
						filter(rho_it>2&rho_it<13) %>%
						# mutate(it=if_else(is.na(it),0,it)) %>%
						filter(parameter=='rho')) %>% 
	ggplot()+
	geom_point(aes(x=N,y=exp(value)*100,colour = factor(rho_it)),alpha=0.5)+
	geom_smooth(aes(x=N,y=exp(value)*100,colour = factor(rho_it)),se=F)+
	ylab('rho')+
	scale_y_log10()+
	scale_x_log10()+
	facet_wrap(~rho_it)+
	theme_bw()

ggsave(here('Plots','expl_3.jpg'))

est_obs_param_wi <- est_obs_param %>%
	filter(grepl('mu',param)) %>%
	mutate(mu_it = as.integer(str_extract(param, "(?<=\\[)\\d+(?=\\])"))) %>%
	transmute(N, species, it,mu_it,mu = mean) %>%
	left_join(.,
	est_obs_param %>%
	filter(param %in% c("alpha", "rho", "sigma")) %>%
	select(N, species, it, param, mean) %>%
	pivot_wider(names_from  = param,values_from = mean),by = c("N", "species", "it")) %>%
	left_join(.,NNI,
						by=c('N','species','it','mu_it'))

est_sim_param_1_wi <- est_sim_param_1 %>%
	filter(grepl('mu',parameter)) %>%
	mutate(mu_it = as.integer(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>%
	transmute(N, rho_it, it,mu_it,mu = value) %>%
	left_join(.,
	est_sim_param_1 %>%
	filter(parameter %in% c("alpha", "rho", "sigma")) %>%
	select(N, rho_it, it, parameter, value) %>%
	pivot_wider(names_from  = parameter,values_from = value),by = c("N", "rho_it", "it")) %>%
	left_join(.,NNI_sim_1,
						by=c('N','rho_it','it','mu_it'))


est_obs_param_wi %>% 
	filter(N==350) %>% 
ggplot()+
	geom_point(aes(x=exp(rho)*100,y=mu,colour = species),size=4)+
	geom_smooth(aes(x=exp(rho)*100,y=mu),se=F,span=2,lty=2,color='grey')+
	scale_colour_manual(values=my_sp_colors)+
	# facet_wrap(~N,ncol=3)+
	# facet_wrap(~species*mu_it,ncol=5,scale='free_y')+
	# facet_wrap(~species,ncol=5)+
	xlab('rho')+
	scale_x_log10()+
	facet_wrap(~depth,ncol=5)+
	theme_bw()
ggsave(here('Plots','expl_4.jpg'),width = 14,height = 6)

est_obs_param_wi %>% filter(N==350) %>% 
ggplot()+
	geom_point(data=est_sim_param_1_wi %>% filter(N==350),
						 aes(x=exp(rho)*100,y=alpha),size=2,alpha=0.04)+
	geom_point(aes(x=exp(rho)*100,y=alpha,color=species),size=4)+
	scale_color_manual(values=my_sp_colors)+
	xlab('rho')+
	scale_x_log10()+
	theme_bw()
ggsave(here('Plots','expl_5.jpg'),width = 9,height = 6)

est_obs_param_wi %>% filter(N==350) %>% 
ggplot()+
	geom_point(data=est_sim_param_1_wi %>% filter(N==350),
						 aes(x=exp(rho)*100,y=sigma),size=2,alpha=0.04)+
	geom_point(aes(x=exp(rho)*100,y=sigma,color=species),size=4)+
	scale_color_manual(values=my_sp_colors)+
	xlab('rho')+
	scale_x_log10()+
	theme_bw()
ggsave(here('Plots','expl_6.jpg'),width = 9,height = 6)


est_obs_param_wi %>% 
	filter(N==350) %>% 
	left_join(.,sp_trophic,by='species') %>% 
ggplot()+
	geom_point(aes(x=depth,y=mu,colour = factor(rho)))+
	geom_line(aes(x=depth,y=mu,colour = factor(rho)))+
	scale_color_manual(values = moma.colors('ustwo', 12))+
	facet_wrap(~trophic)+
	theme_bw()



NNI %>% ggplot()+
	geom_point(aes(x=N,y=NNI,color=species))+
	geom_smooth(aes(x=N,y=NNI,color=species),se=F)+
	facet_wrap(~depth)+
	theme_bw()


pred_summ_1 <- sim_conc_pred_1 %>% 
	rename(rho_it='rho') %>% 
	filter(rho_it>1&rho_it<13) %>% 
	group_by(it,rho_it,depth,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim_df,by=c('rho_it')) %>% 
	mutate(rho=rho*100)


w <- c(range(pred_data$X_utm),range(pred_data$Y_utm))

NNI_sim_1 <- sim_conc_1 %>% 
	rename(rho_it='rho') %>% 
	filter(rho_it>1&rho_it<13) %>% 
	group_by(it,rho_it,depth,N) %>% 
	summarise(NNI = compute_NNI(data   = cur_data(),x= X_utm,y= Y_utm,window = w),.groups = "drop") %>% 
	mutate(
		mu_it = case_when(
			depth == 0   ~ 1,
			depth == 50  ~ 2,
			depth == 150 ~ 3,
			depth == 300 ~ 4,
			depth == 500 ~ 5,
			TRUE ~ NA_real_
		)) %>% 
	mutate(N=as.numeric(N))

pred_summ_1 %>% 
	left_join(.,NNI_sim_1,by=c('it','rho_it','depth','N')) %>% 
	ggplot()+
	geom_point(aes(x=NNI,y=mean_diff,color=(rho)))+
	geom_smooth(aes(x=NNI,y=mean_diff,colow=(rho)),se=F)+
	ylab('epsilon')+
	xlab('Nearest-Neighbor Index')+
	facet_wrap(~N)+
	theme_bw()
ggsave(here('Plots','expl_7.jpg'),width = 14,height = 9)

pred_summ_2 <- sim_conc_pred_2 %>% 
	# rename(rho_it='rho') %>% 
	filter(rho_it>1&rho_it<13) %>% 
	group_by(it,rho_it,depth,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim_df,by=c('rho_it')) %>% 
	mutate(rho=rho*100)


w <- c(range(pred_data$X_utm),range(pred_data$Y_utm))

NNI_sim_2 <- sim_conc_2 %>% 
	rename(rho_it='rho') %>% 
	filter(rho_it>1&rho_it<13) %>% 
	mutate(N=as.numeric(N)) %>% 
	group_by(it,rho_it,depth,N) %>% 
	summarise(NNI = compute_NNI(data   = cur_data(),x= X_utm,y= Y_utm,window = w),.groups = "drop") %>% 
	mutate(
		mu_it = case_when(
			depth == 0   ~ 1,
			depth == 50  ~ 2,
			depth == 150 ~ 3,
			depth == 300 ~ 4,
			depth == 500 ~ 5,
			TRUE ~ NA_real_
		)) %>% 
	mutate(N=as.numeric(N))

pred_summ_2 %>% 
	left_join(.,NNI_sim_2,by=c('it','rho_it','depth','N')) %>% 
	ggplot()+
	geom_point(aes(x=NNI,y=mean_diff,color=(rho)))+
	geom_smooth(aes(x=NNI,y=mean_diff,colow=(rho)),se=F)+
	ylab('epsilon')+
	xlab('Nearest-Neighbor Index')+
	facet_wrap(~N)+
	theme_bw()
ggsave(here('Plots','expl_8.jpg'),width = 14,height = 9)



# Figure 1 ------------------------------------------------------------------------------------


plot_1_data <- 
	est_sim_param_1 %>% 
	filter(parameter=='rho')  %>% 
	left_join(.,sim_param %>% select(-N,-source) %>%rename(sim_val='value'),
						by=c('parameter','rho_it','it')) %>% 
	filter(rho_it>2&rho_it<13) %>% 
	mutate(value=exp(value)*100) %>% 
	mutate(sim_val=sim_val*100) %>% 
	# mutate(lo_q=lo_q*100) %>% 
	# mutate(up_q=up_q*100) %>%
	mutate(
		N_dens = case_when(
			N == 3000 ~ "N==3000*' (s=350 km'^2*')'",
			N == 2000 ~ "N==2000*' (s=500 km'^2*')'",
			N == 1000 ~ "N==1000*' (s=1000 km'^2*')'",
			N == 600  ~ "N==600*' (s=1600 km'^2*')'",
			N == 350  ~ "N==358*' (s=3000 km'^2*')'",
			N == 200  ~ "N==200*' (s=5000 km'^2*')'"),
		N_dens = factor(N_dens,levels = c(
			"N==200*' (s=5000 km'^2*')'",
			"N==358*' (s=3000 km'^2*')'",
			"N==600*' (s=1600 km'^2*')'",
			"N==1000*' (s=1000 km'^2*')'",
			"N==2000*' (s=500 km'^2*')'",
			"N==3000*' (s=350 km'^2*')'")))

p1 <-
plot_1_data %>% filter(!(N==3000)) %>%
	ggplot() +
	# geom_errorbar(aes(x = sim_val, y = value, ymin = lo_q, ymax = up_q),
	# 							color = 'grey', width = 0.02) +
	geom_abline(intercept = 0, slope = 1, color = 'black', lty = 2, linewidth = 0.5) +
	geom_point(aes(x = sim_val, y = value, color=N_dens), alpha = 0.6, size = 1.8) +
	geom_smooth(aes(x = sim_val, y = value, fill=N_dens), level = 0.995,alpha = 0.3, linewidth = 0) +
	# geom_smooth(aes(x = sim_val, y = value, fill=N_dens), level = 0.70,alpha = 0.4, linewidth = 0) +
	# geom_smooth(aes(x = sim_val, y = value, fill=N_dens), level = 0.50,alpha = 0.3, linewidth = 0) +
	# geom_smooth(aes(x = sim_val, y = value, fill=N_dens), level = 0.10,alpha = 0.4, linewidth = 0) +
	scale_color_manual(values = my_colors)+
	scale_fill_manual(values = my_colors) +
	scale_x_log10() +
	scale_y_log10() +
	facet_wrap(~ N_dens, labeller = label_parsed) +
	labs(
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		y = expression("Estimated spatial autocorrelation - " * rho * " (km)")) +
	theme_bw() +
	theme(legend.position = 'none') +
	my_theme()

p1


model <- lm(log(sim_val) ~ log(value), data = plot_1_data)
residuals <- resid(model)

plot_1_data$residuals <- resid(model)

residual_summary <- plot_1_data %>%
	group_by(N,sim_val) %>%
	summarise(
		mean_resid = mean(residuals),
		sd_resid   = sd(residuals),)

residual_summary %>% 	
	filter(!N==3000) %>% 
	ggplot()+
	geom_point(aes(x=sim_val,y=sd_resid,color=factor(N)))+
	geom_line(aes(x=sim_val,y=sd_resid,color=factor(N)))+
	# geom_smooth(aes(x=sim_val,y=sd_resid,color=factor(N)),se=F)+
	scale_color_manual(values = my_colors)+
	scale_fill_manual(values = my_colors) +
	labs(y='Residual Standard Error for estimating rho (ρ)',color='Number of\nsamples\ndeployed',x='Rho')+
	scale_x_log10() +
	theme_bw()+
	my_theme()

ggsave(here('Plots','expl_9.jpg'),width = 14,height = 9)


# Figure 2 ------------------------------------------------------------------------------------


plot_2_data <- est_sim_param_1 %>% 
	left_join(.,sim_param %>% select(-N,-source) %>%rename(sim_val='value'),
						by=c('parameter','rho_it','it')) %>% 
	filter(!parameter%in%c('alpha','rho','sigma')) %>% 
	filter(rho_it>2&rho_it<13) %>% 
	left_join(.,rho_sim_df %>% mutate(rho=rho*100),by=c('rho_it')) %>%
	mutate(
		N_dens = case_when(
			N == 3000 ~ "N==3000*' (s=350 km'^2*')'",
			N == 2000 ~ "N==2000*' (s=500 km'^2*')'",
			N == 1000 ~ "N==1000*' (s=1000 km'^2*')'",
			N == 600  ~ "N==600*' (s=1600 km'^2*')'",
			N == 350  ~ "N==358*' (s=3000 km'^2*')'",
			N == 200  ~ "N==200*' (s=5000 km'^2*')'"),
		N_dens = factor(N_dens,levels = c(
			"N==200*' (s=5000 km'^2*')'",
			"N==358*' (s=3000 km'^2*')'",
			"N==600*' (s=1600 km'^2*')'",
			"N==1000*' (s=1000 km'^2*')'",
			"N==2000*' (s=500 km'^2*')'",
			"N==3000*' (s=350 km'^2*')'"))) %>% 
	mutate(rho_label = paste0("rho==", rho, "*' km'")) %>% 
	mutate(rho_label = factor(
		rho_label,
		levels = paste0("rho==", sort(unique(rho)), "*' km'")
	))


p2 <- plot_2_data %>% filter(!(N==3000)) %>%
	ggplot()+
	geom_abline(intercept = 0, slope = 1, color = 'black', lty = 2, linewidth = 0.5) +
	geom_point(aes(x=sim_val,y=value,color=N_dens),alpha=0.24,size=2)+
	scale_color_manual(values = my_colors)+
	labs(
		x = expression("Simulated global mean - " * mu * " (log copies/L)"),
		y = expression("Estimated global mean - " * mu * " (log copies/L)")
	) +
	facet_wrap(~ N_dens, labeller = label_parsed,ncol = 3) +
	theme_bw()+
	theme(legend.position = 'none') +
	my_theme()

# p2
# ggsave(here('Plots','Figure_2.jpg'),p2,width=16,height = 10)

p2_2 <- plot_2_data %>% filter(!(N==3000)) %>%
	filter(!rho_it%in%c(1,2)) %>%
	ggplot()+
	geom_abline(intercept = 0, slope = 1, color = 'black', lty = 2, linewidth = 0.5) +
	geom_point(aes(x=sim_val,y=value,color=N_dens),alpha=0.24)+
	scale_color_manual(values = my_colors)+
	labs(
		x = expression("Simulated global mean - " * mu * " (log copies/L)"),
		y = expression("Estimated global mean - " * mu * " (log copies/L)")
	) +
	# facet_wrap(~ rho*N_dens, labeller = label_parsed,ncol = 6) +
	facet_wrap(~ rho_label*N_dens, labeller = label_parsed,ncol = 5) +
	theme_bw()+
	theme(legend.position = 'none') +
	my_theme()

# p2_2
# ggsave(here('Plots','Figure_2_extended.jpg'),p2_2,width=16*1,height = 10*3)


model <- lm(sim_val ~ value, data = plot_2_data)
residuals <- resid(model)

plot_2_data$residuals <- resid(model)

residual_summary <- plot_2_data %>%
	group_by(N,rho) %>%
	summarise(
		mean_resid = mean(residuals),
		sd_resid   = sd(residuals),)

residual_summary

residual_summary %>% 	
	filter(!N==3000) %>% 
	ggplot()+
	geom_point(aes(x=rho,y=sd_resid,color=factor(N)))+
	geom_line(aes(x=rho,y=sd_resid,color=factor(N)))+
	# geom_smooth(aes(x=sim_val,y=sd_resid,color=factor(N)),se=F)+
	scale_color_manual(values = my_colors)+
	scale_fill_manual(values = my_colors) +
	scale_x_log10() +
	labs(y='Residual Standard Error for estimating mu (μ)',color='Number of\nsamples\ndeployed',x='Rho')+
	theme_bw()+
	my_theme()
ggsave(here('Plots','expl_10.jpg'),width = 14,height = 9)


# Figure 3 ------------------------------------------------------------------------------------

p3 <- sim_conc_pred_1 %>% 
	rename(rho_it='rho') %>% 
	group_by(it,rho_it,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim_df,by=c('rho_it')) %>% 
	mutate(rho=rho*100) %>% 
	mutate(N=if_else(N==350,358,N)) %>%
	filter(rho_it>2&rho_it<13) %>% 
	# filter(!rho_it%in%c(1,2,13:15)) %>% 
	# mutate(mean_diff=exp(mean_diff)) %>%
	ggplot()+
	geom_point(aes(x=rho,y=mean_diff,color=factor(N)))+
	geom_smooth(aes(x=rho,y=mean_diff,color=factor(N)),se=F)+
	labs(
		# y = expression("Mean prediction error - " * epsilon * " (" * frac(predicted, simulated) * ")"),
		# y = expression("Mean prediction error - " * epsilon * " = | log(" * frac(predicted, simulated) * ") |"),
		y = expression("Mean prediction error - " * epsilon * 
									 	" = |" * hat(Z[i]) - hat(S[i*e]) * "|"),
		# y = expression("Mean prediction error - " * epsilon *
		#              " |" * log(pred) - log(sim) * "|"),
		# y = expression("Mean prediction error - " * epsilon * 
		#              " (" * log(pred/sim) * ")"),
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		color='Number of\nsamples\ndeployed')+
	scale_x_log10()+
	# scale_y_log10(breaks = 1:5, labels = 1:5)+
	scale_color_manual(values=my_colors)+
	theme_bw()+
	my_theme()

p3
# ggsave(here('Plots','Figure_3.jpg'),p3,width=16,height = 10)

# mod <- pred_summ %>%
# 	filter(N == 350, !rho %in% c(1, 2)) %>%
# 	# lm(mean_diff ~ log(rho_sim), data = .) %>% 
# 	lm(mean_diff ~ I((log(rho_sim))^(-1.5)), data = .)
# cowplot::plot_grid(p3,p3_2+scale_x_log10(limits = c(25, 1000)))




# Figure 4 ------------------------------------------------------------------------------------

p4 <-
sim_conc_pred_2 %>%
	group_by(it,rho_it,N,depth) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim_df,by=c('rho_it')) %>% 
	mutate(rho=rho*100) %>% 
	mutate(N=if_else(N==350,358,N)) %>%
	filter(rho_it>2&rho_it<13) %>% 
	ggplot()+
	geom_point(aes(x=rho,y=mean_diff,color=factor(N)))+
	geom_smooth(aes(x=rho,y=mean_diff,color=factor(N)),se=F)+
	labs(
		y = expression("Mean prediction error - " * epsilon * 
									 	" = |" * hat(Z[i]) - hat(S[i*e]) * "|"),
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		color='Number of\nsamples\ndeployed')+
	facet_grid(~depth)+
	scale_x_log10()+
	scale_color_manual(values=my_colors_2)+
	theme_bw()+
	my_theme()

p4
# ggsave(here('Plots','Figure_3.jpg'),p3,width=16,height = 10)



# Figure 4 ------------------------------------------------------------------------------------

# p4 <-
sim_conc_pred_2 %>%
	group_by(it,rho_it,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim_df,by=c('rho_it')) %>% 
	mutate(rho=rho*100) %>% 
	mutate(N=if_else(N==350,358,N)) %>%
	filter(rho_it>2&rho_it<13) %>% 
	mutate(source='Simulation') %>% 
bind_rows(.,
					obs_conc_pred %>%
						group_by(it,species,N) %>% 
						summarise(mean_diff=mean(diff)) %>% 
						left_join(.,est_obs_param %>% filter(param=='rho') %>% filter(N==350) %>% select(species,mean) %>% rename(rho='mean'),by=c('species')) %>% 
						mutate(rho=exp(rho)*100) %>% 
						mutate(N=if_else(N==350,358,N)) %>% 
						mutate(source='Guri et al. (2025)')) %>% 
	ggplot()+
	geom_point(aes(x=rho,y=mean_diff,color=factor(N)))+
	geom_smooth(aes(x=rho,y=mean_diff,color=factor(N)),se=F,span=2)+
	labs(
		y = expression("Mean prediction error - " * epsilon * 
									 	" = |" * hat(Z[i]) - hat(S[i*e]) * "|"),
		x = expression("Simulated spatial autocorrelation - " * rho * " (km)"),
		color='Number of\nsamples\ndeployed')+
	scale_x_log10()+
	scale_color_manual(values=my_colors_2)+
	facet_wrap(~source) +
	theme_bw()+
	my_theme()

