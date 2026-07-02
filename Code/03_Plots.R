
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
library(broom)
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
# my_sp_colors <- c('#6a3d9a','#cab2d6','#ff7f01','#fdbf6f','#e31b1d','#fb9a99','#33a02c','#b2df8a','#1f78b4','#009999','#999900', '#a6cee3') #old version
my_sp_colors <- c('#6a3d9a','#cab2d6','#ff7f01','#fdbf6f','#1f78b4','#fb9a99','#33a02c','#b2df8a','#e31b1d','#999900','#009999', '#a6cee3')
# my_colors <- c("#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d")
my_colors_2 <- c("#0660f1ff","#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d")
N_values <- c(3000, 2000, 1000, 600, 350, 200)
N_values_2 <- c(350, 330, 300, 260, 220, 150, 100)
rho_sim <- c(0.01,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20) # Length scales parameter to be simulated
rho_sim_df <- rho_sim %>% as.data.frame() %>% setNames('rho') %>% mutate(rho_it=row_number(.)) 
edna_data <- readRDS(here('Data','edna_data.rds'))
pred_data <- readRDS(here('Data','pred_data.rds'))
spp <- edna_data %>% distinct(species) %>% slice(-n()) %>% pull(species)
# trophic <- c('Surface','Surface','Midwater','Deepwater','Surface','Surface','Midwater','Midwater','Deepwater','Midwater','Surface','Surface')
trophic <- c('Surface','Surface','Deepwater','Deepwater','Surface','Surface','Midwater','Midwater','Deepwater','Midwater','Surface','Surface')
sp_trophic <- data.frame(species=spp,trophic=trophic)

sim_param <- read.csv(here('Output','Complete_db','sim_GP_param_sim.csv'))

sim_df <- read_rds(here('Output','Complete_db','xy_conc_sim_raw.rds'))

est_sim_param <- read.csv(here('Output','Complete_db','est_GP_param_sim.csv'))
est_obs_param <- read.csv(here('Output','Complete_db','est_obs_param.csv'))

sim_conc_1 <- read_rds(here('Output','Complete_db','xy_conc_sim_1.rds'))
sim_conc_2 <- read_rds(here('Output','Complete_db','xy_conc_sim_2.rds'))
obs_conc <- read_rds(here('Output','Complete_db','xy_conc_obs.rds'))

sim_conc_pred <- read_rds(here('Output','Complete_db','xy_conc_pred_sim_1.rds'))
sim_conc_pred_1_comp <- read_rds(here('Output','Complete_db','xy_conc_pred_sim_1_comp.rds'))
sim_conc_pred_2_comp <- read_rds(here('Output','Complete_db','xy_conc_pred_sim_2_comp.rds'))
# sim_conc_pred_2 <- read_rds(here('Output','Complete_db','xy_conc_pred_sim_2.rds'))
obs_conc_pred <- read_rds(here('Output','Complete_db','xy_conc_pred_obs.rds'))
obs_conc_pred_comp <- read_rds(here('Output','Complete_db','xy_conc_pred_obs_comp.rds'))

NNI <- read_rds(here('Output','Complete_db','NNI_obs.rds'))

est_obs_param_350 <- list.files(here('Output','GP_multifish','Multifish_raw'),pattern = paste0('_350_est_GP_param.csv$')) %>%
	map_df(~ read_csv(here('Output','GP_multifish','Multifish_raw',.x),show_col_types = FALSE) %>%
				 	mutate(file = .x)) %>%
	separate(file,into = c('species','N','type','est','ext1'),sep = '_',remove = FALSE) %>%
	select(-ext1,-est,-type,-`...1`) %>% 
	mutate(N = as.numeric(N))

# Figure 3 --- --------------------------------------------------------------------------------

p1_f3_d <-
est_obs_param_350 %>% 
# est_obs_param %>% filter(N==350)
	filter(param == "rho") %>% 
	select(-file, -N, -param) %>% 
	left_join(sp_trophic, by = "species") %>%
	mutate(rho = exp(mean) * 100,lo  = exp(`25%`) * 100,up  = exp(`75%`) * 100) 

p2_f3_d <- est_obs_param_350 %>% 
	filter(grepl("mu\\[", param)) %>% 
	mutate(param = recode(param,
								 "mu[1]" = "0",
								 "mu[2]" = "50",
								 "mu[3]" = "150",
								 "mu[4]" = "300",
								 "mu[5]" = "500")) %>% 
	rename(depth='param') %>% mutate(depth=as.numeric(depth)) %>% 
	left_join(.,p1_f3_d %>% select(species,rho))


p3_f3_d <- est_obs_param_350 %>% 
	filter(param%in%c('alpha','sigma')) %>% 
	select(species,param,mean,`2.5%`,`97.5%`) %>% 
	pivot_wider(
		names_from  = param,
		values_from = c(mean, `2.5%`, `97.5%`))


p1_f3 <-
	p1_f3_d %>% 
	mutate(trophic = factor(trophic, levels = c("Surface", "Midwater", "Deepwater"))) %>%
	ggplot() +
	geom_point(aes(x = species, y = rho,color=species, size=3)) +
	geom_errorbar(aes(x = species, ymin = lo, ymax = up,color=species),width=0.3,size=1) +
	scale_color_manual(values=my_sp_colors)+
	ylab(expression('Estimated length-scale parameter - '*rho*' (km)'))+
	# xlab('Species')+
	scale_y_log10() +
	facet_grid(~trophic, scales = "free_x") +
	theme_bw() +
	theme(axis.text.x = element_blank(),
				axis.title.x = element_blank(),
				legend.position = 'none')+
	my_theme()
p1_f3

p2_f3 <-
	p2_f3_d %>% 
	ggplot() +
	# geom_point(aes(x = depth, y = mean, color = species, size = rho)) +
	geom_point(aes(x = depth, y = mean, color = species),size=3) +
	geom_line(aes(x = depth, y = mean, color = species),lty=2,alpha=0.6) +
	labs(x = "Depth (m)", y = expression(mu[d])) +
	scale_color_manual(name = "Species", values = my_sp_colors,
		labels = function(x) paste0("italic('", x, "')") |> parse(text = _)) +
	# scale_size_continuous(name = expression(rho*' (km)'),trans = "log10",range = c(1, 10)) +
	theme_bw() +
	my_theme()

p3_f3 <-
	p3_f3_d %>% 
	ggplot() +
	geom_point(aes(x = mean_alpha, y = mean_sigma, color = species),size=3) +
	geom_errorbar(aes(x = mean_alpha, y = mean_sigma, xmin = `2.5%_alpha`, xmax = `97.5%_alpha`, color = species),alpha=0.5,width=0.1) +
	geom_errorbar(aes(x = mean_alpha, y = mean_sigma, ymin = `2.5%_sigma`, ymax = `97.5%_sigma`, color = species),alpha=0.5,width=0.1) +
	scale_color_manual(values = my_sp_colors)+
	theme_bw() +
	labs(x = expression(alpha^2), y = expression(sigma^2)) +
	theme(legend.position = 'none')+
	my_theme()

f3_leg <- cowplot::get_legend(p2_f3)
	
pp1_f3 <- cowplot::plot_grid(p2_f3+theme(legend.position = 'none'),p3_f3,rel_widths = c(5,3), labels = c('B','C'))

pp2_f3 <- cowplot::plot_grid(p1_f3,pp1_f3,ncol = 1, labels = c('A'))

figure_3 <- cowplot::plot_grid(pp2_f3,f3_leg,rel_widths = c(8,2))
figure_3

# ggsave(here('Plots','Figure_3.jpg'),figure_3,width=14,height = 10)
# ggsave(here('Plots','Figure_3.pdf'),figure_3,width=14,height = 10)


# Figure 4 --- --------------------------------------------------------------------------------

# rm(sim_conc_pred_2_comp)
p1_f4_d <- sim_conc_pred_1_comp %>%
	filter(N<3000) %>% 
	rename(rho_it='rho') %>% 
	group_by(it,rho_it,N) %>% 
	summarise(mean_diff=mean(diff),
						bias=mean(conc-Conc_0),
						rmse=sqrt(mean((conc - Conc_0)^2))) %>% 
	left_join(.,rho_sim_df,by=c('rho_it')) %>% 
	mutate(rho=rho*100) %>% 
	mutate(N=if_else(N==350,358,N)) %>%
	filter(rho_it>2&rho_it<13) %>% 
	mutate(eff=mean_diff/N) %>% 
	mutate(dens=(N/5)/(200*1000)*10000) 

# sim_conc_pred_2_comp <- read_rds(here('Output','Complete_db','xy_conc_pred_sim_2_comp.rds'))


p1_f4 <- p1_f4_d %>%
	ggplot() +
	geom_point(aes(x = dens, y = rmse, color = factor(rho)), alpha = 0.3) +
	stat_summary(
		aes(x = dens, y = rmse, color = factor(rho), group = rho),
		fun = mean,
		geom = "line",
		linewidth = 1) +
	labs(x = expression('Sampling density (E = N in 10000 km'^2*')'),
			 y = expression('Prediction error ('*epsilon[2]*')'),
			 color = expression(rho*' (km)')) +
	scale_color_manual(values = pnw_palette("Bay", n = 10)) +
	scale_x_continuous(breaks = c(1, 5, 10, 20, 30)) +
	theme_bw() +
	my_theme()


model_data <-
	p1_f4_d %>% 
	mutate(rmse=log(rmse),
				 dens=log(dens),
				 rho=log(rho)) %>% 
	ungroup()

model <- lm(rmse ~ dens*rho, data = model_data)

new_data <- new_df <- expand.grid(
	dens = seq(1, 20, by = 1),
	rho  = c(25,50,75,100,200,300,400,500,800,1000)) %>% 	
	mutate(dens=log(dens),rho=log(rho))

pred_log <- predict(model, newdata = new_data) %>% as.data.frame() %>% setNames('pred') %>% 
	# pred_log <- predict(model, newdata = model_data %>% select(rho,dens)) %>% as.data.frame() %>% setNames('pred') %>% 
	mutate(rmse=(pred)) %>% 
	# cbind(.,model_data %>% select(rho,dens)) %>% 
	cbind(.,new_data) %>% 
	mutate(rmse=exp(rmse),
				 dens=exp(dens),
				 rho=exp(rho))


p2_f4_d <- pred_log %>% 
	group_by(rho) %>% 
	arrange(dens, .by_group = TRUE) %>%
	mutate(
		dens_next = lead(dens),
		rmse_next = lead(rmse),
		d_rmse_d_dens = (rmse_next - rmse) / (dens_next - dens),
		dens_mid  = (dens + dens_next) / 2,
		dens_mid_2 = paste0(dens, "→", dens_next)
	) %>% 
	ungroup() %>% 
	filter(!is.na(d_rmse_d_dens))

label_df <- p2_f4_d %>% distinct(dens_mid, dens_mid_2)

p2_f4 <-
	p2_f4_d %>% 
	ggplot(aes(x = dens_mid, y = d_rmse_d_dens, color = factor(rho))) +
	geom_point(alpha = 0.4) +
	geom_line() +
	geom_hline(yintercept = -0.05, linetype = "dashed", alpha = 0.4) +
	scale_color_manual(values = pnw_palette("Bay", n = 10)) +
	scale_x_continuous(
		breaks = c(1:10,12,14,16,18,20),
		labels = c(1:10,12,14,16,18,20)
		# breaks = label_df$dens_next,
		# labels = label_df$dens_next
	) +
	labs(x = 'Sampling density (E)', 
			 # x = expression(Delta*'Sampling density (E)'), 
			 y = expression('Marginal change in prediction error ('*Delta * epsilon[2] * ' / ' * Delta*'E)'))+ 
	theme_bw() +
	theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
		panel.grid.minor.x = element_blank(),
		legend.position = "none") +
	my_theme()

f4_leg <- cowplot::get_legend(p1_f4)

figure_4 <- cowplot::plot_grid(p1_f4+theme(legend.position = 'none'),p2_f4,rel_widths = c(5,2.5),labels=c('A','B')) %>% 
	cowplot::plot_grid(.,f4_leg,rel_widths = c(11,1))
figure_4

# ggsave(here('Plots','Figure_4.jpg'),figure_4,width=15,height = 9)
# ggsave(here('Plots','Figure_4.pdf'),figure_4,width=15,height = 9)


# Figure 5 --- --------------------------------------------------------------------------------

p1_f5_d <-
	est_sim_param %>% 
	filter(parameter=='rho') %>% 
	left_join(.,sim_param %>% 
							select(-N,-source) %>% 
							rename(sim_val='value') %>% 
							filter(parameter=='rho') %>% 
							mutate(sim_val=log(sim_val)),
						by=c('parameter','rho_it','it')) %>% 
	filter(rho_it>2&rho_it<13) %>% 
	mutate(dens=(N/5)/(200*1000)*10000) %>% 
	group_by(rho_it,dens) %>% 
	mutate(rmse = sqrt(mean((value - sim_val)^2))) %>%
	ungroup() 

p1_f5 <-
	p1_f5_d %>% 
	filter(N<3000) %>% 
	mutate(sim_val_rho=exp(sim_val)*100) %>% 
	group_by(dens,sim_val_rho) %>% 
	summarise(rmse=first(rmse), .groups = "drop") %>% 
	ggplot() +
	geom_point(aes(x = dens, y = rmse,color=factor(sim_val_rho)))+
	geom_line(aes(x = dens, y = rmse,color=factor(sim_val_rho)))+
	labs(x=expression('Sampling density (N in 100x100 km grid - 10000 km'^2*')'),
			 y=expression(tau),color=expression(rho*' (km)'))+
	theme_bw()+
	annotate("text",x = Inf, y = Inf,label = "rho",parse = TRUE,hjust = 1.5, vjust = 1.5,size = 15) +
	# scale_y_sqrt()+
	scale_color_manual(values = pnw_palette("Bay", n = 10)) +
	theme(axis.title.x = element_blank(),
				axis.title.y = element_blank()) +
	my_theme()

f5_leg <- cowplot::get_legend(p1_f5)

# p1 <- p11+theme(legend.position = 'none')

p2_f5_d <- 
	est_sim_param %>% 
	mutate(
	parameter = recode(parameter,
										 "mu[1]" = "mu_0",
										 "mu[2]" = "mu_50",
										 "mu[3]" = "mu_150",
										 "mu[4]" = "mu_300",
										 "mu[5]" = "mu_500")) %>% 
	left_join(.,sim_param %>% select(-N,-source) %>%rename(sim_val='value'),
						by=c('parameter','rho_it','it')) %>% 
	filter(!parameter%in%c('alpha','rho','sigma')) %>% 
	filter(rho_it>2&rho_it<13) %>% 
	left_join(.,rho_sim_df %>% mutate(rho=rho*100),by=c('rho_it')) %>%
	mutate(dens=(N/5)/(200*1000)*10000) %>% 
	mutate(rho_label = paste0("rho==", rho, "*' km'")) %>% 
	mutate(rho_label = factor(
		rho_label,
		levels = paste0("rho==", sort(unique(rho)), "*' km'")
	))

p2_f5 <-
	p2_f5_d %>% 
	filter(N<3000) %>%
	group_by(rho_it,dens) %>% 
	summarise(rmse = sqrt(mean((value - sim_val)^2))) %>%
	ggplot() +
	geom_point(aes(x = dens, y = rmse,color=factor(rho_it)))+
	geom_line(aes(x = dens, y = rmse,color=factor(rho_it)))+
	scale_color_manual(values = pnw_palette("Bay", n = 10)) +
	labs(x=expression('Sampling density (N in 10000 km'^2*')'),
			 y=expression(tau),color=expression('Simulated '*rho))+
	annotate("text",x = -Inf, y = -Inf,label = "mu",parse = TRUE,hjust = -0.5, vjust = -0.5,size = 15) +
	theme_bw()+
	theme(legend.position = 'none',
				axis.title.x = element_blank()) +
	my_theme()



p3_f5_d <- 
	est_sim_param %>% 
	filter(parameter%in%c('alpha','sigma')) %>% 
	left_join(.,sim_param %>% select(-N,-source) %>%rename(sim_val='value'),
						by=c('parameter','rho_it','it')) %>% 
	filter(rho_it>2&rho_it<13) %>% 
	left_join(.,rho_sim_df %>% mutate(rho=rho*100),by=c('rho_it')) %>%
	mutate(dens=(N/5)/(200*1000)*10000) %>% 
	mutate(rho_label = paste0("rho==", rho, "*' km'"))

p3_f5 <- p3_f5_d %>% 
	filter(parameter=='alpha') %>%
	filter(N<3000) %>%
	group_by(rho_it,dens) %>% 
	summarise(rmse = sqrt(mean((value - sqrt(4))^2))) %>%
	# summarise(rmse = (mean((value - 2)))) %>%
	ggplot() +
	geom_point(aes(x = dens, y = rmse,color=factor(rho_it)))+
	geom_line(aes(x = dens, y = rmse,color=factor(rho_it)))+
	scale_color_manual(values = pnw_palette("Bay", n = 10)) +
	annotate("text",x = -Inf, y = -Inf,label = "alpha^scriptscriptstyle(2)",parse = TRUE,hjust = -0.5, vjust = -0.5,size = 15) +
	labs(x=expression('Sampling density (N in 10000 km'^2*')'),
			 y=expression(tau),color=expression('Simulated '*rho))+
	theme_bw()+
	theme(legend.position = 'none') +
	my_theme()

p4_f5 <- p3_f5_d %>% 
	filter(parameter=='sigma') %>%
	filter(N<3000) %>%
	group_by(rho_it,dens) %>% 
	summarise(rmse = sqrt(mean((value - 3)^2))) %>%
	# summarise(rmse = (mean((value - 3)))) %>%
	ggplot() +
	geom_point(aes(x = dens, y = rmse,color=factor(rho_it)))+
	geom_line(aes(x = dens, y = rmse,color=factor(rho_it)))+
	scale_color_manual(values = pnw_palette("Bay", n = 10)) +
	annotate("text", x = Inf, y = Inf, label = "sigma^scriptscriptstyle(2)", parse = TRUE, hjust = 1.5, vjust = 1.5, size = 15)+
	labs(x=expression('Sampling density (N in 10000 km'^2*')'),
			 y=expression(tau),color=expression('Simulated '*rho))+
	theme_bw()+
	theme(legend.position = 'none',
				axis.title.y = element_blank()) +
	my_theme()


figure_5 <- cowplot::plot_grid(p1_f5+theme(legend.position = 'none'),p4_f5,p2_f5,p3_f5,align = 'v',labels = c('A','B','C','D')) %>% 
	cowplot::plot_grid(.,f5_leg,rel_widths = c(10,1))
figure_5

# ggsave(here('Plots','Figure_5.jpg'),figure_5,width=14,height = 10)
# ggsave(here('Plots','Figure_5.pdf'),figure_5,width=14,height = 10)



# Supplementary Fig 1 -------------------------------------------------------------------------

su_fi_data_1 <- sim_conc_pred_2_comp %>%
	rename(rho_it='rho') %>% 
	group_by(it,rho_it,N) %>% 
	summarise(mean_diff=mean(diff),
						rmse=sqrt(mean((conc - Conc_0)^2))) %>% 
	left_join(.,rho_sim_df,by=c('rho_it'))  %>% 
	mutate(rho=rho*100) %>% 
	mutate(N=if_else(N==350,358,N)) %>%
	filter(rho_it>2&rho_it<13) %>% 
	mutate(source='Simulation') %>% 
	bind_rows(.,
						obs_conc_pred_comp %>%
							group_by(it,species,N) %>% 
							summarise(mean_diff=mean(diff),
												rmse=sqrt(mean((conc - Conc_0)^2))) %>% 
							left_join(.,est_obs_param %>% filter(param=='rho') %>% filter(N==350) %>% select(species,mean) %>% rename(rho='mean'),by=c('species')) %>% 
							mutate(rho=exp(rho)*100) %>% 
							# mutate(N=if_else(N==350,358,N)) %>% 
							mutate(source='Guri et al. (2026)')) %>% 
	mutate(dens=(N/5)/(200*1000)*10000)



model_data <-
	su_fi_data_1 %>% filter(is.na(species)) %>% 
	mutate(rmse=log(rmse),
				 dens=log(dens),
				 rho=log(rho)) %>% 
	ungroup()

model_1 <- lm(rmse ~ dens*rho, data = model_data)

model_data <-
	su_fi_data_1 %>% filter(!is.na(species)) %>% 
	mutate(rmse=log(rmse),
				 dens=log(dens),
				 rho=log(rho)) %>% 
	ungroup()

model_2 <- lm(rmse ~ dens*rho, data = model_data)

new_data <- new_df <- expand.grid(
	dens = seq(1.5, 3.3, by = 0.3),
	rho  = c(25,50,75,100,200,300,400,500,800,1000)) %>% 	
	mutate(dens=log(dens),rho=log(rho))

pred_log_1 <- predict(model_1, newdata = new_data) %>% as.data.frame() %>% setNames('pred') %>% 
	mutate(rmse=(pred)) %>% 
	cbind(.,new_data) %>% 
	mutate(rmse=exp(rmse),
				 dens=exp(dens),
				 rho=exp(rho))

pred_log_2 <- predict(model_2, newdata = new_data) %>% as.data.frame() %>% setNames('pred') %>% 
	mutate(rmse=(pred)) %>% 
	cbind(.,new_data) %>% 
	mutate(rmse=exp(rmse),
				 dens=exp(dens),
				 rho=exp(rho))


c1 <- tidy(model_1) %>% select(term, estimate, std.error) %>%
	rename(est1 = estimate, se1 = std.error)

c2 <- tidy(model_2) %>% select(term, estimate, std.error) %>%
	rename(est2 = estimate, se2 = std.error)

coef_compare <- full_join(c1, c2, by = "term") %>%
	mutate(
		diff   = est1 - est2,
		se_diff = sqrt(se1^2 + se2^2),          # assumes independent fits
		z      = diff / se_diff,
		p      = 2 * pnorm(abs(z), lower.tail = FALSE)
	) %>%
	arrange(desc(abs(z)))

# coef_compare

sf_1 <- coef_compare %>%
	mutate(term = recode(term,
											 "rho" = "θ",
											 "(Intercept)" = "ω",
											 "dens" = "β",
											 "dens:rho" = "γ"
	)) %>%
	ggplot(aes(x = est1, y = est2, label = term)) +
	geom_errorbarh(aes(xmin = est1 - (2*se1), xmax = est1 + (2*se1)),
								 height = 0.03, alpha = 0.6) +
	geom_errorbar(aes(ymin = est2 - (2*se2), ymax = est2 + (2*se2)),
								width = 0.03, alpha = 0.6) +
	geom_hline(yintercept = 0, linewidth = 0.2) +
	geom_vline(xintercept = 0, linewidth = 0.2) +
	geom_abline(slope = 1, intercept = 0, linetype = 2) +
	geom_point(size=2) +
	ggrepel::geom_text_repel(size = 7, max.overlaps = 20) +
	labs(x = "Model coefficients (simulation data)", 
			 y = "Model coefficients (real-world data)") +
	theme_bw() +
	my_theme()

sf_1
ggsave(here('Plots','SI_fig_1.jpg'),sf_1,height = 8,width = 11)


# Supplementary Fig 2 -------------------------------------------------------------------------

sf_2 <-
	su_fi_data_1 %>% 
	ggplot()+
	geom_point(aes(x=dens,y=rmse,color=rho),alpha=0.4)+
	geom_smooth(aes(x=dens,y=rmse,color=rho,group = rho),se=F,span=2)+
	labs(x=expression('Sampling density (N in 100x100 km grid - 10000 km'^2*')'),
			 y='RMSE',color=expression(rho*' (km)'))+
	scale_color_gradientn(colours = pnw_palette("Bay", n = 350,type='continuous'),
												trans='sqrt',limits  = c(25, 1000),oob= scales::squish,
												breaks  = c(25, 100, 250, 500, 750,1000)) +
	facet_wrap(~source) +
	theme_bw()+
	my_theme()

sf_2
ggsave(here('Plots','SI_fig_2.jpg'),sf_2,height = 8,width = 12)


# Supplementary Fig 3 -------------------------------------------------------------------------


sf_3 <-
	p2_f4_d %>% 
	mutate(lower=exp(log(100)-rmse)) %>% 
	mutate(upper=exp(log(100)+rmse)) %>% 
	# mutate(lower=100/exp(-d_rmse_d_dens)) %>% 
	# mutate(upper=100*exp(-d_rmse_d_dens)) %>% 
	ggplot()+
	geom_point(aes(x=dens_mid,y=100,color=factor(rho)))+
	geom_errorbar(aes(x=dens_mid,y=100,ymin = lower,ymax=upper,color=factor(rho)))+
	scale_color_manual(name = expression(rho* ' (km)'),values = pnw_palette("Bay", n = 10)) +
	labs(x = expression('Sampling density (Number of samples / 10000km'^2*')'),
			 y = 'Concentration estimation for true concnetration of 100 copies/L')+
	facet_wrap(~rho,ncol=2)+
	scale_y_log10(breaks=c(10,30,100,300,1000))+
	theme_bw()

sf_3
ggsave(here('Plots','SI_fig_3.jpg'),sf_3,height = 12,width = 10)



# Supplementary Fig 4 -------------------------------------------------------------------------


sf_4 <-
p1_f4_d %>%
	group_by(it,rho,N) %>%
	summarise(mean_bias=mean(bias)) %>%
	ggplot()+
	ylab(expression('Concentration derived from Z - S'[E]))+
	xlab(expression(rho*' (km)'))+
	geom_point(aes(x=rho,y=mean_bias))+
	scale_x_log10()+
	theme_bw()+
	facet_grid(~N)

sf_4
ggsave(here('Plots','SI_fig_4.jpg'),sf_4,height = 6,width = 14)



# Supplementary Fig 5 -------------------------------------------------------------------------

sf_5 <-
p1_f4_d %>% 
	ggplot()+
	geom_line(data=pred_log,aes(x=dens,y=rmse,color=factor(rho)),lty=2)+
	geom_point(aes(x=dens,y=rmse,color=factor(rho),group = rho))+
	# scale_x_log10()+
	labs(x=expression('Sampling density (N in 100x100 km grid - 10000 km'^2*')'),
			 y=expression(epsilon[2]),color=expression(rho*' (km)'))+
	# scale_y_log10(breaks = 1:5, labels = 1:5)+
	scale_color_manual(values = pnw_palette("Bay", n = 10)) +
	scale_x_continuous(breaks=c(1,5,10,20,30))+
	theme_bw()+
	my_theme()

sf_5
ggsave(here('Plots','SI_fig_5.jpg'),sf_5,height = 7,width = 10)
