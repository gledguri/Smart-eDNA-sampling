make_index <- function(data, index_variable, index_name) {
	data %>%
		group_by(across(all_of(index_variable))) %>%
		mutate(!!index_name := cur_group_id()) %>%
		ungroup()
}

combine_index <- function(data, index_variables, combined_index_name) {
	data %>% 
		mutate(across(all_of(index_variables), 
									~ paste0(
										tolower(LETTERS[((. - 1) %/% (26^2)) %% 26 + 1]), # First letter
										tolower(LETTERS[((. - 1) %/% 26) %% 26 + 1]),     # Second letter
										tolower(LETTERS[(. - 1) %% 26 + 1])              # Third letter
									), 
									.names = "letter_{col}")) %>% # Convert to letters
		rowwise() %>% # Row-wise operation for combining letters
		mutate(combined = paste(across(starts_with("letter_")), collapse = "_")) %>%
		ungroup() %>%
		group_by(combined) %>%
		mutate(!!combined_index_name := cur_group_id()) %>% # Create numeric index
		ungroup() %>%
		select(-starts_with("letter_"), -combined) # Remove intermediary letter-based columns
}

join_ext_param <- function(stanmod,par){
	extract_param(stanmod,par) %>% 
		rownames_to_column('param') %>% 
		mutate(g_idx = str_extract(param, "\\d+")) %>% 
		mutate(g_idx = as.numeric(g_idx))
}

library(stringr)
library(tidyr)
library(ggplot2)
library(tibble)
library(rstan);options(mc.cores = parallel::detectCores())
library(dplyr)
library(here)
library(QM)

error <- read.csv(here('Data','error_df.csv'))

df <- error %>% select(-c(N_10, N_30, N_50, N_80)) %>% 
	make_index('N','N_idx') 


df <- df %>% pivot_longer(cols = c(deviation_10,deviation_30,deviation_50,deviation_80),
													names_to = 'thinning',
													values_to = 'dev') %>% 
	make_index('thinning','th_idx') %>% 
	combine_index(index_variables = c('N_idx','th_idx'),combined_index_name = 'N_t_idx') %>% 
	as.data.frame() %>% 
	mutate(th=gsub('deviation_','',thinning)) %>% 
	mutate(th=as.numeric(th))

# df <- error %>% filter(n==30) %>% select(deviation_80,a,d,n) %>% setNames(c('dev','a','d','n'))

stan_data <- list(
	N = nrow(df),
	S = df$dev,
	d = df$d,
	n = df$N,
	th = df$th,
	idx = df$N_t_idx,
	n_idx = length(unique(df$N_t_idx)),
	st_points = df$n,
	sigma_prior=10,
	alpha_prior=1,
	beta_prior=1,
	gamma_prior=1
)

# stan_model <- stan_model(here('Code','Error_relationship_with_d.stan'))
stan_model_2 <- stan_model(here('Code','Error_relationship_with_d_2.stan'))

# fit <- sampling(stan_model,data = stan_data,chains = 4,iter = 2000,warmup = 1000)
fit_2 <- sampling(stan_model_2,data = stan_data,chains = 4,iter = 2000,warmup = 1000)

extract_param(fit_2,c('kappa_0','kappa_1','kappa_2','eta_0','eta_1','eta_2'))

summ_fit <- bind_cols(
join_ext_param(fit_2,'alpha') %>% select(mean,g_idx) %>% rename(alpha='mean'),
join_ext_param(fit_2,'kappa_0') %>% select(mean) %>% rename(kappa_0='mean'),
join_ext_param(fit_2,'kappa_1') %>% select(mean) %>% rename(kappa_1='mean'),
join_ext_param(fit_2,'kappa_2') %>% select(mean) %>% rename(kappa_2='mean'),
join_ext_param(fit_2,'eta_0') %>% select(mean) %>% rename(eta_0='mean'),
join_ext_param(fit_2,'eta_1') %>% select(mean) %>% rename(eta_1='mean'),
join_ext_param(fit_2,'eta_2') %>% select(mean) %>% rename(eta_2='mean'),
join_ext_param(fit_2,'beta') %>% select(mean) %>% rename(beta='mean'),
join_ext_param(fit_2,'mu') %>% select(mean) %>% rename(mu='mean'),
df %>% select(-c(mean_deviation_10,mean_deviation_30,mean_deviation_50,mean_deviation_80)))

summ_fit %>% 
	mutate(pred_dev=exp(log(alpha) + beta * log(d))) %>%
	ggplot()+
	geom_point(aes(x=dev,y=pred_dev))+
	geom_abline(intercept = 0,slope=1)+
	labs(x='Observed ε', y='Predicted ε')+
	scale_x_log10()+
	scale_y_log10()+
	theme_bw()+
	theme(axis.title = element_text(size=15),
				axis.text = element_text(size=14))

summ_fit %>% 
	mutate(pred_dev=exp(log(alpha) + beta * log(d))) %>%
	ggplot()+
	geom_point(aes(x=dev,y=pred_dev))+
	geom_abline(intercept = 0,slope=1)

summ_fit %>% 
	mutate(pred_dev=exp(log(alpha) + beta * log(d))) %>%
	ggplot()+
	geom_point(aes(x=d,y=(dev),color=as.factor(th)))+
	geom_point(aes(x=d+0.2,y=(pred_dev),color=as.factor(th)))+
	facet_grid(~N)
	# geom_abline(intercept = 0,slope=1)

est_prec_2 <- summ_fit %>% 
	mutate(pred_dev=exp(log(alpha) + beta * log(d))) %>%
	ggplot()+
	geom_point(aes(x=d,y=(dev),color=as.factor(th)))+
	geom_line(aes(x=d,y=(pred_dev),color=as.factor(th)),lty=2)+
	facet_grid(~N)+
	scale_y_log10()+
	scale_color_manual(
		name = "Spatial thinning",  # Change legend title here
		values = c("10" = "tomato2", "30" = "orange", "50" = "deepskyblue3","80" = "darkblue"),
		labels = c("10%", "30%", "50%", "80%")
	)+
	xlab('Distance of autocorrelation (ρ)')+
	ylab('Precision term (ε)')+
	theme_bw()+
	theme(legend.position = 'none')

legend_deviation <- cowplot::get_legend(
	data.frame(x = rep(1, 4),y = 1:4,group = factor(c("10%", "30%", "50%", "80%"))) %>% 
		ggplot(aes(x = x, y = y, color = group, group = group)) +
		geom_line(size = 1.2) +
		scale_color_manual(
			values = c("10%" = "tomato2", "30%" = "orange", "50%" = "deepskyblue3", "80%" = "darkblue"),
			name = "Spatial thinning",
			labels = c("10%", "30%", "50%", "80%")) +
		theme_minimal()+
		theme(
			legend.key.size = unit(1.5, "lines"),
			legend.text = element_text(size = 12),
			legend.title = element_text(size = 13))
)

cowplot::plot_grid(est_prec_2,legend_deviation,rel_widths = c(7,1))
cowplot::plot_grid(est_prec,legend_deviation,rel_widths = c(7,1))
cowplot::plot_grid(dev_plot,legend_deviation,rel_widths = c(7,1))
