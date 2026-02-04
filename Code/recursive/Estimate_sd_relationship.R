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
	as.data.frame()

	

# df <- error %>% filter(n==30) %>% select(deviation_80,a,d,n) %>% setNames(c('dev','a','d','n'))

stan_data <- list(
	N = nrow(df),
	S = df$dev,
	d = df$d,
	idx = df$N_t_idx,
	n_idx = length(unique(df$N_t_idx)),
	st_points = df$n,
	sigma_prior=10,
	alpha_prior=0.5,
	intercept_prior=0.5,
	beta_prior=0.4,
	gamma_prior=1
)

df %>% 
	ggplot()+
	geom_point(aes(x=th,y=dev,color=as.factor(d)))+
	# geom_line(data=df %>% group_by(d,N,th) %>% summarise(dev=mean(dev)),
	# 					aes(x=th,y=dev,color=as.factor(d)))+
	geom_smooth(aes(x=th,y=dev,color=as.factor(d)), method = "lm",se=F)+
	facet_wrap(~N)+
	scale_y_log10()+
	theme_bw()+
	theme(legend.position = 'none')

stan_model <- stan_model(here('Code','Error_relationship_with_d.stan'))

fit <- sampling(stan_model,data = stan_data,chains = 4,iter = 2000,warmup = 1000)


# dd_80 <- extract_param(fit) %>% as.data.frame() %>% 
	# slice(1:3)
	
out <- df %>% 
	left_join(.,join_ext_param(fit,'alpha') %>% rename(alpha='mean') %>%
							select(alpha,g_idx),
						by=c('N_t_idx'='g_idx')) %>% 
	left_join(.,join_ext_param(fit,'intercept') %>% rename(intercept='mean') %>%
							select(intercept,g_idx),
						by=c('N_t_idx'='g_idx')) %>%
	left_join(.,join_ext_param(fit,'beta') %>% rename(beta='mean') %>%
							select(beta,g_idx),
						by=c('N_t_idx'='g_idx'))

out %>% 
	ggplot()+
	# geom_line(aes(x=N,y=alpha, colour = thinning))+
	geom_point(aes(x=thinning,y=alpha, colour = as.factor(N)))+
	# geom_point(aes(x=thinning,y=beta, colour = as.factor(N)))+
	# geom_line(aes(x=N,y=beta, colour = thinning))+
	# geom_line(aes(x=N,y=gamma, colour = thinning))+
	theme_bw()

out_s <- out %>% group_by(N,thinning) %>% 
	summarise(alpha=mean(alpha),
						intercept=mean(intercept),
						beta=mean(beta))
out_s

sim <- expand_grid(d=c(1,1.5,3,4,6),
									 N=c(100,400,900),
									 thinning=c('deviation_10','deviation_30','deviation_50','deviation_80')) %>% 
	left_join(.,out_s,by=c('N','thinning')) %>% 
	# mutate(sd=alpha*d^beta+exp(gamma)*d)
	# mutate(sd=alpha*d^beta)
	mutate(sd=intercept+((beta*d)/alpha-d))
	# mutate(sd=exp(alpha)/(exp(beta)+d*log(d+0.1)))
	# mutate(sd=exp(-10.06)/(0.001+d*log(d+0.0000000001)))
	# mutate(sd=alpha*((d*1.05)^beta))
	# mutate(sd=alpha*((d*gamma)^beta))


sim %>% 
	ggplot()+
	geom_line(aes(x=d,y=sd,color=thinning))+
	facet_wrap(~N)+
	# scale_y_log10()+
	theme_bw()+
	theme(legend.position = 'none')

est_prec <- deviation %>% 
	ggplot()+
	# geom_point(data=sim,aes(x=,y=))
	geom_line(data=,aes(x=d,y=mean_deviation_10),color='tomato2',size=1)+
	geom_line(aes(x=d,y=mean_deviation_30),color='orange',size=1)+
	geom_line(aes(x=d,y=mean_deviation_50),color='deepskyblue3',size=1)+
	geom_line(aes(x=d,y=mean_deviation_80),color='darkblue',size=1)+
	geom_line(data=sim,aes(x=d,y=sd,color=as.factor(th)),pch=4,lty=2,alpha=0.6)+
	scale_color_manual(
		name = "Spatial thinning",  # Change legend title here
		values = c("10" = "tomato2", "30" = "orange", "50" = "deepskyblue3","80" = "darkblue"),
		labels = c("10%", "30%", "50%", "80%")
	)+
	# geom_point(aes(x=d,y=deviation_10),color='tomato2',pch=4,size=3,alpha=0.6)+
	# geom_point(aes(x=d,y=deviation_30),color='orange',pch=4,size=3,alpha=0.6)+
	# geom_point(aes(x=d,y=deviation_50),color='deepskyblue3',pch=4,size=3,alpha=0.6)+
	# geom_point(aes(x=d,y=deviation_80),color='darkblue',pch=4,size=3,alpha=0.6)+
	facet_wrap(~N)+
	# xlab('Length-scale correlation (distance of autocorrelation)')+
	# ylab('Magnitude of error term (ε) between the thinned and unthinned GP')+
	xlab('Distance of autocorrelation (ρ)')+
	ylab('Precision term (ε)')+
	scale_x_continuous(breaks=c(1,2,3,4,5,6))+
	scale_y_log10()+
	# scale_x_log10()+
	theme_bw()+
	theme(legend.position = 'none')

cowplot::plot_grid(est_prec,legend_deviation,rel_widths = c(7,1))
