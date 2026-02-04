# library(dplyr)
# library(ggplot2)
# library(here)
# library(fields)
# library(MASS)
# library(PNWColors)
# library(purrr)
# library(tidyr)
# library(tibble)
# library(rstan);options(mc.cores = parallel::detectCores()); rstan_options(auto_write = TRUE)

# select <- dplyr::select

library(purrr)
library(ggplot2)
library(dplyr)
library(readr)
library(here)
library(tidyr)

# combine_csv_files <- function(csv_files) {
# 	map_dfr(csv_files, function(file) {
# 		filename <- basename(file)
# 		
# 		# Extract numbers from filename
# 		numbers <- as.numeric(unlist(regmatches(filename, gregexpr("\\d+", filename))))
# 		
# 		# Extract thinning parameter (th followed by number with optional decimal)
# 		thinning_match <- regmatches(filename, regexpr("th\\d+(?:\\.\\d+)?", filename))
# 		thinning <- if(length(thinning_match) > 0) {
# 			as.numeric(sub("th", "", thinning_match))
# 		} else {
# 			NA  # or some default value if thinning not found
# 		}
# 		
# 		# Read CSV and add index columns
# 		read_csv(file, show_col_types = FALSE) %>%
# 			select(-1) %>%
# 			mutate(
# 				scenario = numbers[1],
# 				iteration = numbers[2],
# 				thinning = thinning
# 			)
# 	})
# }

combine_csv_files <- function(csv_files) {
	map_dfr(csv_files, function(file) {
		filename <- basename(file)
		
		# Extract numbers for scenario & iteration (first two numeric chunks in the filename)
		numbers <- as.numeric(unlist(regmatches(filename, gregexpr("\\d+", filename))))
		
		# Extract thinning parameter: th<digits>[.<digits>] e.g., th0.7
		thinning_match <- regmatches(filename, regexpr("th\\d+(?:\\.\\d+)?", filename))
		thinning <- if (length(thinning_match) > 0) as.numeric(sub("^th", "", thinning_match)) else NA_real_
		
		# Safe read: empty file -> empty tibble
		df <- tryCatch(
			read_csv(file, show_col_types = FALSE),
			error = function(e) tibble()
		)
		
		# Skip truly empty (0-column) files
		if (ncol(df) == 0) {
			warning(sprintf("Skipping empty file (0 columns): %s", file))
			return(tibble())
		}
		
		# If first column exists and is just an unnamed index (e.g., ...1), drop it.
		# (Adjust the predicate if your index column has a specific name.)
		first_col <- names(df)[1]
		if (!is.na(first_col) && first_col %in% c("", "...1", "X1")) {
			df <- df %>% select(-all_of(first_col))
		}
		
		df %>%
			mutate(
				scenario  = numbers[1],
				iteration = numbers[2],
				thinning  = thinning,
				.before = 1
			)
	})
}

combine_rds_files <- function(rds_files) {
	map_dfr(rds_files, function(file) {
		filename <- basename(file)
		
		# Extract numbers for scenario & iteration (first two numeric chunks in the filename)
		numbers <- suppressWarnings(as.numeric(unlist(regmatches(filename, gregexpr("\\d+", filename)))))
		
		# Extract thinning parameter: th<digits>[.<digits>] e.g., th0.7
		# (use a standard capturing group; no perl=TRUE needed)
		thinning_match <- regmatches(filename, regexpr("th\\d+(\\.\\d+)?", filename))
		thinning <- if (length(thinning_match) > 0) as.numeric(sub("^th", "", thinning_match)) else NA_real_
		
		# Safe read
		obj <- try(readRDS(file), silent = TRUE)
		if (inherits(obj, "try-error")) {
			warning(sprintf("Skipping unreadable RDS: %s (%s)", file, as.character(obj)))
			return(tibble())
		}
		
		# Coerce to tibble if needed
		df <-
			if (inherits(obj, "data.frame")) {
				as_tibble(obj)
			} else if (is.list(obj)) {
				# try to turn named list into columns; unnamed list becomes one-column tibble
				if (!is.null(names(obj)) && any(nzchar(names(obj)))) {
					as_tibble(obj, .name_repair = "unique")
				} else {
					tibble(value = I(obj))
				}
			} else if (length(obj) == 0) {
				tibble()
			} else if (is.atomic(obj)) {
				tibble(value = obj)
			} else {
				tibble(value = I(list(obj)))
			}
		
		# Skip truly empty (0-column) tibbles
		if (ncol(df) == 0) {
			warning(sprintf("Skipping empty object (0 columns): %s", file))
			return(tibble())
		}
		
		# If first column exists and is just an unnamed index (e.g., ...1), drop it
		first_col <- names(df)[1]
		if (!is.na(first_col) && first_col %in% c("", "...1", "X1")) {
			df <- df %>% select(-all_of(first_col))
		}
		
		df %>%
			mutate(
				scenario  = numbers[1],
				iteration = numbers[2],
				thinning  = thinning,
				.before = 1
			)
	})
}
colors <- c("black", "#ff7f01", "#fdbf6f", "#33a02c", "#1f78b4", "#e31b1d")

# Get all est_GP_param.csv files
# est_files <- list.files(here('Output'), pattern = "*_est_GP_param\\.csv$", full.names = TRUE)
# sim_files <- list.files(here('Output'), pattern = "*_sim_GP_param\\.csv$", full.names = TRUE)
# est_files <- list.files(here('Output','GP_6.10x0.3'), pattern = "*_est_GP_param\\.csv$", full.names = TRUE)
# sim_files <- list.files(here('Output','GP_6.10x0.3'), pattern = "*_sim_GP_param\\.csv$", full.names = TRUE)
# est_files <- list.files(here('Output','GP_7.2'), pattern = "*_est_GP_param\\.csv$", full.names = TRUE)
# sim_files <- list.files(here('Output','GP_7.2'), pattern = "*_sim_GP_param\\.csv$", full.names = TRUE)
# est_files <- list.files(here('Output','GP_7.1'), pattern = "*_est_GP_param\\.csv$", full.names = TRUE)
# sim_files <- list.files(here('Output','GP_7.1'), pattern = "*_sim_GP_param\\.csv$", full.names = TRUE)
# est_files <- list.files(here('Output','GP_7.2.1'), pattern = "*_est_GP_param\\.csv$", full.names = TRUE)
# sim_files <- list.files(here('Output','GP_7.2.1'), pattern = "*_sim_GP_param\\.csv$", full.names = TRUE)
est_files <- list.files(here('Output','GP_8.2'), pattern = "*_est_GP_param\\.csv$", full.names = TRUE)
sim_files <- list.files(here('Output','GP_8.2'), pattern = "*_sim_GP_param\\.csv$", full.names = TRUE)
pred_files <- list.files(here('Output','GP_8.2'), pattern = "*_pred_GP\\.rds$", full.names = TRUE)

pred_df <- combine_rds_files(pred_files)

est_df <- combine_csv_files(est_files) %>% 
	mutate(param = case_when(
		param == "mu[1]" ~ "mu_0",
		param == "mu[2]" ~ "mu_50",
		param == "mu[3]" ~ "mu_150",
		param == "mu[4]" ~ "mu_300",
		param == "mu[5]" ~ "mu_500",
		TRUE ~ param
	))

sim_df <- combine_csv_files(sim_files) %>% 
	mutate(parameter = case_when(
		parameter == "0" ~ "mu_0",
		parameter == "50" ~ "mu_50",
		parameter == "150" ~ "mu_150",
		parameter == "300" ~ "mu_300",
		parameter == "500" ~ "mu_500",
		parameter == "sigma_sim" ~ "sigma",
		TRUE ~ parameter
	)) %>% 
	rename(param='parameter',
				 sim_val='value')


# sum_sim <- est_df %>% 
# 	left_join(.,sim_df,by=c('param','scenario','iteration','thinning')) %>% 
# 	select(-thinning)

sum_sim <- est_df %>% 
	left_join(.,sim_df,by=c('param','scenario','iteration','thinning'))

# saveRDS(sum_sim,here('Output','GP_8.2_sum_sim.rds'))
# saveRDS(est_df,here('Output','GP_8.2_sim_df.rds'))
# saveRDS(pred_df,here('Output','GP_8.3_pred_df.rds'))




# Figure 1 ------------------------------------------------------------------------------------

sum_sim_2000 <- read_rds(here('Output','GP_8.3_sum_sim.rds'))
sum_sim_600 <- read_rds(here('Output','GP_8.1_sum_sim.rds'))
sum_sim_200 <- read_rds(here('Output','GP_8.2_sum_sim.rds'))

sim <- bind_rows(sum_sim_2000 %>% mutate(N=2000),
										 sum_sim_600 %>% mutate(N=600),
										 sum_sim_200 %>% mutate(N=200))

sum_sim_rho <- sim %>% 
	filter(param=='rho') %>%
	mutate(sim_rho_km=sim_val*100) %>% 
	mutate(N_effective=N*thinning) %>% 
	mutate(est_rho_km=exp(mean)*100) %>% 
	mutate(N = factor(N,
										levels = c(200, 600, 2000),
										labels = c("N = 200", "N = 600", "N = 2000"))) %>% 
	rename(sim_rho='sim_val') %>% 
	rename(est_rho='mean')

pred_df <- read_rds(here('Output','GP_8.2_pred_df.rds'))

pred_df_summ_200 <- pred_df %>% 
	group_by(scenario, iteration, thinning, depth) %>% 
	summarise(mean_conc = mean(conc)) %>% 
	mutate(param = case_when(
		depth == 0   ~ "Depth = 0",
		depth == 50  ~ "Depth = 50", 
		depth == 150 ~ "Depth = 150",
		depth == 300 ~ "Depth = 300",
		depth == 500 ~ "Depth = 500",
		TRUE         ~ as.character(depth)
	)) %>% 
	mutate(N= 'N = 200')

pred_df <- read_rds(here('Output','GP_8.1_pred_df.rds'))
pred_df_summ_600 <- pred_df %>% 
	group_by(scenario, iteration, thinning, depth) %>% 
	summarise(mean_conc = mean(conc)) %>% 
	mutate(param = case_when(
		depth == 0   ~ "Depth = 0",
		depth == 50  ~ "Depth = 50", 
		depth == 150 ~ "Depth = 150",
		depth == 300 ~ "Depth = 300",
		depth == 500 ~ "Depth = 500",
		TRUE         ~ as.character(depth)
	)) %>% 
	mutate(N= 'N = 600')

pred_df <- read_rds(here('Output','GP_8.3_pred_df.rds'))
pred_df_summ_2000 <- pred_df %>% 
	group_by(scenario, iteration, thinning, depth) %>% 
	summarise(mean_conc = mean(conc)) %>% 
	mutate(param = case_when(
		depth == 0   ~ "Depth = 0",
		depth == 50  ~ "Depth = 50", 
		depth == 150 ~ "Depth = 150",
		depth == 300 ~ "Depth = 300",
		depth == 500 ~ "Depth = 500",
		TRUE         ~ as.character(depth)
	)) %>% 
	mutate(N= 'N = 2000')


pred_df_summ <- 
	bind_rows(pred_df_summ_200,
						pred_df_summ_600,
						pred_df_summ_2000)

# fig_1 <-
sum_sim_rho %>% 
	ggplot()+
	# geom_errorbar(aes(x=sim_val,ymin=`2.5%`,ymax=`97.5%`),color='grey',width=0.1,alpha=0.4)+
	labs(x='Simulated ρ (km)',y='Estimated ρ (km)',
			 colour = "Thinning %")+
	geom_point(aes(x=sim_rho_km,y=est_rho_km),size=0.6,alpha=0.2)+
	scale_x_log10()+
	scale_y_log10()+
	geom_smooth(data= sum_sim_rho %>% 
								group_by(sim_rho_km,thinning,N) %>%
								summarise(est_rho_km_mean=mean(est_rho_km)),
						 aes(x=sim_rho_km,y=est_rho_km_mean,color=factor(thinning)),size=1,alpha=1,se=F)+
	scale_color_manual(values = c('#e31b1d','#1f78b4','#33a02c','#fdbf6f','#ff7f01','black')) +
	geom_abline(intercept = 0,slope=1,lty=2)+
	facet_wrap(~N, scales = "free")+
	theme_bw()

# fig_1 <-
sum_sim_rho %>% 
	mutate(N_effective=as.factor(N_effective)) %>%
	ggplot()+
	labs(x='Simulated ρ (km)',y='Estimated ρ (km)',
			 colour = "Number of samples")+
	geom_point(aes(x=sim_rho_km,y=est_rho_km),size=0.6,alpha=0.2)+
	geom_smooth(aes(x=sim_rho_km,y=est_rho_km,color=N_effective),se=F)+
	scale_x_log10()+
	scale_y_log10()+
	scale_color_manual(
		# values = pnw_palette("Bay", n = length(unique(sum_sim_rho$N_effective)), type = "continuous")
		values = moma.colors('ustwo',length(unique(sum_sim_rho$N_effective)),
												 type = 'continuous',direction = -1)) +
	geom_abline(intercept = 0,slope=1,lty=2)+
	theme_bw()

fig_1
# ggsave(here('Plots','Paper_figs','Figure_1_3.jpg'), fig_1,width = 12,height = 5)
ggsave(here('Plots','Paper_figs','Figure_1_3.jpg'), fig_1,width = 12,height = 8)

fig_2 <-
	sim %>%
	filter(!(param%in%c('alpha','rho','sigma'))) %>%
	mutate(N = factor(N, levels = c(200, 600, 2000), 
										labels = c("N = 200", "N = 600", "N = 2000")),
				 param = case_when(
				 	param == "mu_0" ~ "Depth = 0",
				 	param == "mu_50" ~ "Depth = 50", 
				 	param == "mu_150" ~ "Depth = 150",
				 	param == "mu_300" ~ "Depth = 300",
				 	param == "mu_500" ~ "Depth = 500",
				 	TRUE ~ param
				 )) %>%
		left_join(.,pred_df_summ,by=c('scenario','iteration','thinning','param','N')) %>% 
		# filter(N=='N = 200') %>% 
		filter(thinning==1) %>% 
	ggplot()+
	geom_errorbar(aes(x=mean_conc,ymin=`25%`,ymax=`75%`),color='grey',width=0.1,alpha=0.4)+
	xlab('Simulated μ')+
	ylab('Estimated μ')+
	geom_point(aes(x=mean_conc,y=mean),size=1,alpha=0.6)+
	geom_abline(intercept = 0,slope=1,lty=2)+
	facet_wrap(~param*N, scales = "free", nrow = 5, ncol = 3)+ 
	theme_bw()

fig_2
ggsave(here('Plots','Paper_figs','Figure_2_2.jpg'), fig_2,width = 12,height = 14)


# fig_3 <- 
# sum_sim %>% 
# 	filter(!is.na(sim_val)) %>% 
# 	mutate(sim_val=if_else(param=='rho',log(sim_val),sim_val))
	
	

# Recursive

# sum_sim %>% 
# 	# filter(!(param=='alpha'|param=='sigma')) %>%
# 	# filter(thinning==1) %>%
# 	filter(param=='rho') %>%
# 	# filter((scenario%in%sc&iteration%in%it)) %>%
# 	# mutate(mean=if_else(param=='rho',exp(mean),mean)) %>%
# 	mutate(sim_val=if_else(param=='rho',log(sim_val),sim_val)) %>%
# 	# mutate(`2.5%`=if_else(param=='rho',exp(`2.5%`),`2.5%`)) %>%
# 	# mutate(`97.5%`=if_else(param=='rho',exp(`97.5%`),`97.5%`)) %>%
# 	# filter(!(param=='real_rho'&mean>30)) %>%
# 	ggplot()+
# 	# geom_errorbar(aes(x=sim_val,ymin=`2.5%`,ymax=`97.5%`),color='grey')+
# 	# geom_point(aes(x=sim_val,y=mean),size=2,alpha=0.4)+
# 	geom_point(aes(x=sim_val,y=mean),size=2,alpha=0.6)+
# 	# geom_point(aes(x=sim_val,y=mean,color= as.factor(scenario)),size=2,alpha=1.2)+
# 	# geom_point(aes(x=sim_val,y=mean,color= as.factor(scenario),shape=as.factor(thinning)),size=2,alpha=0.9)+
# 	geom_abline(intercept = 0,slope=1,lty=2)+
# 	facet_wrap(~N, scales = "free")+
# 	# facet_wrap(~thinning, scales = "free")+
# 	theme_bw()

# END Recursive

# Figure 2 ------------------------------------------------------------------------------------

# dd <- expand_grid(iteration=paste0(c(1:20)),
# 									scenario=c(1:12)) %>%
# 									# scenario=c(1:6,11:12)) %>%
# 									# scenario=c(1,2,5)) %>%
# 	mutate(epsilon_10=NA,
# 				 epsilon_20=NA,
# 				 epsilon_30=NA,
# 				 epsilon_40=NA,
# 				 epsilon_50=NA)
# 
# post_list <- vector("list", nrow(dd))
# 
# for (j in 1:nrow(dd)) {
# 	it <- dd$iteration[j]
# 	sc <- dd$scenario[j]
# 	
# 	pred_GP <- lapply(seq(1, 0.5, by = -0.1), function(i) {
# 		readRDS(here('Output', 'GP_8.2',paste0(sc,'_',it,'_th',i,'_pred_GP.rds')))})
# 	
# 	conc_matrix <- sapply(pred_GP, function(df) df$conc)
# 	post_list[[j]] <- conc_matrix
# 	differences <- abs(conc_matrix[, 1] - conc_matrix[, 2:ncol(conc_matrix)])
# 	temp <- differences %>% colMeans()
# 	
# 	dd[dd$iteration==it&dd$scenario==sc,3:7] <- as.list(temp)
# }

# saveRDS(dd,here('Output','GP_8.2_pred_thinn.rds'))

dd <- bind_rows(read_rds(here('Output','GP_8.1_pred_thinn.rds')) %>% mutate(N=600),
					read_rds(here('Output','GP_8.3_pred_thinn.rds')) %>% mutate(N=2000),
					read_rds(here('Output','GP_8.2_pred_thinn.rds')) %>% mutate(N=200)
					)


sim_prec_df <- dd %>% 
	mutate(iteration=as.numeric(iteration)) %>% 
	left_join(.,sim_df %>% group_by(scenario,iteration,param) %>% 
							summarise(sim_val=mean(sim_val)) %>% 
							filter(param=='rho'),
						by=c('scenario','iteration')) %>% 
	left_join(.,sim %>% filter(param=='rho') %>% filter(thinning==0.5) %>% select(scenario,iteration,mean),
						by=c('scenario','iteration')) %>%
	rename(est_val='mean') %>%
	# filter(!is.na(epsilon_10)) %>%
	pivot_longer(cols = c('epsilon_10','epsilon_20','epsilon_30','epsilon_40','epsilon_50'),
							 names_to = 'epsilon',
							 values_to = 'value') %>%
	mutate(epsilon=gsub('epsilon_','',epsilon)) %>% 
	mutate(est_val=exp(est_val))

fig_3 <-
	sim_prec_df %>% 
		ggplot()+
	geom_point(aes(x=epsilon,y=value,colour = as.factor(sim_val),group = sim_val),alpha=0.4,size=2)+
	geom_smooth(aes(x=epsilon,y=value,colour = as.factor(sim_val),group = sim_val),se=F,span=1,lty=2,size=0.8)+
	# scale_x_log10()+
	labs(x='Estimated (length-scale parameter) ρ',y='Error ε',
			 colour = "ρ - Length scale parameter")+
	# scale_color_manual(values = colors[-1]) +
	scale_color_manual(values = moma.colors('ustwo',length(unique(sim_prec_df$sim_val)))) +
 	# ylim(0,1.5)+
 	# xlim(0,30)+
	facet_wrap(~N)+
	theme_bw()

fig_3
ggsave(here('Plots','Paper_figs','Figure_3_2.jpg'), fig_3,width = 14,height = 8)



# Figure 4 ------------------------------------------------------------------------------------
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

est_gp_param_sp_list
est_gp_param_sp <- bind_rows(est_gp_param_sp_list)

est_gp_param_sp %>% 
	filter(param=='rho') %>%
	ggplot()+
	geom_point(aes(x=thinning,y=exp(mean),color=sp))+
	geom_smooth(aes(x=thinning,y=exp(mean),color=sp),se=F)

fig_supp_1 <- est_gp_param_sp %>% 
	filter(!param%in%c('alpha','rho','sigma')) %>%  
	# filter(sp=='Sebastes entomelas') %>% 
	ggplot()+
	geom_point(aes(x=thinning,y=mean,color=as.factor(iteration)))+
	geom_smooth(aes(x=thinning,y=mean,color=as.factor(iteration)),se=F,size=0.5)+
	facet_wrap(~sp+param,scales = "free_y",
						 ncol = 5, 
						 nrow = 12)+
	theme_bw()

ggsave(here('Plots','Paper_figs','Figure_supp_1.jpg'), fig_supp_1,width = 16,height = 28)


dd_sp <- expand_grid(iteration=paste0(c(1:30)),
									scenario=spp) %>%
	mutate(epsilon_10=NA,
				 epsilon_20=NA,
				 epsilon_30=NA,
				 epsilon_40=NA,
				 epsilon_50=NA)

post_list <- vector("list", nrow(dd_sp))

for (j in 1:nrow(dd_sp)) {
	it <- dd_sp$iteration[j]
	sc <- dd_sp$scenario[j]

	pred_GP <- lapply(seq(1, 0.5, by = -0.1), function(i) {
		readRDS(here('Output', 'multifish_spp',paste0(sc,'_',it,'_th',i,'_pred_GP.rds')))})

	conc_matrix <- sapply(pred_GP, function(df) df$conc)
	post_list[[j]] <- conc_matrix
	differences <- abs(conc_matrix[, 1] - conc_matrix[, 2:ncol(conc_matrix)])
	temp <- differences %>% colMeans()

	dd_sp[dd_sp$iteration==it&dd_sp$scenario==sc,3:7] <- as.list(temp)
}

rho_sp <- est_gp_param_sp %>% filter(thinning==1) %>% filter(param=='rho') %>% group_by(sp) %>% summarise(rho=mean(mean))

dat_prec_df <- 
	dd_sp %>% 
	mutate(iteration=as.numeric(iteration)) %>% 
	# filter(!scenario=='Thaleichthys pacificus') %>% 
	pivot_longer(cols = c('epsilon_10','epsilon_20','epsilon_30','epsilon_40','epsilon_50'),
							 names_to = 'epsilon',
							 values_to = 'value') %>%
	mutate(epsilon=gsub('epsilon_','',epsilon)) %>% 
	rename(sp='scenario') %>% 
	left_join(.,rho_sp,by='sp') 
fig_4 <-
	dat_prec_df %>% 
	ggplot()+
	geom_point(aes(x=exp(rho),y=value,colour = as.factor(epsilon),group = epsilon),alpha=0.4,size=2)+
	geom_smooth(aes(x=exp(rho),y=value,colour = as.factor(epsilon),group = epsilon),se=F,span=1,lty=2,size=0.8)+
	scale_x_log10()+
	labs(x='Estimated (length-scale parameter) ρ',y='Error ε',
			 colour = "Thinning %")+
	scale_color_manual(values = colors[-1]) +
	ylim(0,1.5)+
	# xlim(0,30)+
	# facet_wrap(~N)+
	theme_bw()
	
fig_4
ggsave(here('Plots','Paper_figs','Figure_4.jpg'), fig_4,width = 6,height = 8)



# Figure 5 ------------------------------------------------------------------------------------

fig_5 <- 
sim_prec_df %>% 
	ggplot()+
	geom_smooth(aes(x=sim_val,y=value,colour = as.factor(epsilon),group = epsilon),se=F,span=1,lty=2,size=0.8)+
	geom_point(data=dat_prec_df %>% 
							filter(!sp=='Thaleichthys pacificus') %>% 
							group_by(sp,epsilon) %>% 
						 	summarise(mean_est_epsilon=mean(value),
						 						mean_est_rho=mean(rho)) %>% 
						 	mutate(N=200) %>% 
						 	ungroup() %>% 
						 	select(N,epsilon,mean_est_rho,mean_est_epsilon),
		aes(x=exp(mean_est_rho),y=mean_est_epsilon,colour = as.factor(epsilon)))+
	# scale_x_log10()+
	labs(x='Estimated (length-scale parameter) ρ',y='Error ε',
			 colour = "Thinning %")+
	scale_color_manual(values = colors[-1]) +
	ylim(0,1.5)+
	# xlim(0,30)+
	facet_wrap(~N)+
	theme_bw()

fig_5
ggsave(here('Plots','Paper_figs','Figure_5.jpg'), fig_5,width = 14,height = 8)




