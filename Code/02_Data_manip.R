
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
# my_colors <- c("#0660f1ff","#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d") 
N_values <- c(4000, 2000, 1000, 600, 350)
N_values_2 <- c(350, 330, 300, 260, 220, 150)
rho_sim <- c(0.01,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20) # Length scales parameter to be simulated
rho_sim_df <- rho_sim %>% as.data.frame() %>% setNames('rho') %>% mutate(rho_it=row_number(.)) 
edna_data <- readRDS(here('Data','edna_data.rds'))
spp <- edna_data %>% distinct(species) %>% slice(-n()) %>% pull(species)


## Raw GP simulated fields ---------------------------------------------------------------------------------
sim_files <- list.files(here('Output','Raw_GP_fileds_simulated'),pattern = '.rds',full.names = TRUE)

sim_df <- map_df(sim_files, function(f) {
	nm <- basename(f)
	m <- str_match(nm, "^(\\d+?)_(\\d+?)\\.rds$")
	
	readRDS(f) %>% as_tibble() %>%
		mutate(rho = as.integer(m[, 2]),it  = as.integer(m[, 3]))
})

saveRDS(sim_df,here('Output','Complete_db','xy_conc_sim_raw.rds'))


## Import simulated parameters ---------------------------------------------------------------------------------
sim_param_files <- list.files(here('Output','simulated_parameters'),pattern = 'sim_GP_param.csv')

sim <- here("Output", "simulated_parameters", sim_param_files) %>% setNames(basename(.)) %>% 
	lapply(read.csv) %>%  
	bind_rows(.id = "source") %>%   
	select(-X) %>% 
	mutate(rho_it = sub("^(\\d+)_.*", "\\1", source),
				 it = sub("^\\d+_(\\d+)_.*", "\\1", source)) %>% 
	mutate(rho_it=as.numeric(rho_it),
				 it=as.numeric(it)) %>% 
	arrange(rho_it,it) %>% 
	mutate(N=9999) %>% 
	mutate(
		parameter = recode(parameter,
											 "0" = "mu_0",
											 "50" = "mu_50",
											 "150" = "mu_150",
											 "300" = "mu_300",
											 "500" = "mu_500",
											 "sigma_sim" = "sigma")) %>% as_tibble()

write.csv(sim,here('Output','Complete_db','sim_GP_param_sim.csv'),row.names = F)


## Import estimated parameters ---------------------------------------------------------------------------------

est_sim_param_1 <- map_dfr(N_values, function(N) {
	list.files(here("Output", paste0("GP_", N)), pattern = "est_GP_param.csv", full.names = TRUE) %>%
		setNames(basename(.)) %>%
		map_dfr(read_csv, .id = "source", show_col_types = FALSE) %>%
		mutate(rho_it = as.numeric(sub("^(\\d+)_.*", "\\1", source)),
					 it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", source)),
					 N = N)}) %>% 
	select(-`...2`) %>% 
	rename(parameter = param, value = mean) %>%
	select(-matches("%"))  %>% 
	as_tibble()

write.csv(est_sim_param_1,here('Output','Complete_db','est_GP_param_sim.csv'),row.names = F)




# est_sim_param_2 <- map_dfr(N_values_2[-1], function(N) {
# 	list.files(here("Output", paste0("GP_", N)), pattern = "est_GP_param.csv", full.names = TRUE) %>%
# 		setNames(basename(.)) %>%
# 		map_dfr(read_csv, .id = "source", show_col_types = FALSE) %>%
# 		mutate(rho_it = as.numeric(sub("^(\\d+)_.*", "\\1", source)),
# 					 it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", source)),
# 					 N = N)}) %>% 
# 	select(-`...2`) %>%
# 	rename(parameter = param, value = mean) %>%
# 	select(-matches("%"))  %>% 
# 	as_tibble()
# 
# write.csv(est_sim_param_2,here('Output','Complete_db','est_GP_param_sim_2.csv'),row.names = F)


## Import species-specific prediction results ---------------------------------------------------------------------------------

est_obs_param_350 <- list.files(here('Output','GP_multifish','Multifish_raw'),pattern = paste0('_350_est_GP_param.csv$')) %>%
	map_df(~ read_csv(here('Output','GP_multifish','Multifish_raw',.x),show_col_types = FALSE) %>%
				 	mutate(file = .x)) %>%
	separate(file,into = c('species','N','type','est','ext1'),sep = '_',remove = FALSE) %>%
	select(-ext1,-est,-type,-`...1`) %>% 
	mutate(N = as.numeric(N))


est_th_obs_param_list <- map(
	set_names(N_values_2[-1]),
	~ {
		folder <- here("Output", "GP_multifish", paste0("Multifish_", .x))
		
		list.files(folder, pattern = "est_GP_param\\.csv$", full.names = TRUE) %>%
			map_df(~ {
				species_name <- spp[str_detect(.x, spp)]
				iteration <- str_extract(basename(.x),"(?<=_)\\d+(?=_est_GP_param)")
				read_csv(.x, show_col_types = FALSE) %>%
					select(-...1) %>% 
					mutate(source  = basename(.x),species = species_name,it = as.integer(iteration))
			}) %>%
			select(-source)
	}
)


est_th_obs_param <- bind_rows(est_th_obs_param_list, .id = "N") %>%
	mutate(N = as.integer(N))

est_obs_param <- bind_rows(est_th_obs_param,
													 est_obs_param_350 %>% mutate(it=0) %>% select(param,mean,species,N,it))

write.csv(est_obs_param,here('Output','Complete_db','est_obs_param.csv'),row.names = F)




## Import simulated spatial data ---------------------------------------------------------------------------------
sim_conc_1_files <- list.files(here('Output','GP_4000'),pattern = 'simulated_data.rds')


sim_conc_1 <- map_df(N_values, function(N) {
	sim_conc_1_files %>%
		file.path(here("Output", paste0("GP_", N)), .) %>%
		map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>%
		mutate(
			rho = sub("^(\\d+)_.*", "\\1", source),
			it  = sub("^\\d+_(\\d+)_.*", "\\1", source),
			rho = as.numeric(rho),
			it  = as.numeric(it),
			N   = N        # add identifier for N
		) %>%
		select(-source) %>%
		as_tibble()
})

saveRDS(sim_conc_1,here('Output','Complete_db','xy_conc_sim_1.rds'))



## Import simulated spatial data ---------------------------------------------------------------------------------
sim_conc_2_list <- map(
	set_names(N_values_2),
	~ {N <- .x  # current N value

		folder <- here("Output", paste0("GP_", N))

		list.files(folder, pattern = "simulated_data\\.rds$", full.names = TRUE) %>%
			map_df(~ {
				file_name <- basename(.x)

				readRDS(.x) %>%
					as_tibble() %>%
					mutate(rho = sub("^(\\d+)_.*", "\\1", file_name),it  = sub("^\\d+_(\\d+)_.*", "\\1", file_name),
						rho = as.numeric(rho),it = as.numeric(it),N=N)})})

sim_conc_2_350 <- list.files(here("Output", "GP_350"), pattern = "simulated_data\\.rds$", full.names = TRUE) %>%
	map_df(~ {
		file_name <- basename(.x)

		readRDS(.x) %>%
			as_tibble() %>%
			mutate(
				rho = as.numeric(sub("^(\\d+)_.*", "\\1", file_name)),
				it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", file_name)),
				N   = 350
			)
	})

sim_conc_2_list[[1]] <- sim_conc_2_350

sim_conc_2 <- bind_rows(sim_conc_2_list,.id='N')

saveRDS(sim_conc_2,here('Output','Complete_db','xy_conc_sim_2.rds'))



# Check that all points belong to the data frame previous to that (insuring thinning)
map2_lgl(
	head(sort(N_values, decreasing = TRUE), -1),
	tail(sort(N_values, decreasing = TRUE), -1),
	~ expand_grid(rho = 1:15, it = 1:20) %>%
		mutate(valid = map2_lgl(rho, it, \(r, t)
														nrow(sim_all %>% filter(N == .y, rho == r, it == t) %>% select(X_utm, Y_utm) %>%
																 	anti_join(sim_all %>% filter(N == .x, rho == r, it == t) %>% select(X_utm, Y_utm),
																 						by = c("X_utm", "Y_utm"))) == 0)) %>%
		summarise(all(valid)) %>% pull()) %>%
	set_names(paste0(head(sort(N_values, decreasing = TRUE), -1),
									 "→",tail(sort(N_values, decreasing = TRUE), -1)))


## Import observation spatial data ---------------------------------------------------------------------------------


obs_conc_list <- map(
	set_names(N_values_2),
	~ {
		folder <- here("Output", "GP_multifish", paste0("Multifish_", .x))
		
		list.files(folder, pattern = "thinned_data\\.rds$", full.names = TRUE) %>%
			map_df(~ {
				file_name <- basename(.x)
				
				# species is everything before "_<number>_thinned_data.rds"
				species <- str_remove(file_name, "_\\d+_thinned_data\\.rds$")
				
				# iteration is the number between the last "_" and "_thinned_data"
				iteration <- str_extract(file_name, "(?<=_)\\d+(?=_thinned_data)")
				
				readRDS(.x) %>%
					# if readRDS returns a tibble/data.frame, this will add the columns
					mutate(
						species = species,
						it      = as.integer(iteration)
					)
			})
	}
)

obs_conc <- bind_rows(obs_conc_list,.id = "N") %>% 
	mutate(N=as.numeric(N)) %>% 
	bind_rows(.,edna_data %>% 
							filter(!species=='Zz_Merluccius productus') %>% 
							mutate(N=350) %>% 
							mutate(it=0) %>% 
							arrange(species,depth))

write_rds(obs_conc,here('Output','Complete_db','xy_conc_obs.rds'))




# NNI compute ---------------------------------------------------------------------------------

w <- c(range(pred_data$X_utm),range(pred_data$Y_utm))
NNI_obs <- obs_conc %>% 
	# filter(species=='Clupea pallasii',it==1,N==330,depth==0) %>% 
	group_by(species, it, N, depth) %>%
	summarise(NNI = compute_NNI(data   = cur_data(),x= X_utm,y= Y_utm,window = w),.groups = "drop") %>% 
	mutate(
		mu_it = case_when(
			depth == 0   ~ 1,
			depth == 50  ~ 2,
			depth == 150 ~ 3,
			depth == 300 ~ 4,
			depth == 500 ~ 5,
			TRUE ~ NA_real_
		))


saveRDS(NNI_obs,here('Output','Complete_db','NNI_obs.rds'))



## Import predicted spatial data ---------------------------------------------------------------------------------
N_values <- c(4000, 2000, 1000, 600, 350)

# est_rds_files <- list.files(here('Output','GP_4000'),pattern = 'pred_GP.rds')

rho_it_name <- rho_sim_df %>% filter(rho_it>1&rho_it<13) %>% pull(rho_it)
it_name <- c(2:20)

est_rds_files <- paste0(rep(rho_it_name,each=length(it_name)),'_',
												rep(it_name,length(rho_it_name)),
												'_pred_GP.rds')

sim_conc_pred_1_list <- map(
	set_names(N_values),
	~ est_rds_files %>%
		file.path(here("Output", paste0("GP_", .x)), .) %>%
		map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>%
		mutate(
			rho = as.numeric(sub("^(\\d+)_.*", "\\1", source)),
			it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", source)),
			N   = .x) %>%
		select(-source))

sim_conc_pred_1 <- map_df(sim_conc_pred_1_list, ~ .x)

sim_conc_pred_1_comp <- map_dfr(
	names(sim_conc_pred_1_list)[names(sim_conc_pred_1_list) != "4000"],
	~ sim_conc_pred_1_list[[.x]] %>%
		mutate(
			Conc_0 = sim_conc_pred_1_list[["4000"]]$conc,
			diff   = abs(Conc_0 - conc),
			comparison = paste0(.x, "→4000")))


write_rds(sim_conc_pred_1,here('Output','Complete_db','xy_conc_pred_sim_1.rds'))
write_rds(sim_conc_pred_1_comp,here('Output','Complete_db','xy_conc_pred_sim_1_comp.rds'))


sim_conc_pred_list <- c(sim_conc_pred_1_list,sim_conc_pred_2_list)

sim_conc_pred <- map_dfr(
	names(sim_conc_pred_list)[names(sim_conc_pred_list) != "4000"],
	~ sim_conc_pred_list[[.x]] %>%
		mutate(
			Conc_0 = sim_conc_pred_list[["4000"]]$conc,
			diff   = abs(Conc_0 - conc),
			comparison = paste0(.x, "→4000")))

## Import prediction results ---------------------------------------------------------------------------------
rho_it_name <- rho_sim_df %>% filter(rho_it>1&rho_it<13) %>% pull(rho_it)
it_name <- c(2:20)
N_values_2 <- c(330, 300, 260, 220, 150)

est_rds_files <- paste0(rep(rho_it_name,each=length(it_name)),'_',
												rep(it_name,length(rho_it_name)),
												'_pred_GP.rds')

sim_conc_pred_2_list <- map(
	set_names(N_values_2),
	~ {
		folder <- if (.x == 350) {here("Output", paste0("GP_", .x))} else {
			here("Output", "GP_sth", paste0("GP_", .x))}
		
		list.files(folder, pattern = "pred_GP.rds", full.names = TRUE) %>%
			map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>%
			mutate(
				rho = as.numeric(sub("^(\\d+)_.*", "\\1", source)),
				it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", source)),
				N   = .x) %>%
			select(-source)})

sim_conc_pred_2_list <- map(
	set_names(N_values_2),
	~ est_rds_files %>%
		file.path(here("Output",'GP_sth', paste0("GP_", .x)), .) %>%
		map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>%
		mutate(
			rho = as.numeric(sub("^(\\d+)_.*", "\\1", source)),
			it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", source)),
			N   = .x) %>%
		select(-source))

sim_conc_pred_2 <- map_df(sim_conc_pred_2_list, ~ .x)

sim_conc_pred_2_comp <- map_dfr(
	names(sim_conc_pred_2_list)[names(sim_conc_pred_2_list) != "350"],
	~ sim_conc_pred_2_list[[.x]] %>%
		left_join(
			sim_conc_pred_2_list[["350"]] %>%
				select(rho, it, X_utm, Y_utm, depth, Conc_0 = conc),
			by = c("rho", "it", "X_utm", "Y_utm", "depth")
		) %>%
		mutate(
			diff = abs(Conc_0 - conc),
			comparison = paste0(.x, "→350")
		)
)

write_rds(sim_conc_pred_2,here('Output','Complete_db','xy_conc_pred_sim_2.rds'))
write_rds(sim_conc_pred_2_comp,here('Output','Complete_db','xy_conc_pred_sim_2_comp.rds'))




## Import prediction results form observations ---------------------------------------------------------------------------------

obs_conc_pred_list <- map(
	set_names(N_values_2), ~{
		folder <- if (.x == 350) {
			here("Output", "GP_multifish", "Multifish_raw")
		} else {
			here("Output", "GP_multifish", paste0("Multifish_", .x))
		}
		
		list.files(folder, pattern = "pred_GP.rds$", full.names = TRUE) %>%
			map_df(~{
				species_name <- spp[str_detect(.x, spp)]
				iteration <- str_extract(basename(.x), "(?<=_)[0-9]+(?=_pred_GP)")
				
				readRDS(.x) %>%
					mutate(source = basename(.x),
								 species = species_name,
								 it = as.integer(iteration))
			}) %>%
			select(-source)
	})

obs_conc_pred <- map_df(obs_conc_pred_list, ~ .x,.id = "N") %>%
	mutate(N = as.numeric(N))

obs_conc_pred_comp <- map_dfr(
	names(obs_conc_pred_list)[names(obs_conc_pred_list) != "350"],
	~ obs_conc_pred_list[[.x]] %>%
		left_join(
			obs_conc_pred_list[["350"]] %>%
				select(species, X_utm, Y_utm, depth, Conc_0 = conc),
			by = c("species", "X_utm", "Y_utm", "depth")
		) %>%
		mutate(
			diff = abs(Conc_0 - conc),
			N = as.numeric(.x),
			comparison = paste0(.x, "→350")
		)
)


write_rds(obs_conc_pred,here('Output','Complete_db','xy_conc_pred_obs.rds'))
write_rds(obs_conc_pred_comp,here('Output','Complete_db','xy_conc_pred_obs_comp.rds'))

