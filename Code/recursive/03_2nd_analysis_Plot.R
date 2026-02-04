
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
library(stringr)
library(tibble)
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















# Import data ---------------------------------------------------------------------------------

## Declare variables ---------------------------------------------------------------------------------
my_colors <- c("#0660f1ff","#1f78b4","#33a02c","#fdbf6f", "#ff7f01", "#e31b1d")
my_sp_colors <- c('#6a3d9a','#cab2d6','#ff7f01','#fdbf6f','#e31b1d','#fb9a99','#33a02c','#b2df8a','#1f78b4','#009999','#999900', '#a6cee3')
N_values <- c(350, 330, 300, 260, 220, 150, 100)
rho_sim <- c(0.01,0.1,0.25,0.5,0.75,1,2,3,4,5,8,10,12,15,20) # Length scales parameter to be simulated
rho_sim_df <- rho_sim %>% as.data.frame() %>% setNames('rho') %>% mutate(rho_it=row_number(.)) 
edna_data <- readRDS(here('Data','edna_data.rds'))
spp <- edna_data %>% distinct(species) %>% slice(-n()) %>% pull(species)


## Import prediction results ---------------------------------------------------------------------------------
pred_sim_list <- map(
  set_names(N_values),
  ~ {
    folder <- if (.x == 350) {here("Output", paste0("GP_", .x))} else {
      here("Output", "GP_sth", paste0("GP_", .x))}

    list.files(folder, pattern = "pred_GP.rds", full.names = TRUE) %>%
      map_df(~ readRDS(.x) %>% mutate(source = basename(.x))) %>%
      mutate(
        rho_it = as.numeric(sub("^(\\d+)_.*", "\\1", source)),
        it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", source)),
        N   = .x) %>%
      select(-source)})

pred_sim <- map_dfr(
  names(pred_sim_list)[names(pred_sim_list) != "350"],
  ~ pred_sim_list[[.x]] %>%
      left_join(
        pred_sim_list[["350"]] %>%
          select(rho_it, it, X_utm, Y_utm, depth, Conc_0 = conc),
        by = c("rho_it", "it", "X_utm", "Y_utm", "depth")
      ) %>%
      mutate(
        diff = abs(Conc_0 - conc),
        comparison = paste0(.x, "→350")
      )
)

pred_sim_summ <- pred_sim %>% 
	group_by(it,rho_it,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,rho_sim_df,by=c('rho_it')) %>% 
  mutate(rho=rho*100)

## Import species-specific prediction results ---------------------------------------------------------------------------------

est_rho_sp <- list.files(here('Output','GP_multifish','Multifish_raw'),pattern = paste0('_350_est_GP_param.csv$')) %>%
  map_df(~ read_csv(here('Output','GP_multifish','Multifish_raw',.x),show_col_types = FALSE) %>%
           mutate(file = .x)) %>%
  separate(file,into = c('species','N','type','est','ext1'),sep = '_',remove = FALSE) %>%
  select(-ext1,-est,-type,-`...1`) %>% 
  mutate(N = as.numeric(N)) %>% 
  filter(param=='rho') %>% 
  group_by(species,param) %>%
  summarise(
    est_rho = mean(mean,na.rm=TRUE)) %>% 
  select(-param)


# pred_obs_list <- map(
#   set_names(N_values),~ {
#     folder <- if (.x == 350) {here("Output", "GP_multifish", "Multifish_raw")} else {
#       here("Output", "GP_multifish", paste0("Multifish_", .x))}

#     list.files(folder, pattern = "pred_GP.rds$", full.names = TRUE) %>%
#       map_df(~ {
#         species_name <- spp[str_detect(.x, spp)] 
#         readRDS(.x) %>%
#           mutate(source = basename(.x),
#                  species = species_name)}) %>%
#       select(-source)})

pred_obs_list <- map(
  set_names(N_values), ~{
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

pred_obs <- map_dfr(
  names(pred_obs_list)[names(pred_obs_list) != "350"],
  ~ pred_obs_list[[.x]] %>%
      left_join(
        pred_obs_list[["350"]] %>%
          select(species, X_utm, Y_utm, depth, Conc_0 = conc),
        by = c("species", "X_utm", "Y_utm", "depth")
      ) %>%
      mutate(
        diff = abs(Conc_0 - conc),
								N = as.numeric(.x),
        comparison = paste0(.x, "→350")
      )
)

pred_obs_summ <- pred_obs %>% 
	group_by(species,it,N) %>% 
	summarise(mean_diff=mean(diff)) %>% 
	left_join(.,est_rho_sp,by=c('species')) %>% 
  mutate(rho=exp(est_rho)*100) 
  # filter(!species=='Leuroglossus stilbius') %>% 
  # filter(!species=='Microstomus pacificus') 
  # filter(!species=='Thaleichthys pacificus') %>% 
  # filter(!species=='Trachurus symmetricus') %>% 
  # filter(!species=='Tarletonbeania crenularis')

pred_summ <- bind_rows(
	pred_obs_summ %>% mutate(source='Guri et al. (2025)'), 
	pred_sim_summ %>% mutate(source='Simulation'))





















p4 <- pred_summ %>%
  filter(!rho_it%in%c(1,2,13:15)) %>% 
  # filter(N==100) %>% 
  # filter(!species=='Microstomus pacificus') %>% 
  ggplot() +
  geom_point(aes(x = rho, y = mean_diff, color = factor(N)),alpha=0.25) +
  geom_smooth(aes(x = rho, y = mean_diff, color = factor(N)),se = FALSE,span=2) +
  scale_x_log10() +
  # scale_y_log10() +
  scale_color_manual(values = my_colors) +
  labs(
    # y = expression("Absolute mean ratio of error - " * epsilon * " (" * frac(predicted, simulated) * ")"),
    	y = expression("Mean prediction error - " * epsilon * 
               " = |" * hat(Z[i]) - hat(S[i*e]) * "|"),
    x = expression("Spatial autocorrelation - " * rho * " (km)"),
    color = "Number of\nsamples\ndeployed",
    shape = "Data source"
  ) +
  facet_wrap(~source) +
  theme_bw() +
  my_theme()
# p4
# ggsave(here('Plots','Figure_4.jpg'),p4,width=16,height = 10)


















# Take a scenario where we know the rho but let's pretend we don't know 
# and estimate it from N=350 batch to see how wrong we could be on Guri et al. (2025)

est <- map_dfr(N_values, function(N) {
  list.files(here("Output", paste0("GP_", N)), pattern = "est_GP_param.csv", full.names = TRUE) %>%
    setNames(basename(.)) %>%
    map_dfr(read_csv, .id = "source", show_col_types = FALSE) %>%
    mutate(rho_it = as.numeric(sub("^(\\d+)_.*", "\\1", source)),
											it  = as.numeric(sub("^\\d+_(\\d+)_.*", "\\1", source)),
      					N = N)}) %>% 
  select(-`...2`) %>% 
  rename(parameter = param, value = mean) %>%
  mutate(
    parameter = recode(parameter,
                       "mu[1]" = "mu_0",
                       "mu[2]" = "mu_50",
                       "mu[3]" = "mu_150",
                       "mu[4]" = "mu_300",
                       "mu[5]" = "mu_500"),
    value = if_else(parameter == "rho", exp(value), value),
    lo_q  = if_else(parameter == "rho", exp(`2.5%`), `2.5%`),
    up_q  = if_else(parameter == "rho", exp(`97.5%`), `97.5%`)) %>%
  select(-matches("%"))  %>% 
  mutate(it=as.numeric(it)) %>%
  mutate(rho_it=as.numeric(rho_it)) %>%
  as_tibble()


pred_sim_summ_2 <- pred_sim_summ %>% left_join(est %>% filter(parameter=='rho')  %>% 
  filter(!rho_it==1)  %>%  
  filter(N==350)  %>% 
  rename(rho_est_sim='value') %>% 
  select(rho_est_sim,it,rho_it),by=c('it','rho_it'))  %>% 
  mutate(rho_est_sim=as.numeric(rho_est_sim))  %>% 
  mutate(rho=100*(rho_est_sim)) %>% select(-rho_est_sim)



pred_summ <- bind_rows(
	pred_obs_summ %>% mutate(source='Guri et al. (2025)'), 
	pred_sim_summ %>% mutate(source='Simulation'),
	pred_sim_summ_2 %>% mutate(source='Simulation_est'))

p5 <- 
pred_summ %>%
  filter(!rho_it%in%c(1,2,13:15)) %>% 
  # filter(!rho_it%in%c(1)) %>% 
  # filter(N==100) %>% 
  # filter(!species=='Microstomus pacificus') %>% 
  ggplot() +
  geom_point(aes(x = rho, y = mean_diff, color = factor(N)),alpha=0.25) +
  geom_smooth(aes(x = rho, y = mean_diff, color = factor(N)),se = FALSE,span=2) +
  scale_x_log10() +
  # scale_y_log10() +
  scale_color_manual(values = my_colors) +
  labs(
    # y = expression("Absolute mean ratio of error - " * epsilon * " (" * frac(predicted, simulated) * ")"),
    	y = expression("Mean prediction error - " * epsilon * 
               " = |" * hat(Z[i]) - hat(S[i*e]) * "|"),
    x = expression("Spatial autocorrelation - " * rho * " (km)"),
    color = "Number of\nsamples\ndeployed",
    shape = "Data source"
  ) +
  facet_wrap(~source) +
  theme_bw() +
  my_theme()
# ggsave(here('Plots','Figure_5.jpg'),p5,width=16,height = 10)

