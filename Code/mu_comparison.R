comp <- vector('list');i=0
for (h in 1:20) {
 for (j in 1:15) {
  i=i+1
  name_file <- paste0(j,'_',h)
  sim_data_raw <- readRDS(here('Output','Raw_GP_fileds_simulated',paste0(name_file,'.rds')))
  comp[[i]] <- sim_data_raw %>% group_by(depth_cat) %>% 
   summarise(mu_est=mean(z),
  mu_sim=mean(mu_sim))
 }
}

comp %>% bind_rows() %>%
 ggplot(aes(x=mu_sim,y=mu_est))+
 geom_point(alpha=0.5)+
 geom_abline(slope=1,intercept=0,colour='red')+
 facet_wrap(~depth_cat)+
 labs(x='Simulated mean',y='Estimated mean')+
 theme_minimal()