library(PNWColors)
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
# r <- 'Round 3'
# # # sp <- spp[1]
# sp <- 'Scomber japonicus'
# 
# readRDS(here('Plots','Thinning_maps',paste0(r,' 4_',sp),
# 						 paste0(sp,'4_','stan_param_1','.rds')))

dd <- expand_grid(Round=paste0('Round ',c(1:20)),
								 species=spp) %>% 
	mutate(rho=NA,
				 epsilon_10=NA,
				 epsilon_20=NA,
				 epsilon_30=NA,
				 epsilon_40=NA,
				 epsilon_50=NA)

post_list <- vector("list", nrow(dd))


for (j in 1:nrow(dd)) {
	r <- dd$Round[j]
	sp <- dd$species[j]
temp <- readRDS(here('Plots','Thinning_maps',paste0(r,' 4_',sp),
										 paste0(sp,'4_','stan_param_1','.rds'))) %>%
	filter(param=='real_rho') %>% pull(mean)
dd$rho[dd$Round==r&dd$species==sp] <- temp
	
# temp <- readRDS(here('Plots','Thinning_maps',paste0(r,' 4_',sp),
# 										 paste0(sp,'4_','stan_param_1','.rds')))
# dd$rho[dd$Round==r&dd$species==sp] <- temp$mean[temp$param=='rho']*10^temp$mean[temp$param=='mag_rho']

pred_data_stan_4 <- lapply(1:6, function(i) {
	readRDS(here('Plots','Thinning_maps', paste0(r,' 4_', sp),
							 paste0(sp, '4_', 'stan_pred_', i, '.rds')))
})

conc_matrix <- sapply(pred_data_stan_4, function(df) df$conc)
post_list[[j]] <- conc_matrix
differences <- abs(conc_matrix[, 1] - conc_matrix[, 2:ncol(conc_matrix)])
# differences <- abs(exp(conc_matrix[, 1]) - exp(conc_matrix[, 2:ncol(conc_matrix)]))
# differences <- log(abs(exp(conc_matrix[, 1]) - exp(conc_matrix[, 2:ncol(conc_matrix)])))
temp <- differences %>% colMeans()

dd[dd$Round==r&dd$species==sp,4:8] <- as.list(temp)
}

dd %>%
	ggplot()+
	geom_point(aes(x=rho,y=species))+
	scale_x_log10()+
	theme_bw()

prec <- dd %>% 
	group_by(species) %>% 
	mutate(rho_mean=round(mean(rho),1)) %>% 
	# filter(!(species=='Scomber japonicus')) %>%
	# filter(!(species=='Stenobrachius leucopsarus')) %>%
	pivot_longer(cols = epsilon_10:epsilon_50,
										values_to = 'epsilon',
										names_to = 'thinning') %>% 
	mutate(th=if_else(thinning=='epsilon_10',0.1,NA)) %>% 
	mutate(th=if_else(thinning=='epsilon_20',0.2,th)) %>% 
	mutate(th=if_else(thinning=='epsilon_30',0.3,th)) %>% 
	mutate(th=if_else(thinning=='epsilon_40',0.4,th)) %>% 
	mutate(th=if_else(thinning=='epsilon_50',0.5,th)) 

prec_summ <- prec %>% 
	group_by(species,th,rho_mean) %>% 
	summarise(epsilon_mean=mean(epsilon))


prec %>% 
	ggplot()+
	geom_point(aes(x=th,y=epsilon,color=log(rho_mean)),pch=3)+
	geom_point(data=prec_summ,aes(x=th,y=epsilon_mean,color=log(rho_mean)))+
	geom_line(data=prec_summ,aes(x=th,y=epsilon_mean,group = species, color=log(rho_mean)))+
	scale_y_log10()+
	# geom_smooth(aes(x=th,y=epsilon_mean,color=log(rho_mean)),method='lm',se=F)+
	scale_color_gradientn(colors = rev(pnw_palette("Bay", 11, type = "continuous")))+
	theme_bw()


dd %>%
	filter(!(species=='Scomber japonicus')) %>%
	# filter(!(species=='Stenobrachius leucopsarus')) %>%
	pivot_longer(cols = epsilon_10:epsilon_50,
										values_to = 'epsilon',
										names_to = 'thinning') %>%
	mutate(th=if_else(thinning=='epsilon_10',0.1,NA)) %>%
	mutate(th=if_else(thinning=='epsilon_20',0.2,th)) %>%
	mutate(th=if_else(thinning=='epsilon_30',0.3,th)) %>%
	mutate(th=if_else(thinning=='epsilon_40',0.4,th)) %>%
	mutate(th=if_else(thinning=='epsilon_50',0.5,th)) %>%
	ggplot()+
	geom_point(aes(x=rho,y=epsilon,color=as.factor(th)))+
	# geom_smooth(aes(x=rho,y=epsilon,color=as.factor(th)), se=F)+
	# scale_y_log10()+
	# scale_x_log10()+
	theme_bw()


colors <- c("black", "#ff7f01", "#fdbf6f", "#33a02c", "#1f78b4", "#e31b1d")
column_names <- c('Initial_filed',paste("Thinning", c('10%','20%','30%','40%','50%')))


post_data <- post_list[[13]]
# Create the plot
hist(post_data[,1], col = colors[1], border = NA, 
		 main = "Overlapping Histograms of All Columns", 
		 xlab = "Values", ylab = "Frequency",
		 xlim = range(post_data[,1]),
		 breaks = 300)

# Add histograms for remaining columns
for(i in c(2,2)) {
	hist(post_data[,i], col = colors[i], xlim = range(post_data[,1]),
			 breaks = 300, border = NA, add = TRUE)
}

# Add legend
legend("topright", legend = column_names, fill = colors, bty = "n")
