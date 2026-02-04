#
# # # Gaussian Process ----------------------------------------------------------------------------
library(spatstat.geom)
library(spatstat)
library(dplyr)
library(ggspatial)
library(ggplot2)
library(tibble)
library(tidyr)
library(purrr)
library(sf)
library(sp)
library(MASS)     # For mvrnorm
library(fields)   # For rdist for distance matrix
library(mgcv)
library(here)
library(spdep)

cond_raw <- expand.grid(a = c(1), d = c(10,15,30,40), n = c(400,1600))

# cond_raw <- data.frame(a=1,d=c(5,15,30,40),n=1600)
cond <- bind_rows(replicate(20, cond_raw, simplify = FALSE))

dim(cond)

ripley_i <- list()
moran_i <- list()
for (j in 1:nrow(cond)) {

# Simulate coordinates
# set.seed(12345)
n_points <- cond$n[j]
x <- rep(seq(0,100,length.out=sqrt(n_points)),sqrt(n_points))
y <- rep(seq(0,100,length.out=sqrt(n_points)),each=sqrt(n_points))

coords <- cbind(x, y) %>% as.data.frame()


# Create spatial covariance matrix using exponential or Gaussian kernel
# Choose kernel - Gaussian: exp(- (d^2) / (2 * l^2))
alpha <- cond$a[j]        # Variance
length_scale <- cond$d[j]  # Controls smoothness
# length_scale <- 45  # Controls smoothness

# set.seed(12345)
# thining <- sample(seq(1,cond$n[j]),cond$n[j]*0.5,replace = F)
# coords <- coords[thining,] %>% arrange(y,x)
# n_points <- n_points-length(thining)


# Compute distance matrix
dist_mat <- rdist(coords)
cov_mat <- alpha * exp(- (dist_mat^2) / (2 * length_scale^2))

# center_point <- c(50, 50)
# center_effect <- function(point, center, spread = 15) {
# 	dist_to_center <- sqrt(sum((point - center)^2))
# 	return(3 * exp(-(dist_to_center^2) / (2 * spread^2)))
# }
# 
# # Generate mean values for each point based on distance to center
# mu <- apply(coords, 1, function(point) center_effect(point, center_point))
# z <- mvrnorm(n = 1,mu = mu, Sigma = cov_mat)

# Simulate GP with the specified mean and covariance
z <- mvrnorm(mu = rep(0,n_points), Sigma = cov_mat)



# Ripley's K function -------------------------------------------------------------------------

study_area <- owin(c(min(coords$x), max(coords$x)),
									 c(min(coords$y), max(coords$y)))

k <- seq(min(z),max(z),length.out=11)

fish_data_list <- list()
k_ripley_list <- list()

for (i in 1:length(k)) {
	C <- k[i]
	fish_data_list[[i]] <- cbind(coords,z) %>%
		filter(z>C) %>%
		mutate(Conc=C)

	fish_ppp <- ppp(fish_data_list[[i]]$x, fish_data_list[[i]]$y, window = study_area,marks = exp(fish_data_list[[i]]$z))
	k_ripley_list[[i]] <- Kest(fish_ppp, correction = "Ripley",rmax = 60) %>%
		as.data.frame() %>% mutate(diff=iso-theo) %>%
		mutate(Conc=C)
}

k_ripley_list <- k_ripley_list[-c(10:11)]
k_ripley <- bind_rows(k_ripley_list, .id = "source") %>%
	rename(id='source') %>%
	mutate(id=as.numeric(id))

r_new <- data.frame(r = unique(k_ripley$r))

pred <- k_ripley %>%
	group_by(id) %>%
	group_split() %>%
	map_dfr(~{
		model <- gam(iso ~ s(r, k = 3), data = .x)
		preds <- predict(model, newdata = r_new)
		tibble(id = unique(.x$id), prediction = list(preds))
	}) %>%
	unnest(prediction) %>%
	group_by(id) %>%
	mutate(row = row_number()) %>%
	ungroup()


k_ripley$iso_pred <- pred$prediction

k_ripley <- k_ripley %>% mutate(diff_pred=iso_pred-theo)

r_pred <- k_ripley %>%
	group_by(id) %>%
	mutate(min_range = quantile(diff_pred, 0.85)) %>%
	ungroup() %>%
	filter(diff_pred > min_range) %>%
	group_by(Conc) %>%
	summarise(min_r=min(r),
						mean_r=mean(r)) %>% 
	mutate(C_p=Conc/sum(Conc)) %>%
	mutate(w_min_r=weighted.mean(min_r, C_p))

ripley_i[[j]] <- r_pred %>% distinct(w_min_r) %>% as.data.frame() %>% setNames('R') %>% 
	mutate(alpha=alpha,
				 distance=length_scale,
				 N=n_points)


# Moran's Index -------------------------------------------------------------------------------

window_range_x <- coords %>% distinct(x) %>% arrange(x) %>% pull(x) %>% diff() %>% min()
window_range_y <- coords %>% distinct(y) %>% arrange(y) %>% pull(y) %>% diff() %>% min()
window_range <- (window_range_x+window_range_y)/2
# coords %>% mutate(diff_x=c(0,diff(x)))%>% mutate(diff_y=c(0,diff(y)))
# window_range_x <- coords %>% distinct(x) %>% arrange(x) %>% pull(x) %>% diff() %>% min()
# window_range_y <- coords %>% distinct(y) %>% arrange(y) %>% pull(y) %>% diff() %>% min()
# window_range <- (window_range_x+window_range_y)/2

l <- seq(10, 70, by = floor(window_range))
m_i <- vector(length = length(l),mode = 'numeric')
for (i in 2:(length(l))) {
	nb <- dnearneigh(x=coords, l[i-1],l[i])
	nb_list <- nb2listw(nb, style = "W",zero.policy=F)
	m_test <- moran.test(z, nb_list,alternative = 'greater')
	m_i[i] <- m_test$estimate[1]%>% as.vector()	
}

moran_i[[j]] <- m_i %>% diff(.) %>% 
	as.data.frame() %>% 
	cbind(l[-1]) %>% 
	set_names(c('M','r')) %>%
	slice(-1) %>% 
	mutate(min_M=min(M)) %>% filter(M==min_M) %>% 
	dplyr::select(r) %>% as.data.frame() %>% setNames('M') %>% 
	mutate(alpha=alpha,
				 distance=length_scale,
				 N=n_points)

cat(paste0(100*j/nrow(cond),'%'));cat('\n')
}

output <- bind_rows(ripley_i,.id='id') %>% 
	left_join(bind_rows(moran_i,.id='id'),by=c('id','alpha','distance','N'))

output %>% 
	filter(R<100) %>%
	group_by(alpha,distance,N) %>% 
	mutate(mean_R=mean(R),
						mean_M=mean(M)) %>% 
	ggplot()+
	# geom_point(aes(x=distance,y=M),color='deepskyblue2',size=2,alpha=0.5)+
	# geom_jitter(aes(x=distance,y=R),color='tomato2',size=2,alpha=0.5,width = 0.1)+
	# geom_smooth(aes(x=distance,y=mean_R),color='tomato2',size=1,se=F)+
	# geom_line(aes(x=distance,y=mean_R),color='tomato2',size=1)+
	geom_jitter(aes(x=distance,y=M),color='deepskyblue',size=2,pch=4,alpha=0.5,width = 0.1)+
	geom_point(aes(x=distance,y=mean_M),color='deepskyblue',size=2)+
	geom_line(aes(x=distance,y=mean_M),color='deepskyblue',size=1)+
	# geom_smooth(aes(x=distance,y=M),color='deepskyblue',size=0.5,se=F,span=1)+
	geom_abline(intercept = 0,slope=1,lty=2)+
	# ylim(15,60)+
	facet_wrap(~N)+
	labs(x='Distance simulated',y='Distance predicted')+
	theme_bw()

output %>% 
	ggplot()+
	# geom_point(aes(x=distance,y=M),color='deepskyblue2',size=2,alpha=0.5)+
	geom_jitter(aes(x=distance,y=R),color='tomato2',size=2,alpha=0.5,width = 0.1)+
	geom_smooth(aes(x=distance,y=R),color='tomato2',size=1,se=F,span=1)+
	# geom_abline(intercept = 0,slope=1,lty=2)+
	ylim(15,60)+
	labs(x='Distance simulated',y='Distance predicted')+
	theme_bw()

# 
# # Plot GP
# p1 <- data.frame(x = coords$x, y = coords$y, abundance = z) %>%
# 	# slice(-slice) %>%
# 	ggplot(aes(x = x, y = y, fill = abundance)) +
# 	geom_raster() +  # or use geom_tile()
# 	scale_fill_viridis_c() +
# 	coord_equal() +
# 	theme_minimal() +
# 	labs(title = paste0('Gaussian Process, μ=3; ',"α=",alpha,'; d=',length_scale),
# 			 fill = "Abundance")
# 
# # # Add hotspots manually (optional but useful to force 2-3 abundance peaks)
# # hotspot_centers <- matrix(c(25, 25, 75, 75, 50, 50), ncol=2, byrow=TRUE)
# # hotspot_sd <- 4
# # hotspot_intensity <- rep(1, 3)  # adjust strength
# # 
# # for (i in 1:nrow(hotspot_centers)) {
# # 	d2 <- rowSums((coords - hotspot_centers[i, ])^2)
# # 	z <- z + hotspot_intensity[i] * exp(-d2 / (2 * hotspot_sd^2))
# # }
# # 
# # # Make sure values are positive (as abundance can't be negative)
# # abundance <- exp(z - mean(z))  # Exponentiate to make abundances
# 
# # # Plot result
# # library(ggplot2)
# # data.frame(x = x, y = y, abundance = z) %>%
# # 	ggplot(aes(x = x, y = y, color = z)) +
# # 	geom_point(size = 4) +
# # 	scale_color_viridis_c() +
# # 	coord_equal() +
# # 	theme_minimal() +
# # 	labs(title = "Simulated Fish Abundance (Gaussian Process with Hotspots)",
# # 			 color = "Abundance")
# 
# 
# 
# # Ripley's K function -------------------------------------------------------------------------
# fish_data <- data.frame(x = x, y = y, abundance = z)
# 
# # fish_data %>% 
# # ggplot(aes(x = x, y = y, fill = abundance)) +
# # 	geom_raster() +
# # 	scale_fill_viridis_c() +
# # 	coord_equal() +
# # 	theme_minimal() +
# # 	labs(title = "Simulated Fish Abundance (Gaussian Process)",
# # 			 color = "Abundance")
# 
# k <- seq(min(fish_data$abundance),max(fish_data$abundance),length.out=11)
# 
# # plot_list <- list()
# fish_data_list <- list()
# k_ripley_list <- list() 
# 
# study_area <- owin(c(min(fish_data$x), max(fish_data$x)),
# 									 c(min(fish_data$y), max(fish_data$y)))
# 
# for (i in 1:length(k)) {
# 	C <- k[i]
# 	fish_data_list[[i]] <- fish_data %>% 
# 		filter(abundance>C) %>% 
# 		mutate(Conc=C)
# 	
# 	fish_ppp <- ppp(fish_data_list[[i]]$x, fish_data_list[[i]]$y, window = study_area,marks = exp(fish_data_list[[i]]$abundance))
# 	k_ripley_list[[i]] <- Kest(fish_ppp, correction = "Ripley",rmax = 60) %>%
# 		as.data.frame() %>% mutate(diff=iso-theo) %>% 
# 		mutate(Conc=C)
# }
# 
# k_ripley_list <- k_ripley_list[-c(10:11)]
# df <- bind_rows(k_ripley_list, .id = "source") %>% 
# 	rename(id='source') %>%
# 	mutate(id=as.numeric(id))
# 
# # df %>% filter(diff>1) %>% 
# # 	mutate(diff=log(diff)) %>% 
# # 	pull(diff) %>% hist(,breaks = 300)
# 
# # # p2 <- 
# # df %>% 
# # 	# filter(id==8) %>%
# # 	ggplot()+
# # 	# geom_line(aes(x=r,y=iso,colour=Conc))+
# # 	geom_line(aes(x=r,y=theo,colour=Conc))+
# # 	# geom_line(aes(x=r,y=theo,colour=Conc))+
# # 	geom_smooth(aes(x=r,y=iso,colour=Conc))+
# # 	ggtitle('Ripley\'s K function')+
# # 	# geom_smooth(aes(x=r,y=diff,colour=C),se=F,span = 0.2)+
# # 	# scale_color_gradient(low = "blue", high = "red") +
# # 	theme_bw()
# 
# # fig <- cowplot::plot_grid(p1,p2,align = 'h')
# # ggsave(here('Plots',paste0('n',n_points,'_a=',alpha,'_d',length_scale,'.jpg')),fig,width=12,height=5)
# 
# # library(dplyr)
# # library(purrr)
# # library(mgcv)  # for gam()
# 
# # Assume r_new is a dataframe containing 'r' values for prediction
# # and you want predictions for each id using their own model
# 
# # r_new <- unique(df$r)
# 
# r_new <- data.frame(r = unique(df$r))
# 
# # Split data by id, fit a model for each, then predict
# pred <- df %>%
# 	group_by(id) %>%
# 	group_split() %>%
# 	map_dfr(~{
# 		model <- gam(iso ~ s(r, k = 3), data = .x)
# 		preds <- predict(model, newdata = r_new)
# 		tibble(id = unique(.x$id), prediction = list(preds))
# 	}) %>%
# 	unnest(prediction) %>%
# 	group_by(id) %>%
# 	mutate(row = row_number()) %>%
# 	ungroup()
# 
# df$iso_pred <- pred$prediction
# 
# df <- df %>% mutate(diff_pred=iso_pred-theo)
# 
# r_pred <- df %>%
# 	group_by(id) %>%
# 	mutate(
# 		min_range = quantile(diff_pred, 0.85),
# 		# max_range = quantile(diff_pred, 0.95)
# 	) %>%
# 	ungroup() %>%
# 	filter(diff_pred > min_range) %>%
# 	group_by(Conc) %>%
# 	summarise(min_r=min(r),
# 						mean_r=mean(r)) %>%
# 	mutate(Conc_exp=exp(Conc)) %>%
# 	mutate(C_p=Conc_exp/sum(Conc_exp)) %>%
# 	mutate(w_mean_r=weighted.mean(mean_r, C_p)) %>%
# 	mutate(w_min_r=weighted.mean(min_r, C_p)) %>%
# 	mutate(sweet_spot=(w_mean_r+w_min_r)/2)
# 
# output_data[[j]] <- r_pred %>%
# 	mutate(alpha=alpha,
# 				 distance=length_scale,
# 				 N=n_points)
# 
# print(j)

# # p2 <-
# df %>%
# 	# filter(id==8) %>%
# 	ggplot()+
# 	# geom_line(aes(x=r,y=iso,colour=Conc))+
# 	geom_line(aes(x=r,y=diff_pred,colour=Conc))+
# 	# geom_line(aes(x=r,y=theo,colour=Conc))+
# 	# geom_smooth(aes(x=r,y=iso,colour=Conc))+
# 	ggtitle('Ripley\'s K function')+
# 	# geom_smooth(aes(x=r,y=diff,colour=C),se=F,span = 0.2)+
# 	# scale_color_gradient(low = "blue", high = "red") +
# 	theme_bw()

#
# predict(gam(iso ~ s(r,k = 3), data = dd), newdata = r_new)
#
# dd <- df %>% filter(id==8)
# model <- gam(iso ~ s(r,k = 3), data = dd)
#
# # Predict on new or existing values
# r_new <- data.frame(r = seq(min(dd$r), max(dd$r), length.out = 100))
# r_new <- data.frame(r = dd$r)
#
# iso_pred <- predict(model, newdata = r_new)
# iso_pred_df <- iso_pred %>% as.data.frame() %>% setNames('iso') %>%
# 	cbind(r_new %>% as.data.frame())
}

oo <- bind_rows(output_data)

oo_sum <- oo %>%
	mutate(N=as.factor(N)) %>%
	mutate(alpha=as.factor(alpha)) %>%
	group_by(N,alpha,distance) %>%
	mutate(R=mean(mean_r))

# f <-
# oo_sum %>%
# 	mutate(N=as.factor(N)) %>%
# 	mutate(alpha=as.factor(alpha)) %>%
# 	ggplot()+
# 	geom_point(aes(x=distance,y=min_r,color=alpha),size=1,alpha=0.4)+
# 	# geom_point(aes(x=distance,y=mean_r,color=alpha),size=4)+
# 	geom_line(aes(x=distance,y=mean_r,color=alpha))+
# 	# geom_smooth(aes(x=distance,y=r,color=alpha),se=F)+
# 	geom_abline(intercept = 0,slope=1,lty=2)+
# 	facet_wrap(~N)+
# 	labs(x='Distance simulated',y='Distance predicted')+
# 	theme_bw()
#
ggsave(here('Plots','Sim_pred_GP_2.jpg'),f,width = 10,height=5)
# library(ggridges)
# library(ggforce)
#
# oo_sum %>%
# 	mutate(distance=as.factor(distance)) %>%
# ggplot(aes(x = min_r, y = `alpha`, fill = `N`)) +
# 	geom_density_ridges(scale = 0.6, alpha = 0.5) +
# 	facet_wrap(~ distance)+
# 	geom_abline(slope = c(0, 0, 0), intercept = c(10, 20, 40),
# 							color = "red", linetype = "dashed")+
# 	theme_ridges()


f <- oo_sum %>%
	mutate(N=as.factor(N)) %>%
	mutate(alpha=as.factor(alpha)) %>%
	ggplot()+
	geom_point(aes(x=distance,y=mean_r,color=alpha),size=1,alpha=0.4)+
	geom_point(aes(x=distance,y=R,color=alpha),size=4)+
	geom_line(aes(x=distance,y=R,color=alpha))+
	# geom_smooth(aes(x=distance,y=r,color=alpha),se=F)+
	geom_abline(intercept = 0,slope=1,lty=2)+
	facet_wrap(~N)+
	labs(x='Distance simulated',y='Distance predicted')+
	theme_bw()
# 
# # Moran's I -----------------------------------------------------------------------------------
# 
# # Example: define neighbors within 0.1 units
# nb <- dnearneigh(coords, 0, 32)
# 
# # Check for isolated points (no neighbors)
# # summary(nb)
# 
# lw <- nb2listw(nb, style = "W")  # row-standardized weights
# moran.test(fish_data$abundance, lw)

out <- list()
for (k in 1:5) {
z <- mvrnorm(mu = rep(3, n_points), Sigma = cov_mat)

l <- seq(5, 70, by = 3)
m_i <- vector(length = length(l),mode = 'numeric')
for (i in 2:(length(l))) {
	nb <- dnearneigh(x=coords, l[i-1],l[i])
	nb_list <- nb2listw(nb, style = "W",zero.policy=F)
	m_test <- moran.test(z, nb_list,alternative = 'greater')
	m_i[i] <- m_test$estimate[1]%>% as.vector()	
}

out[[k]] <- m_i %>% diff(.) %>% 
	as.data.frame() %>% 
	cbind(l[-1]) %>% 
	set_names(c('M','r')) %>%
	slice(-1) 
print(k)
}

out_c <- bind_rows(out,.id = "id")

out_c %>% 
ggplot()+
	labs(y='Derivative of Moran\'s index',x = 'Distance window',
			 # title=paste0('d=',length_scale,'_α=',alpha,'_N=',n_points))+
			 title=paste0('d=',length_scale))+
	# geom_smooth(aes(x=r,y=M,color=id),se=F)+
	geom_point(aes(x=r,y=M),size = 2,alpha=0.6)+
	geom_smooth(aes(x=r,y=M),se=F)+
	theme_bw()

# l <- c(15,20,25,30,35,40,45,50,60,70)
l <- seq(5, 70, by = 5)
m_i <- vector(length = length(l),mode = 'numeric')
# m_i_var <- vector(length = length(l),mode = 'numeric')
for (i in 2:(length(l))) {
nb <- dnearneigh(x=coords, l[i-1],l[i])
# for (i in 1:(length(l))) {
# 	nb <- dnearneigh(coords, 0, l[i])
	nb_list <- nb2listw(nb, style = "W",zero.policy=F)
	# m_test <- moran.test(z, nb_list,alternative = 'two.sided')
	m_test <- moran.test(z, nb_list,alternative = 'greater')
	m_i[i] <- m_test$estimate[1]%>% as.vector()	
	m_i_var[i] <- m_test$estimate[3]%>% as.vector()	
}

# f2 <- 
	m_i %>% diff(.) %>% 
	as.data.frame() %>% 
	cbind(l[-1]) %>% 
	# cbind(m_i_var) %>% 
	set_names(c('M','r')) %>%
	# set_names(c('M','r','m_i_var')) %>%
	# mutate(alpha=alpha,
	# 			 distance=length_scale,
	# 			 N=n_points) %>% 
		slice(-1) %>% 
	ggplot()+
	labs(y='Moran\'s index',x = 'Distance window',
			 title=paste0('d=',length_scale,'_α=',alpha,'_N=',n_points))+
	geom_smooth(aes(x=r,y=M),se=F)+
	geom_point(aes(x=r,y=M),size = 3)+
	theme_bw()

ggsave(here('Plots','Morans_I.jpg'),f2,width = 8,height = 6)

df <-
m_i %>% 
	as.data.frame() %>% 
	cbind(l) %>% 
	cbind(m_i_var) %>% 
	set_names(c('M','r','m_i_var')) %>%
	mutate(alpha=alpha,
				 distance=length_scale,
				 N=n_points) %>% slice(-1) %>% 
	mutate(M=M+(-min(M))) %>% 
	mutate(log_M=log(M+0.00001))

# # Sample data: x and y from Gaussian kernel
# rho_raw <- 20
# a <- seq(10, 70, by = 0.1)
# b <- exp(-a^2 / (2 * rho_raw^2))  # simulate with known rho = 40
# df <- data.frame(a = a, b = b)

# df$logb <- log(df$b)
lin_fit <- lm(log_M ~ I(r^2), data = df)
slope <- coef(lin_fit)[2]
rho_est <- sqrt(-1 / (2 * slope))
rho_est

model <- nls(b ~ exp(-a^2 / (2 * rho^2)),
						 data = df,
						 start = list(rho = 10))


model

output_data[[j]] <- m_i %>% 
	as.data.frame() %>% 
	cbind(l) %>% 
	cbind(m_i_ex) %>% 
	set_names(c('M','r','m_i_ex')) %>%
	mutate(alpha=alpha,
				 distance=length_scale,
				 N=n_points)

print(j)

}
# moran.mc(fish_data$abundance, lw, nsim = 999)

oo <- bind_rows(output_data)

oo_sum <- oo %>%
	mutate(N=as.factor(N)) %>%
	mutate(alpha=as.factor(alpha)) %>%
	group_by(N,alpha,distance) %>%
	mutate(max_M=max(M)) %>% 
	ungroup() %>% 
	filter(M==max_M)

# f <-
oo_sum %>%
	mutate(N=as.factor(N)) %>%
	mutate(alpha=as.factor(alpha)) %>%
	ggplot()+
	geom_point(aes(x=distance,y=r,color=alpha),size=1,alpha=0.4,size=4)+
	# geom_point(aes(x=distance,y=mean_r,color=alpha),size=4)+
	geom_line(aes(x=distance,y=r,color=alpha))+
	# geom_smooth(aes(x=distance,y=r,color=alpha),se=F)+
	# geom_abline(intercept = 0,slope=1,lty=2)+
	facet_wrap(~N)+
	labs(x='Distance simulated',y='Distance predicted')+
	theme_bw()
#
ggsave(here('Plots','Sim_pred_GP_MoransI.jpg'),f,width = 10,height=5)
# library(ggridges)
# library(ggforce)
#
# oo_sum %>%
# 	mutate(distance=as.factor(distance)) %>%
# ggplot(aes(x = min_r, y = `alpha`, fill = `N`)) +
# 	geom_density_ridges(scale = 0.6, alpha = 0.5) +
# 	facet_wrap(~ distance)+
# 	geom_abline(slope = c(0, 0, 0), intercept = c(10, 20, 40),
# 							color = "red", linetype = "dashed")+
# 	theme_ridges()

# Plot GP
data.frame(x = coords$x, y = coords$y, abundance = z) %>%
	# slice(-slice) %>%
	ggplot(aes(x = x, y = y, fill = abundance)) +
	geom_raster() +  # or use geom_tile()
	scale_fill_viridis_c() +
	coord_equal() +
	theme_minimal() +
	labs(title = paste0('Gaussian Process, μ=3; ',"α=",alpha,'; d=',length_scale),
			 fill = "Abundance")
