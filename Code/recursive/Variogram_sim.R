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
# install.packages('gstat')
library(gstat)
library(sp)

cond_raw <- expand.grid(a = c(0.1,1), d = c(10,15,20,40), n = c(100,900))

# cond_raw <- data.frame(a=1,d=c(5,15,30,40),n=1600)
cond <- bind_rows(replicate(10, cond_raw, simplify = FALSE))

dim(cond)

variogram_list <- list()
z_list <- list()
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
	
	# Compute distance matrix
	dist_mat <- rdist(coords)
	cov_mat <- alpha * exp(- (dist_mat^2) / (2 * length_scale^2))

# Simulate GP with the specified mean and covariance
z <- mvrnorm(mu = rep(0,n_points), Sigma = cov_mat)
z_list[[j]] <- as.data.frame(z)
fish_density <- cbind(coords,z)

coordinates(fish_density) <- ~x+y
variogram_obj <- variogram(z ~ 1, width=1,fish_density)
variogram_list[[j]] <- variogram_obj %>% 
	mutate(alpha=alpha,
				 distance=length_scale,
				 N=n_points)

cat(j);cat('\n')
}



variogram_comb <- bind_rows(variogram_list,.id='id')
variogram_comb <- variogram_comb %>% mutate(diff=c(0,diff(variogram_comb$gamma)))

bind_rows(variogram_list,.id = 'id') %>% 
	group_by(id) %>% 
	mutate(gamma_adj=gamma/max(gamma)) %>%
	ungroup() %>% 
	group_by(N,distance,dist,alpha) %>% 
	mutate(gamma_adj_mean=mean(gamma_adj)) %>% 
	filter(N==900) %>% 
	# filter(distance==15) %>% 
# variogram_obj %>% 
	ggplot()+
	# geom_line(aes(x=dist,y=gamma))+
	geom_point(aes(x=dist,y=gamma_adj_mean,color=as.factor(distance)),alpha=0.5)+
	# geom_point(aes(x=dist,y=gamma_adj,color=as.factor(distance)),alpha=0.5)+
	# geom_line(aes(x=dist,y=gamma_adj,color=as.factor(distance)))+
	geom_smooth(aes(x=dist,y=gamma_adj_mean,color=as.factor(distance)),se=F)+
	# facet_wrap(~N)+
	# geom_hline(yintercept=0.02)+
	# geom_vline(xintercept=length_scale)+
	theme_bw()


bind_rows(z_list,.id='id') %>% group_by(id) %>% summarize(sd2=sqrt(sd(z)^2)) %>% pull(sd2) %>% hist(,breaks = 100)


dat <- readRDS('/Users/gledguri/Library/CloudStorage/OneDrive-UW/UW/QM-qPCR-joint_clone/data/Zenodo/Log_D_est.rds')
glimpse(dat)
dat <- dat %>% mutate(species=ifelse(species=="Zz_Merluccius productus","Merluccius productus",species))


dat_2 <- dat %>% dplyr::select(mean,lat,lon,species,depth_cat)

dat_2 %>% distinct(depth_cat)
dat_2 %>% distinct(species)

fish_density <- dat_2 %>% 
	# filter(species=='Engraulis mordax') %>%
	filter(species=='Merluccius productus') %>%
	# filter(species=='Sebastes entomelas') %>% 
	filter(depth_cat==50) %>%
	# filter(depth_cat==50) %>%
	# filter(depth_cat==150) %>% 
	rename(x='lon',y='lat',z=mean) %>% 
	dplyr::select(x,y,z)

coordinates(fish_density) <- ~x+y
variogram_obj <- variogram(z ~ 1, width=0.005,fish_density)

variogram_obj %>% 
mutate(gamma_adj=gamma/max(gamma)) %>%
	# filter(distance==15) %>% 
	# variogram_obj %>% 
	ggplot()+
	# geom_line(aes(x=dist,y=gamma))+
	geom_point(aes(x=dist,y=gamma_adj),alpha=0.5)+
	# geom_point(aes(x=dist,y=gamma),alpha=0.5)+
	# geom_point(aes(x=dist,y=gamma_adj,color=as.factor(distance)),alpha=0.5)+
	# geom_line(aes(x=dist,y=gamma_adj,color=as.factor(distance)))+
	geom_smooth(aes(x=dist,y=gamma_adj),se=F)+
	# facet_wrap(~N)+
	# geom_hline(yintercept=0.02)+
	# geom_vline(xintercept=length_scale)+
	theme_bw()



l <- seq(0, 6, by = 1.0)
m_i <- vector(length = length(l),mode = 'numeric')

l <- seq(0, 6, by = 0.5)
m_i <- vector(length = length(l),mode = 'numeric')
dnearneigh(x=fish_density %>% dplyr::select(-z), l[i-1],l[i])

for (i in 2:(length(l))) {
	nb <- dnearneigh(x=fish_density %>% dplyr::select(-z), l[i-1],l[i])
	nb_list <- nb2listw(nb, style = "W",zero.policy=F)
	m_test <- moran.test(fish_density$z, nb_list,alternative = 'greater')
	m_i[i] <- m_test$estimate[1]%>% as.vector()	
}

m_i %>% diff(.) %>% 
	as.data.frame() %>% 
	cbind(l[-1]) %>% 
	set_names(c('M','r')) %>%
	slice(-1) 
