srun --cpus-per-task=4 --mem=40G --time=30-00:00:00 --pty bash

.libPaths(c('/mnt/nfs/home/KellyCEG/ggu/Rpack','/mnt/nfs/home/KellyCEG/R/x86_64-pc-linux-gnu-library/4.3',.libPaths()))
library(here)


prop_samp_thinn <- 200
ww <- '1'
sim_param_1 <- data.frame(alpha=4,rho=1,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '2'
sim_param_1 <- data.frame(alpha=4,rho=2,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '3'
sim_param_1 <- data.frame(alpha=4,rho=3,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '4'
sim_param_1 <- data.frame(alpha=4,rho=4,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '5'
sim_param_1 <- data.frame(alpha=4,rho=5,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '6'
sim_param_1 <- data.frame(alpha=4,rho=10,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '7'
sim_param_1 <- data.frame(alpha=4,rho=15,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '8'
sim_param_1 <- data.frame(alpha=4,rho=20,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '9'
sim_param_1 <- data.frame(alpha=4,rho=25,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '10'
sim_param_1 <- data.frame(alpha=4,rho=30,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '11'
sim_param_1 <- data.frame(alpha=4,rho=8,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 200
ww <- '12'
sim_param_1 <- data.frame(alpha=4,rho=12,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))




prop_samp_thinn <- 0.3
ww <- '13'
sim_param_1 <- data.frame(alpha=4,rho=5,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 0.3
ww <- '14'
sim_param_1 <- data.frame(alpha=4,rho=10,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))

prop_samp_thinn <- 0.3
ww <- '15'
sim_param_1 <- data.frame(alpha=4,rho=1,sigma_sim=3)
source(here('Code','GP_simulation_Jul.R'))