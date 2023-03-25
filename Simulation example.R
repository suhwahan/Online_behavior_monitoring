### Example simulation

source("SimulationDriver_Mar2023_Git.R")
library(foreach)
library(doParallel)

##--| For parallel computing -----------------------
n.cores <- detectCores() - 2
cl <- makeCluster(n.cores)
registerDoParallel(cl)

##--| Simulation -----------------------
sim <- foreach(r=c(1:2), .packages = c("LNIRT", "purrr", "plyr"), .errorhandling = "pass") %dopar% {
  
  Sim_driver (nitem=40, nexaminee=1000, rho_item=.25, 
              rho_person=.3, c_f=0, p_ab_p=.1, c.pt=21,
              mu_a=.91, var_a=.79^2, mu_alp=1.49, var_alp=.39^2,
              ref_v_W_j = c(.25, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4), 
              lambda_Wj = seq(0.005, 0.995, by = 0.005),
              ref_v_V_j = c(0.5, 1, 2.5, 5, 7.5, 10, 12.5, 15), 
              lambda_Vj = seq(0.005, 0.995, by = 0.005),
              lambda_G = seq(0.005, 0.995, by = 0.005),
              lambda_G_excl = seq(0.005, 0.995, by = 0.005) )
  
}


