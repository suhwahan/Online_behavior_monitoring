### ------------------------------------------------------------ ###
### Function to calculate monitoring statistics for CUSUM chart  ###
### monitoring indicator variables based on a fixed training set ###
### ------------------------------------------------------------ ###

calc_W_j_ref_fixed <- function(resp, rt, n_rf_s, ir_est_mat, rt_est_mat,
                               mu_p_t, cov_p_t, sampling_th_1){ 
  
  ## Input values for the function -------------------
  # resp: a vector of item responses
  # rt: a vector of response times
  # n_rf_s: Number of items in the training set
  # ir_est_mat: item parameters across the items
  # rt_est_mat: response time parameters across the items
  # mu_p_t: prior for the person parameters (mean ability and speed)
  # cov_p_t: prior for the covariance matrix of the person parameters
  # sampling_th_1: A character variable that determines the sampling 
  #                strategy for the evaluation set
  # Note. The function was written to take the entire vector of item
  #       responses and response times for the computational efficiency
  #       Below the observations are evaluated sequentially 
  #       within this function
  
  ## Set up -------------------------------
  Nitem <- dim(resp)[2] # number of items
  Nexaminee <- dim(resp)[1] # number of examinees
  st_eval <- n_rf_s+1 # start of the evaluation point
  neval <- Nitem-st_eval+1 # number of evaluations
  
  ## Expected theta and tau ---------------
  # Estimate the person parameters using the fixed training set
  # The map function below is from Kang (2016)
  th_tau_null <- map(resp=resp[,1:n_rf_s], RT=rt[,1:n_rf_s], 
                     ipar = cbind(ir_est_mat[1:n_rf_s,], 
                                  rt_est_mat[1:n_rf_s,]),
                     ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t),
                     SE=FALSE, D=1.702) 
  
  th_0 <- th_tau_null$est[,"th"] # ability estimate based on the training set
  tau_0 <- th_tau_null$est[,"tau"] # speed estimate based on the training set
  
  ## Obtain sequential person parameter estimates ----------
  # This is to calculate p in Equation 4
  th_mat <- matrix(NA, nrow = Nexaminee, ncol = neval) 
  for (j in 1:neval){ 
    # If calculating p using all observations until j
    if (sampling_th_1 == "union_S_R"){ 
      tmp <- map(resp=resp[,1:(n_rf_s+j)], RT=rt[,1:(n_rf_s+j)], 
                 ipar=cbind(ir_est_mat[1:(n_rf_s+j),], 
                            par_rt_t=rt_est_mat[1:(n_rf_s+j),]),
                 ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t))
      
      # if calculating using observations after the training until j 
    } else if (sampling_th_1 == "S_only"){ 
      if (j < n_rf_s) {
        tmp <- map(resp=resp[,(j+1):(n_rf_s+j)], RT=rt[,(j+1):(n_rf_s+j)], 
                   ipar=cbind(ir_est_mat[(j+1):(n_rf_s+j),], 
                              par_rt_t=rt_est_mat[(j+1):(n_rf_s+j),]),
                   ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t))
      } else{
        tmp <- map(resp=resp[,(n_rf_s+1):(n_rf_s+j)], 
                   RT=rt[,(n_rf_s+1):(n_rf_s+j)], 
                   ipar=cbind(ir_est_mat[(n_rf_s+1):(n_rf_s+j),], 
                              par_rt_t=rt_est_mat[(n_rf_s+1):(n_rf_s+j),]),
                   ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t))
      }
    }
    th_mat[,j] <- tmp$est[,"th"]
  }
  
  ## Calculate expected item response values--------
  exp_X_null_mat <- matrix(NA, nrow=Nexaminee, ncol=neval, byrow=FALSE) 
  for (j in 1:neval){ 
    exp_X_null_mat[,j] <- irf_3pl(th=th_0, ipars=ir_est_mat[n_rf_s+j,], D=1.702)
  }
  
  ## Calculate expected log RT -----------
  # time intensity parameters
  bet_mat <- matrix(rt_est_mat[st_eval:Nitem,"beta"], 
                    nrow=Nexaminee, ncol=neval, byrow=TRUE) 
  # person speed estimate from the training set
  tau_mat <- matrix(tau_0, nrow=Nexaminee, ncol=neval, byrow=FALSE) 
  exp_logRT_null_mat <- bet_mat-tau_mat  
  
  ## Calculate the multivariate monitoring statistic ------------
  Wj_mat <- matrix(NA, nrow=Nexaminee, ncol=neval)   
  for (p in 1:Nexaminee){ 
    cov_pi <- matrix(0, 2, 2) # Equation 4
    for (j in 1:neval){ 
      pr <- irf_3pl(th=th_mat[p,j], ipars=ir_est_mat[n_rf_s+j,], D=1.702)
      cov_pi[1,1] <- pr*(1-pr)
      cov_pi[1,2] <- cov_pi[2,1] <- 0
      cov_pi[2,2] <- 1/(rt_est_mat[n_rf_s+j,"alpha"]^2)
      U_i <- c(resp[p,n_rf_s+j], log(rt[p,n_rf_s+j]))
      names(U_i) <- c("X_i", "logT_i") 
      U_0 <- c(exp_X_null_mat[p,j], exp_logRT_null_mat[p,j])
      names(U_0) <- c("X_0", "logT_0")
      Wj <- t(as.matrix(U_i-U_0)) %*% MASS::ginv(cov_pi) %*% as.matrix(U_i-U_0)
      Wj_mat[p,j] <- Wj # Equation 3
    }
  }
  return(list(Wj_mat=Wj_mat, th_0=th_0, tau_0=tau_0, th_mat=th_mat))
}

### ------------------------------------------------------------ ###
### Function to calculate monitoring statistics                  ###
### for the trait-based CUSUM control chart                      ###
### based on a fixed training set and batches of items           ###
### ------------------------------------------------------------ ###

calc_V_j_CUSUM_ref_fixed <- function(resp, rt, n_rf_s, ir_est_mat, rt_est_mat, 
                                     mu_p_t, cov_p_t){
  
  ## Input values for the function -----------------
  # resp: a vector of item responses
  # rt: a vector of response times
  # n_rf_s: Number of items in the training set
  # ir_est_mat: item parameters across the items
  # rt_est_mat: response time parameters across the items
  # mu_p_t: prior for the person parameters (mean ability and speed)
  # cov_p_t: prior for the covariance matrix of the person parameters
  # Note. The function was written to take the entire vector of item
  #       responses and response times for the computational efficiency
  #       Below the observations are evaluated sequentially 
  #       within this function
  
  ## Set up --------------------
  Nitem <- dim(resp)[2] # number of items
  Nexaminee <- dim(resp)[1] # number of examinees
  
  ## Expected theta and tau --------------
  # Estimate the person parameters using the fixed training set
  # The map function below is from Kang (2016)
  exp_th_tau <- map(resp=resp[,1:n_rf_s], RT=rt[,1:n_rf_s], 
                    ipar = cbind(ir_est_mat[1:n_rf_s,], 
                                 rt_est_mat[1:n_rf_s,]),
                    ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t),
                    SE=TRUE, D=1.702)
  
  th_0 <- exp_th_tau$est[,"th"]  # ability estimate based on the training set
  tau_0 <- exp_th_tau$est[,"tau"]   # speed estimate based on the training set
  eta_0 <- cbind(th_0, tau_0) # a vector of latent trait estimates from the training set
  S_0 <- array(NA, c(2,2,length(th_0))) # a var-cov matrix from the null estimates
  for (i in 1:nrow(eta_0)){ 
    S_0[1,2,i] <- S_0[2,1,i] <- 0
    S_0[1,1,i] <- (exp_th_tau$se^2)[i,1] # standard error from the map function
    S_0[2,2,i] <- (exp_th_tau$se^2)[i,2]
  }
  
  ## Specify indices of the item sets ---------------------
  idx_sample <- seq(2, (Nitem/n_rf_s))  # indices of evaluation item sets
  # The size of moving samples is equal to `n_rf_s` throughout
  # length(idx_sample): # of evaluation points
  st_eval_pts <- idx_sample*n_rf_s-(n_rf_s-1) # First item in each evaluation set
  end_eval_pts <- idx_sample*n_rf_s # Last item in each evaluation set
  
  ## Latent trait estimates from each evaluation set ------------- 
  th_1_mat <- matrix(NA, nrow=Nexaminee, ncol=length(idx_sample)) # ability estimate at each j
  tau_1_mat <- matrix(NA, nrow=Nexaminee, ncol=length(idx_sample)) # speed estimate at each j
  V_j_mat <- matrix(NA, nrow=Nexaminee, ncol=length(idx_sample))
  
  for (j in 1:length(idx_sample)){ 
    th_tau_1 <- map(resp=resp[ ,st_eval_pts[j]:end_eval_pts[j]], 
                    RT=rt[ ,st_eval_pts[j]:end_eval_pts[j]], 
                    ipar = cbind(ir_est_mat[st_eval_pts[j]:end_eval_pts[j], ], 
                                 rt_est_mat[st_eval_pts[j]:end_eval_pts[j], ]),
                    ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t),
                    SE=TRUE, D=1.702)
    
    th_1_mat[,j] <- th_tau_1$est[ ,"th"]
    tau_1_mat[,j] <- th_tau_1$est[ ,"tau"]
    
    eta_1 <- cbind(th_tau_1$est[ ,"th"],  
                   th_tau_1$est[ ,"tau"]) # a vector of latent trait estimates from evaluation set j
    colnames(eta_1) <- c("th_1", "tau_1")
    S_1 <- array(NA, c(2,2,length(th_0))) # a var-cov matrix from the trait estimates at j
    
    for (i in 1:nrow(eta_1)){ 
      S_1[1,2,i] <- S_1[2,1,i] <- 0
      S_1[1,1,i] <- (th_tau_1$se^2)[i,1] # standard error from the map function
      S_1[2,2,i] <- (th_tau_1$se^2)[i,2]
    }
    
    ## Calculate the multivariate monitoring statistic --------------
    Htl_T2 <- rep(NA, Nexaminee) # Hotelling's T^2
    for (p in 1:Nexaminee){ 
        pooled_S <- ((n_rf_s-1)*(S_0[,,p])+(n_rf_s-1)*(S_1[,,p])) / (2*n_rf_s-2)
        Htl_T2[p] <- ((n_rf_s^2)/(2*n_rf_s))*(eta_1[p,]-eta_0[p,]) %*% solve(pooled_S) %*% as.matrix(eta_1[p,]-eta_0[p,])
      }
    V_j_mat[,j] <- Htl_T2
  } # end of j
  return(list(V_j_mat=V_j_mat, th_0=th_0, tau_0=tau_0, 
              th_1_mat=th_1_mat, tau_1_mat=tau_1_mat))
}  

### ------------------------------------------------------------ ###
### Function to implement the SGLRT procedure                    ###
### based on a fixed training set                                ###
### ------------------------------------------------------------ ###

calc_SGLRT <- function(resp, rt, n_rf_s, ir_est_mat, rt_est_mat, mu_p_t, cov_p_t,
                         sampling){
  
  ## Input values for the function -----------------
  # resp: a vector of item responses
  # rt: a vector of response times
  # n_rf_s: Number of items in the training set
  # ir_est_mat: item parameters across the items
  # rt_est_mat: response time parameters across the items
  # mu_p_t: prior for the person parameters (mean ability and speed)
  # cov_p_t: prior for the covariance matrix of the person parameters
  # sampling: a character variable specifying sampling strategies
  #     "ref_mov_fixed": Sliding window sampling with a fixed training set
  #     "exclusive": Sampling with exclusive item batches with a fixed training set
  # Note. The function was written to take the entire vector of item
  #       responses and response times for the computational efficiency
  #       Below the observations are evaluated sequentially 
  #       within this function
  
  ## Set up ---------------------------
  Nitem <- ncol(resp)
  Nexaminee <- nrow(resp)
  
  if (sampling %in% c("ref_mov_fixed")){
    st_eval <- n_rf_s+1; neval <- Nitem-st_eval+1
  } else if (sampling == "exclusive"){
    st_eval <- 2*n_rf_s; neval <- Nitem/n_rf_s-1
  }
  
  th_1_mat <- matrix(NA, nrow=Nexaminee, ncol=neval) 
  tau_1_mat <- matrix(NA, nrow=Nexaminee, ncol=neval)
  
  ## Obtain the trait estimates for different sampling strategies
  # If sliding window sampling with a fixed training set
  if (sampling == "ref_mov_fixed"){ 
    
    ## Expected theta and tau from a fixed training set -----
    # The map function below is from Kang (2016)
    th_tau_0 <- map(resp=resp[,1:n_rf_s], RT=rt[,1:n_rf_s], 
                    ipar = cbind(ir_est_mat[1:n_rf_s,], rt_est_mat[1:n_rf_s,]),
                    ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t),
                    SE=FALSE, D=1.702)
 
    th_0_mat <- matrix(th_tau_0$est[,"th"], nrow=Nexaminee, ncol=neval) 
    tau_0_mat <- matrix(th_tau_0$est[,"tau"], nrow=Nexaminee, ncol=neval) 
    
    ## Sequential theta and tau at j ------------------
    for (j in 1:neval){ 
      th_tau_1 <- map(resp=resp[,(j+1):(j+n_rf_s)], RT=rt[,(j+1):(j+n_rf_s)], 
                      ipar = cbind(ir_est_mat[(j+1):(j+n_rf_s),], 
                                   rt_est_mat[(j+1):(j+n_rf_s),]),
                      ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t),
                      SE=FALSE, D=1.702)
      
      th_1_mat[,j] <- th_tau_1$est[,"th"]
      tau_1_mat[,j] <- th_tau_1$est[,"tau"]
    }
    
  # If sampling with item batches with a fixed training set  
  } else if (sampling == "exclusive"){
  
    ## Expected theta and tau from a fixed training set -----
    th_tau_0 <- map(resp=resp[,1:n_rf_s], RT=rt[,1:n_rf_s],
                    ipar=cbind(ir_est_mat[1:n_rf_s,], rt_est_mat[1:n_rf_s,]),
                    ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t),
                    SE=FALSE, D=1.702) 
    th_0_mat <- matrix(th_tau_0$est[,"th"], nrow=Nexaminee, ncol=neval) 
    tau_0_mat <- matrix(th_tau_0$est[,"tau"], nrow=Nexaminee, ncol=neval) 
    
    ## Sequential theta and tau at j ------------------
    for (j in 1:neval){ 
      th_tau_1 <- map(resp=resp[,(n_rf_s*j+1):(n_rf_s*j+n_rf_s)], 
                      RT=rt[,(n_rf_s*j+1):(n_rf_s*j+n_rf_s)], 
                      ipar=cbind(ir_est_mat[(n_rf_s*j+1):(n_rf_s*j+n_rf_s),], 
                                 rt_est_mat[(n_rf_s*j+1):(n_rf_s*j+n_rf_s),]),
                      ppar_prior=list(mu_p=mu_p_t, cov_p=cov_p_t), 
                      SE=FALSE, D=1.702)
      
      th_1_mat[,j] <- th_tau_1$est[,"th"]
      tau_1_mat[,j] <- th_tau_1$est[,"tau"]
    }
  }
    ## Calculate the multivariate monitoring statistic ----------------
    LR_mat <- matrix(NA, Nexaminee, neval)
    for (p in 1:Nexaminee){ 
      LR <- vector(mode="numeric", length=neval)
      for (j in 1:neval){ 
        if (sampling == "ref_mov_fixed"){
          RG <- c(j+1, j+n_rf_s) # range of item indices within the item set at j
        } else if (sampling == "exclusive"){
          RG <- c((n_rf_s*j+1),(n_rf_s*j+n_rf_s)) # range of item indices within the item set at j
        }
        f_rt_num <- rt_dfn(rt=rt[p,RG[1]:RG[2]], # Calculate the lognormal RT density given tau at j
                           Alp=rt_est_mat[RG[1]:RG[2], "alpha"], 
                           Bet=rt_est_mat[RG[1]:RG[2], "beta"], 
                           Tau=tau_1_mat[p,j])
        x <- resp[p, RG[1]:RG[2]]
        pr <- irf_3pl_v(xi_ir=ir_est_mat[RG[1]:RG[2],], th=th_1_mat[p,j]) # IRF at the interim theta at j
        f_x_num <- pr^x*(1-pr)^(1-x)
        f_num <- prod(f_rt_num)*prod(f_x_num)
        f_rt_den <- rt_dfn(rt=rt[p,RG[1]:RG[2]], # Calculate the lognormal RT density given the null tau
                           Alp=rt_est_mat[RG[1]:RG[2], "alpha"], 
                           Bet=rt_est_mat[RG[1]:RG[2], "beta"], 
                           Tau=tau_0_mat[p,j])
        pr <- irf_3pl_v(xi_ir=ir_est_mat[RG[1]:RG[2],], th=th_0_mat[p,j]) # IRF at the null theta
        f_x_den <- pr^x*(1-pr)^(1-x)
        f_den <- prod(f_rt_den)*prod(f_x_den)
        LR[j]<- log(f_num/f_den)
      }
      LR_mat[p,] <- LR
    } 
  return(list(LR_mat=LR_mat))
}

### ------------------------------------------------------------ ###
### Calculate CUSUM charting statistics                          ###
### ------------------------------------------------------------ ###

cusumUp <- function(stat, ref_v = 0){
  
  ## Input values for the function -----------------
  # stat: monitoring statistics
  # ref_v: Effect size for the CUSUM chart (a.k.a. reference value)
  
  nexaminee <- length(stat[,1])
  nitem <- length(stat[1,])
  cusum_up <- matrix(NA, nexaminee, nitem)
  for (j in 1:nexaminee){
    cusum <- numeric(1)
    for (i in 1:nitem ){
      cusum <- cusum + stat[j,i] - ref_v
      if (cusum > 0){
        cusum <- cusum
        cusum_up[j,i] <- cusum 
      }
      else{
        cusum <- 0
        cusum_up[j,i] <- 0
      }
    }
  }
  return(cusum_up)
}


