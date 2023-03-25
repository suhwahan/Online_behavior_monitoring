##--| Author: Suhwa Han (suhwa@utexas.edu)

setwd("/Users/suhwahan/Git_project/JEM_Mar2023_Git")

##--| Sources & libraries ---
source('Data generation.R', echo=TRUE)
source('Calculate_sequential_stats.R', echo=TRUE)
source('Evaluation.R', echo=TRUE)
source('mmap_cfix.R', echo=TRUE)
source('mmap_fsb_cfix.R', echo=TRUE)
source("map.R")

library(LNIRT)
library(MASS)
library(plyr)
library(tidyverse)

##--| Simulation driver
Sim_driver <- function(nitem, nexaminee, rho_item, rho_person, c_f, p_ab_p, c.pt,
                       mu_a, var_a, mu_alp, var_alp,
                       ref_v_W_j, lambda_Wj=NULL,
                       ref_v_V_j, lambda_Vj=NULL,
                       lambda_G=NULL,
                       lambda_G_excl=NULL){
  
  
  rn <- sample(seq(1, 10000), 1) # a random number to set seed
  set.seed(rn)
  
  ##--| Generate null data -----------------
  null_Data <- nullDataGen(nitem=nitem, nexaminee=nexaminee, rho_item=rho_item, rho_person=rho_person, c_f=c_f,
                           mu_a=mu_a, var_a=var_a, mu_alp=mu_alp, var_alp=var_alp)
  
  #--| Item parameter calibration (MMLE) ----------------
  start_time <- Sys.time()
  mod_mbme <- mmap_cfix(null_Data$response, null_Data$rt,
                        ppar_prior=list(mu_p=null_Data$mu_p, cov_p=null_Data$cov_p),
                        iparst_prior = list(mu_is=null_Data$mu_item_str, cov_is=null_Data$cov_item_str),
                        cfix=c_f, tol=list(max_em=1000, tol_em=0.005, loglike=1e-4, max_nr=500, tol_nr=0.005))
  end_time <- Sys.time()
  elps_time_mbme <- end_time-start_time

  # Item parameter matrix
  ir_est_mat <- cbind(mod_mbme$iest[,"a"], mod_mbme$iest[,"b"], mod_mbme$iest[,"c"]); colnames(ir_est_mat) <- c("a", "b", "c")
  # Response time parameter matrix
  rt_est_mat <- cbind(mod_mbme$iest[,"alp"], mod_mbme$iest[,"bet"]); colnames(rt_est_mat) <- c("alpha", "beta")
  
  ##--| Bayesian item parameter calibration (LNIRT R package) ----------------
  mod_LNIRT <- LNIRT::LNIRT(RT=log(null_Data$rt), Y=null_Data$response,
                            data=list(as.matrix(null_Data$rt), as.matrix(null_Data$response)),
                            XG=3000, burnin=10,
                            guess=FALSE, par1=TRUE, td=FALSE, WL=TRUE, residual = F)
  

  ##--| Design variables -----------------
  
  # Change size
  delta <- data.frame(c(-1, -1.5, -2), c(1, 1.5, 2)); colnames(delta) <- c("th", "tau")
  # Number of items for training and evaluation set
  n_rf_s <- c(5, 10) 
  # Data frame of the simulation conditions
  sim_fts <- data.frame(merge(n_rf_s, delta)); colnames(sim_fts) <- c("n_rf_s", "delta_th", "delta_tau")

  ##--| An empty vector of lists to store results from the simulation conditions
  Res_w_fts <- vector(mode="list", length=nrow(sim_fts))
  
  ##--| A simulation instance for simulation condition f -------------
  for (f in 1:nrow(sim_fts)){ # f: index of a simulation condition
    
    # Specify simulation condition
    n_rf_s <- sim_fts[f, "n_rf_s"] # number of items within training or evaluation set
    delta_th <- sim_fts[f,"delta_th"] # change in the ability
    delta_tau <- sim_fts[f,"delta_tau"] # change in the speed
    
    ##--| Generate data with aberrancy according to f -------------
    Ab_Data <- genAbData(p_ab_p=p_ab_p, c.pt=c.pt, 
                         delta_th=delta_th, delta_tau=delta_tau, ipars=null_Data$ipars, ppars=null_Data$ppars)
    
    # Store aberrant data
    Res_w_fts[[f]]$Ab_Data <- Ab_Data
    
    ##--| Calculate monitoring statistics for indicator-based control chart -----------------
    
    # Calculate monitoring statistics with data with aberrancy
    Wj_mat_ab <- calc_W_j_ref_fixed(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s, 
                                     ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                     mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                                     sampling_th_1 = "union_S_R") 
    
    # Calculate monitoring statistics with null data
    Wj_mat_null <- calc_W_j_ref_fixed(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s, 
                                       ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                       mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                                       sampling_th_1 ="union_S_R") 
    # Extract monitoring statistics
    Wj_ab <- Wj_mat_ab$Wj_mat
    Wj_null <- Wj_mat_null$Wj_mat  

    ##--| An empty vector of lists to store results across the different effect-size parameters
    res_CX_list <- vector(mod="list", length=length(ref_v_W_j))
    names(res_CX_list) <- paste0("ref_v=", as.character(ref_v_W_j))
    
    ##--| A simulation instance for effect size l -------------
    for (l in 1:length(ref_v_W_j)){ 
      
      ##--| CUSUM control chart -----------
      # Construct a CUSUM control chart with the null data
      CUSUM_Wj_null <- cusumUp(Wj_null, ref_v=ref_v_W_j[l]) 
      # Construct a CUSUM control chart with the data with aberrancy
      CUSUM_Wj_ab <- cusumUp(Wj_ab, ref_v=ref_v_W_j[l]) 

      ##--| Obtain decision limit (i.e., critical value) -----------
      CV <- calc_CV(pct=.05, stats=CUSUM_Wj_null)
      
      ##--| Evaluate the performance (power, Type I, ARL)
      res_CX_list[[l]] <- Eval_performance(Proc="CUSUM_obs", Res=CUSUM_Wj_ab, CV=CV, 
                                           c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)
    }
    
    # Collapse the result lists and store the results
    res_CX_list <- plyr::ldply(res_CX_list, data.frame)
    Res_w_fts[[f]]$CX <- res_CX_list

    ##--| Weighted control chart that monitors indicator variables  ----------
    tmp_a <- Wj_ab
    tmp_n <- Wj_null 

    ##--| An empty vector of lists to store results across the smoothing parameter ------
    res_CX_W_list <- vector(mode = "list", length = length(lambda_Wj));
    names(res_CX_W_list) <- paste0("lambda=",as.character(lambda_Wj))

    ##--| Exponentially weighted control chart with smoothing parameter l------
    if (!is.null(lambda_Wj)) {
      for (l in 1:length(lambda_Wj)){ # l=1

        Wj_ab[,1] <- (1-lambda_Wj[l])*0 + lambda_Wj[l]*tmp_a[,1]
        Wj_null[,1] <- (1-lambda_Wj[l])*0 + lambda_Wj[l]*tmp_n[,1]

        for (k in 2:ncol(tmp_a)){
          Wj_ab[,k] <- (1-lambda_Wj[l])*Wj_ab[,k-1] + lambda_Wj[l]*tmp_a[,k]
          Wj_null[,k] <- (1-lambda_Wj[l])*Wj_null[,k-1] + lambda_Wj[l]*tmp_n[,k]
        }

        ##--| Obtain decision limit (i.e., critical value) -----------
        CV <- calc_CV(pct=.05, stats=Wj_null)
        
        ##--| Evaluate the performance (power, Type I, ARL) ---------
        res_CX_W_list[[l]] <- Eval_performance(Proc="CUSUM_obs", Res=Wj_ab, CV=CV, 
                                               c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)
      }
    }

    ##--| Arrange and store the results ------------------------------------------------
    res_CX_W_list <- plyr::ldply(res_CX_W_list, data.frame)
    Res_w_fts[[f]]$CX_W <- res_CX_W_list

    ##--| Calculate monitoring statistics for trait-based control chart -----------------
    
    # Calculate monitoring statistics with data with aberrancy
    tmp_ab <- calc_V_j_CUSUM_ref_fixed(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s,
                                       ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                       mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p) 

    V_j_mat_ab <- tmp_ab$V_j_mat 
    th_0_ab <- tmp_ab$th_0; tau_0_ab <- tmp_ab$tau_0
    th_1_mat_ab <- tmp_ab$th_1_mat; tau_1_mat_ab <- tmp_ab$tau_1_mat

    # Calculate monitoring statistics with null data
    tmp_null <- calc_V_j_CUSUM_ref_fixed(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s,
                                         ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                         mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p)
    
    V_j_mat_null <- tmp_null$V_j_mat 
    th_0_null <- tmp_null$th_0; tau_0_null <- tmp_null$tau_0
    th_1_mat_null <- tmp_null$th_1_mat; tau_1_mat_null <- tmp_ab$tau_1_mat

    ##--| Earliest possible detection point for the batch-wise sampling------------
    e_dt_pt  <- c.pt/n_rf_s
    if (e_dt_pt - floor(e_dt_pt) == 0 && e_dt_pt > 1){
      e_dt_pt  <- e_dt_pt -1
    } else if (e_dt_pt-floor(e_dt_pt) != 0 && e_dt_pt > 1){
      e_dt_pt <- floor(e_dt_pt)
    } else {
      e_dt_pt  <- 1
    }
    
    ##--| An empty vector of lists to store results across the different effect-size parameters
    res_CE_list <- vector(mode="list", length=length(ref_v_V_j))
    names(res_CE_list) <- paste0("ref_v=", as.character(ref_v_V_j))
  
    ##--| A simulation instance for effect size l -------------
    for (l in 1:length(ref_v_V_j)){ 
      
      ##--| CUSUM control chart -----------
      # Construct a CUSUM control chart with the null data
      CUSUM_Vj_null <- cusumUp(V_j_mat_null, ref_v=ref_v_V_j[l]) 
      # Construct a CUSUM control chart with the data with aberrancy
      CUSUM_Vj_ab <- cusumUp(V_j_mat_ab, ref_v=ref_v_V_j[l])

      ##--| Obtain decision limit (i.e., critical value) -----------
      CV <- calc_CV(pct=.05, stats=CUSUM_Vj_null)
      
      ##--| Evaluate the performance (power, Type I, ARL)
      res_CE_list[[l]] <- Eval_performance(Proc="CUSUM_est", Res=CUSUM_Vj_ab, CV=CV, 
                                           e_dt_pt=e_dt_pt, c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)
    }

    # Arrange and store the results ----------------------
    res_CE_list <- plyr::ldply(res_CE_list, data.frame)
    Res_w_fts[[f]]$CE <- res_CE_list

    ##--| Weighted control chart that monitors latent trait variables ----------
    tmp_a <- V_j_mat_ab
    tmp_n <- V_j_mat_null
    
    ##--| An empty vector of lists to store results across the smoothing parameter ------
    res_CE_W_list <- vector(mode="list", length=length(lambda_Vj))
    names(res_CE_W_list) <- paste0("lambda=", as.character(lambda_Vj))
    
    ##--| Exponentially weighted control chart with smoothing parameter l------
    if (!is.null(lambda_Vj)){
      for (l in 1:length(lambda_Vj)){
        V_j_mat_null[,1] <- (1-lambda_Vj[l])*0 + lambda_Vj[l]*tmp_n[,1]
        V_j_mat_ab[,1] <- (1-lambda_Vj[l])*0 + lambda_Vj[l]*tmp_a[,1]

        for (k in 2:ncol(tmp_a)){
          V_j_mat_null[,k] <- (1-lambda_Vj[l])*V_j_mat_null[,k-1] + lambda_Vj[l]*tmp_n[,k]
          V_j_mat_ab[,k] <- (1-lambda_Vj[l])*V_j_mat_ab[,k-1] + lambda_Vj[l]*tmp_a[,k]
        }

        ##--| Obtain decision limit (i.e., critical value) -----------
        CV <- calc_CV(pct=.05, stats=V_j_mat_null)
        
        ##--| Evaluate the performance (power, Type I, ARL) ---------
        res_CE_W_list[[l]] <- Eval_performance(Proc="CUSUM_est", Res=V_j_mat_ab, CV=CV, 
                                               e_dt_pt=e_dt_pt, c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)
      }
    }

    ##--| Arrange and store the results ----------------------------
    res_CE_W_list <- plyr::ldply(res_CE_W_list, data.frame)
    Res_w_fts[[f]]$CE_W <- res_CE_W_list
    
    ##--| SGLRT with the sliding door sampling----------------------
    sampling <- "ref_mov_fixed" # sliding door sampling with a fixed training data
    
    # Calculate likelihood-ratio statistics with null data
    LR_null <- calc_SGLRT(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s,
                             ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                             sampling=sampling)

    # Calculate likelihood-ratio statistics with data with aberrancy
    LR_ab <- calc_SGLRT(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s,
                           ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                           sampling=sampling)

    # Extract the likelihood-ratio statistics
    LR_null_mat <- tmp_n <- LR_null$LR_mat 
    LR_ab_mat <- tmp_a <- LR_ab$LR_mat

    ##--| Obtain decision limit (i.e., critical value) -----------
    CV <- calc_CV(pct=.05, stats=LR_null_mat)
    
    ##--| Evaluate the performance (power, Type I, ARL) -----------
    res_G <- Eval_performance(Proc="SGLRT", Res=LR_ab_mat, CV=CV, sampling=sampling,
                              c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)
    
    ##--| Arrange and store the results ---------------------------
    res_G <- cbind(" ", as.data.frame(res_G)); colnames(res_G) <- c(".id", names(res_G[2:ncol(res_G)]))
    Res_w_fts[[f]]$G <- res_G


    ##--| Weighted SGLRT with the sliding door sampling-----------------------
    
    ##--| An empty vector of lists to store results across the smoothing parameter ------
    res_G_W_list <- vector(mode="list", length=length(lambda_G))
    names(res_G_W_list) <- paste0("lambda=",as.character(lambda_G))

    ##--| Exponentially weighted control chart with smoothing parameter l------
    if (!is.null(lambda_G)){
      for (l in 1:length(lambda_G)){
        LR_null_mat[,1] <- (1-lambda_G[l])*0 + lambda_G[l]*tmp_n[,1]
        LR_ab_mat[,1] <- (1-lambda_G[l])*0 + lambda_G[l]*tmp_a[,1]

        for (k in 2:ncol(tmp_n)){
          LR_null_mat[,k] <- (1-lambda_G[l])*LR_null_mat[,k-1] + lambda_G[l]*tmp_n[,k]
          LR_ab_mat[,k] <- (1-lambda_G[l])*LR_ab_mat[,k-1] + lambda_G[l]*tmp_a[,k]
        }

        ##--| Obtain decision limit (i.e., critical value) -----------
        CV <- calc_CV(pct=.05, stats=LR_null_mat)
        
        ##--| Evaluate the performance (power, Type I, ARL) ---------
        res_G_W_list[[l]] <- Eval_performance(Proc="SGLRT", Res=LR_ab_mat, CV=CV, sampling=sampling,
                                              c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)
      }
    }
    
    ##--| Arrange and store the results
    res_G_W_list <- plyr::ldply(res_G_W_list, data.frame)
    Res_w_fts[[f]]$G_W <- res_G_W_list

    ##--| G statistic with the exclusive sampling----------------------
    sampling <- "exclusive" # exclusive item-batch sampling with a fixed training data
    
    # Calculate likelihood-ratio statistics with null data
    LR_null <- calc_SGLRT(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s,
                             ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                             sampling=sampling)
    
    # Calculate likelihood-ratio statistics with data with aberrancy
    LR_ab <- calc_SGLRT(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s,
                           ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                           sampling=sampling)

    # Extract the likelihood-ratio statistics
    LR_null_mat <- tmp_n <- LR_null$LR_mat 
    LR_ab_mat <- tmp_a <- LR_ab$LR_mat
    
    ##--| Obtain decision limit (i.e., critical value) -----------
    CV <- calc_CV(pct=.05, stats=LR_null_mat)
    
    ##--| Evaluate the performance (power, Type I, ARL) -----------
    res_G_excl <- Eval_performance(Proc="SGLRT", Res=LR_ab_mat, CV=CV, e_dt_pt=e_dt_pt, sampling=sampling,
                                      c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)

    ##--| Arrange and store the results ---------------------------
    res_G_excl <- cbind(" ", as.data.frame(res_G_excl));
    colnames(res_G_excl) <- c(".id", names(res_G_excl[2:ncol(res_G_excl)]))
    Res_w_fts[[f]]$G_excl <- res_G_excl

    ##--| Weighted SGLRT with the exclusive sampling ----------------------
    res_G_W_excl_list <- vector(mode="list", length=length(lambda_G_excl))
    names(res_G_W_excl_list) <- paste0("lambda=",as.character(lambda_G_excl))

    ##--| Exponentially weighted control chart with smoothing parameter l------
    if (!is.null(lambda_G_excl)){
      for (l in 1:length(lambda_G_excl)){
        LR_null_mat[,1] <- (1-lambda_G_excl[l])*0 + lambda_G_excl[l]*tmp_n[,1]
        LR_ab_mat[,1] <- (1-lambda_G_excl[l])*0 + lambda_G_excl[l]*tmp_a[,1]

        for (k in 2:ncol(tmp_n)){
          LR_null_mat[,k] <- (1-lambda_G_excl[l])*LR_null_mat[,k-1] + lambda_G_excl[l]*tmp_n[,k]
          LR_ab_mat[,k] <- (1-lambda_G_excl[l])*LR_ab_mat[,k-1] + lambda_G_excl[l]*tmp_a[,k]
        }
        ##--| Obtain decision limit (i.e., critical value) -----------
        CV <- calc_CV(pct=.05, stats=LR_null_mat)
        
        ##--| Evaluate the performance (power, Type I, ARL) ---------
        res_G_W_excl_list[[l]] <- Eval_performance(Proc="SGLRT", Res=LR_ab_mat, CV=CV, e_dt_pt=e_dt_pt, sampling=sampling,
                                                   c.pt=c.pt, n_rf_s=n_rf_s, idx_ab_p=Ab_Data$idx_ab_p)
      }
    }

    ##--| Arrange and store the results
    res_G_W_excl_list <- plyr::ldply(res_G_W_excl_list, data.frame)
    Res_w_fts[[f]]$G_W_excl <- res_G_W_excl_list
    
    ##--| Fox & Marianti (2017) ---------------------------
    mod_hf <- LNIRT::LNIRT(RT=log(Ab_Data$Ab_RT), Y=Ab_Data$Ab_response, data=list(log(Ab_Data$Ab_RT), Ab_Data$Ab_response),
                           alpha=mod_LNIRT$Post.Means$Item.Discrimination,
                           beta=mod_LNIRT$Post.Means$Item.Difficulty,
                           phi=mod_LNIRT$Post.Means$Sigma2,
                           lambda=mod_LNIRT$Post.Means$Time.Intensity,
                           XG=3000, burnin=10,
                           guess=FALSE, par1=TRUE, td=TRUE, WL=TRUE, residual=TRUE)
    
    ##--| Arrange and store the results ------------------------
    FM <- list()
    FM$RT <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p, estdAb=which(mod_hf$EAPCP1 > .95))) # Based on RT
    FM$RA <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p, estdAb=which(mod_hf$EAPCP2 > .95))) # Based on RA
    FM$Both <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p, estdAb=which(mod_hf$EAPCP3 > .95))) # Based on Both
    FM$Union <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p,
                                          estdAb=unique(c(which(mod_hf$EAPCP1 > .95), which(mod_hf$EAPCP2 > .95))) ))
    Res_w_fts[[f]]$FM <- plyr::ldply(FM, data.frame)
    
  } # end of f
   
  return(list(null_Data=null_Data, mod_mbme=mod_mbme,
              sim_fts=sim_fts, rn=rn,
              Res_w_fts=Res_w_fts))
  
} # end of simulation driver


