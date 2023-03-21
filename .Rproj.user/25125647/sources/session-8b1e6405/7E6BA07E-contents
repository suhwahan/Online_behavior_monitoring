##--| Author: Suhwa Han

setwd("/Users/suhwahan/Git_project/JEM_Mar2023_Git")

##--| Sources & libraries ---
source('source_Apr2022.R', echo=TRUE)
source('CUSUM_observed_source_Apr2022.R')
source('CUSUM_estimated_source_Apr2022.R')
source('SGLRT_source_Apr2022.R')
source('/Users/suhwahan/Library/CloudStorage/Box-Box/SuHwa Han/Response Time Simulation/RTM R Program/mmap_cfix.R', echo=TRUE)
source('/Users/suhwahan/Library/CloudStorage/Box-Box/SuHwa Han/Response Time Simulation/RTM R Program/mmap_fsb_cfix.R', echo=TRUE)
source("/Users/suhwahan/Library/CloudStorage/Box-Box/SuHwa Han/Response Time Simulation/RTM R Program/map.R")

library(LNIRT)
library(foreach)
library(doParallel)
library(MASS)
library(plyr)
library(tidyverse)

# # ##--| For debugging ---------------
# nitem=40; nexaminee=1000; rho_item=0.25; rho_person=0.3; c_f=0; p_ab_p=0.1; c.pt=21;
# mu_a=0.91; var_a=0.79^2; mu_alp=1.49; var_alp=0.39^2;
# rn_vec <- get(load("randN_r100.RData")); rm(rn)
# rn_vec=rn_vec;
# ref_v_W_j=c(.25, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4);
# lambda_Wj=seq(0.005, 0.995, by = 0.005);
# ref_v_V_k=c(0.5, 1, 2.5, 5, 7.5, 10, 12.5, 15);
# lambda_Vk=seq(0.005, 0.995, by = 0.005);
# lambda_G=seq(0.005, 0.995, by = 0.005);
# lambda_G_excl=seq(0.005, 0.995, by = 0.005)
# dat_grp_ind="grp"

##--| Simulation driver
Sim_driver <- function(nitem, nexaminee, rho_item, rho_person, c_f, p_ab_p, c.pt,
                       rn_vec,
                       mu_a, var_a, mu_alp, var_alp,
                       ref_v_W_j, lambda_Wj=NULL,
                       test_type = "two_samples", ref_v_V_k, lambda_Vk=NULL,
                       lambda_G=NULL, dat_grp_ind="grp",
                       lambda_G_excl=NULL){
  
  # r <- 1
  rn <- rn_vec[r] #rn <- sample(seq(1,10000), size=1)
  set.seed(rn)
  
  ##--| Generate null data -----------------
  null_Data <- nullDataGen(nitem=nitem, nexaminee=nexaminee, rho_item=rho_item, rho_person=rho_person, c_f=c_f,
                           mu_a=mu_a, var_a=var_a, mu_alp=mu_alp, var_alp=var_alp)
  
  #--| Calibration ----------------
  start_time <- Sys.time()
  mod_mbme <- mmap_cfix(null_Data$response, null_Data$rt,
                        ppar_prior=list(mu_p=null_Data$mu_p, cov_p=null_Data$cov_p),
                        iparst_prior = list(mu_is=null_Data$mu_item_str, cov_is=null_Data$cov_item_str),
                        cfix=c_f, tol=list(max_em=1000, tol_em=0.005, loglike=1e-4, max_nr=500, tol_nr=0.005))
  end_time <- Sys.time()
  elps_time_mbme <- end_time-start_time

  ir_est_mat <- cbind(mod_mbme$iest[,"a"], mod_mbme$iest[,"b"], mod_mbme$iest[,"c"]); colnames(ir_est_mat) <- c("a", "b", "c")
  rt_est_mat <- cbind(mod_mbme$iest[,"alp"], mod_mbme$iest[,"bet"]); colnames(rt_est_mat) <- c("alpha", "beta")
  
  # plot(null_Data$ipars[,"a"], mod_mbme$iest[,"a"]); abline(0,1)
  # plot(null_Data$ipars[,"b"], mod_mbme$iest[,"b"]); abline(0,1)
  # plot(null_Data$ipars[,"alp"], mod_mbme$iest[,"alp"]); abline(0,1)
  # plot(null_Data$ipars[,"bet"], mod_mbme$iest[,"bet"]); abline(0,1)
  
  mod_LNIRT <- LNIRT::LNIRT(RT=log(null_Data$rt), Y=null_Data$response,
                            data=list(as.matrix(null_Data$rt), as.matrix(null_Data$response)),
                            XG=3000, burnin=10,
                            guess=FALSE, par1=TRUE, td=FALSE, WL=TRUE, residual = F)
  
  # plot(null_Data$ipars[,"a"], mod_LNIRT$Post.Means$Item.Discrimination); abline(0,1)
  # plot(null_Data$ipars[,"b"], mod_LNIRT$Post.Means$Item.Difficulty); abline(0,1)
  
  
  ##--| Implementing different effect size -----------------
  
  delta <- data.frame(c(-1, -1.5, -2), c(1, 1.5, 2)); colnames(delta) <- c("th", "tau")
  ##delta <- data.frame(c(-0.5), c(0.5)); colnames(delta) <- c("th", "tau")
  n_rf_s <- c(5, 10) 
  sim_fts <- data.frame(merge(n_rf_s, delta)); colnames(sim_fts) <- c("n_rf_s", "delta_th", "delta_tau")

  # prob_th <- c(0.3, 0.2, 0.1)
  # n_rf_s <- c(5, 10) 
  # sim_fts <- data.frame(merge(n_rf_s, prob_th)); colnames(sim_fts) <- c("n_rf_s", "prob_th")
  # 
  Res_w_fts <- vector(mode="list", length=nrow(sim_fts))
  
  for (f in 1:nrow(sim_fts)){ ## f=2; f=1; f=3
    
    ##--| Generate aberrant data -----------------
    n_rf_s <- sim_fts[f, "n_rf_s"]
    
    # i) Change in traits
    delta_th <- sim_fts[f,"delta_th"]
    delta_tau <- sim_fts[f,"delta_tau"]
    Ab_Data <- genAbData(p_ab_p=p_ab_p, c.pt=c.pt, # p_ab_p=.1; f=1 # p_ab_p=.05; f=2
                         delta_th=delta_th, delta_tau=delta_tau, ipars=null_Data$ipars, ppars=null_Data$ppars)
    
    # ii) Change in observed data
    # prob_th <- sim_fts[f,"prob_th"]  
    # Ab_Data <- genAbData_obs(p_ab_p=p_ab_p, c.pt=c.pt, ipars=null_Data$ipars, ppars=null_Data$ppars, prob_th=prob_th)
    # barplot(Ab_Data$Ab_response[Ab_Data$idx_ab_p[2],]) 
    # barplot(Ab_Data$Ab_RT[Ab_Data$idx_ab_p[7],]) 
    # barplot(Ab_Data$Ab_RT[Ab_Data$idx_nab_p[12],]) 
    
    # ##--| Calibration ----------------
    # start_time <- Sys.time()
    # mod_mbme <- mmap_cfix(Ab_Data$Ab_response, Ab_Data$Ab_RT,
    #                       ppar_prior=list(mu_p=null_Data$mu_p, cov_p=null_Data$cov_p),
    #                       iparst_prior = list(mu_is=null_Data$mu_item_str, cov_is=null_Data$cov_item_str),
    #                       cfix=c_f, tol=list(max_em=1000, tol_em=0.005, loglike=1e-4, max_nr=500, tol_nr=0.005))
    # end_time <- Sys.time()
    # elps_time_mbme <- end_time-start_time
    # 
    # ir_est_mat <- cbind(mod_mbme$iest[,"a"], mod_mbme$iest[,"b"], mod_mbme$iest[,"c"]); colnames(ir_est_mat) <- c("a", "b", "c")
    # rt_est_mat <- cbind(mod_mbme$iest[,"alp"], mod_mbme$iest[,"bet"]); colnames(rt_est_mat) <- c("alpha", "beta")
    # 
    
    ##--| Store aberrant data
    Res_w_fts[[f]]$Ab_Data <- Ab_Data
    
    ##--| C^X statistic -----------------
    # Wj_mat_ab <- calc_W_j_ref_inc(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s, # hist(log(Ab_Data$Ab_RT))
    #                                 ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
    #                                 mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
    #                                 sampling_th_1 = "union_S_R") # hist(Wj_mat_ab)
    
    # Wj_mat_null <- calc_W_j_ref_inc(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s, # hist(log(Ab_Data$Ab_RT))
    #                                   ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
    #                                   mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
    #                                   sampling_th_1 ="union_S_R")
    
    Wj_mat_ab <- calc_W_j_ref_fixed(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s, # hist(log(Ab_Data$Ab_RT))
                                     ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                     mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                                     sampling_th_1 = "union_S_R") # hist(Wj_mat_ab)

    Wj_mat_null <- calc_W_j_ref_fixed(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s, # hist(log(Ab_Data$Ab_RT))
                                       ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                       mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                                       sampling_th_1 ="union_S_R") # hist(Wj_mat_null)

    Wj_ab <- Wj_mat_ab$Wj_mat
    Wj_null <- Wj_mat_null$Wj_mat  # hist(Wj_mat_null$Wj_mat)

    res_CX_list <- vector(mod="list", length=length(ref_v_W_j))
    names(res_CX_list) <- paste0("ref_v=", as.character(ref_v_W_j))

    for (l in 1:length(ref_v_W_j)){ # l=3
      CUSUM_Wj_null <- cusumUp(Wj_null, ref_v=ref_v_W_j[l]) # hist(CUSUM_Wj_null)
      CUSUM_Wj_ab <- cusumUp(Wj_ab, ref_v=ref_v_W_j[l]) # hist(Wj_mat_ab$Wj_mat)

      CV <- calc_CV(pct=.05, stats=CUSUM_Wj_null)
      res_CX_list[[l]] <- calcPw_CUSUM_observed(Res=CUSUM_Wj_ab, CV=CV, c.pt=c.pt, n_rf_s=n_rf_s, st_eval=n_rf_s+1, idx_ab_p=Ab_Data$idx_ab_p)

    }

    res_CX_list <- plyr::ldply(res_CX_list, data.frame)
    Res_w_fts[[f]]$CX <- res_CX_list

    ##--| C^O statistic with weight ------------------
    tmp_a <- Wj_ab
    tmp_n <- Wj_null # View(Wj_mat)

    res_CX_W_list <- vector(mode = "list", length = length(lambda_Wj));
    names(res_CX_W_list) <- paste0("lambda=",as.character(lambda_Wj))

    if (!is.null(lambda_Wj)) {
      for (l in 1:length(lambda_Wj)){ # l=1

        Wj_ab[,1] <- (1-lambda_Wj[l])*0 + lambda_Wj[l]*tmp_a[,1]
        Wj_null[,1] <- (1-lambda_Wj[l])*0 + lambda_Wj[l]*tmp_n[,1]

        for (k in 2:ncol(tmp_a)){
          Wj_ab[,k] <- (1-lambda_Wj[l])*Wj_ab[,k-1] + lambda_Wj[l]*tmp_a[,k]
          Wj_null[,k] <- (1-lambda_Wj[l])*Wj_null[,k-1] + lambda_Wj[l]*tmp_n[,k]
        }

        CV <- calc_CV(pct=.05, stats=Wj_null)
        res_CX_W_list[[l]] <- calcPw_CUSUM_observed(Res=Wj_ab, CV=CV, c.pt=c.pt, n_rf_s=n_rf_s, st_eval=n_rf_s+1, idx_ab_p=Ab_Data$idx_ab_p)

      }
    }

    res_CX_W_list <- plyr::ldply(res_CX_W_list, data.frame)
    Res_w_fts[[f]]$CX_W <- res_CX_W_list

    ##--| C^E statistic -----------------
    # tmp_ab <- calc_V_k_CUSUM_ref_inc(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s,
    #                                    ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
    #                                    mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
    #                                    test_type="two_samples") # View(V_k_mat_ab)
    
    tmp_ab <- calc_V_k_CUSUM_ref_fixed(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s,
                                       ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                       mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                                       test_type="two_samples") # View(V_k_mat_ab)

    V_k_mat_ab <- tmp_ab$V_k_mat #mean(V_k_mat_ab)
    th_0_ab <- tmp_ab$th_0; tau_0_ab <- tmp_ab$tau_0
    th_1_mat_ab <- tmp_ab$th_1_mat; tau_1_mat_ab <- tmp_ab$tau_1_mat

    # tmp_null <- calc_V_k_CUSUM_ref_inc(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s,
    #                                      ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
    #                                      mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
    #                                      test_type = "two_samples")
    
    tmp_null <- calc_V_k_CUSUM_ref_fixed(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s,
                                         ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat,
                                         mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                                         test_type = "two_samples")
    
    V_k_mat_null <- tmp_null$V_k_mat # hist(V_k_mat_null)  # sd(V_k_mat_null)
    th_0_null <- tmp_null$th_0; tau_0_null <- tmp_null$tau_0
    th_1_mat_null <- tmp_null$th_1_mat; tau_1_mat_null <- tmp_ab$tau_1_mat

    x <- c.pt/n_rf_s
    if (x-floor(x) == 0 && x > 1){
      x <- x-1
    } else if (x-floor(x) != 0 && x > 1){
      x <- floor(x)
    } else {
      x <- 1
    }
    e_dt_pt <- x  # earliest possible detection point
    
    res_CE_list <- vector(mode="list", length=length(ref_v_V_k))
    names(res_CE_list) <- paste0("ref_v=", as.character(ref_v_V_k))

    for (l in 1:length(ref_v_V_k)){ # l=4
      CUSUM_Vk_null <- cusumUp(V_k_mat_null, ref_v=ref_v_V_k[l]) # hist(CUSUM_Vk_null)
      CUSUM_Vk_ab <- cusumUp(V_k_mat_ab, ref_v=ref_v_V_k[l])

      CV <- calc_CV(pct=.05, stats=CUSUM_Vk_null)
      res_CE_list[[l]] <- calcPw_estimated_CUSUM(Res=CUSUM_Vk_ab, CV=CV, c.pt=c.pt, e_dt_pt=e_dt_pt, idx_ab_p=Ab_Data$idx_ab_p)
    }

    res_CE_list <- plyr::ldply(res_CE_list, data.frame)
    Res_w_fts[[f]]$CE <- res_CE_list


    ##--| C^eta statistic with weight ------------------------
    tmp_a <- V_k_mat_ab
    tmp_n <- V_k_mat_null

    res_CE_W_list <- vector(mode="list", length=length(lambda_Vk))
    names(res_CE_W_list) <- paste0("lambda=", as.character(lambda_Vk))

    if (!is.null(lambda_Vk)){
      for (l in 1:length(lambda_Vk)){
        V_k_mat_null[,1] <- (1-lambda_Vk[l])*0 + lambda_Vk[l]*tmp_n[,1]
        V_k_mat_ab[,1] <- (1-lambda_Vk[l])*0 + lambda_Vk[l]*tmp_a[,1]

        for (k in 2:ncol(tmp_a)){
          V_k_mat_null[,k] <- (1-lambda_Vk[l])*V_k_mat_null[,k-1] + lambda_Vk[l]*tmp_n[,k]
          V_k_mat_ab[,k] <- (1-lambda_Vk[l])*V_k_mat_ab[,k-1] + lambda_Vk[l]*tmp_a[,k]
        }

        CV <- calc_CV(pct=.05, stats=V_k_mat_null)
        res_CE_W_list[[l]] <- calcPw_estimated_CUSUM(Res=V_k_mat_ab, CV=CV, c.pt=c.pt, e_dt_pt=e_dt_pt, idx_ab_p=Ab_Data$idx_ab_p)
      }
    }

    res_CE_W_list <- plyr::ldply(res_CE_W_list, data.frame)
    Res_w_fts[[f]]$CE_W <- res_CE_W_list
    
    ##--| G statistic with the sliding door sampling----------------------
    sampling <- "ref_mov_fixed"
    LR_null <- calcLR_SGLRT(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s,
                             ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                             sampling=sampling, dat_grp_ind=dat_grp_ind, null_ppars=null_Data$ppars, c.pt=c.pt,
                             idx_ab_p=Ab_Data$idx_ab_p, idx_nab_p=Ab_Data$idx_nab_p, delta_th=delta_th, delta_tau=delta_tau)

    LR_ab <- calcLR_SGLRT(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s,
                           ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                           sampling=sampling, dat_grp_ind=dat_grp_ind, null_ppars=null_Data$ppars, c.pt=c.pt,
                           idx_ab_p=Ab_Data$idx_ab_p, idx_nab_p=Ab_Data$idx_nab_p, delta_th=delta_th, delta_tau=delta_tau)

    LR_null_mat <- tmp_n <- LR_null$LR_mat # View(LR_mat)
    LR_ab_mat <- tmp_a <- LR_ab$LR_mat

    CV <- calc_CV(pct=.05, stats=LR_null_mat)
    res_G <- calcPw_SGLRT(Res=LR_ab_mat, CV=CV, e_dt_pt_excl=e_dt_pt, sampling=sampling,
                                 c.pt=c.pt, st_eval=LR_ab$st_eval, idx_ab_p=Ab_Data$idx_ab_p)

    res_G <- cbind(" ", as.data.frame(res_G)); colnames(res_G) <- c(".id", names(res_G[2:ncol(res_G)]))
    Res_w_fts[[f]]$G <- res_G


    ##--| G statistic with weight with the sliding door sampling-----------------------
    res_G_W_list <- vector(mode="list", length=length(lambda_G))
    names(res_G_W_list) <- paste0("lambda=",as.character(lambda_G))

    if (!is.null(lambda_G)){
      for (l in 1:length(lambda_G)){
        LR_null_mat[,1] <- (1-lambda_G[l])*0 + lambda_G[l]*tmp_n[,1]
        LR_ab_mat[,1] <- (1-lambda_G[l])*0 + lambda_G[l]*tmp_a[,1]

        for (k in 2:ncol(tmp_n)){
          LR_null_mat[,k] <- (1-lambda_G[l])*LR_null_mat[,k-1] + lambda_G[l]*tmp_n[,k]
          LR_ab_mat[,k] <- (1-lambda_G[l])*LR_ab_mat[,k-1] + lambda_G[l]*tmp_a[,k]
        }

        CV <- calc_CV(pct=.05, stats=LR_null_mat)
        res_G_W_list[[l]] <- calcPw_SGLRT(Res=LR_ab_mat, CV=CV, e_dt_pt_excl=e_dt_pt, sampling=sampling,
                                          c.pt=c.pt, st_eval=LR_ab$st_eval, idx_ab_p=Ab_Data$idx_ab_p)
      }
    }

    res_G_W_list <- plyr::ldply(res_G_W_list, data.frame)
    Res_w_fts[[f]]$G_W <- res_G_W_list

    ##--| G statistic with the exclusive sampling----------------------
    sampling <- "exclusive"
    LR_null <- calcLR_SGLRT(resp=null_Data$response, rt=null_Data$rt, n_rf_s=n_rf_s,
                             ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                             sampling=sampling, dat_grp_ind=dat_grp_ind, null_ppars=null_Data$ppars, c.pt=c.pt,
                             idx_ab_p=Ab_Data$idx_ab_p, idx_nab_p=Ab_Data$idx_nab_p, delta_th=delta_th, delta_tau=delta_tau)

    LR_ab <- calcLR_SGLRT(resp=Ab_Data$Ab_response, rt=Ab_Data$Ab_RT, n_rf_s=n_rf_s,
                           ir_est_mat=ir_est_mat, rt_est_mat=rt_est_mat, mu_p_t=null_Data$mu_p, cov_p_t=null_Data$cov_p,
                           sampling=sampling, dat_grp_ind=dat_grp_ind, null_ppars=null_Data$ppars, c.pt=c.pt,
                           idx_ab_p=Ab_Data$idx_ab_p, idx_nab_p=Ab_Data$idx_nab_p, delta_th=delta_th, delta_tau=delta_tau)


    LR_null_mat <- tmp_n <- LR_null$LR_mat # View(LR_mat)
    LR_ab_mat <- tmp_a <- LR_ab$LR_mat

    CV <- calc_CV(pct=.05, stats=LR_null_mat)
    res_G_excl <- calcPw_SGLRT(Res=LR_ab_mat, CV=CV, e_dt_pt_excl=e_dt_pt, sampling=sampling,
                                      c.pt=c.pt, st_eval=LR_ab$st_eval, idx_ab_p=Ab_Data$idx_ab_p)

    res_G_excl <- cbind(" ", as.data.frame(res_G_excl));
    colnames(res_G_excl) <- c(".id", names(res_G_excl[2:ncol(res_G_excl)]))
    Res_w_fts[[f]]$G_excl <- res_G_excl

    ##--| G statistic with weight with the exclusive sampling ----------------------
    res_G_W_excl_list <- vector(mode="list", length=length(lambda_G_excl))
    names(res_G_W_excl_list) <- paste0("lambda=",as.character(lambda_G_excl))

    if (!is.null(lambda_G_excl)){
      for (l in 1:length(lambda_G_excl)){
        LR_null_mat[,1] <- (1-lambda_G_excl[l])*0 + lambda_G_excl[l]*tmp_n[,1]
        LR_ab_mat[,1] <- (1-lambda_G_excl[l])*0 + lambda_G_excl[l]*tmp_a[,1]

        for (k in 2:ncol(tmp_n)){
          LR_null_mat[,k] <- (1-lambda_G_excl[l])*LR_null_mat[,k-1] + lambda_G_excl[l]*tmp_n[,k]
          LR_ab_mat[,k] <- (1-lambda_G_excl[l])*LR_ab_mat[,k-1] + lambda_G_excl[l]*tmp_a[,k]
        }

        CV <- calc_CV(pct=.05, stats=LR_null_mat)
        res_G_W_excl_list[[l]] <- calcPw_SGLRT(Res=LR_ab_mat, CV=CV, e_dt_pt_excl=e_dt_pt, sampling=sampling,
                                               c.pt=c.pt, st_eval=LR_ab$st_eval, idx_ab_p=Ab_Data$idx_ab_p)
      }
    }

    res_G_W_excl_list <- plyr::ldply(res_G_W_excl_list, data.frame)
    Res_w_fts[[f]]$G_W_excl <- res_G_W_excl_list
    
    ##--| Fox & Marianti (2017)
    mod_hf <- LNIRT::LNIRT(RT=log(Ab_Data$Ab_RT), Y=Ab_Data$Ab_response, data=list(log(Ab_Data$Ab_RT), Ab_Data$Ab_response),
                           alpha=mod_LNIRT$Post.Means$Item.Discrimination,
                           beta=mod_LNIRT$Post.Means$Item.Difficulty,
                           phi=mod_LNIRT$Post.Means$Sigma2,
                           lambda=mod_LNIRT$Post.Means$Time.Intensity,
                           XG=3000, burnin=10,
                           guess=FALSE, par1=TRUE, td=TRUE, WL=TRUE, residual=TRUE)
    
    # mod_hf <- LNIRT::LNIRT(RT=log(Ab_Data$Ab_RT), Y=Ab_Data$Ab_response, data=list(log(Ab_Data$Ab_RT), Ab_Data$Ab_response),
    #                        #alpha=mod_LNIRT$Post.Means$Item.Discrimination,
    #                        #beta=mod_LNIRT$Post.Means$Item.Difficulty,
    #                        #phi=mod_LNIRT$Post.Means$Sigma2,
    #                        #lambda=mod_LNIRT$Post.Means$Time.Intensity,
    #                        XG=3000, burnin=10,
    #                        guess=FALSE, par1=TRUE, td=TRUE, WL=TRUE, residual=TRUE)
    # 
    # # #summary(mod_hf)
    # # 
    FM <- list()
    FM$RT <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p, estdAb=which(mod_hf$EAPCP1 > .95))) # Based on RT
    FM$RA <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p, estdAb=which(mod_hf$EAPCP2 > .95))) # Based on RA
    FM$Both <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p, estdAb=which(mod_hf$EAPCP3 > .95))) # Based on Both
    FM$Union <- as.data.frame(calc_pfcase(nexaminee=nexaminee, trueAb=Ab_Data$idx_ab_p,
                                          estdAb=unique(c(which(mod_hf$EAPCP1 > .95), which(mod_hf$EAPCP2 > .95))) ))
    Res_w_fts[[f]]$FM <- plyr::ldply(FM, data.frame)
    
    #Res_w_fts[[f]]$mod_hf <- mod_hf

  } # end of f
   
  return(list(null_Data=null_Data, mod_mbme=mod_mbme,
              sim_fts=sim_fts, rn=rn,
              Res_w_fts=Res_w_fts))
  
} # end of simulation driver


