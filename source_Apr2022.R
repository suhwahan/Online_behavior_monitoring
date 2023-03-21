### Update history ###
# Apr/25/2022: Migrate to a new wd to be in line with Sim_All_Apr2022.R

# test1: practice commit --amend
# test2: practice switch branch to oldies
# test3: another commit under oldies
# test: emptyplaylist
# test: just deleted the line that was here and replaced with this line
# test: compare git diff vs git diff HEAD
# test: Make changes to odd numbers
# Test: Make changes to the rainbow branch
# Test: Make another change to the rainbow branch

### Item response function
irf_3pl <- function(th, ipars, D){
  a=ipars[1]; b=ipars[2]; c=ipars[3]
  p <- c + (1 - c)/(1 + exp(-D * (a * (th - b))))
  return(p)
}

### Response time function

# for one RT observation
rt_fn <- function(tau, time, xi_rt){
  alpha <- xi_rt[1]; beta <- xi_rt[2]
  p_rt <- alpha/(time * sqrt(2*pi)) * exp(- 1/2 * alpha^2 * (log(time) - (beta - tau))^2)
  #p_rt <- dlnorm(time, meanlog = beta - tau, 1/alpha)
  return(p_rt)
}

# for an RT vector
rt_dfn <- function(rt, Alp, Bet, Tau){ 
  n_dt_p <- length(rt)  
  density <- vector("numeric", length = n_dt_p)
  for (i in 1:n_dt_p){ # i=1
    density[i] <- (Alp[i]/(rt[i]*sqrt(2*pi)))*exp(-0.5*((Alp[i]*(log(rt[i])-(Bet[i]-Tau)))^2))
  } # dlnorm(rt, Bet-Tau,1/Alp)
  return(density)
}

### Generating responses based on 3PL---------
genResp3PL <- function(a, b, c, theta, D = 1.702){ # a, b, c, theta as vectors
  Nexaminee <- length(theta)
  Nitem <- length(a)
  prob <- matrix(NA, Nexaminee, Nitem)
  for (j in 1:Nexaminee){
    prob[j,] <- c + (1 - c)/(1 + exp(-D*a * (theta[j] - b)))
  }
  response <- (matrix(runif(Nexaminee*Nitem), Nexaminee, Nitem) < prob) * 1
  colnames(response) <- paste0("Item", seq(1,Nitem))
  return(list(response = response, prob = prob))
}

### Generating response times based on the log-normal model--------
genLNRTs <- function (alp, bet, tau){
  Nexaminee <- length(tau) # length(ppars[,"tau"])
  Nitem <- length(alp) # length(ipars[,"alpha"])
  RT <- matrix(NA, Nexaminee, Nitem)
  for (j in 1:Nexaminee){
    RT[j,] <- rlnorm(Nitem, bet-tau[j], 1/alp) # pdf function; variance was much larger than the actual item parameter values
  }
  return(list(RT=RT))
}

### Null data generation function--------
nullDataGen <- function(nitem, nexaminee, rho_item, rho_person, c_f, 
                        mu_a, var_a, mu_alp, var_alp){
  
  ### Hyperparameters ---
  
  ### Kang (2017, 2020): item parameters are assumed to follow a multivariate normal distribution
  mu_a = mu_a; var_a = var_a; mu_alp = mu_alp; var_alp = var_alp; 
  mu_b = 0; var_b = 1; mu_bet = 0; var_bet = 1; 
  mu_cstr = -1.386; var_cstr = 0.04; # correspond to c=0.2
  rho_item = rho_item; 
  rho_person = rho_person; mu_person <- c(0,0)
  c_f <- c_f # constant value for the guessing parameter
  trunc_item = matrix(c(0.3, 3, -4, 4, 0.0, 0.3), byrow = T, ncol = 2, nrow = 3)
  
  
  ### Generating item parameters for responses & RTs based on multivariate normal ---
  
  ### transform a, alpha, c parameters to log/logit scale for appropriate domain for priors
  mu_astr <- log(mu_a^2/sqrt(var_a + mu_a^2))
  var_astr <- log(var_a/mu_a^2 + 1)
  mu_alpstr <- log(mu_alp^2/sqrt(var_alp + mu_alp^2))
  var_alpstr <- log(var_alp/mu_alp^2 + 1)
  
  mu_item_str <- c(mu_astr, mu_b, mu_cstr, mu_alpstr, mu_bet); 
  names(mu_item_str) <- c("mu_astr", "mu_b", "mu_cstr", "mu_alpstr", "mu_bet")
  var_item_str <- c(var_astr, var_b, var_cstr, var_alpstr, var_bet)
  names(var_item_str) <- c("var_astr", "var_b", "var_cstr", "var_alpstr", "var_bet")
  
  ## var-cov matrix for the item parameters 
  rhomat <- matrix(rho_item, length(mu_item_str), length(mu_item_str)); diag(rhomat) <- 1
  cov_item_str <- sqrt(as.matrix(var_item_str) %*% var_item_str) * rhomat
  
  ## Fix c to a constant
  if(!is.null(c_f)){
    mu_item_str <- mu_item_str[-3]
    cov_item_str <- cov_item_str[-3,]; cov_item_str <- cov_item_str[,-3]
  }
  
  while (1 > 0) {
    
    if (!is.null(c_f)){
      partemp <- mvtnorm::rmvnorm(nitem, mean = mu_item_str, sigma = cov_item_str)
      ipars <-  cbind(exp(partemp[,1]), partemp[,2], c_f, # transform back to the original scale
                      exp(partemp[,3]), partemp[,4])
    } else {
      partemp <- mvtnorm::rmvnorm(nitem, mean = mu_item_str, sigma = cov_item_str)
      ipars <-  cbind(exp(partemp[,1]), partemp[,2], 1/(1 + exp(-partemp[,3])), # transform back to the original scale
                      exp(partemp[,4]), partemp[,5])
    }
    colnames(ipars) <- c("a", "b", "c", "alp", "bet")
    if(min(ipars[,"a"]) >= trunc_item[1,1] && max(ipars[,"a"]) <= trunc_item[1,2] &&
       min(ipars[,"alp"]) >= trunc_item[1,1] && max(ipars[,"alp"]) <= trunc_item[1,2] &&
       max(abs(ipars[,"b"])) <= trunc_item[2,2] && max(abs(ipars[,"bet"])) <= trunc_item[2,2] &&
       if (!is.null(c_f)){
         min(ipars[,"c"]) >= trunc_item[3,1] && max(ipars[,"c"]) <= trunc_item[3,2] 
       }
    ) break;
  }
  
  rm(partemp)
  
  ### Generating person parameters for responses & RTs ---

  cov_person <- matrix(c(1, rho_person, rho_person, 1), 2, 2)
  
  while (1 > 0) {
    ppars <- mvtnorm::rmvnorm(nexaminee, mean = mu_person, sigma = cov_person)
    if(min(ppars) >= -4 && max(ppars) <= 4) break
  }
  
  colnames(ppars) <- c("th", "tau") 
  
  
  ### Generating null responses ---
  resp <- genResp3PL(a=ipars[,"a"], b=ipars[,"b"], c=ipars[,"c"], theta=ppars[,"th"])
  
  ### Generating null RTs ---
  
  RT <- genLNRTs(alp=ipars[,"alp"], bet=ipars[,"bet"], tau=ppars[,"tau"])
  
  #dat <- gendata(ppars, ipars, seednum=12345)
  
  return( list(ipars = ipars, ppars = ppars, 
               response = resp$response, rt = RT$RT, mu_p = mu_person, cov_p = cov_person,
               mu_item_str = mu_item_str, cov_item_str = cov_item_str) )
}

### Generating aberrant data from a c.pt until the end
genAbData <- function(p_ab_p, c.pt, delta_th, delta_tau, ipars, ppars){
  
  Nexaminee <- dim(ppars)[1]
  Nitem <- dim(ipars)[1]
  
  idx_ab_p <- sort(sample(x = 1:Nexaminee, size = Nexaminee*p_ab_p, replace = FALSE)) # View(idx_ab_p)
  idx_nab_p <- setdiff(seq(1:Nexaminee), idx_ab_p)
  
  Ab_response <- matrix(NA, Nexaminee, Nitem)
  Ab_RT <- matrix(NA, Nexaminee, Nitem)
  Ab_prob <- matrix(NA, Nexaminee, Nitem)
  
  ### Fill in for aberrant examinees
  temp <- genResp3PL(a=ipars[,"a"][1:(c.pt-1)], b=ipars[,"b"][1:(c.pt-1)], c=ipars[,"c"][1:(c.pt-1)], theta=ppars[,"th"][idx_ab_p], D=1.702) 
  Ab_response[idx_ab_p,1:(c.pt-1)] <- temp$response
  Ab_prob[idx_ab_p,1:(c.pt-1)] <- temp$prob
  temp <- genResp3PL(a=ipars[,"a"][c.pt:Nitem], b=ipars[,"b"][c.pt:Nitem], c=ipars[,"c"][c.pt:Nitem], theta=ppars[,"th"][idx_ab_p]+delta_th, D=1.702) 
  Ab_response[idx_ab_p,c.pt:Nitem] <- temp$response
  Ab_prob[idx_ab_p,c.pt:Nitem] <- temp$prob
  
  temp <- genLNRTs(alp=ipars[,"alp"][1:(c.pt-1)], bet=ipars[,"bet"][1:(c.pt-1)], tau=ppars[,"tau"][idx_ab_p])
  Ab_RT[idx_ab_p,1:(c.pt-1)] <- temp$RT
  temp <- genLNRTs(alp=ipars[,"alp"][c.pt:Nitem], bet=ipars[,"bet"][c.pt:Nitem], tau=ppars[,"tau"][idx_ab_p]+delta_tau)
  Ab_RT[idx_ab_p,c.pt:Nitem] <- temp$RT
  
  ### Fill in for non-aberrant examinees
  temp <- genResp3PL(a=ipars[,"a"], b=ipars[,"b"], c=ipars[,"c"], theta=ppars[,"th"][idx_nab_p], D=1.702) 
  Ab_response[idx_nab_p,] <- temp$response
  Ab_prob[idx_nab_p,] <- temp$prob
  
  temp <- genLNRTs(alp=ipars[,"alp"], bet=ipars[,"bet"], tau=ppars[,"tau"][idx_nab_p])
  Ab_RT[idx_nab_p,] <- temp$RT
  
  return(list(Ab_response=Ab_response, Ab_RT=Ab_RT, Ab_prob=Ab_prob,
              idx_ab_p=idx_ab_p, idx_nab_p=idx_nab_p))
}

### Generating responses based on uniform probability---------
genResp3PL_prob_th <- function(a, b, c, theta, D = 1.702, prob_th){ # a, b, c, theta as vectors # prob_threshold
  Nexaminee <- length(theta)
  Nitem <- length(a)
  # prob <- matrix(NA, Nexaminee, Nitem)
  # for (j in 1:Nexaminee){
  #   prob[j,] <- c + (1 - c)/(1 + exp(-D*a * (theta[j] - b)))
  # }
  response <- (matrix(runif(Nexaminee*Nitem), Nexaminee, Nitem) < prob_th) * 1
  colnames(response) <- paste0("Item", seq(1,Nitem))
  return(list(response = response))
}

# response <- (matrix(runif(Nexaminee*Nitem), Nexaminee, Nitem) < 0.3) * 1
# mean(response)
# response <- (matrix(runif(Nexaminee*Nitem), Nexaminee, Nitem) < 0.2) * 1
# mean(response)
# response <- (matrix(runif(Nexaminee*Nitem), Nexaminee, Nitem) < 0.1) * 1
# mean(response)

### Generating aberrant data by directly manipulating observed values
genAbData_obs <- function(p_ab_p, c.pt, ipars, ppars, prob_th){
  
  Nexaminee <- dim(ppars)[1]
  Nitem <- dim(ipars)[1]
  
  idx_ab_p <- sort(sample(x = 1:Nexaminee, size = Nexaminee*p_ab_p, replace = FALSE)) # View(idx_ab_p)
  idx_nab_p <- setdiff(seq(1:Nexaminee), idx_ab_p)
  
  Ab_response <- matrix(NA, Nexaminee, Nitem)
  Ab_RT <- matrix(NA, Nexaminee, Nitem)
  Ab_prob <- matrix(NA, Nexaminee, Nitem)
  
  ### Responses
  # aberrant examinee
  temp <- genResp3PL(a=ipars[,"a"][1:(c.pt-1)], b=ipars[,"b"][1:(c.pt-1)], c=ipars[,"c"][1:(c.pt-1)], theta=ppars[,"th"][idx_ab_p], D=1.702) 
  Ab_response[idx_ab_p,1:(c.pt-1)] <- temp$response
  #Ab_prob[idx_ab_p,1:(c.pt-1)] <- temp$prob
  temp <- genResp3PL_prob_th(a=ipars[,"a"][c.pt:Nitem], b=ipars[,"b"][c.pt:Nitem], c=ipars[,"c"][c.pt:Nitem], theta=ppars[,"th"][idx_ab_p], 
                             prob_th=prob_th, D=1.702) 
  Ab_response[idx_ab_p,c.pt:Nitem] <- temp$response
  #Ab_prob[idx_ab_p,c.pt:Nitem] <- temp$prob
  
  # normal examinee
  temp <- genResp3PL(a=ipars[,"a"], b=ipars[,"b"], c=ipars[,"c"], theta=ppars[,"th"][idx_nab_p], D=1.702) 
  Ab_response[idx_nab_p,] <- temp$response
  Ab_prob[idx_nab_p,] <- temp$prob
  
  ### Response times
  temp <- genLNRTs(alp=ipars[,"alp"], bet=ipars[,"bet"], tau=ppars[,"tau"]) # response times for all examinees for all items
  Ab_RT <- temp$RT
  Ab_RT[idx_ab_p, c.pt:Nitem] <-  exp(log(Ab_RT[idx_ab_p, c.pt:Nitem]) - apply(log(temp$RT[,c.pt:Nitem]), 2, sd))
  
  #barplot(Ab_RT[idx_ab_p[2],])
  return(list(Ab_response=Ab_response, Ab_RT=Ab_RT, 
              idx_ab_p=idx_ab_p, idx_nab_p=idx_nab_p))
}


### Generating aberrant data (specify aberrant items)
genAbData2 <- function(p_ab_p, ab_item_idx, delta_th, delta_tau, ipars, ppars){ # ab_item_idx = seq(21,30)
  
  Nexaminee <- dim(ppars)[1]
  Nitem <- dim(ipars)[1]
  
  idx_ab_p <- sort(sample(x = 1:Nexaminee, size = Nexaminee*p_ab_p, replace = FALSE)) # View(idx_ab_p)
  idx_nab_p <- setdiff(seq(1:Nexaminee), idx_ab_p)
  
  Ab_response <- matrix(NA, Nexaminee, Nitem)
  Ab_RT <- matrix(NA, Nexaminee, Nitem)
  Ab_prob <- matrix(NA, Nexaminee, Nitem)
  
  nab_item_idx <- setdiff(seq(1,Nitem), ab_item_idx)
  
  ### Fill in for aberrant examinees
  temp <- genResp3PL(a=ipars[,"a"][ab_item_idx], b=ipars[,"b"][ab_item_idx], c=ipars[,"c"][ab_item_idx], theta=ppars[,"th"][idx_ab_p], D=1.702) 
  Ab_response[idx_ab_p,ab_item_idx] <- temp$response
  Ab_prob[idx_ab_p,ab_item_idx] <- temp$prob
  temp <- genResp3PL(a=ipars[,"a"][nab_item_idx], b=ipars[,"b"][nab_item_idx], c=ipars[,"c"][nab_item_idx], theta=ppars[,"th"][idx_ab_p]+delta_th, D=1.702) 
  Ab_response[idx_ab_p,nab_item_idx] <- temp$response
  Ab_prob[idx_ab_p,nab_item_idx] <- temp$prob
  
  temp <- genLNRTs(alp=ipars[,"alp"][ab_item_idx], bet=ipars[,"bet"][ab_item_idx], tau=ppars[,"tau"][idx_ab_p])
  Ab_RT[idx_ab_p,ab_item_idx] <- temp$RT
  temp <- genLNRTs(alp=ipars[,"alp"][nab_item_idx], bet=ipars[,"bet"][nab_item_idx], tau=ppars[,"tau"][idx_ab_p]+delta_tau)
  Ab_RT[idx_ab_p,nab_item_idx] <- temp$RT
  
  ### Fill in for non-aberrant examinees
  temp <- genResp3PL(a=ipars[,"a"], b=ipars[,"b"], c=ipars[,"c"], theta=ppars[,"th"][idx_nab_p], D=1.702) 
  Ab_response[idx_nab_p,] <- temp$response
  Ab_prob[idx_nab_p,] <- temp$prob
  
  temp <- genLNRTs(alp=ipars[,"alp"], bet=ipars[,"bet"], tau=ppars[,"tau"][idx_nab_p])
  Ab_RT[idx_nab_p,] <- temp$RT
  
  return(list(Ab_response=Ab_response, Ab_RT=Ab_RT, Ab_prob=Ab_prob,
              idx_ab_p=idx_ab_p, idx_nab_p=idx_nab_p))
}

### Generating aberrant data according to linear change of theta/tau
genAbData_Linear <- function(p_ab_p, c.pt, delta_th, delta_tau, ipars, ppars){ # ab_item_idx = seq(21,30)
  
  Nexaminee <- dim(ppars)[1]
  Nitem <- dim(ipars)[1]
  abNitem <- Nitem-c.pt+1
  
  idx_ab_p <- sort(sample(x = 1:Nexaminee, size = Nexaminee*p_ab_p, replace = FALSE)) # View(idx_ab_p)
  idx_nab_p <- setdiff(seq(1:Nexaminee), idx_ab_p)
  
  Ab_response <- matrix(NA, Nexaminee, Nitem)
  Ab_RT <- matrix(NA, Nexaminee, Nitem)
  
  ### Specify gradual theta/tau change
  
  delta_th_itrm <- seq(0, delta_th, length.out = abNitem+1)
  delta_th_itrm <- delta_th_itrm[2:(abNitem+1)]
  delta_tau_itrm <- seq(0, delta_tau, length.out = abNitem+1)
  delta_tau_itrm <- delta_tau_itrm[2:(abNitem+1)]
  
  th_mat <- matrix(ppars[,"th"], nrow=Nexaminee, ncol=Nitem)
  tau_mat <- matrix(ppars[,"tau"], nrow=Nexaminee, ncol=Nitem)
        
  ### Specify th and tau mat
  for (i in 1:Nexaminee){ # i=idx_ab_p[1]
    
    if (i %in% idx_ab_p){
      print(i)
      for (j in 1:abNitem){ # j=1
        th_mat[i,c.pt+j-1] <- th_mat[i,c.pt+j-1]+delta_th_itrm[j]
        tau_mat[i,c.pt+j-1] <- tau_mat[i,c.pt+j-1]+delta_tau_itrm[j]
      }
    }
    
  } # end of i
  
  ### Generate aberrant data based on the updated th/tau mat
  for (j in 1:Nitem){
    tmp <- genResp3PL(a=ipars[,"a"][j], b=ipars[,"b"][j], c=ipars[,"c"][j], 
                       theta=th_mat[,j], D=1.702) 
    Ab_response[,j] <- tmp$response
    tmp <- genLNRTs(alp=ipars[,"alp"][j], bet=ipars[,"bet"][j], 
                          tau=tau_mat[,j])
    Ab_RT[,j] <- tmp$RT
  }

  # checking data generation
  # mod_mbme <- mmap_cfix(Ab_response[idx_nab_p,], Ab_RT[idx_nab_p,],
  #                       ppar_prior=list(mu_p=null_Data$mu_p, cov_p=null_Data$cov_p),
  #                       iparst_prior = list(mu_is=null_Data$mu_item_str, cov_is=null_Data$cov_item_str),
  #                       cfix=c_f, tol=list(max_em=1000, tol_em=0.005, loglike=1e-4, max_nr=500, tol_nr=0.005))
  # 
  # ir_est_mat <- cbind(mod_mbme$iest[,"a"], mod_mbme$iest[,"b"], mod_mbme$iest[,"c"]); colnames(ir_est_mat) <- c("a", "b", "c")
  # rt_est_mat <- cbind(mod_mbme$iest[,"alp"], mod_mbme$iest[,"bet"]); colnames(rt_est_mat) <- c("alpha", "beta")
  # 
  # plot(null_Data$ipars[,"a"], mod_mbme$iest[,"a"])
  # plot(null_Data$ipars[,"b"], mod_mbme$iest[,"b"])
  # plot(null_Data$ipars[,"alp"], mod_mbme$iest[,"alp"])
  # plot(null_Data$ipars[,"bet"], mod_mbme$iest[,"bet"])
  # 
  # pest <- map(resp=Ab_response[idx_nab_p,], RT=Ab_RT[idx_nab_p,], 
  #             ipar = cbind(ir_est_mat, rt_est_mat),
  #             SE=FALSE, D=1.702)
  # 
  # pest <- map(resp=Ab_response[idx_ab_p,], RT=Ab_RT[idx_ab_p,], 
  #             ipar = cbind(ir_est_mat, rt_est_mat),
  #             SE=FALSE, D=1.702)
  # 
  # plot(th_mat[idx_ab_p,1], pest$est[,"th"]); abline(0,1)
  # plot(tau_mat[idx_ab_p,1], pest$est[,"tau"]); abline(0,1)
  # 

  return(list(Ab_response=Ab_response, Ab_RT=Ab_RT, 
              idx_ab_p=idx_ab_p, idx_nab_p=idx_nab_p))
  
}


### Calculate CUSUM charting statistics-------
cusumUp <- function(dev_RT, ref_v = 0){
  nexaminee <- length(dev_RT[,1])
  nitem <- length(dev_RT[1,])
  cusum_up <- matrix(NA, nexaminee, nitem)
  for (j in 1:nexaminee){
    cusum <- numeric(1)
    for (i in 1:nitem ){
      cusum <- cusum + dev_RT[j,i] - ref_v
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


### Calculate CVs--------------------
calc_CV <- function(pct, stats){
  CV <- apply(stats, 1, max)
  CV <- sort.int(CV, decreasing = TRUE, index.return=T)
  CV <- CV$x[floor(nrow(stats)*pct)]
  return(CV)
}


### Calculate densities for an RT vector
rt_dfn <- function(rt, Alp, Bet, Tau){ 
  n_dt_p <- length(rt)  
  density <- vector("numeric", length = n_dt_p)
  for (i in 1:n_dt_p){ # i=1
    density[i] <- (Alp[i]/(rt[i]*sqrt(2*pi)))*exp(-0.5*((Alp[i]*(log(rt[i])-(Bet[i]-Tau)))^2))
  } # dlnorm(rt, Bet-Tau,1/Alp)
  return(density)
}

### Calculate PF cases
calc_pfcase <- function (nexaminee, trueAb, estdAb){
  pfcase <- matrix(NA, nexaminee, 2, dimnames=list(NULL, c("true", "estd")))
  pfcase[trueAb,"true"] <- 1; pfcase[-trueAb,"true"] <- 0
  pfcase[estdAb, "estd"] <- 2; pfcase[-estdAb, "estd"] <- 0
  pfcaseSum <- rowSums(pfcase)
  pfprop <- rep(0, 4); names(pfprop) <- c("TN", "FN", "TypeI", "Power")  
  for (c in 0:3){ # c=0
    pfprop[c+1] <- sum(1*(pfcaseSum==c)) 
  }
  power <- pfprop["Power"]/length(trueAb)
  typeI <- pfprop["TypeI"]/(nexaminee-length(trueAb))
  
  return(list(power=power, typeI=typeI))
}

### Item response function vectorized
irf_3pl_v <- function(th, xi_ir, D = 1.702){
  a <- xi_ir[,1]
  b <- xi_ir[,2]
  c <- xi_ir[,3]
  p <- c + (1 - c)/(1 + exp(-D * (a * (th - b))))
  return(p)
}
