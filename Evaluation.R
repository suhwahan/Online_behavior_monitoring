### Functions to conduct evaluations ----

###   Calculate power/type I error for SGLRT -------------------
calcPw_SGLRT <- function (Proc, Res, CV, e_dt_pt=NULL, sampling, c.pt, n_rf_s, idx_ab_p){
  
  neval <- length(Res[1,])
  Nexaminee <- length(Res[,1])
  
  # Find a cell which exceeded CV ------------
  dt_cht <- abs(Res) > abs(CV)  # LR has been set to be higher when the statistic favors H1
  
  # matrix containing idx of people who were "detected" anywhere along the items (power+type1) and their detection point; 
  # 1st column=index; 2nd column=detection point  
  
  dt_pt <- matrix(NA, Nexaminee, 2, dimnames=list(NULL, c("idx", "dt.pt_eval")))  
  
  k <- 0  # counting variable for examinee
  for (i in 1:Nexaminee){ 
    found <- FALSE
    j <- 1  # counting variable for item on the eval_point scale
    while (j <= neval && found == FALSE){
      if (dt_cht[i,j]){
        found <- TRUE
        k <- k+1
        dt_pt[k,"idx"] <- i
        dt_pt[k,"dt.pt_eval"] <- j
      }
      j <- j+1
    }
  }
  
  tmp <- which(is.na(dt_pt[,1]))  # row number for NAs in `dt_pt`
  dt_pt <- dt_pt[-tmp,]  # entire detected ppl
  pw_tent_mat <- dt_pt[which(dt_pt[,"idx"] %in% idx_ab_p),, drop = F] # Correctly detected ppl (including early detection)
  f_dt_mat <- dt_pt[-which(dt_pt[,"idx"] %in% idx_ab_p),, drop = F] # Incorrectly detected ppl (false detection t1 error)  # earliest detection point in the detection chart
  
  # Define premature detection points for each procedure
  if (Proc == "CUSUM_obs"){ # CUSUM control chart with manifest indicators
    st_eval = n_rf_s+1
    e_dt_pt = c.pt-st_eval+1
  } else if (Proc == "CUSUM_est" && !is.null(e_dt_pt)){
    st_eval = n_rf_s*2
    e_dt_pt = e_dt_pt
  } else if (Proc == "SGLRT" && sampling == "exclusive" && !is.null(e_dt_pt)){
    st_eval = n_rf_s*2
    e_dt_pt = e_dt_pt
  } else if (Proc == "SGLRT" && sampling == "ref_mov_fixed"){
    st_eval = n_rf_+1
    e_dt_pt = c.pt-st_eval+1
  }
  
  # Find prematurely detected cases --------
  if (c.pt >= st_eval){
    e_dt_mat <- pw_tent_mat[which(pw_tent_mat[,"dt.pt_eval"] < e_dt_pt),, drop = F]  # index and dt point of early detected t1 error
    pw_mat <- pw_tent_mat[which(pw_tent_mat[,"dt.pt_eval"] >= e_dt_pt),, drop = F]  # index and dt point of true power ppl
  } else {
    e_dt_mat <- NULL
    pw_mat <- pw_tent_mat
  }
  t1_mat <- rbind(f_dt_mat, e_dt_mat)
  
  # Power, Type I, early detection rates ---------
  e_dt_rat <- length(e_dt_mat[,"idx"])/length(pw_tent_mat[,"idx"])
  pw_rat <- length(pw_mat[,"idx"])/length(idx_ab_p)
  t1_rat <- length(t1_mat[,"idx"])/(Nexaminee-length(idx_ab_p))
  if (c.pt >= st_eval){
    ARL1 <- mean(pw_mat[,"dt.pt_eval"]-e_dt_pt) # run length from the earliest possible point
    ARL0 <- mean(f_dt_mat[,"dt.pt_eval"])  
  } else {
    ARL1 <- mean(pw_mat[,"dt.pt_eval"]+st_eval-c.pt) # has to add (st_eval-c.pt); assuming no detection happened before the eval point
    ARL0 <- mean(f_dt_mat[,"dt.pt_eval"]) # ARL0 is calculated from the onset of evaluation point
  }
  
  return (list(
    pw_rat=pw_rat, t1_rat=t1_rat, ARL1=ARL1,
    ARL0=ARL0, e_dt_rat=e_dt_rat, CV=CV, neval=neval
  )
  )
}

### Calculate CVs -----------------------
calc_CV <- function(pct, stats){
  CV <- apply(stats, 1, max)
  CV <- sort.int(CV, decreasing = TRUE, index.return=T)
  CV <- CV$x[floor(nrow(stats)*pct)]
  return(CV)
}

### Calculate PF cases -------------------
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