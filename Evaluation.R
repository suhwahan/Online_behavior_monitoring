### Functions to conduct evaluations ----

###   Calculate power/type I error for SGLRT -------------------
calcPw_SGLRT <- function (Proc,Res, CV, e_dt_pt_excl, sampling, c.pt, st_eval, idx_ab_p){
  
  neval <- length(Res[1,])
  Nexaminee <- length(Res[,1])
  dt_cht <- abs(Res) > abs(CV) # View(dt_cht) # LR has been set to be higher when the statistic favors H1
  dt_pt <- matrix(NA, Nexaminee, 2, dimnames=list(NULL, c("idx", "dt.pt_eval")))  # matrix containing idx of people who were "detected" anywhere along the items (power+type1) and their detection point; 1st column=index; 2nd column=detection point  # View(dt_pt)
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
  if (sampling == "exclusive"){
    e_dt_pt <- e_dt_pt_excl
  } else {
    e_dt_pt <- c.pt-st_eval+1
  }
  if (c.pt >= st_eval){
    e_dt_mat <- pw_tent_mat[which(pw_tent_mat[,"dt.pt_eval"] < e_dt_pt),, drop = F]  # index and dt point of early detected t1 error
    pw_mat <- pw_tent_mat[which(pw_tent_mat[,"dt.pt_eval"] >= e_dt_pt),, drop = F]  # index and dt point of true power ppl
  } else {
    e_dt_mat <- NULL
    pw_mat <- pw_tent_mat
  }
  t1_mat <- rbind(f_dt_mat, e_dt_mat)
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