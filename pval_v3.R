
##code for pvalue for toroidal shift and variance correction with different kernel choices
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/stats.R", echo = FALSE)
wts.kernel <- function(x, bw, type = "epanechnikov"){
  if(type == "epanechnikov"){
    wts <- 3/4 *(1/bw)* (1 - (x/bw)^2) * (abs(x) <= bw)
  } else if(type == "gaussian"){
    wts <- dnorm(x, mean = 0, sd = bw)
  } else if(type == "uniform"){
    wts <- 0.5 * (abs(x) <= bw)
  } else {
    stop("Unsupported kernel type. Choose 'epanechnikov', 'gaussian', or 'uniform'.")
  }
  return(wts)
}
trimmed_mean_by_reference <- function(Kcross.simulated, nsimu, trim = 0.05) {
  # Kcross.simulated: rlen × (n.perm + 1) matrix (e.g., 50 × 200)
  #                   Each row corresponds to a distance r
  #                   Each column corresponds to a simulation (last column is observed)
  # nsimu: 2 × (n.perm + 1) matrix
  #        Row 1: counts of base taxa for each simulation
  #        Row 2: counts of shift taxa for each simulation
  # trim: proportion to trim from each tail of each row of nsimu
  #
  # Returns: rlen × 1 vector of trimmed means (one per distance r)
  
  # Check dimensions
  if (ncol(Kcross.simulated) != ncol(nsimu)) {
    stop("Kcross.simulated and nsimu must have the same number of columns")
  }
  
  n_sims <- ncol(Kcross.simulated)  # n.perm + 1
  n_distances <- nrow(Kcross.simulated)  # rlen (e.g., 50)
  
  # Calculate how many simulations to trim from each tail
  n_trim <- floor(n_sims * trim)
  
  # If trim proportion is too large, issue warning
  if (n_trim * 4 >= n_sims) {
    warning(paste0("Trim proportion may be too large. Trimming ", 
                   4 * n_trim, " out of ", n_sims, " simulations."))
  }
  
  # Find extreme indices based on nsimu
  # For each row of nsimu, identify the indices of extreme values
  extreme_indices <- unique(c(
    order(nsimu[1, ])[1:n_trim],                    # lowest n1 values
    order(nsimu[1, ], decreasing = TRUE)[1:n_trim], # highest n1 values
    order(nsimu[2, ])[1:n_trim],                    # lowest n2 values
    order(nsimu[2, ], decreasing = TRUE)[1:n_trim]  # highest n2 values
  ))
  
  # Keep only non-extreme simulations
  keep_indices <- setdiff(1:n_sims, extreme_indices)
  
  if (length(keep_indices) == 0) {
    stop("All simulations were trimmed. Reduce trim proportion.")
  }
  
  # Calculate row means for the kept simulations
  # Result is an rlen × 1 vector
  Kmean.trim <- rowMeans(Kcross.simulated[, keep_indices, drop = FALSE])
  
  return(Kmean.trim)
}

nw_predict <- function(X, Y, X0 = X, bw) {
  n <- nrow(X)
  m <- nrow(X0)
  
  # Handle Y as either matrix or vector
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  p <- ncol(Y)     # number of response variables
  out <- matrix(0, nrow = m, ncol = p)
  
  # Gaussian kernel function (isotropic)
  gaussian_kernel <- function(d2, bw) {
    exp(-0.5 * d2 / (bw^2))
  }
  
  for (j in seq_len(m)) {
    # squared Euclidean distances from query point X0[j,] to all X
    d2 <- rowSums((X - matrix(X0[j, ], nrow = n, ncol = ncol(X), byrow = TRUE))^2)
    w <- gaussian_kernel(d2, bw)
    sw <- sum(w)
    if (sw > 0) {
      out[j, ] <- crossprod(w, Y) / sw   # weighted average of the rows of Y
    } else {
      out[j, ] <- NA
    }
  }
  
  # Return as vector if p = 1
  if (p == 1) {
    out <- as.vector(out)
  }
  
  out
}

pval.assoc <- function(data, base.taxa = 1, shift.taxa = 2, r = r, n.perm = 199, bw = "silverman", ind = 25){
  rlen = length(r)
  freq_marks <- table(data$marks)
  original_window <- data$window
  x_range <- diff(original_window$xrange)
  y_range <- diff(original_window$yrange)
  n1 <- freq_marks[as.character(base.taxa)]
  n2 <- freq_marks[as.character(shift.taxa)]
  
  ##Compute stats##
  
  obj_K <- compute_K(data, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = n1, lambda2 = n2, r = r)
  
  Kcross_stat <- obj_K$Kcross_stat
  Kstar <- obj_K$Kstar
  Kcor <- obj_K$Kcor[-1]
  rcor <- r[-1]
  
  auc_Kcross <- trapz(r, Kcross_stat)
  med_Kcross <- Kcross_stat[ind]
  
  
  
  ##rshift data and compute stats##
  Kcross_stat.tor <- matrix(NA, nrow = n.perm, ncol = rlen)
  Kcross_stat.vc <- matrix(NA, nrow = n.perm, ncol = rlen)
  
  Kstar.tor <- matrix(NA, nrow = n.perm, ncol = rlen)
  Kstar.vc <- matrix(NA, nrow = n.perm, ncol = rlen)
  
  
  Kcor.tor <- matrix(NA, nrow = n.perm, ncol = rlen-1)
  Kcor.vc <- matrix(NA, nrow = n.perm, ncol = rlen-1)
  
  
  auc.tor <- numeric(n.perm)
  auc.vc <- numeric(n.perm)
  medr.tor <- numeric(n.perm)
  medr.vc <- numeric(n.perm)
  
  
  v_perm <- matrix(nrow = 2, ncol = n.perm)
  nstar_reduced <- matrix(NA, nrow = 2, ncol = n.perm)
  perm <- 1
  
  data_tor <- data
  ##check if window is rectangular
  if(!is.rectangle(data$win)){
    Window(data_tor) <- boundingbox(data$win)
  }
  
  
  while(perm <= n.perm){
    #toroidal shift
    data.shifted <- rshift(data_tor, which = as.character(shift.taxa), edge = "torus")
    freq_marks <- table(data.shifted$marks)
    lambda1 <- freq_marks[as.character(base.taxa)]
    lambda2 <- freq_marks[as.character(shift.taxa)]
    
    if(lambda1 == 0 || lambda2 == 0){
      next  # Skip this iteration, don't increment perm
    }
    
    obj_K_shifted <- compute_K(data.shifted, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = lambda1, lambda2 = lambda2, r = r)
    Kcross_stat.tor[perm, ] <- obj_K_shifted$Kcross_stat
    Kstar.tor[perm, ] <- obj_K_shifted$Kstar
    Kcor.tor[perm, ] <- obj_K_shifted$Kcor[-1]
    auc.tor[perm] <- trapz(r, obj_K_shifted$Kcross_stat)
    medr.tor[perm] <- obj_K_shifted$Kcross_stat[ind]
    
    #non-toroidal shift
    jump_rad <- incircle(data$window)$r
    shift_vector <- runifdisc(1, radius = jump_rad)
    shift_x <- shift_vector$x
    shift_y <- shift_vector$y
    points_to_shift <- subset(data, marks == shift.taxa)
    other_points <- subset(data, marks != shift.taxa)
    shifted_points <- shift(points_to_shift, vec = c(shift_x, shift_y))
    W_shifted <- shift(original_window, vec = c(shift_x, shift_y))
    W.reduced <- intersect.owin(original_window, W_shifted)
    if(is.null(W.reduced) || area.owin(W.reduced) == 0){
      next
    }
    pp_original_reduced <- data[W.reduced]
    pp_shifted_full <- superimpose(other_points, shifted_points)
    pp_shifted_reduced <- pp_shifted_full[W.reduced]
    freq_marks <- table(pp_shifted_reduced$marks)
    lambda1_new <- freq_marks[as.character(base.taxa)]
    lambda2_new <- freq_marks[as.character(shift.taxa)]
    
    if(lambda1_new == 0 || lambda2_new == 0){
      next
    }
    obj_K_minus <- compute_K(pp_shifted_reduced, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = lambda1_new, lambda2 = lambda2_new, r = r)
    Kcross_stat.vc[perm, ] <- obj_K_minus$Kcross_stat
    Kstar.vc[perm, ] <- obj_K_minus$Kstar
    Kcor.vc[perm, ] <- obj_K_minus$Kcor[-1]
    auc.vc[perm] <- trapz(r, obj_K_minus$Kcross_stat)
    medr.vc[perm] <- obj_K_minus$Kcross_stat[ind]
    v_perm[,perm] <- c(shift_x, shift_y)
    pp_prob <- table(pp_shifted_reduced$marks)
    nstar_reduced[,perm] <- c(pp_prob[base.taxa] , pp_prob[shift.taxa])
    # Only increment counter if we pass all checks
    perm <- perm + 1
  }
  
  #pval for auc.tor and medr.tor
  pval_auc <- (1+ sum(auc.tor >= auc_Kcross))/(1+n.perm)
  pval_medr <- (1+ sum(medr.tor >= med_Kcross))/(1+n.perm)
  
  #Global envelop test for Kcross, Kstar and Kcor
  df <- data.frame(r = r, Kcross_stat = Kcross_stat, t(Kcross_stat.tor))
  clean_indices <- apply(df,1, function(x)all(is.finite(x)))
  df <- df[clean_indices,]
  if(nrow(df)==0){
    pval_Kcross <- NA
  } else {
    r1 <- df$r
    obs <- df$Kcross_stat
    sim <- as.matrix(df[,-c(1,2)])
    cset <- curve_set(obs = obs, sim = sim, r = r1)
    GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank",alpha = 0.05, alternative = "greater")
    pval_Kcross <- attr(GET_erl, "p")
  }
  
  df_star <- data.frame(r = r, Kstar = Kstar, t(Kstar.tor))
  clean_indices_star <- apply(df_star,1, function(x)all(is.finite(x)))
  df_star <- df_star[clean_indices_star,]
  if(nrow(df_star)==0){
    pval_Kstar <- NA
  } else {
    r1_star <- df_star$r
    obs_star <- df_star$Kstar
    sim_star <- as.matrix(df_star[,-c(1,2)])
    cset_star <- curve_set(obs = obs_star, sim = sim_star, r = r1_star)
    GET_erl_star <- global_envelope_test(cset_star, typeone = "fwer", type = "rank",alpha = 0.05, alternative = "greater")
    pval_Kstar <- attr(GET_erl_star, "p")
  }
  
  df_cor <- data.frame(r = rcor, Kcor = Kcor, t(Kcor.tor))
  clean_indices_cor <- apply(df_cor,1, function(x)all(is.finite(x)))
  df_cor <- df_cor[clean_indices_cor,]
  if(nrow(df_cor)==0){
    pval_Kcor <- NA
  } else {
    r1_cor <- df_cor$r
    obs_cor <- df_cor$Kcor
    sim_cor <- as.matrix(df_cor[,-c(1,2)])
    cset_cor <- curve_set(obs = obs_cor, sim = sim_cor, r = r1_cor)
    GET_erl_cor <- global_envelope_test(cset_cor, typeone = "fwer", type = "rank",alpha = 0.05, alternative = "greater")
    pval_Kcor <- attr(GET_erl_cor, "p")
  }
  
 
##variance correction. weights from three different kernels.
  vsimu <- cbind(v_perm, c(0,0))
  nsimu <- cbind(nstar_reduced, c(n1, n2))
  #aij <- pairdist.default(t(vsimu))
  #nij <- pairdist.default(t(nsimu))
  if(bw == "silverman" || bw == "scott"){
    n_sims <- ncol(vsimu)
    d <- nrow(vsimu)
    
    if(bw == "silverman"){
      # Silverman's rule
      sd_v <- apply(vsimu, 1, sd)
      sigma_v <- exp(mean(log(sd_v)))
      bw <- sigma_v * (n_sims * (d + 2) / 4)^(-1/(d+4))
      
      sd_n <- apply(nsimu, 1, sd)
      sigma_n <- exp(mean(log(sd_n)))
      bw_n <- sigma_n * (n_sims * (d + 2) / 4)^(-1/(d+4))
    } else {
      # Scott's rule: often better for multivariate data
      sd_v <- apply(vsimu, 1, sd)
      sigma_v <- exp(mean(log(sd_v)))
      bw <- sigma_v * n_sims^(-1/(d+4))
      
      sd_n <- apply(nsimu, 1, sd)
      sigma_n <- exp(mean(log(sd_n)))
      bw_n <- sigma_n * n_sims^(-1/(d+4))
    }
  }
 
  
 
  
  #Kcross - vc, evc and minus
  Kcross.simulated <- cbind(t(Kcross_stat.vc), Kcross_stat)
  Kmean <- apply(Kcross.simulated, 1, mean)
  Kmean1 <- trimmed_mean_by_reference(Kcross.simulated, nsimu, trim = 0.05)
  Si <- sweep(Kcross.simulated, 1, Kmean, "-")
  si2_gauss <- t(nw_predict(X = t(vsimu), Y = t(Si^2), X0 = t(vsimu), bw = bw))
  #S_gauss <- (si2_gauss)^(-0.5) * Si
  si2_gauss[si2_gauss <= 1e-10] <- 1e-10  # Floor the variance
  S_gauss <- Si / sqrt(si2_gauss)
  S_gauss[!is.finite(S_gauss)] <- 0 
  CS_gauss <- create_curve_set(list(r = r, obs = S_gauss[ , n.perm + 1], sim_m = S_gauss[ , 1:n.perm]))
  pval.Kcross.vc.gauss <- attr(rank_envelope(CS_gauss, type = "erl"), "p")
  
  si2_n <- t(nw_predict(X = t(nsimu), Y = t(Si^2), X0 = t(nsimu), bw = bw_n))
  S_n <- (si2_n)^(-0.5) * Si
  S_n[is.nan(S_n)] <- 0
  CS_n <- create_curve_set(list(r = r, obs = S_n[ , n.perm + 1], sim_m = S_n[ , 1:n.perm]))
  pval.Kcross.vc.n <- attr(rank_envelope(CS_n, type = "erl"), "p")
  
 
  EKcross <- t(nw_predict(X = t(vsimu), Y = t(Kcross.simulated), X0 = t(vsimu), bw = bw))
  E2Kcross <- t(nw_predict(X = t(vsimu), Y = t(Kcross.simulated^2), X0 = t(vsimu), bw = bw))
  var_Kcross <- E2Kcross - EKcross^2
  S_Kcross <- (Kcross.simulated - EKcross) / sqrt(var_Kcross)
  S_Kcross[is.nan(S_Kcross)] <- 0
  CS <- create_curve_set(list(r = r, obs = S_Kcross[ , n.perm + 1], sim_m = S_Kcross[ , 1:n.perm]))
  pval.Kcross.evc <- attr(rank_envelope(CS, type = "erl"), "p")
  
  EKcross_n <- t(nw_predict(X = t(nsimu), Y = t(Kcross.simulated), X0 = t(nsimu), bw = bw_n))
  E2Kcross_n <- t(nw_predict(X = t(nsimu), Y = t(Kcross.simulated^2), X0 = t(nsimu), bw = bw_n))
  var_Kcross_n <- E2Kcross_n - EKcross_n^2
  S_Kcross_n <- (Kcross.simulated - EKcross) / sqrt(var_Kcross_n)
  S_Kcross_n[is.nan(S_Kcross_n)] <- 0
  CS_n <- create_curve_set(list(r = r, obs = S_Kcross_n[ , n.perm + 1], sim_m = S_Kcross_n[ , 1:n.perm]))
  pval.Kcross.evc.n <- attr(rank_envelope(CS_n, type = "erl"), "p")
  
  CS_minus <- create_curve_set(list(r = r, obs = Kcross.simulated[,n.perm+1], sim_m = Kcross.simulated[ , 1:n.perm]))
  pval.Kcross.minus <- attr(rank_envelope(CS_minus, type = "erl"), "p")
  
  ##AUC and medr - vc and minus correction
  auc.simulated.vc <- c(auc.vc, auc_Kcross)
  aucmean.vc <- mean(auc.simulated.vc)
  Si_auc.vc <- auc.simulated.vc - aucmean.vc
  Si2_auc.vc <- nw_predict(X = t(vsimu), Y = Si_auc.vc^2, X0 = t(vsimu), bw = bw)
  S_auc.vc <- Si_auc.vc / sqrt(Si2_auc.vc)
  pval.auc.vc <- (1 + sum(S_auc.vc[1:n.perm] >= S_auc.vc[n.perm + 1]))/(1+n.perm)
  
  medr.simulated.vc <- c(medr.vc, med_Kcross)
  medrmean.vc <- mean(medr.simulated.vc)
  Si_medr.vc <- medr.simulated.vc - medrmean.vc
  Si2_medr.vc <- nw_predict(X = t(vsimu), Y = Si_medr.vc^2, X0 = t(vsimu), bw = bw)
  S_medr.vc <- Si_medr.vc / sqrt(Si2_medr.vc)
  pval.medr.vc <- (1 + sum(S_medr.vc[1:n.perm] >= S_medr.vc[n.perm + 1]))/(1+n.perm)
  
  pval.auc.minus <- (1 + sum(auc.vc >= auc_Kcross))/(1 + n.perm)
  pval.medr.minus <- (1 + sum(medr.vc >= med_Kcross))/(1 + n.perm)
  
  return(list(pval_Kcross = pval_Kcross, pval_Kstar = pval_Kstar, pval_Kcor = pval_Kcor, 
             pval.Kcross.vc.gauss = pval.Kcross.vc.gauss, pval.Kcross.evc = pval.Kcross.evc, 
             pval.Kcross.vc.n = pval.Kcross.vc.n, pval.Kcross.evc.n = pval.Kcross.evc.n,
             pval.Kcross.minus = pval.Kcross.minus, pval.auc.vc = pval.auc.vc, pval.medr.vc = pval.medr.vc, 
             pval.auc.tor = pval_auc, pval.medr.tor = pval_medr, pval.auc.minus = pval.auc.minus, pval.medr.minus = pval.medr.minus))
}


#plot(nsimu[2,], (Si^2)[ind,])
#plot(nsimu[1,], (Si^2)[ind,])
#plot(vsimu[1,], (Si^2)[ind,])
#plot(vsimu[2,], (Si^2)[ind,])
