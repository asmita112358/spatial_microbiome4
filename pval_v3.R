
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
nw_predict <- function(X, Y, X0 = X, bw) {
  n <- nrow(X)
  m <- nrow(X0)
  p <- ncol(Y)     # 50
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
  out
}

#wts.kernel(x = c(-0.1, 0, 0.12), bw = 0.15, type = "epanechnikov")
#epanechnikov(x = c(-0.1, 0, 0.12), bw = 0.15)


pval.assoc <- function(data, base.taxa = 1, shift.taxa = 2, r = r, n.perm = 199, bw = "silverman"){
  rlen = length(r)
  freq_marks <- table(data$marks)
  original_window <- data$window
  x_range <- diff(original_window$xrange)
  y_range <- diff(original_window$yrange)
  lambda1 <- freq_marks[as.character(base.taxa)]
  lambda2 <- freq_marks[as.character(shift.taxa)]
  
  ##Compute stats##
  
  obj_K <- compute_K(data, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = lambda1, lambda2 = lambda2, r = r)
  
  Kcross_stat <- obj_K$Kcross_stat
  Kstar <- obj_K$Kstar
  Kcor <- obj_K$Kcor[-1]
  rcor <- r[-1]
  
  Nij <- compute_NN(data, base.taxa, shift.taxa)
  
  ##rshift data and compute stats##
  Kcross_stat.tor <- matrix(NA, nrow = n.perm, ncol = rlen)
  Kcross_stat.vc <- matrix(NA, nrow = n.perm, ncol = rlen)
  
  
  Kstar.tor <- matrix(NA, nrow = n.perm, ncol = rlen)
  Kstar.vc <- matrix(NA, nrow = n.perm, ncol = rlen)
  
  
  Kcor.tor <- matrix(NA, nrow = n.perm, ncol = rlen-1)
  Kcor.vc <- matrix(NA, nrow = n.perm, ncol = rlen-1)
  
  
  Nij.tor <- c()
  Nij.vc <- c()
  
  
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
    Nij.tor[perm] <- compute_NN(data.shifted, base.taxa, shift.taxa)
    
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
    Nij.vc[perm] <- compute_NN(pp_shifted_reduced, base.taxa, shift.taxa)
    v_perm[,perm] <- c(shift_x, shift_y)
    pp_prob <- table(pp_shifted_reduced$marks)
    nstar_reduced[,perm] <- c(pp_prob[base.taxa] , pp_prob[shift.taxa])
    # Only increment counter if we pass all checks
    perm <- perm + 1
  }
  
  
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
  
  ## pvalue for NN
  
  pval_NN <- (1 + sum(Nij.tor >= Nij))/(1 + n.perm)
  
  
  ##variance correction. weights from three different kernels.
  vsimu <- cbind(v_perm, c(0,0))
  nsimu <- cbind(nstar_reduced, c(lambda1, lambda2))
  aij <- pairdist.default(t(vsimu))
  nij <- pairdist.default(t(nsimu))
  if(bw == "silverman"){
    bw <- 1.06 * mean(apply(vsimu, 1, sd)) * nrow(vsimu)^(-1/5)
    bw_n <- 1.06 * mean(apply(nsimu, 1, sd)) * nrow(nsimu)^(-1/5)
  }
  #wts_epanechnikov <- wts.kernel(aij, bw = bw, type = "epanechnikov")
  #wts_epanechnikov <- wts_epanechnikov/sum(wts_epanechnikov)
  #wts_gaussian <- wts.kernel(aij, bw = bw, type = "gaussian")
  #wts_gaussian <- wts_gaussian/sum(wts_gaussian)
  #wts_uniform <- wts.kernel(aij, bw = bw, type = "uniform")
  #wts_uniform <- wts_uniform/sum(wts_uniform)
  
  ##Kcross, Kstar, Kcor and NN in combination with the different kernels
  
  #Kcross - vc, evc and minus
  Kcross.simulated <- cbind(t(Kcross_stat.vc), Kcross_stat)
  Kmean <- apply(Kcross.simulated, 1, mean)
  Si <- sweep(Kcross.simulated, 1, Kmean, "-")
  si2_gauss <- t(nw_predict(X = t(vsimu), Y = t(Si^2), X0 = t(vsimu), bw = bw))
  S_gauss <- (si2_gauss)^(-0.5) * Si
  S_gauss[is.nan(S_gauss)] <- 0
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
  S_Kcross_n <- (Kcross.simulated - EKcross_n) / sqrt(var_Kcross_n)
  S_Kcross_n[is.nan(S_Kcross_n)] <- 0
  CS_n <- create_curve_set(list(r = r, obs = S_Kcross_n[ , n.perm + 1], sim_m = S_Kcross_n[ , 1:n.perm]))
  pval.Kcross.evc.n <- attr(rank_envelope(CS_n, type = "erl"), "p")
  
  CS_minus <- create_curve_set(list(r = r, obs = Kcross.simulated[,n.perm+1], sim_m = Kcross.simulated[ , 1:n.perm]))
  pval.Kcross.minus <- attr(rank_envelope(CS_minus, type = "erl"), "p")
  
  return(list(pval_Kcross = pval_Kcross, pval_Kstar = pval_Kstar, pval_Kcor = pval_Kcor, pval_NN = pval_NN,
             pval.Kcross.vc.gauss = pval.Kcross.vc.gauss, pval.Kcross.evc = pval.Kcross.evc, 
             pval.Kcross.vc.n = pval.Kcross.vc.n, pval.Kcross.evc.n = pval.Kcross.evc.n,
             pval.Kcross.minus = pval.Kcross.minus))
  }
