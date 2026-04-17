
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
  aij <- pairdist.default(t(vsimu))
  if(bw == "silverman"){
    bw <- 1.06 * min(apply(vsimu, 2, sd)) * nrow(vsimu)^(-1/5)
  }
  wts_epanechnikov <- wts.kernel(aij, bw = bw, type = "epanechnikov")
  wts_epanechnikov <- wts_epanechnikov/sum(wts_epanechnikov)
  wts_gaussian <- wts.kernel(aij, bw = bw, type = "gaussian")
  wts_gaussian <- wts_gaussian/sum(wts_gaussian)
  wts_uniform <- wts.kernel(aij, bw = bw, type = "uniform")
  wts_uniform <- wts_uniform/sum(wts_uniform)
  
  ##Kcross, Kstar, Kcor and NN in combination with the different kernels
  
  #Kcross
  Kcross.simulated <- cbind(t(Kcross_stat.vc), Kcross_stat)
  Kmean <- apply(Kcross.simulated, 1, mean)
  Si <- sweep(Kcross.simulated, 1, Kmean, "-")
  si2_ep <- (Si^2)%*% wts_epanechnikov
  si2_gauss <- (Si^2)%*% wts_gaussian
  si2_uniform <- (Si^2)%*% wts_uniform
  
  S_ep <- (si2_ep)^(-0.5) * Si
  S_ep[is.nan(S_ep)] <- 0
  CS <- create_curve_set(list(r = r, obs = S_ep[ , n.perm + 1], sim_m = S_ep[ , 1:n.perm]))
  pval.Kcross.vc.ep <- attr(rank_envelope(CS, type = "erl"), "p")
  
  S_gauss <- (si2_gauss)^(-0.5) * Si
  S_gauss[is.nan(S_gauss)] <- 0
  CS_gauss <- create_curve_set(list(r = r, obs = S_gauss[ , n.perm + 1], sim_m = S_gauss[ , 1:n.perm]))
  pval.Kcross.vc.gauss <- attr(rank_envelope(CS_gauss, type = "erl"), "p")
  
  S_uniform <- (si2_uniform)^(-0.5) * Si
  S_uniform[is.nan(S_uniform)] <- 0
  CS_uniform <- create_curve_set(list(r = r, obs = S_uniform[ , n.perm + 1], sim_m = S_uniform[ , 1:n.perm]))
  pval.Kcross.vc.uniform <- attr(rank_envelope(CS_uniform, type = "erl"), "p")
  
  #Kstar
  Kstar.simulated <- cbind(t(Kstar.vc), Kstar)
  Kstar.mean <- apply(Kstar.simulated, 1, mean)
  Si_star <- sweep(Kstar.simulated, 1, Kstar.mean, "-")
  si2_star_ep <- (Si_star^2)%*% wts_epanechnikov
  si2_star_gauss <- (Si_star^2)%*% wts_gaussian
  si2_star_uniform <- (Si_star^2)%*% wts_uniform
  
  S_star_ep <- (si2_star_ep)^(-0.5) * Si_star
  S_star_ep[is.nan(S_star_ep)] <- 0
  CS_star_ep <- create_curve_set(list(r = r, obs = S_star_ep[ , n.perm + 1], sim_m = S_star_ep[ , 1:n.perm]))
  pval.Kstar.vc.ep <- attr(rank_envelope(CS_star_ep, type = "erl"), "p")
  
  S_star_gauss <- (si2_star_gauss)^(-0.5) * Si_star
  S_star_gauss[is.nan(S_star_gauss)] <- 0
  CS_star_gauss <- create_curve_set(list(r = r, obs = S_star_gauss[ , n.perm + 1], sim_m = S_star_gauss[ , 1:n.perm]))
  pval.Kstar.vc.gauss <- attr(rank_envelope(CS_star_gauss, type = "erl"), "p")
  
  S_star_uniform <- (si2_star_uniform)^(-0.5) * Si_star
  S_star_uniform[is.nan(S_star_uniform)] <- 0
  CS_star_uniform <- create_curve_set(list(r = r, obs = S_star_uniform[ , n.perm + 1], sim_m = S_star_uniform[ , 1:n.perm]))
  pval.Kstar.vc.uniform <- attr(rank_envelope(CS_star_uniform, type = "erl"), "p")
  
  #Kcor
  Kcor.simulated <- cbind(t(Kcor.vc), Kcor)
  Kcor.mean <- apply(Kcor.simulated, 1, mean)
  Si_cor <- sweep(Kcor.simulated, 1, Kcor.mean, "-")
  si2_cor_ep <- (Si_cor^2)%*% wts_epanechnikov
  si2_cor_gauss <- (Si_cor^2)%*% wts_gaussian
  si2_cor_uniform <- (Si_cor^2)%*% wts_uniform
  
  S_cor_ep <- (si2_cor_ep)^(-0.5) * Si_cor
  S_cor_ep[is.nan(S_cor_ep)] <- 0
  CS_cor_ep <- create_curve_set(list(r = rcor, obs = S_cor_ep[ , n.perm + 1], sim_m = S_cor_ep[ , 1:n.perm]))
  pval.Kcor.vc.ep <- attr(rank_envelope(CS_cor_ep, type = "erl"), "p")
  
  S_cor_gauss <- (si2_cor_gauss)^(-0.5) * Si_cor
  S_cor_gauss[is.nan(S_cor_gauss)] <- 0
  CS_cor_gauss <- create_curve_set(list(r = rcor, obs = S_cor_gauss[ , n.perm + 1], sim_m = S_cor_gauss[ , 1:n.perm]))
  pval.Kcor.vc.gauss <- attr(rank_envelope(CS_cor_gauss, type = "erl"), "p")
  
  S_cor_uniform <- (si2_cor_uniform)^(-0.5) * Si_cor
  S_cor_uniform[is.nan(S_cor_uniform)] <- 0
  CS_cor_uniform <- create_curve_set(list(r = rcor, obs = S_cor_uniform[ , n.perm + 1], sim_m = S_cor_uniform[ , 1:n.perm]))
  pval.Kcor.vc.uniform <- attr(rank_envelope(CS_cor_uniform, type = "erl"), "p")
  
  #NN
  Nij.simulated <- c(Nij.vc, Nij)
  Nij.mean <- mean(Nij.simulated)
  Si_NN <- Nij.simulated - Nij.mean
  
  si2_NN_ep <- sum(Si_NN^2 * wts_epanechnikov)  
  si2_NN_gauss <- sum(Si_NN^2 * wts_gaussian)
  si2_NN_uniform <- sum(Si_NN^2 * wts_uniform)
  
  S_NN_ep <- (si2_NN_ep)^(-0.5) * Si_NN
  S_NN_ep[is.nan(S_NN_ep)] <- 0
  pval.NN.vc.ep <- (1 + sum(S_NN_ep[1:n.perm] >= S_NN_ep[n.perm + 1]))/(1 + n.perm)

  S_NN_gauss <- (si2_NN_gauss)^(-0.5) * Si_NN
  S_NN_gauss[is.nan(S_NN_gauss)] <- 0
  pval.NN.vc.gauss <- (1 + sum(S_NN_gauss[1:n.perm] >= S_NN_gauss[n.perm + 1]))/(1 + n.perm)
  
  S_NN_uniform <- (si2_NN_uniform)^(-0.5) * Si_NN
  S_NN_uniform[is.nan(S_NN_uniform)] <- 0
  pval.NN.vc.uniform <- (1 + sum(S_NN_uniform[1:n.perm] >= S_NN_uniform[n.perm + 1]))/(1 + n.perm)
  
  return(list(pval_Kcross = pval_Kcross, pval_Kstar = pval_Kstar, pval_Kcor = pval_Kcor, pval_NN = pval_NN,
              pval.Kcross.vc.ep = pval.Kcross.vc.ep, pval.Kcross.vc.gauss = pval.Kcross.vc.gauss, pval.Kcross.vc.uniform = pval.Kcross.vc.uniform,
              pval.Kstar.vc.ep = pval.Kstar.vc.ep, pval.Kstar.vc.gauss = pval.Kstar.vc.gauss, pval.Kstar.vc.uniform = pval.Kstar.vc.uniform,
              pval.Kcor.vc.ep = pval.Kcor.vc.ep, pval.Kcor.vc.gauss = pval.Kcor.vc.gauss, pval.Kcor.vc.uniform = pval.Kcor.vc.uniform,
              pval.NN.vc.ep = pval.NN.vc.ep, pval.NN.vc.gauss = pval.NN.vc.gauss, pval.NN.vc.uniform = pval.NN.vc.uniform))
}
