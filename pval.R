##computing p value using M1 = toroidal shift, M2 = minus correction, M3 = variance correction.


source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/stats.R", echo = FALSE)
bw <- 0.10

# Epanechnikov kernel used for determining the weights in the regression
epanechnikov <- function(x,bw){
  return(pmax(0,(1/bw)*3/4*(1-(x/bw)^2)))
}
pval.assoc <- function(data, base.taxa = 1, shift.taxa = 2,r = r, n.perm = 199){
  rlen = length(r)
  freq_marks <- table(data$marks)
  original_window <- data$window
  x_range <- diff(original_window$xrange)
  y_range <- diff(original_window$yrange)
  lambda1 <- freq_marks[as.character(base.taxa)]
  lambda2 <- freq_marks[as.character(shift.taxa)]
  obj_K <- compute_K(data, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = lambda1, lambda2 = lambda2, r = r)
  
  #Extract all observed stats from obj_K
  
  K1 <- obj_K$K1
  K2 <- obj_K$K2
  K3 <- obj_K$K3
  K4 <- obj_K$K4
  K5 <- obj_K$K5
  K6 <- obj_K$K6
  K7 <- obj_K$K7
  K8 <- obj_K$K8
  K9 <- obj_K$K9
  
  
  # Permutation test using toroidal shift, minus correction and variance correction
  K1M1.perm <- K1M2.perm <- K1M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K2M1.perm <- K2M2.perm <- K2M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K3M1.perm <- K3M2.perm <- K3M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K4M1.perm <- K4M2.perm <- K4M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K5M1.perm <- K5M2.perm <- K5M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K6M1.perm <- K6M2.perm <- K6M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K7M1.perm <- K7M2.perm <- K7M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K8M1.perm <- K8M2.perm <- K8M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  K9M1.perm <- K9M2.perm <- K9M3.perm <- matrix(NA, nrow = n.perm, ncol = rlen)
  
  v_perm <- matrix(nrow = 2, ncol = n.perm)
  perm <- 0  # Counter for valid permutations
  
  while(perm < n.perm){
    
    #toroidal shift
    data.shifted <- rshift(data, which = as.character(shift.taxa), edge = "torus")
    data.shifted <- rshift(data.shifted, which = as.character(base.taxa), edge = "torus")
    freq_marks <- table(data.shifted$marks)
    lambda1 <- freq_marks[as.character(base.taxa)]
    lambda2 <- freq_marks[as.character(shift.taxa)]
    
    if(lambda1 == 0 || lambda2 == 0){
      next  # Skip this iteration, don't increment perm
    }
    
    # Only increment counter if we pass the check
    perm <- perm + 1
    
    obj_K_shifted <- compute_K(data.shifted, base.taxa = base.taxa, shift.taxa = shift.taxa, 
                               lambda1 = lambda1, lambda2 = lambda2, r = r)
    K1M1.perm[perm, ] <- obj_K_shifted$K1
    K2M1.perm[perm, ] <- obj_K_shifted$K2
    K3M1.perm[perm, ] <- obj_K_shifted$K3
    K4M1.perm[perm, ] <- obj_K_shifted$K4
    K5M1.perm[perm, ] <- obj_K_shifted$K5
    K6M1.perm[perm, ] <- obj_K_shifted$K6
    K7M1.perm[perm, ] <- obj_K_shifted$K7
    K8M1.perm[perm, ] <- obj_K_shifted$K8
    K9M1.perm[perm, ] <- obj_K_shifted$K9
    
    
    #minus correction
    random_x_shift <- runif(1, -x_range/2, x_range/2)
    random_y_shift <- runif(1, -y_range/2, y_range/2)
    shifted.points <- shift(subset(data, marks == as.character(shift.taxa)),
                            vec = c(random_x_shift, random_y_shift))  
    common_window <- intersect.owin(data$window, shifted.points$window)
    shifted.points <- shifted.points[common_window]
    data.minus <- subset(data, marks != as.character(shift.taxa))
    data.minus <- data.minus[common_window]
    data.new <- superimpose(shifted.points, data.minus, W = common_window)
    freq_marks <- table(data.new$marks)
    lambda1_new <- freq_marks[as.character(base.taxa)]
    lambda2_new <- freq_marks[as.character(shift.taxa)]
    
    if(lambda1_new == 0 || lambda2_new == 0){
      perm <- perm - 1  # Decrement because this iteration failed
      next
    }
    
    obj_K_minus <- compute_K(data.new, base.taxa = base.taxa, shift.taxa = shift.taxa, 
                             lambda1 = lambda1_new, lambda2 = lambda2_new, r = r)
    K1M2.perm[perm, ] <- obj_K_minus$K1
    K2M2.perm[perm, ] <- obj_K_minus$K2
    K3M2.perm[perm, ] <- obj_K_minus$K3
    K4M2.perm[perm, ] <- obj_K_minus$K4
    K5M2.perm[perm, ] <- obj_K_minus$K5
    K6M2.perm[perm, ] <- obj_K_minus$K6
    K7M2.perm[perm, ] <- obj_K_minus$K7
    K8M2.perm[perm, ] <- obj_K_minus$K8
    K9M2.perm[perm, ] <- obj_K_minus$K9
    
    #variance correction
    v_perm[, perm] <- c(random_x_shift, random_y_shift)
  }
  ##Clean K1, M1 and run GET
  df <- data.frame(r = r, K1 = K1, t(K1M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  if(nrow(df) == 0){
    pvalK1M1 <- NA}else{
  r1 <- df$r
  K1M1 <- df$K1
  sim_K1 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K1M1, sim = sim_K1, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK1M1 <- attr(GET_erl, "p")
    }
  
  ##Clean K1, M2 and run GET
  df <- data.frame(r = r, K1 = K1, t(K1M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K1M2 <- df$K1
  sim_K1 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K1M2, sim = sim_K1, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK1M2 <- attr(GET_erl, "p")
    
  
  ##Clean K2, M1 and run GET
  df <- data.frame(r = r, K2 = K2, t(K2M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K2M1 <- df$K2
  sim_K2 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K2M1, sim = sim_K2, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK2M1 <- attr(GET_erl, "p")
  
  ##Clean K2, M2 and run GET
  df <- data.frame(r = r, K2 = K2, t(K2M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K2M2 <- df$K2
  sim_K2 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K2M2, sim = sim_K2, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK2M2 <- attr(GET_erl, "p")
  
  ##Clean K3, M1 and run GET
  df <- data.frame(r = r, K3 = K3, t(K3M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  if(nrow(df) == 0){
    pvalK3M1 <- NA}else{
  r1 <- df$r
  K3M1 <- df$K3
  sim_K3 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K3M1, sim = sim_K3, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK3M1 <- attr(GET_erl, "p")
    }
  
  
  ##Clean K3, M2 and run GET
  df <- data.frame(r = r, K3 = K3, t(K3M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  if(nrow(df) == 0){
    pvalK3M2 <- NA
  }else{
    r1 <- df$r
    K3M2 <- df$K3
    sim_K3 <- as.matrix(df[,-c(1,2)])
    cset <- curve_set(obs = K3M2, sim = sim_K3, r = r1)
    plot(cset)
    GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05,
                                    alternative = "greater")
    pvalK3M2 <- attr(GET_erl, "p")
  }
  
  
  ##Clean K4, M1 and run GET
  df <- data.frame(r = r, K4 = K4, t(K4M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K4M1 <- df$K4
  sim_K4 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K4M1, sim = sim_K4, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK4M1 <- attr(GET_erl, "p")
  
  ##Clean K4, M2 and run GET
  df <- data.frame(r = r, K4 = K4, t(K4M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K4M2 <- df$K4
  sim_K4 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K4M2, sim = sim_K4, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK4M2 <- attr(GET_erl, "p")
  
  
  ##Clean K5, M1 and run GET
  df <- data.frame(r = r, K5 = K5, t(K5M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K5M1 <- df$K5
  sim_K5 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K5M1, sim = sim_K5, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK5M1 <- attr(GET_erl, "p")
  
  ##Clean K5, M2 and run GET
  df <- data.frame(r = r, K5 = K5, t(K5M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K5M2 <- df$K5
  sim_K5 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K5M2, sim = sim_K5, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK5M2 <- attr(GET_erl, "p")
  
  ##Clean K6, M1 and run GET
  df <- data.frame(r = r, K6 = K6, t(K6M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  if(nrow(df) == 0){
    pvalK6M1 <- NA
  }else{
    r1 <- df$r
    K6M1 <- df$K6
    sim_K6 <- as.matrix(df[,-c(1,2)])
    cset <- curve_set(obs = K6M1, sim = sim_K6, r = r1)
    plot(cset)
    GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                    alternative = "greater")
    pvalK6M1 <- attr(GET_erl, "p")
  }
  
  
  ##Clean K6, M2 and run GET
  df <- data.frame(r = r, K6 = K6, t(K6M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  if(nrow(df) == 0){
    pvalK6M2 <- NA
  }else{
  r1 <- df$r
  K6M2 <- df$K6
  sim_K6 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K6M2, sim = sim_K6, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK6M2 <- attr(GET_erl, "p")
  }
  
  
  ##Clean K7, M1 and run GET
  df <- data.frame(r = r, K7 = K7, t(K7M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K7M1 <- df$K7
  sim_K7 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K7M1, sim = sim_K7, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK7M1 <- attr(GET_erl, "p")
  
  ##Clean K7, M2 and run GET
  df <- data.frame(r = r, K7 = K7, t(K7M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K7M2 <- df$K7
  sim_K7 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K7M2, sim = sim_K7, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK7M2 <- attr(GET_erl, "p")
  
  
  ##Clean K8, M1 and run GET
  df <- data.frame(r = r, K8 = K8, t(K8M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K8M1 <- df$K8
  sim_K8 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K8M1, sim = sim_K8, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK8M1 <- attr(GET_erl, "p")
  
  ##Clean K8, M2 and run GET
  df <- data.frame(r = r, K8 = K8, t(K8M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  r1 <- df$r
  K8M2 <- df$K8
  sim_K8 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K8M2, sim = sim_K8, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK8M2 <- attr(GET_erl, "p")

  ##Clean K9, M1 and run GET
  df <- data.frame(r = r, K9 = K9, t(K9M1.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  if(nrow(df) == 0){
    pvalK9M1 <- NA
  }else{
    r1 <- df$r
    K9M1 <- df$K9
    sim_K9 <- as.matrix(df[,-c(1,2)])
    cset <- curve_set(obs = K9M1, sim = sim_K9, r = r1)
    plot(cset)
    GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                    alternative = "greater")
    pvalK9M1 <- attr(GET_erl, "p")
  }
  

  ##Clean K9, M2 and run GET
  df <- data.frame(r = r, K9 = K9, t(K9M2.perm))
  df <- df[apply(df, 1, function(x) all(is.finite(x))), ]
  if(nrow(df) == 0){
    pvalK9M2 <- NA}else{
  r1 <- df$r
  K9M2 <- df$K9
  sim_K9 <- as.matrix(df[,-c(1,2)])
  cset <- curve_set(obs = K9M2, sim = sim_K9, r = r1)
  plot(cset)
  GET_erl <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                  alternative = "greater")
  pvalK9M2 <- attr(GET_erl, "p")
    }
  
  
  ##M3. Calculate weights seperately as they are same for all stats
  vsimu <- cbind(v_perm, c(0,0))
  aij <- pairdist.default(t(vsimu))
  KK <- epanechnikov(aij,bw=bw)
  KK <- matrix(KK, nrow = n.perm + 1)
  wij <- KK/sum(KK)
  
  #K1
  K.simulated <- cbind(t(K1M2.perm), K1)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  
  pvalK1M3 <- p.value.kernel.K
  
  #K2
  K.simulated <- cbind(t(K2M2.perm), K2)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK2M3 <- p.value.kernel.K
  
  #K3
  K.simulated <- cbind(t(K3M2.perm), K3)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK3M3 <- p.value.kernel.K
  
  
  #K4
  K.simulated <- cbind(t(K4M2.perm), K4)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK4M3 <- p.value.kernel.K
  
  
  #K5
  K.simulated <- cbind(t(K5M2.perm), K5)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK5M3 <- p.value.kernel.K
  
  #K6
  K.simulated <- cbind(t(K6M2.perm), K6)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK6M3 <- p.value.kernel.K
  
  
  #K7
  K.simulated <- cbind(t(K7M2.perm), K7)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK7M3 <- p.value.kernel.K
  
  #K8
  K.simulated <- cbind(t(K8M2.perm), K8)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK8M3 <- p.value.kernel.K
  
  #K9
  K.simulated <- cbind(t(K9M2.perm), K9)
  Kmean <- apply(K.simulated, 1, mean)
  Si <- sweep(K.simulated, 1, Kmean, FUN = "-")
  si2 <- (Si ^ 2) %*% wij
  S <- (si2 ^ (-1/2)) * (Si)
  S[is.nan(S)] <- 0
  CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
  p.value.kernel.K <- attr(rank_envelope(CS, type = "erl"), "p")
  pvalK9M3 <- p.value.kernel.K
  
  
 return(list(pvalK1M1 = pvalK1M1, pvalK1M2 = pvalK1M2, pvalK1M3= pvalK1M3,
             pvalK2M1 = pvalK2M1, pvalK2M2 = pvalK2M2, pvalK2M3 = pvalK2M3,
             pvalK3M1 = pvalK3M1, pvalK3M2 = pvalK3M2, pvalK3M3 = pvalK3M3,
             pvalK4M1 = pvalK4M1, pvalK4M2 = pvalK4M2, pvalK4M3 = pvalK4M3,
             pvalK5M1 = pvalK5M1, pvalK5M2 = pvalK5M2, pvalK5M3 = pvalK5M3,
             pvalK6M1 = pvalK6M1, pvalK6M2 = pvalK6M2, pvalK6M3 = pvalK6M3,
             pvalK7M1 = pvalK7M1, pvalK7M2 = pvalK7M2, pvalK7M3 = pvalK7M3,
             pvalK8M1 = pvalK8M1, pvalK8M2 = pvalK8M2, pvalK8M3 =pvalK8M3,
             pvalK9M1 =pvalK9M1,  pvalK9M2 =pvalK9M2 ,pvalK9M3=pvalK9M3))
}






