##principled code for variance correction, count correction and modified count correction
##deprecated as of apr 15, 2026, as it did not perform well in terms of size. Will explore other methods for variance correction.
library(spatstat)
library(ggplot2)
library(psych)
library(parallel)
library(GET)
library(dplyr)
library(tidyr)
bw <- 10
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
# Epanechnikov kernel used for determining the weights in the regression
epanechnikov <- function(x,bw){
  return(pmax(0,(1/bw)*3/4*(1-(x/bw)^2)))
}
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/pval.R", echo = FALSE)

#load window
path <- "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv"
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/windows/", path, ".RData"))

W_new <- W
W_new$bdry[[1]]$x <- scale(W$bdry[[1]]$x, center = 0, scale = 6000)
W_new$bdry[[1]]$y <- scale(W$bdry[[1]]$y, center = 0, scale = 6000)
W_new$xrange <- range(W_new$bdry[[1]]$x)
W_new$yrange <- range(W_new$bdry[[1]]$y)

#simulation parameters
M = 4
p.cells = c(0.3, 0.4, 0.2, 0.1)
win = W_new

rmax = 0.15#incircle(win)$r/4
r = seq(0, rmax, length.out = 50)
n.sim = 500
n.perm = 199
cluster_sigma = 0.05
base.taxa = 1
shift.taxa = 2

#generate data


var.corr.test <- function(data, base.taxa, shift.taxa, r = NULL, n.perm, max.attempts = 10 * n.perm, bw =  "silverman"){
  if(is.null(r)){
    rmax <- incircle(data$window)$r/6
    r <- seq(0, rmax, length.out = 50)
  }
  
  freq_marks <- table(data$marks)
  lambda1 <- freq_marks[as.character(base.taxa)]
  lambda2 <- freq_marks[as.character(shift.taxa)]
  
  K0 <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa), r = r, correction = "border")[[3]]
  K_perm <- matrix(NA, nrow = length(r), ncol = n.perm)
  Kstar_perm <- matrix(NA, nrow = length(r), ncol = n.perm)
  n_reduced <- matrix(NA, nrow = 2, ncol = n.perm)
  nstar_reduced <- matrix(NA, nrow = 2, ncol = n.perm)
  W.area.reduced <- numeric(n.perm)
  
  perm <- 0
  attempts <- 0
  
  while(perm < n.perm){
    attempts <- attempts + 1
    
    # Safety check to prevent infinite loop
    if(attempts > max.attempts){
      stop(paste("Maximum number of attempts (", max.attempts, ") exceeded. Only", perm, 
                 "valid permutations found out of", n.perm, "requested."))
    }
    
    jump_rad <- incircle(data$window)$r
    shift_vector <- runifdisc(1, radius = jump_rad)
    shift_x <- shift_vector$x
    shift_y <- shift_vector$y
    
    # Step 2: Extract points of the type to shift
    points_to_shift <- subset(data, marks == shift.taxa)
    other_points <- subset(data, marks != shift.taxa)
    
    # Step 3: Shift the selected points
    shifted_points <- shift(points_to_shift, vec = c(shift_x, shift_y))
    
    # Step 4: Create shifted window
    W_shifted <- shift(W_new, vec = c(shift_x, shift_y))
    
    # Step 5: Define W.reduced as intersection
    W.reduced <- intersect.owin(W_new, W_shifted)
    
    # Check if W.reduced is valid
    if(is.null(W.reduced) || area.owin(W.reduced) == 0){
      next
    }
    W.area.reduced[perm] <- area.owin(W.reduced)
    # Step 6: Extract original pattern on W.reduced
    pp_original_reduced <- data[W.reduced]
    
    # Step 7: Create shifted version and extract on W.reduced
    pp_shifted_full <- superimpose(other_points, shifted_points)
    pp_shifted_reduced <- pp_shifted_full[W.reduced]
    
    # Check that both taxa are present in pp_original_reduced
    marks_orig <- table(pp_original_reduced$marks)
    has_base_orig <- as.character(base.taxa) %in% names(marks_orig) && marks_orig[as.character(base.taxa)] > 0
    has_shift_orig <- as.character(shift.taxa) %in% names(marks_orig) && marks_orig[as.character(shift.taxa)] > 0
    
    # Check that both taxa are present in pp_shifted_reduced
    marks_shifted <- table(pp_shifted_reduced$marks)
    has_base_shifted <- as.character(base.taxa) %in% names(marks_shifted) && marks_shifted[as.character(base.taxa)] > 0
    has_shift_shifted <- as.character(shift.taxa) %in% names(marks_shifted) && marks_shifted[as.character(shift.taxa)] > 0
    
    # If any required taxa is missing, skip this iteration
    if(!has_base_orig || !has_shift_orig || !has_base_shifted || !has_shift_shifted){
      next
    }
    
    # Step 8: Compute Kcross for both patterns on W.reduced
    # This permutation is valid, so increment counter
    perm <- perm + 1
    
    K_perm[, perm] <- Kcross(pp_original_reduced, i = as.character(base.taxa), j = as.character(shift.taxa), r = r, correction = "border")[[3]]
    ppred_prob <- table(pp_original_reduced$marks)
    n_reduced[,perm] <- c(ppred_prob[base.taxa] , ppred_prob[shift.taxa])
    Kstar_perm[, perm] <- Kcross(pp_shifted_reduced, i = as.character(base.taxa), j = as.character(shift.taxa), r = r, correction = "border")[[3]]
    pp_prob <- table(pp_shifted_reduced$marks)
    nstar_reduced[,perm] <- c(pp_prob[base.taxa] , pp_prob[shift.taxa])
  }
  
  #cset <- curve_set(obs = K0, sim = K_perm, r = r)
  #plot(cset)
  #cset_star <- curve_set(obs = K0, sim = Kstar_perm, r = r)
  #plot(cset_star)
  vsimu <- cbind(nstar_reduced, c(lambda1, lambda2))
  aij <- pairdist.default(t(vsimu))
  if(bw == "silverman"){
    bw <- 1.06 * mean(apply(vsimu, 1, sd)) * nrow(vsimu)^(-1/5)
  }
  print(bw)
  wts_gaussian <- wts.kernel(aij, bw = bw, type = "gaussian")
  wts_gaussian <- wts_gaussian/sum(wts_gaussian)
  
  Kcross.simulated <- cbind(Kstar_perm, K0)
  Kmean <-  apply(Kcross.simulated, 1, mean)
  Si <- sweep(Kcross.simulated, 1, Kmean, "-")
  si2_gauss <- (Si^2)%*% wts_gaussian
  S_gauss <- (si2_gauss)^(-0.5) * Si
  S_gauss[is.nan(S_gauss)] <- 0
  CS_gauss <- create_curve_set(list(r = r, obs = S_gauss[ , n.perm + 1], sim_m = S_gauss[ , 1:n.perm]))
  obj <- global_envelope_test(CS_gauss, typeone = "fwer", alternative = "greater", type = "erl")
  pval_n <- attr(obj, "p")
  #wsimu <- c(area.owin(data$window), W.area.reduced)
  #aij <- pairdist.default(matrix(wsimu, ncol = 1))
  #bw <- 1.06 * sd(wsimu) * length(wsimu)^(-1/5)
  #wts_gaussian <- wts.kernel(aij, bw = bw, type = "gaussian")
  #wts_gaussian <- wts_gaussian/sum(wts_gaussian)
  
  #si2_gauss  <- (Si^2)%*% wts_gaussian
  #S_gauss <- (si2_gauss)^(-0.5) * Si
  #S_gauss[is.nan(S_gauss)] <- 0
  #CS_gauss <- create_curve_set(list(r = r, obs = S_gauss[ , n.perm + 1], sim_m = S_gauss[ , 1:n.perm]))
 # pval_w <- attr(rank_envelope(CS_gauss, type = "erl"), "p")
 
  return(pval_n)
}
win = square(1)#W_new
#Run simulations
#pval <- matrix(NA, nrow = M, ncol = M)
type1_error <- matrix(NA, nrow = M, ncol = M)
for(base.taxa in 1:(M-1)){
  for(shift.taxa in (base.taxa+1):M){
    results <- mclapply(1:n.sim, function(i){
    data <- rcluster_marked_ppp_dependant(win = win, n_parent = 100, M = M, p.cells = p.cells, mu_offspring = 50, offspring_dist = "nbinom", sigma = 0.1, pair_types = c(1,3), pair_distance = 0.05)
    pval <- var.corr.test(data, base.taxa, shift.taxa, n.perm = 199)
    return(pval)
    }, mc.cores = detectCores() - 1)
    type1_error[base.taxa, shift.taxa] <- mean(unlist(results) < 0.05)
    print(paste("Base taxa:", base.taxa, "Shift taxa:", shift.taxa, "type1 err: ", type1_error[base.taxa, shift.taxa]))
  }
}


