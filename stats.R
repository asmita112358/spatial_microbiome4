##call libraries
library(spatstat)
library(ggplot2)
library(psych)
library(parallel)
library(GET)
library(dplyr)
library(tidyr)
library(RANN)
library(FNN)
library(pracma)
##Compute statistics


##All K functions with border correction

compute_K <- function(data, base.taxa, shift.taxa, lambda1, lambda2, r = NULL){
  obj12 <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa),
                 r =r, correction = "border")
  obj21 <- Kcross(data, i = as.character(shift.taxa), j = as.character(base.taxa),
                 r = r, correction = "border")
  
  
  #stat1: Kcross
  Kcross_stat <- obj12$border
  
  #stat2: Kstar
  Kstar <- (lambda2*obj12$border + lambda1*obj21$border)/(lambda1+lambda2)
  
  #stat3: Kcor
  obj11 <- Kest(data[data$marks == as.character(base.taxa)], r = r, correction = "border")
  obj22 <- Kest(data[data$marks == as.character(shift.taxa)], r = r, correction = "border")
  Kcor <- (Kstar )/(sqrt(obj11$border*obj22$border))
  
  return(list(Kcross_stat = Kcross_stat, Kstar = Kstar, Kcor = Kcor, r = obj12$r))
}

##NN function

compute_NN <- function(data,base.taxa, shift.taxa){
  base.taxa.points <- subset(data, marks == as.character(base.taxa))
  all.coords <- cbind(data$x, data$y)
  base.taxa.coords <- cbind(base.taxa.points$x, base.taxa.points$y)
  nn <- get.knnx(all.coords, base.taxa.coords, k = 2, algorithm = "kd_tree")
  nn.types <- data$marks[nn$nn.index[,2]]
  Nij <- sum(nn.types == as.character(shift.taxa))/sum(data$marks == as.character(shift.taxa))
  
  return(Nij)
}

##AUC under K
compute_AUC <- function(K_stat, r){
  # Assuming K_stat is a vector of K values corresponding to the r values
  # We will compute the AUC using the trapezoidal rule
  auc <- trapz(r, K_stat)
  return(auc)
}

compute_NNX <- function(data, base.taxa, shift.taxa){
  # Compute the proportion of base.taxa points whose nearest cross-type neighbor 
  # (i.e., nearest neighbor that is NOT base.taxa) is of type shift.taxa,
  # normalized by the total number of shift.taxa points
  #
  # Returns: (# of base.taxa points with shift.taxa as cross-type NN) / (# of shift.taxa points)
  
  base.taxa.points <- subset(data, marks == as.character(base.taxa))
  other.taxa.points <- subset(data, marks != as.character(base.taxa))
  
  n_base <- npoints(base.taxa.points)
  n_shift <- sum(data$marks == as.character(shift.taxa))
  n_other <- npoints(other.taxa.points)
  
  if (n_base == 0 || n_shift == 0 || n_other == 0) {
    warning("One or more taxa groups have zero points")
    return(NA)
  }
  
  base.taxa.coords <- cbind(base.taxa.points$x, base.taxa.points$y)
  other.taxa.coords <- cbind(other.taxa.points$x, other.taxa.points$y)
  
  # Find nearest neighbor among non-base.taxa points for each base.taxa point
  nn <- get.knnx(other.taxa.coords, base.taxa.coords, k = 1, algorithm = "kd_tree")
  
  # Get the types of the nearest cross-type neighbors
  nn.types <- other.taxa.points$marks[nn$nn.index[, 1]]
  
  # Count how many base.taxa points have shift.taxa as their nearest cross-type neighbor
  n_cross_nn <- sum(nn.types == as.character(shift.taxa))
  
  # Normalize by number of shift.taxa points
  NNX <- n_cross_nn / n_shift
  
  return(NNX)
}

compute_NNX_d_all <- function(data, base.taxa, shift.taxa){
  # Compute the average distance from ALL base.taxa points to their nearest shift.taxa point
  # (considering only shift.taxa points, not all cross-type points)
  #
  # This differs from compute_NNX_d in that it doesn't filter by whether shift.taxa
  # is the nearest cross-type neighbor overall
  #
  # Returns: mean distance from base.taxa points to nearest shift.taxa point
  
  base.taxa.points <- subset(data, marks == as.character(base.taxa))
  shift.taxa.points <- subset(data, marks == as.character(shift.taxa))
  
  n_base <- npoints(base.taxa.points)
  n_shift <- npoints(shift.taxa.points)
  
  if (n_base == 0 || n_shift == 0) {
    warning("One or both taxa have zero points")
    return(NA)
  }
  
  base.coords <- cbind(base.taxa.points$x, base.taxa.points$y)
  shift.coords <- cbind(shift.taxa.points$x, shift.taxa.points$y)
  
  # Find nearest shift.taxa neighbor for each base.taxa point
  nn <- get.knnx(shift.coords, base.coords, k = 1, algorithm = "kd_tree")
  
  # Get distances to nearest shift.taxa neighbors
  nn.distances <- nn$nn.dist[, 1]
  
  # Compute mean distance across all base.taxa points
  mean_distance <- mean(nn.distances)
  
  return(mean_distance)
}

