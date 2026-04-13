##call libraries
library(spatstat)
library(ggplot2)
library(psych)
library(parallel)
library(GET)
library(dplyr)
library(tidyr)
##Compute statistics

compute_K <- function(data, base.taxa, shift.taxa, lambda1, lambda2, r = NULL){
  obj12_b <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa),
                 r =r, correction = "border")
  obj21_b <- Kcross(data, i = as.character(shift.taxa), j = as.character(base.taxa),
                 r = r, correction = "border")
  obj12_t <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa), r = r, correction = "translate")
  obj21_t <- Kcross(data, i = as.character(shift.taxa), j = as.character(base.taxa), r = r, correction = "translate")
  
  obj12_i <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa), r = r, correction = "isotropic")
  obj21_i <- Kcross(data, i = as.character(shift.taxa), j = as.character(base.taxa), r = r, correction = "isotropic")
  
  obj11_b <- Kest(data[data$marks == as.character(base.taxa)], r = r, correction = "border")
  obj22_b <- Kest(data[data$marks == as.character(shift.taxa)], r = r, correction = "border")
  
  obj11_t <- Kest(data[data$marks == as.character(base.taxa)], r = r, correction = "translate")
  obj22_t <- Kest(data[data$marks == as.character(shift.taxa)], r = r, correction = "translate")
  
  obj11_i <- Kest(data[data$marks == as.character(base.taxa)], r = r, correction = "isotropic")
  obj22_i <- Kest(data[data$marks == as.character(shift.taxa)], r = r, correction = "isotropic")
  
  
  #Stat1: K12
  stat1_border <- obj12_b$border
  stat1_translate <- obj12_t$trans
  stat1_isotropic <- obj12_i$iso
  
  #Stat2: K_star
  Kstar_b <- (lambda2*obj12_b$border + lambda1*obj21_b$border)/(lambda1+lambda2)
  Kstar_t <- (lambda2*obj12_t$trans + lambda1*obj21_t$trans)/(lambda1 + lambda2)
  Kstar_i <- (lambda2*obj12_i$iso + lambda1*obj21_i$iso)/(lambda1 + lambda2)
  
  #Stat3: Kcor
  Kcor_b <- (Kstar_b )/(sqrt(obj11_b$border*obj22_b$border))
  Kcor_t <- (Kstar_t )/(sqrt(obj11_t$trans*obj22_t$trans))
  Kcor_i <- (Kstar_i )/(sqrt(obj11_i$iso*obj22_i$iso))
  
  return(list(K1 = stat1_border, K2 = Kstar_b, K3 = Kcor_b,
              K4 = stat1_translate, K5 = Kstar_t, K6 = Kcor_t,
              K7 = stat1_isotropic, K8 = Kstar_i, K9 = Kcor_i, r = obj12_b$r))
  
  
}

## After we get simulation results and decide on a statistic, will put in the function below that computes the shortlisted statistics

compute_K_edge <- function(data, base.taxa, shift.taxa, lambda1, lambda2, r = NULL, correction = "border"){
  
  if(correction == "border"){
    obj12 <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa),
                    r =r, correction = "border")
    obj21 <- Kcross(data, i = as.character(shift.taxa), j = as.character(base.taxa),
                    r = r, correction = "border")
    
    obj11 <- Kest(data[data$marks == as.character(base.taxa)], r = r, correction = "border")
    obj22 <- Kest(data[data$marks == as.character(shift.taxa)], r = r, correction = "border")
    
  } else if(correction == "translate"){
    obj12 <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa),
                    r =r, correction = "translate")
    obj21 <- Kcross(data, i = as.character(shift.taxa), j = as.character(base.taxa),
                    r = r, correction = "translate")
    
    obj11 <- Kest(data[data$marks == as.character(base.taxa)], r = r, correction = "translate")
    obj22 <- Kest(data[data$marks == as.character(shift.taxa)], r = r, correction = "translate")
    
  } else if(correction == "isotropic"){
    obj12 <- Kcross(data, i = as.character(base.taxa), j = as.character(shift.taxa),
                    r =r, correction = "isotropic")
    obj21 <- Kcross(data, i = as.character(shift.taxa), j = as.character(base.taxa),
                    r = r, correction = "isotropic")
    
    obj11 <- Kest(data[data$marks == as.character(base.taxa)], r = r, correction = "isotropic")
    obj22 <- Kest(data[data$marks == as.character(shift.taxa)], r = r, correction = "isotropic")
  }
  
  
  #Stat1: K12
  stat1 <- obj12[[3]]
  #Stat2: K_star
  Kstar <- (lambda2*obj12[[3]] + lambda1*obj21[[3]])/(lambda1+lambda2)
  #Stat3: Kcor
  Kcor <- (Kstar)/(sqrt(obj11[[3]]*obj22[[3]]))
  
  return(list(K12 = stat1, Kstar = Kstar, Kcor = Kcor, r = obj12$r))
  
}
