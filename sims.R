##R script to run simulations

##Source scripts
rm(list = ls())
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/pval.R", echo = FALSE)

M = 4
p.cells = c(0.3, 0.4, 0.2, 0.1)
win = square(1)
rmax = 0.3
r = seq(0, rmax, length.out = 50)
n.sim = 500
n.perm = 199
cluster_sigma = 0.03





##Type1 error matrices

K1M1 <- K1M2 <- K1M3 <- K2M1 <- K2M2 <- K2M3 <- K3M1 <- K3M2 <- K3M3 <- matrix(NA, nrow = M, ncol = M)
K4M1 <- K4M2 <- K4M3 <- K5M1 <- K5M2 <- K5M3 <- K6M1 <- K6M2 <- K6M3 <- matrix(NA, nrow = M, ncol = M)
K7M1 <- K7M2 <- K7M3 <- K8M1 <- K8M2 <- K8M3 <- K9M1 <- K9M2 <- K9M3 <- matrix(NA, nrow = M, ncol = M)

for(base.taxa in 1:(M-1)){
  for(shift.taxa in (base.taxa+1):M){
    results <- mclapply(1:n.sim, function(i){
 
      data <- rcluster_marked_ppp(win = win, n_parent = 100, M = M, p.cells = p.cells, mu_offspring = 50, offspring_dist = "nbinom", sigma = cluster_sigma)
      #freq_marks <- table(data$marks)
      #lambda1 <- freq_marks[as.character(base.taxa)]
      #lambda2 <- freq_marks[as.character(shift.taxa)]
      pvals <- pval.assoc(data, base.taxa = base.taxa, shift.taxa = shift.taxa,r = r, n.perm = 199)
      return(pvals)
      
    }, mc.cores = detectCores() - 1)
    # Process results and fill matrices
    results <- Filter(Negate(is.null), results)
    K1M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK1M1 < 0.05))
    K1M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK1M2 < 0.05))
    K1M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK1M3 < 0.05))
    K2M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK2M1 < 0.05))
    K2M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK2M2 < 0.05))
    K2M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK2M3 < 0.05))
    K3M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK3M1 < 0.05))
    K3M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK3M2 < 0.05))
    K3M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK3M3 < 0.05))
    K4M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK4M1 < 0.05))
    K4M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK4M2 < 0.05))
     K4M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK4M3 < 0.05))
    K5M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK5M1 < 0.05))
    K5M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK5M2 < 0.05))
     K5M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK5M3 < 0.05))
    K6M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK6M1 < 0.05))
    K6M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK6M2 < 0.05))
     K6M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK6M3 < 0.05))
    K7M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK7M1 < 0.05))
    K7M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK7M2 < 0.05))
     K7M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK7M3 < 0.05))
    K8M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK8M1 < 0.05))
    K8M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK8M2 < 0.05))
     K8M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK8M3 < 0.05))
    K9M1[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK9M1 < 0.05))
    K9M2[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK9M2 < 0.05))
     K9M3[base.taxa, shift.taxa] <- mean(sapply(results, function(res) res$pvalK9M3 < 0.05))
    print(K1M1[base.taxa, shift.taxa])
    print(sprintf("Completed base taxa %d and shift taxa %d", base.taxa, shift.taxa))
  }
}
saveRDS(list(K1M1 = K1M1, K1M2 = K1M2, K1M3 = K1M3, K2M1 = K2M1, K2M2 = K2M2, K2M3 = K2M3, K3M1 = K3M1, K3M2 = K3M2, K3M3 = K3M3,
         K4M1 = K4M1, K4M2 = K4M2, K4M3 = K4M3, K5M1 = K5M1, K5M2 = K5M2, K5M3 = K5M3, K6M1 = K6M1, K6M2 = K6M2, K6M3 = K6M3,
         K7M1 = K7M1, K7M2 = K7M2, K7M3 = K7M3, K8M1 = K8M1, K8M2 = K8M2, 	K8M3 = 	K8M3,	K9M1 =	K9M1,	K9M2 =	K9M2,	K9M3 =	K9M3), 
        file = "size3.rds")
beepr::beep(4)