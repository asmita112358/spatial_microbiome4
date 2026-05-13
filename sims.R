##R script to run simulations

##Source scripts
rm(list = ls())
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/pval_v3.R", echo = FALSE)

#Irregular window if needed
path <- "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv"
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/windows/", path, ".RData"))

W_new <- W
W_new$bdry[[1]]$x <- scale(W$bdry[[1]]$x, center = 0, scale = 6000)
W_new$bdry[[1]]$y <- scale(W$bdry[[1]]$y, center = 0, scale = 6000)
W_new$xrange <- range(W_new$bdry[[1]]$x)
W_new$yrange <- range(W_new$bdry[[1]]$y)

M = 4
p.cells = c(0.3, 0.4, 0.2, 0.1)
win = square(1)
rmax = 0.15
r = seq(0, rmax, length.out = 50)
n.sim = 500
n.perm = 499
cluster_sigma = 0.05





##Type1 error matrices/rejection matrices

Kcross_torshift <- matrix(NA, nrow = M, ncol = M)
Kstar_torshift <- matrix(NA, nrow = M, ncol = M)
Kcor_torshift <- matrix(NA, nrow = M, ncol = M)


Kcross.vc.gauss <- matrix(NA, nrow = M, ncol = M)



Kcross.evc <- matrix(NA, nrow = M, ncol = M)
Kcross.minus <- matrix(NA, nrow = M, ncol = M)
Kcross.vc.n <- matrix(NA, nrow = M, ncol = M)
Kcross.evc.n <- matrix(NA, nrow = M, ncol = M)

auc.torshift <- matrix(NA, nrow = M, ncol = M)
auc.vc <- matrix(NA, nrow = M, ncol = M)
auc.minus <- matrix(NA, nrow = M, ncol = M)

median.torshift <- matrix(NA, nrow = M, ncol = M)
median.vc <- matrix(NA, nrow = M, ncol = M)
median.minus <- matrix(NA, nrow = M, ncol = M)

for(base.taxa in 1:(M-1)){
  for(shift.taxa in (base.taxa+1):M){
    results <- mclapply(1:n.sim, function(i){
 
      data <- rcluster_marked_ppp_dependant(win = win, n_parent = 100, M = M, p.cells = p.cells, mu_offspring = 50, offspring_dist = "nbinom", sigma = cluster_sigma, pair_types = c(1,3), pair_distance = 0.1)
      pvals <- pval.assoc(data, base.taxa = base.taxa, shift.taxa = shift.taxa,r = r, n.perm = 199, bw = "silverman")
      return(pvals)
      
    }, mc.cores = detectCores() - 1)
    # Process results and fill matrices
    results <- Filter(Negate(is.null), results)
    results <- results[!sapply(results, is.character)]
    Kcross_torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval_Kcross <= 0.05))
    Kstar_torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval_Kstar <= 0.05))
    Kcor_torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval_Kcor <= 0.05))
       
    Kcross.vc.gauss[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.vc.gauss <= 0.05))
    Kcross.vc.n[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.vc.n <= 0.05))
    Kcross.evc[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.evc <= 0.05))
    Kcross.evc.n[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.evc.n <= 0.05))
    Kcross.minus[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.minus <= 0.05))
    
     auc.torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.auc.tor <= 0.05))
     auc.vc[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.auc.vc <= 0.05))
     auc.minus[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.auc.minus <= 0.05))
     
     median.torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.medr.tor <= 0.05))
     median.vc[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.medr.vc <= 0.05))
     median.minus[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.medr.minus <= 0.05))
    
    print(median.vc[base.taxa, shift.taxa])
    print(sprintf("Completed base taxa %d and shift taxa %d", base.taxa, shift.taxa))
  }
}
saveRDS(list(Kcross_torshift = Kcross_torshift, Kstar_torshift = Kstar_torshift, Kcor_torshift = Kcor_torshift,
             Kcross.vc.gauss = Kcross.vc.gauss, Kcross.vc.n = Kcross.vc.n, Kcross.evc = Kcross.evc, Kcross.evc.n = Kcross.evc.n, Kcross.minus = Kcross.minus,
             auc.torshift = auc.torshift, median.torshift = median.torshift, auc.vc = auc.vc, median.vc = median.vc), file = "sim_5.rds")