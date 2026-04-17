##R script to run simulations

##Source scripts
rm(list = ls())
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/pval_v2.R", echo = TRUE)

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
win = W_new#square(1)
rmax = 0.15
r = seq(0, rmax, length.out = 50)
n.sim = 500
n.perm = 199
cluster_sigma = 0.05





##Type1 error matrices/rejection matrices

Kcross_torshift <- matrix(NA, nrow = M, ncol = M)
Kstar_torshift <- matrix(NA, nrow = M, ncol = M)
Kcor_torshift <- matrix(NA, nrow = M, ncol = M)
NN_torshift <- matrix(NA, nrow = M, ncol = M)

Kcross.vc.ep <- matrix(NA, nrow = M, ncol = M)
Kcross.vc.gauss <- matrix(NA, nrow = M, ncol = M)
Kcross.vc.uniform <- matrix(NA, nrow = M, ncol = M)

Kstar.vc.ep <- matrix(NA, nrow = M, ncol = M)
Kstar.vc.gauss <- matrix(NA, nrow = M, ncol = M)
Kstar.vc.uniform <- matrix(NA, nrow = M, ncol = M)

Kcor.vc.ep <- matrix(NA, nrow = M, ncol = M)
Kcor.vc.gauss <- matrix(NA, nrow = M, ncol = M)
Kcor.vc.uniform <- matrix(NA, nrow = M, ncol = M)

NN.vc.ep <- matrix(NA, nrow = M, ncol = M)
NN.vc.gauss <- matrix(NA, nrow = M, ncol = M)
NN.vc.uniform <- matrix(NA, nrow = M, ncol = M)

for(base.taxa in 1:(M-1)){
  for(shift.taxa in (base.taxa+1):M){
    results <- mclapply(1:n.sim, function(i){
 
      data <- rcluster_marked_ppp(win = win, n_parent = 100, M = M, p.cells = p.cells, mu_offspring = 50, offspring_dist = "nbinom", sigma = cluster_sigma)
      pvals <- pval.assoc(data, base.taxa = base.taxa, shift.taxa = shift.taxa,r = r, n.perm = 199, bw = 0.15)
      return(pvals)
      
    }, mc.cores = detectCores() - 1)
    # Process results and fill matrices
    results <- Filter(Negate(is.null), results)
    Kcross_torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval_Kcross <= 0.05))
    Kstar_torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval_Kstar <= 0.05))
    Kcor_torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval_Kcor <= 0.05))
    NN_torshift[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval_NN <= 0.05))
    
    Kcross.vc.ep[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.vc.ep <= 0.05))
    Kcross.vc.gauss[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.vc.gauss <= 0.05))
    Kcross.vc.uniform[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcross.vc.uniform <= 0.05))
    
    Kstar.vc.ep[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kstar.vc.ep <= 0.05))
    Kstar.vc.gauss[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kstar.vc.gauss <= 0.05))
    Kstar.vc.uniform[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kstar.vc.uniform <= 0.05))
    
    Kcor.vc.ep[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcor.vc.ep <= 0.05))
    Kcor.vc.gauss[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcor.vc.gauss <= 0.05))
    Kcor.vc.uniform[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.Kcor.vc.uniform <= 0.05))
    
    NN.vc.ep[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.NN.vc.ep <= 0.05))
    NN.vc.gauss[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.NN.vc.gauss <= 0.05))
    NN.vc.uniform[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.NN.vc.uniform <= 0.05))
    
    print(Kstar_torshift[base.taxa, shift.taxa])
    print(sprintf("Completed base taxa %d and shift taxa %d", base.taxa, shift.taxa))
  }
}
saveRDS(list(Kcross_torshift = Kcross_torshift, Kstar_torshift = Kstar_torshift, Kcor_torshift = Kcor_torshift, NN_torshift = NN_torshift,
         Kcross.vc.ep = Kcross.vc.ep, Kcross.vc.gauss = Kcross.vc.gauss, Kcross.vc.uniform = Kcross.vc.uniform,
         Kstar.vc.ep = Kstar.vc.ep, Kstar.vc.gauss = Kstar.vc.gauss, Kstar.vc.uniform = Kstar.vc.uniform,
         Kcor.vc.ep = Kcor.vc.ep, Kcor.vc.gauss = Kcor.vc.gauss, Kcor.vc.uniform = Kcor.vc.uniform,
         NN.vc.ep = NN.vc.ep, NN.vc.gauss = NN.vc.gauss, NN.vc.uniform = NN.vc.uniform), file = "simulation_results_irregW_5.rds")
beepr::beep(4)