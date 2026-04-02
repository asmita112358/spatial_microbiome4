##R code for both sided test. 

library(spatstat)
library(ggplot2)
library(psych)
library(parallel)
library(GET)

f <- function(m) m^2
rcluster_marked_ppp <- function(win,
                                n_parent,
                                M,
                                p.cells = NULL,
                                mu_offspring = 50,
                                offspring_dist = c("poisson", "fixed", "nbinom"),
                                size_nbinom = 10,
                                sigma = 0.1,
                                keep_parents = FALSE) {
  offspring_dist <- match.arg(offspring_dist)
  
  if(is.null(p.cells)) p.cells <- rep(1/M, M)
  stopifnot(length(p.cells) == M, all(p.cells >= 0), abs(sum(p.cells) - 1) < 1e-8)
  
  # parents iid uniform in window
  parents <- runifpoint(n_parent, win = win)
  
  # assign parent types
  parent_type <- factor(sample(seq_len(M), size = n_parent, replace = TRUE, prob = p.cells),
                        levels = seq_len(M))
  
  # offspring counts per parent
  n_off <- switch(offspring_dist,
                  poisson = rpois(n_parent, lambda = mu_offspring),
                  fixed   = rep(as.integer(mu_offspring), n_parent),
                  nbinom  = rnbinom(n_parent, mu = mu_offspring, size = size_nbinom))
  
  if(sum(n_off) == 0L) {
    X <- ppp(numeric(0), numeric(0), window = win,
             marks = factor(integer(0), levels = seq_len(M)))
    if(keep_parents) attr(X, "parents") <- data.frame(x=parents$x, y=parents$y, type=parent_type)
    return(X)
  }
  
  # repeat parents according to offspring count
  px <- rep(parents$x, times = n_off)
  py <- rep(parents$y, times = n_off)
  off_type <- rep(parent_type, times = n_off)
  
  # jitter
  x <- px + rnorm(sum(n_off), , sigma)
  y <- py + rnorm(sum(n_off), , sigma)
  
  # drop outside
  keep <- inside.owin(x, y, win)
  
  X <- ppp(x[keep], y[keep], window = win,
           marks = off_type[keep],
           check = FALSE)
  
  if(keep_parents) {
    attr(X, "parents") <- data.frame(x = parents$x, y = parents$y, type = parent_type)
  }
  X
}







generate_csr <- function(datatype = "nocluster", lambda = 1, M, p.cells = NULL,
                         mean_cluster_size = 100, cluster_spread = 0.1, win = square(10)) {
  if(is.null(p.cells)) {
    p.cells <- runif(M)
    p.cells <- p.cells/sum(p.cells)
  }
  
  if(datatype == "nocluster") {
    pp <- rpoispp(lambda * mean_cluster_size, win = win)
    n_pp <- pp$n
    pp$marks <- as.factor(sample(1:M, size = n_pp, replace = TRUE, prob = p.cells))
    
  } else if(datatype == "eqcluster") {
    pp <- rpoispp(lambda, win = win)
    n_cluster <- pp$n
    cluster_size <- rnbinom(n_cluster, mu = mean_cluster_size, size = 15)
    pp$marks <- as.factor(sample(1:M, size = n_cluster, replace = TRUE, prob = p.cells))
    x <- rep(pp$x, times = cluster_size) + rnorm(sum(cluster_size), sd = cluster_spread)
    y <- rep(pp$y, times = cluster_size) + rnorm(sum(cluster_size), sd = cluster_spread)
    type <- rep(pp$marks, times = cluster_size)
    pp <- suppressWarnings(ppp(x, y, window = win, marks = factor(type)))
    
  } else if(datatype == "uneqcluster") {
    pp <- rpoispp(lambda, win = win)
    n_cluster <- pp$n
    pp$marks <- as.factor(sample(1:M, size = n_cluster, replace = TRUE, prob = p.cells))
    cluster_size <- rnbinom(n_cluster, mu = mean_cluster_size, size = 15) + 
      sapply(as.numeric(pp$marks), f)
    x <- rep(pp$x, times = cluster_size) + rcauchy(sum(cluster_size), scale = cluster_spread)
    y <- rep(pp$y, times = cluster_size) + rcauchy(sum(cluster_size), scale = cluster_spread)
    type <- rep(pp$marks, times = cluster_size)
    pp <- suppressWarnings(ppp(x, y, window = win, marks = factor(type)))
  }
  return(pp)
}

p.CCT <- function(x) {
  CCT <- sum(tan((0.5 - x) * pi))
  pval <- 0.5 - atan(CCT) / pi
  return(pval)
}

compute_K_sym <- function(data, base.taxa, shift.taxa, lambda1, lambda2, 
                          lambda1.overall, lambda2.overall) {
  obj1 <- Kcross.inhom(data, i = as.character(base.taxa), j = as.character(shift.taxa), 
                       lambdaI = lambda1, lambdaJ = lambda2, 
                       correction = "isotropic")
  obj2 <- Kcross.inhom(data, i = as.character(shift.taxa), j = as.character(base.taxa), 
                       lambdaI = lambda2, lambdaJ = lambda1, 
                      correction = "isotropic")
  K_sym <- (lambda2.overall * obj1$iso + lambda1.overall * obj2$iso) /  (lambda1.overall + lambda2.overall)
  #auc_K <- integrate(K_sym, lower = min(obj1$r), upper = max(obj1$r))$value
  return(list( K_sym = K_sym, r = obj1$r))
}

# Optimized single simulation function

run_single_simulation <- function(sim, W_new, p.cells, M, n.perm, 
                                  shift.taxa, base.taxa, cluster_sigma = 0.1) {
  
  # Generate data
  #data <- generate_csr(datatype = "uneqcluster", lambda = 0.1, M = M, 
   #                    p.cells = p.cells, mean_cluster_size = 100, 
    #                   cluster_spread = 0.2, win = W_new)
  data <- rcluster_marked_ppp(win = W_new, n_parent = 50, M = M, p.cells = p.cells, 
                        mu_offspring = 50, offspring_dist = "nbinom", size_nbinom = 20, 
                        sigma = cluster_sigma, keep_parents = FALSE)
  freq_marks <- table(data$marks)
  
  # Check if both taxa exist - FIXED to use as.character
  if(!all(c(as.character(base.taxa), as.character(shift.taxa)) %in% names(freq_marks))) {
    return(NULL)
  }
  
  lambda1.overall <- as.numeric(freq_marks[as.character(base.taxa)])
  lambda2.overall <- as.numeric(freq_marks[as.character(shift.taxa)])
  
  # Compute bandwidths ONCE for original data - FIXED to use as.character for factor subsetting
  data1 <- data[data$marks == as.character(base.taxa)]
  data2 <- data[data$marks == as.character(shift.taxa)]
  
  bw1 <- bw.ppl(data1)
  bw2 <- bw.ppl(data2)
  
  lambda1 <- density.ppp(data1, sigma = bw1)
  lambda2 <- density.ppp(data2, sigma = bw2)
  
  # Observed statistic - FIXED to pass base.taxa and shift.taxa
  obs_result <- compute_K_sym(data, base.taxa, shift.taxa, 
                              lambda1, lambda2, lambda1.overall, lambda2.overall)
  K_sym <- obs_result$K_sym
  r <- obs_result$r
  #auc_K <- obs_result$auc_K
  
  medr <- which(r==median(r))
  
  # Prepare for toroidal shift
  bb <- boundingbox(W_new)
  R <- owin(xrange = bb$xrange, yrange = bb$yrange)
  data_tor <- data
  Window(data_tor) <- R
  
  # Pre-compute window boundary info
  Wlen <- length(W_new$bdry[[1]]$x)
  W_x <- W_new$bdry[[1]]$x
  W_y <- W_new$bdry[[1]]$y
  
  # Permutation loop
  K_sym.shifted <- vector("list", n.perm)
  #auc_K.shifted <- numeric(n.perm)
  
  for(perm in 1:n.perm) {
    # FIXED: pass marks as character to rshift
    data.shifted <- rshift(data_tor, which = as.character(shift.taxa), edge = "torus")
    data.shifted <- rshift(data.shifted, which = as.character(base.taxa), edge = "torus")
    
    # Bring points inside after toroidal shift
    inside <- inside.owin(data.shifted$x, data.shifted$y, W_new)
    outside <- !inside
    max_iter <- 100
    iter <- 0
    
    while(sum(outside) > 0.01 * min(freq_marks[as.character(shift.taxa)], 
                                    freq_marks[as.character(base.taxa)]) && iter < max_iter) {
      ##shift taxa 1 inside - FIXED to use as.character
      outside.1 <- outside & (data.shifted$marks == as.character(base.taxa))
      alpha <- runif(Wlen)
      alpha <- alpha / sum(alpha)
      x_in <- sum(alpha * W_x)
      y_in <- sum(alpha * W_y)
      centroid <- c(mean(data.shifted$x[outside.1]), mean(data.shifted$y[outside.1]))
      shift_vec <- centroid - c(x_in, y_in)
      data.shifted$x[outside.1] <- data.shifted$x[outside.1] - shift_vec[1]
      data.shifted$y[outside.1] <- data.shifted$y[outside.1] - shift_vec[2]
      
      ##shift taxa 2 inside independently - FIXED to use as.character
      outside.2 <- outside & (data.shifted$marks == as.character(shift.taxa))
      alpha <- runif(Wlen)
      alpha <- alpha / sum(alpha)
      x_in <- sum(alpha * W_x)
      y_in <- sum(alpha * W_y)
      centroid <- c(mean(data.shifted$x[outside.2]), mean(data.shifted$y[outside.2]))
      shift_vec <- centroid - c(x_in, y_in)
      data.shifted$x[outside.2] <- data.shifted$x[outside.2] - shift_vec[1]
      data.shifted$y[outside.2] <- data.shifted$y[outside.2] - shift_vec[2]
      
      inside <- inside.owin(data.shifted$x, data.shifted$y, W_new)
      outside <- !inside
      iter <- iter + 1
    }
    
    # Remove remaining outside points if any
    if(sum(outside) > 0) {
      data.shifted$x <- data.shifted$x[inside]
      data.shifted$y <- data.shifted$y[inside]
      data.shifted$marks <- data.shifted$marks[inside]
    }
    
    data.shifted.clean <- ppp(x = data.shifted$x,
                              y = data.shifted$y,
                              window = W_new,
                              marks = marks(data.shifted),
                              check = FALSE)
    #do this untill theres atleast 10 marks of shift.taxa in shifted process
    #data.shifted.clean = rshift(data_tor, which = shift.taxa,edge = "torus")
    #if(sum(data.shifted.clean$marks == as.character(shift.taxa)) < 10) {
     # next
    #}
    # Use FIXED bandwidth from original data - FIXED to use as.character
    #bw1 <- bw.ppl(data.shifted.clean[data.shifted.clean$marks == as.character(base.taxa)])
    #bw2 <- bw.ppl(data.shifted.clean[data.shifted.clean$marks == as.character(shift.taxa)])
    lambda1.shifted <- density.ppp(data.shifted.clean[data.shifted.clean$marks == as.character(base.taxa)], 
                                   sigma = bw1)
    lambda2.shifted <- density.ppp(data.shifted.clean[data.shifted.clean$marks == as.character(shift.taxa)], 
                                   sigma = bw2)
    lambda1.overall <- sum(data.shifted.clean$marks == as.character(base.taxa))
    lambda2.overall <- sum(data.shifted.clean$marks == as.character(shift.taxa))
    
    # FIXED: pass base.taxa and shift.taxa
    shifted_result <- compute_K_sym(data.shifted.clean, base.taxa, shift.taxa,
                                    lambda1.shifted, lambda2.shifted,
                                    lambda1.overall, lambda2.overall)
    K_sym.shifted[[perm]] <- shifted_result$K_sym
    #auc_K.shifted[perm] <- shifted_result$auc_K
  }
  
  
  # Compute p-values
  tmp <- do.call(rbind, lapply(K_sym.shifted, function(x) x >= K_sym))
  pval_torshift <- (1 + colSums(tmp)) / (1 + n.perm)
  
  hmp.pval <- harmonic.mean(pval_torshift)
  #cct.pval <- p.CCT(pval_torshift)
  
  #auc_pval <- (1 + sum(auc_K.shifted >= auc_K)) / (1 + n.perm)
  
  cset <- curve_set(obs = K_sym, sim = do.call(cbind, K_sym.shifted), r = r)
     
  test_result_fwer <- global_envelope_test(cset,typeone = "fwer", type = "rank", alpha = 0.05, 
                                            alternative = "greater", ties = "conservative")
  pval.fwer <- attr(test_result_fwer, "p")
  
  pval.med <- pval_torshift[medr]
  
  GET_random <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                 alternative = "greater", ties = "random")
  GET_midrank <- global_envelope_test(cset, typeone = "fwer", type = "rank", alpha = 0.05, 
                                 alternative = "greater", ties = "midrank")
  pval.random <- attr(GET_random, "p")
  pval.midrank <- attr(GET_midrank, "p")
  
  return(list(pval_torshift = pval_torshift, 
              pval.random = pval.random, 
              pval.fwer = pval.fwer,
              pval.midrank = pval.midrank,
              pval.med = pval.med,
              pval.hmp = hmp.pval))
}

# Main execution
#set.seed(123)  # For reproducibility

M <- 4
p.cells <- runif(M)
p.cells <- p.cells / sum(p.cells)
n.perm <- 199
n.sim <- 100


# Load and scale window
path <- "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv"
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome3/windows/", path, ".RData"))

W_new <- W
W_new$bdry[[1]]$x <- scale(W$bdry[[1]]$x, center = 0, scale = 600)
W_new$bdry[[1]]$y <- scale(W$bdry[[1]]$y, center = 0, scale = 600)
W_new$xrange <- range(W_new$bdry[[1]]$x)
W_new$yrange <- range(W_new$bdry[[1]]$y)

# Run simulations
GET.random <- matrix(NA, nrow = M, ncol = M)
GET.midrank <- matrix(NA, nrow = M, ncol = M)
GET.cons <- matrix(NA, nrow = M, ncol = M)
med_R_test <- matrix(NA, nrow = M, ncol = M)
hmp.pval <- matrix(NA, nrow = M, ncol = M)
pval_torshift_all <- list()

for(base.taxa in 1:(M - 1)){
  for(shift.taxa in (base.taxa+1):M){
    results <- mclapply(1:n.sim, function(sim) {
      tryCatch({
        run_single_simulation(sim, W_new, p.cells, M, n.perm, shift.taxa, base.taxa, cluster_sigma = 0.1)
      }, error = function(e) {
        message(sprintf("Simulation %d failed: %s", sim, e$message))
        return(NULL)
      })
    }, mc.cores = 10)
    # Remove failed simulations
    results <- Filter(Negate(is.null), results)
    GET.random[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.random < 0.05))
    GET.midrank[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.midrank < 0.05))
    GET.cons[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.fwer < 0.05))
    med_R_test[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.med < 0.05))
    hmp.pval[base.taxa, shift.taxa] <- mean(sapply(results, function(x) x$pval.hmp < 0.05))
    pval_torshift_all[[paste(base.taxa, shift.taxa, sep = "_")]] <- sapply(results, function(x) x$pval_torshift)
    print(paste("GET random:", GET.random[base.taxa, shift.taxa], 
                "GET midrank:", GET.midrank[base.taxa, shift.taxa], 
                "GET cons:", GET.cons[base.taxa, shift.taxa], 
                "med_R_test:", med_R_test[base.taxa, shift.taxa],
                "hmp.pval:", hmp.pval[base.taxa, shift.taxa]))
    print(paste("Completed simulations for base taxa", base.taxa, "and shift taxa", shift.taxa))
  }
}


##save results in RDS
saveRDS(list(hmp = hmp.pval, GET.random = GET.random, GET.midrank = GET.midrank, GET.cons = GET.cons,med_R_test = med_R_test,p.cells = p.cells), 
        file = "simulation_results_bothsidedM4_s1.RDS")
beepr::beep(4)

# Analyze Type I error
alpha <- 0.05
type1_hmp <- mean(sapply(results, function(x) x$hmp.pval < alpha))
type1_cct <- mean(sapply(results, function(x) x$cct.pval < alpha))
type1_fwer <- mean(sapply(results, function(x) x$pval.fwer < alpha))
type1_fdr <- mean(sapply(results, function(x) x$pval.fdr < alpha))

pval_fwer <- sapply(results, function(x) x$pval.fwer)
pval_fdr <- sapply(results, function(x) x$pval.fdr)


