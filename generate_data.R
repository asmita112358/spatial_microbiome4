#generate all types of data
library(spatstat)
library(ggplot2)

##Fucntion 1: generate parent clusters from csr, generate child clusters from parent
#sigma controls the spread of the cluster

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

##Function 2: LGCP

multi_type_lcgp <- function(win, M, mu, var, scale, p.cells){
  ppp_list <- list()
  for(i in 1:M){
    ppp_list[[i]] <- rLGCP("exp", mu = 3+10*p.cells[i], var = var, scale = scale, win = win)
  }
  X <- superimpose.ppp(ppp_list[[1]],ppp_list[[2]], ppp_list[[3]], W = win)
  marks(X) <- factor(rep(seq_len(M), sapply(ppp_list, npoints)), levels = seq_len(M))
  X
  
}

rcluster_marked_ppp_dependant <- function(win,
                                n_parent,
                                M,
                                p.cells = NULL,
                                mu_offspring = 50,
                                offspring_dist = c("poisson", "fixed", "nbinom"),
                                size_nbinom = 10,
                                sigma = 0.1,
                                keep_parents = FALSE,
                                pair_types = c(1, 2),  # types to pair together
                                pair_distance = 0.05)  # distance between paired parents
{
  offspring_dist <- match.arg(offspring_dist)
  
  if(is.null(p.cells)) p.cells <- rep(1/M, M)
  stopifnot(length(p.cells) == M, all(p.cells >= 0), abs(sum(p.cells) - 1) < 1e-8)
  
  # Generate parent locations and types
  parents <- runifpoint(n_parent, win = win)
  parent_type <- factor(sample(seq_len(M), size = n_parent, replace = TRUE, prob = p.cells),
                        levels = seq_len(M))
  
  # Pair up type 1 and type 2 parents
  type1_idx <- which(parent_type == pair_types[1])
  type2_idx <- which(parent_type == pair_types[2])
  
  # Determine how many pairs we can make
  n_pairs <- min(length(type1_idx), length(type2_idx))
  
  if(n_pairs > 0) {
    # Randomly pair them
    type1_paired <- sample(type1_idx, n_pairs)
    type2_paired <- sample(type2_idx, n_pairs)
    
    # Move type 2 parents close to their type 1 partners
    for(i in seq_len(n_pairs)) {
      idx1 <- type1_paired[i]
      idx2 <- type2_paired[i]
      
      # Place type 2 parent near type 1 parent
      angle <- runif(1, 0, 2*pi)
      parents$x[idx2] <- parents$x[idx1] + pair_distance * cos(angle)
      parents$y[idx2] <- parents$y[idx1] + pair_distance * sin(angle)
    }
    
    # Keep paired parents inside window
    for(idx in type2_paired) {
      if(!inside.owin(parents$x[idx], parents$y[idx], win)) {
        # If outside, try opposite direction or re-sample angle
        partner_idx <- type1_paired[which(type2_paired == idx)]
        for(attempt in 1:10) {
          angle <- runif(1, 0, 2*pi)
          new_x <- parents$x[partner_idx] + pair_distance * cos(angle)
          new_y <- parents$y[partner_idx] + pair_distance * sin(angle)
          if(inside.owin(new_x, new_y, win)) {
            parents$x[idx] <- new_x
            parents$y[idx] <- new_y
            break
          }
        }
      }
    }
  }
  
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

