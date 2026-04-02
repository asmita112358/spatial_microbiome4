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
data <- multi_type_lcgp(square(1), M = 3, mu = 10, var = 1, scale = 0.5, p.cells = c(0.5, 0.3, 0.2))
data <- rcluster_marked_ppp(square(1), n_parent = 100, M = 4, p.cells = c(0.3, 0.4, 0.2, 0.1), mu_offspring = 50, offspring_dist = "poisson", sigma = 0.03)
df <- data.frame(data)
ggplot(df, aes(x = x, y = y, color = marks)) +
  geom_point(size = 1) +
  theme_minimal() +
  labs(title = "Marked PP", color = "Cell Type")
