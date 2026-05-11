##exploratory code to check behaviour of sd versus sample sizes and shift vectors
library(keras3)
library(reticulate)
library(bspline)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/pval.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/stats.R", echo = FALSE)
predict_S_curve_nn <- function(X, Y, newX,
                               hidden = c(128, 128),
                               dropout = 0.1,
                               lr = 1e-3,
                               batch_size = 32,
                               epochs = 100,
                               patience = 15,
                               val_split = 0.2,
                               scale_X = TRUE,
                               scale_Y = TRUE,
                               seed = 1,
                               verbose = 0) {
  stopifnot(is.matrix(X), is.matrix(Y))
  if (!is.matrix(newX)) newX <- as.matrix(newX)
  stopifnot(ncol(newX) == ncol(X))
  
  suppressPackageStartupMessages(library(keras3))
  
  set.seed(seed)
  
  # --- Scale X ---
  if (scale_X) {
    X_center <- colMeans(X)
    X_scale  <- apply(X, 2, sd)
    X_scale[X_scale == 0] <- 1
    
    Xs <- scale(X, center = X_center, scale = X_scale)
    newXs <- scale(newX, center = X_center, scale = X_scale)
  } else {
    X_center <- rep(0, ncol(X)); X_scale <- rep(1, ncol(X))
    Xs <- X; newXs <- newX
  }
  
  # --- Scale Y ---
  if (scale_Y) {
    Y_center <- colMeans(Y)
    Y_scale  <- apply(Y, 2, sd)
    Y_scale[Y_scale == 0] <- 1
    
    Ys <- scale(Y, center = Y_center, scale = Y_scale)
  } else {
    Y_center <- rep(0, ncol(Y)); Y_scale <- rep(1, ncol(Y))
    Ys <- Y
  }
  
  # --- Build model (keras3 style) ---
  inputs <- layer_input(shape = ncol(Xs))
  
  x <- inputs |>
    layer_dense(units = hidden[1], activation = "relu")
  
  if (!is.null(dropout) && dropout > 0) {
    x <- x |> layer_dropout(rate = dropout)
  }
  
  if (length(hidden) > 1) {
    for (u in hidden[-1]) {
      x <- x |> layer_dense(units = u, activation = "relu")
    }
  }
  
  outputs <- x |> layer_dense(units = ncol(Ys), activation = "linear")
  
  model <- keras_model(inputs = inputs, outputs = outputs)
  
  model |> compile(
    optimizer = optimizer_adam(learning_rate = lr),
    loss = "mse",
    metrics = "mae"
  )
  
  # --- Fit ---
  model |> fit(
    x = Xs, y = Ys,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = val_split,
    callbacks = list(callback_early_stopping(
      patience = patience,
      restore_best_weights = TRUE
    )),
    verbose = verbose
  )
  
  # --- Predict ---
  pred_scaled <- model |> predict(newXs, verbose = 0)
  
  # --- Unscale Y ---
  pred <- sweep(pred_scaled, 2, Y_scale, `*`)
  pred <- sweep(pred,       2, Y_center, `+`)
  
  pred
}
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
n.perm = 999
cluster_sigma = 0.05

base.taxa = 1
shift.taxa = 3

data <- rcluster_marked_ppp_dependant(win = win, n_parent = 100, M = M, p.cells = p.cells, mu_offspring = 50, offspring_dist = "nbinom", sigma = cluster_sigma, pair_types = c(1,3), pair_distance = 0.05)
rlen = length(r)
freq_marks <- table(data$marks)
original_window <- data$window
x_range <- diff(original_window$xrange)
y_range <- diff(original_window$yrange)
lambda1 <- freq_marks[as.character(base.taxa)]
lambda2 <- freq_marks[as.character(shift.taxa)]

##Compute stats##

obj_K <- compute_K(data, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = lambda1, lambda2 = lambda2, r = r)

Kcross_stat <- obj_K$Kcross_stat
Kstar <- obj_K$Kstar
Kcor <- obj_K$Kcor[-1]
rcor <- r[-1]
##rshift data and compute stats##
Kcross_stat.tor <- matrix(NA, nrow = n.perm, ncol = rlen)
Kcross_stat.vc <- matrix(NA, nrow = n.perm, ncol = rlen)


Kstar.tor <- matrix(NA, nrow = n.perm, ncol = rlen)
Kstar.vc <- matrix(NA, nrow = n.perm, ncol = rlen)


Kcor.tor <- matrix(NA, nrow = n.perm, ncol = rlen-1)
Kcor.vc <- matrix(NA, nrow = n.perm, ncol = rlen-1)



v_perm <- matrix(nrow = 2, ncol = n.perm)
nstar_reduced <- matrix(NA, nrow = 2, ncol = n.perm)
perm <- 1
data_tor <- data
##check if window is rectangular
if(!is.rectangle(data$win)){
  Window(data_tor) <- boundingbox(data$win)
}
while(perm <= n.perm){
  #toroidal shift
  data.shifted <- rshift(data_tor, which = as.character(shift.taxa), edge = "torus")
  freq_marks <- table(data.shifted$marks)
  lambda1 <- freq_marks[as.character(base.taxa)]
  lambda2 <- freq_marks[as.character(shift.taxa)]
  
  if(lambda1 == 0 || lambda2 == 0){
    next  # Skip this iteration, don't increment perm
  }
  
  obj_K_shifted <- compute_K(data.shifted, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = lambda1, lambda2 = lambda2, r = r)
  Kcross_stat.tor[perm, ] <- obj_K_shifted$Kcross_stat
  Kstar.tor[perm, ] <- obj_K_shifted$Kstar
  Kcor.tor[perm, ] <- obj_K_shifted$Kcor[-1]
  #non-toroidal shift
  jump_rad <- incircle(data$window)$r
  shift_vector <- runifdisc(1, radius = jump_rad)
  shift_x <- shift_vector$x
  shift_y <- shift_vector$y
  points_to_shift <- subset(data, marks == shift.taxa)
  other_points <- subset(data, marks != shift.taxa)
  shifted_points <- shift(points_to_shift, vec = c(shift_x, shift_y))
  W_shifted <- shift(original_window, vec = c(shift_x, shift_y))
  W.reduced <- intersect.owin(original_window, W_shifted)
  if(is.null(W.reduced) || area.owin(W.reduced) == 0){
    next
  }
  pp_original_reduced <- data[W.reduced]
  pp_shifted_full <- superimpose(other_points, shifted_points)
  pp_shifted_reduced <- pp_shifted_full[W.reduced]
  freq_marks <- table(pp_shifted_reduced$marks)
  lambda1_new <- freq_marks[as.character(base.taxa)]
  lambda2_new <- freq_marks[as.character(shift.taxa)]
  
  if(lambda1_new == 0 || lambda2_new == 0){
    next
  }
  obj_K_minus <- compute_K(pp_shifted_reduced, base.taxa = base.taxa, shift.taxa = shift.taxa, lambda1 = lambda1_new, lambda2 = lambda2_new, r = r)
  Kcross_stat.vc[perm, ] <- obj_K_minus$Kcross_stat
  Kstar.vc[perm, ] <- obj_K_minus$Kstar
  Kcor.vc[perm, ] <- obj_K_minus$Kcor[-1]
  v_perm[,perm] <- c(shift_x, shift_y)
  pp_prob <- table(pp_shifted_reduced$marks)
  nstar_reduced[,perm] <- c(pp_prob[base.taxa] , pp_prob[shift.taxa])
  # Only increment counter if we pass all checks
  perm <- perm + 1
}

ind <- 10 #median value of r
#K0 <- Kcross_stat
Kcross.simulated <- cbind(t(Kcross_stat.vc), Kcross_stat)
Kmean <- apply(Kcross.simulated, 1, mean)
Si <- sweep(t(Kcross_stat.vc), 1, Kmean, "-")
Si <- Si^2

#extract values at ind
Si_ind <- Si[ ind, ]
plot(nstar_reduced[2,], Si_ind)
plot(nstar_reduced[1,], Si_ind)
plot(v_perm[1,], Si_ind)
plot(v_perm[2,], Si_ind)

X <- cbind(nstar_reduced[1,], nstar_reduced[2,], v_perm[1,], v_perm[2,])
sp_X1 <- bs(X[,1], df = 5, degree = 3)
sp_X2 <- bs(X[,2], df = 5, degree = 3)
sp_X3 <- bs(X[,3], df = 5, degree = 3)
sp_X4 <- bs(X[,4], df = 5, degree = 3)
sp_X <- cbind(sp_X1, sp_X2, sp_X3, sp_X4)
obj <- lm(t(Si) ~ sp_X)
Si_hat <- obj$fitted.values

Si.all <- t(newYhat)
S <- Si.all^(-0.5) * Kcross.simulated
S[is.nan(S)] <- 0
CS <- create_curve_set(list(r = r, obs = S[ , n.perm + 1], sim_m = S[ , 1:n.perm]))
pval.Kcross.RNN <- attr(rank_envelope(CS, type = "erl"), "p")

j_points <- subset(data, marks == shift.taxa)
i_points <- subset(data, marks == base.taxa)
#mean(nncross(i_points, j_points, k = 1)$dist)
