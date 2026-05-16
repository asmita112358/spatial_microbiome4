source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/stats.R", echo = FALSE)

# Kernel weight function
kernel_weights <- function(x, bw, type = "epanechnikov") {
  if (type == "epanechnikov") {
    wts <- (3/4) * (1/bw) * (1 - (x/bw)^2) * (abs(x) <= bw)
  } else if (type == "gaussian") {
    wts <- dnorm(x, mean = 0, sd = bw)
  } else if (type == "uniform") {
    wts <- 0.5 * (1/bw) * (abs(x) <= bw)
  } else {
    stop("Unsupported kernel type. Choose 'epanechnikov', 'gaussian', or 'uniform'.")
  }
  return(wts)
}

# Trimmed mean excluding extreme simulations based on sample sizes
trimmed_mean_by_sample_size <- function(K_matrix, sample_sizes, trim = 0.05) {
  # K_matrix: r_length × n_sim matrix of K-function values
  # sample_sizes: 2 × n_sim matrix (row 1: base taxa counts, row 2: shift taxa counts)
  # trim: proportion to trim from each tail
  
  if (ncol(K_matrix) != ncol(sample_sizes)) {
    stop("K_matrix and sample_sizes must have the same number of columns")
  }
  
  n_sim <- ncol(K_matrix)
  n_trim <- floor(n_sim * trim)
  
  if (n_trim * 4 >= n_sim) {
    warning(sprintf("Trim proportion may be too large. Trimming %d out of %d simulations.",
                    4 * n_trim, n_sim))
  }
  
  # Identify extreme indices based on sample sizes
  extreme_idx <- unique(c(
    order(sample_sizes[1, ])[1:n_trim],                      # lowest base taxa counts
    order(sample_sizes[1, ], decreasing = TRUE)[1:n_trim],   # highest base taxa counts
    order(sample_sizes[2, ])[1:n_trim],                      # lowest shift taxa counts
    order(sample_sizes[2, ], decreasing = TRUE)[1:n_trim]    # highest shift taxa counts
  ))
  
  keep_idx <- setdiff(1:n_sim, extreme_idx)
  
  if (length(keep_idx) == 0) {
    stop("All simulations were trimmed. Reduce trim proportion.")
  }
  
  return(rowMeans(K_matrix[, keep_idx, drop = FALSE]))
}

# Nadaraya-Watson kernel regression
nadaraya_watson <- function(X, Y, X_pred = X, bw) {
  n <- nrow(X)
  m <- nrow(X_pred)
  
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  n_response <- ncol(Y)
  predictions <- matrix(0, nrow = m, ncol = n_response)
  
  # Gaussian kernel
  gaussian_kernel <- function(dist_sq, bandwidth) {
    exp(-0.5 * dist_sq / bandwidth^2)
  }
  
  for (i in seq_len(m)) {
    # Squared Euclidean distances
    dist_sq <- rowSums((X - matrix(X_pred[i, ], nrow = n, ncol = ncol(X), byrow = TRUE))^2)
    weights <- gaussian_kernel(dist_sq, bw)
    sum_weights <- sum(weights)
    
    if (sum_weights > 1e-10) {
      predictions[i, ] <- crossprod(weights, Y) / sum_weights
    } else {
      predictions[i, ] <- NA
    }
  }
  
  if (n_response == 1) {
    predictions <- as.vector(predictions)
  }
  
  return(predictions)
}

# Main function for association testing
test_spatial_association <- function(data, base_taxa = 1, shift_taxa = 2, r, 
                                     n_perm = 199, bw = "silverman", 
                                     median_index = 25) {
  
  n_dist <- length(r)
  freq_marks <- table(data$marks)
  original_window <- data$window
  n_base <- freq_marks[as.character(base_taxa)]
  n_shift <- freq_marks[as.character(shift_taxa)]
  
  # Compute observed statistics
  obs_K <- compute_K(data, base.taxa = base_taxa, shift.taxa = shift_taxa, 
                     lambda1 = n_base, lambda2 = n_shift, r = r)
  
  obs_Kcross <- obs_K$Kcross_stat
  obs_Kstar <- obs_K$Kstar
  obs_Kcor <- obs_K$Kcor[-1]
  r_cor <- r[-1]
  
  obs_auc <- trapz(r, obs_Kcross)
  obs_median <- obs_Kcross[median_index]
  
  # Initialize storage for permutation results
  Kcross_toroidal <- matrix(NA, nrow = n_perm, ncol = n_dist)
  Kcross_variance_corrected <- matrix(NA, nrow = n_perm, ncol = n_dist)
  
  Kstar_toroidal <- matrix(NA, nrow = n_perm, ncol = n_dist)
  Kstar_variance_corrected <- matrix(NA, nrow = n_perm, ncol = n_dist)
  
  Kcor_toroidal <- matrix(NA, nrow = n_perm, ncol = n_dist - 1)
  Kcor_variance_corrected <- matrix(NA, nrow = n_perm, ncol = n_dist - 1)
  
  auc_toroidal <- numeric(n_perm)
  auc_vc <- numeric(n_perm)
  median_toroidal <- numeric(n_perm)
  median_vc <- numeric(n_perm)
  
  shift_vectors <- matrix(nrow = 2, ncol = n_perm)
  sample_sizes_vc <- matrix(NA, nrow = 2, ncol = n_perm)
  
  # Prepare toroidal shift data
  data_toroidal <- data
  if (!is.rectangle(data$win)) {
    Window(data_toroidal) <- boundingbox(data$win)
  }
  
  perm_idx <- 1
  
  while (perm_idx <= n_perm) {
    # Toroidal shift
    data_shifted_tor <- rshift(data_toroidal, which = as.character(shift_taxa), edge = "torus")
    freq_tor <- table(data_shifted_tor$marks)
    n_base_tor <- freq_tor[as.character(base_taxa)]
    n_shift_tor <- freq_tor[as.character(shift_taxa)]
    
    if (n_base_tor == 0 || n_shift_tor == 0) next
    
    K_tor <- compute_K(data_shifted_tor, base.taxa = base_taxa, shift.taxa = shift_taxa,
                       lambda1 = n_base_tor, lambda2 = n_shift_tor, r = r)
    
    Kcross_toroidal[perm_idx, ] <- K_tor$Kcross_stat
    Kstar_toroidal[perm_idx, ] <- K_tor$Kstar
    Kcor_toroidal[perm_idx, ] <- K_tor$Kcor[-1]
    auc_toroidal[perm_idx] <- trapz(r, K_tor$Kcross_stat)
    median_toroidal[perm_idx] <- K_tor$Kcross_stat[median_index]
    
    # Non-toroidal shift (variance correction)
    jump_radius <- incircle(data$window)$r
    shift_vector <- runifdisc(1, radius = jump_radius)
    shift_x <- shift_vector$x
    shift_y <- shift_vector$y
    
    points_shift <- subset(data, marks == shift_taxa)
    points_other <- subset(data, marks != shift_taxa)
    points_shifted <- shift(points_shift, vec = c(shift_x, shift_y))
    
    window_shifted <- shift(original_window, vec = c(shift_x, shift_y))
    window_reduced <- intersect.owin(original_window, window_shifted)
    
    if (is.null(window_reduced) || area.owin(window_reduced) == 0) next
    
    pp_combined <- superimpose(points_other, points_shifted)
    pp_reduced <- pp_combined[window_reduced]
    
    freq_vc <- table(pp_reduced$marks)
    n_base_vc <- freq_vc[as.character(base_taxa)]
    n_shift_vc <- freq_vc[as.character(shift_taxa)]
    
    if (n_base_vc == 0 || n_shift_vc == 0) next
    
    K_vc <- compute_K(pp_reduced, base.taxa = base_taxa, shift.taxa = shift_taxa,
                      lambda1 = n_base_vc, lambda2 = n_shift_vc, r = r)
    
    Kcross_variance_corrected[perm_idx, ] <- K_vc$Kcross_stat
    Kstar_variance_corrected[perm_idx, ] <- K_vc$Kstar
    Kcor_variance_corrected[perm_idx, ] <- K_vc$Kcor[-1]
    auc_vc[perm_idx] <- trapz(r, K_vc$Kcross_stat)
    median_vc[perm_idx] <- K_vc$Kcross_stat[median_index]
    
    shift_vectors[, perm_idx] <- c(shift_x, shift_y)
    sample_sizes_vc[, perm_idx] <- c(n_base_vc, n_shift_vc)
    
    perm_idx <- perm_idx + 1
  }
  
  # P-values for toroidal shift (AUC and median)
  pval_auc_tor <- (1 + sum(auc_toroidal >= obs_auc)) / (1 + n_perm)
  pval_median_tor <- (1 + sum(median_toroidal >= obs_median)) / (1 + n_perm)
  
  # Global envelope tests for toroidal shift
  pval_Kcross_tor <- compute_envelope_pval(r, obs_Kcross, Kcross_toroidal)
  pval_Kstar_tor <- compute_envelope_pval(r, obs_Kstar, Kstar_toroidal)
  pval_Kcor_tor <- compute_envelope_pval(r_cor, obs_Kcor, Kcor_toroidal)
  
  # Variance correction with bandwidth selection
  shift_vectors_full <- cbind(shift_vectors, c(0, 0))
  sample_sizes_full <- cbind(sample_sizes_vc, c(n_base, n_shift))
  
  bw_shift <- select_bandwidth(shift_vectors_full, method = bw)
  bw_sample <- select_bandwidth(sample_sizes_full, method = bw)
  
  # Variance-corrected tests
  vc_results <- compute_variance_corrected_tests(
    Kcross_variance_corrected, obs_Kcross, shift_vectors_full, 
    sample_sizes_full, bw_shift, bw_sample, r, n_perm
  )
  
  # AUC and median variance-corrected tests
  auc_vc_results <- compute_scalar_vc_test(auc_vc, obs_auc, shift_vectors_full, bw_shift, n_perm)
  median_vc_results <- compute_scalar_vc_test(median_vc, obs_median, shift_vectors_full, bw_shift, n_perm)
  
  # Simple minus tests (no variance correction)
  pval_auc_minus <- (1 + sum(auc_vc >= obs_auc)) / (1 + n_perm)
  pval_median_minus <- (1 + sum(median_vc >= obs_median)) / (1 + n_perm)
  
  Kcross_minus_matrix <- cbind(t(Kcross_variance_corrected), obs_Kcross)
  CS_minus <- create_curve_set(list(r = r, obs = Kcross_minus_matrix[, n_perm + 1],
                                    sim_m = Kcross_minus_matrix[, 1:n_perm]))
  pval_Kcross_minus <- attr(rank_envelope(CS_minus, type = "erl"), "p")
  
  return(list(
    # Toroidal shift p-values
    pval_Kcross_tor = pval_Kcross_tor,
    pval_Kstar_tor = pval_Kstar_tor,
    pval_Kcor_tor = pval_Kcor_tor,
    pval_auc_tor = pval_auc_tor,
    pval_median_tor = pval_median_tor,
    
    # Variance-corrected p-values (shift-based)
    pval_Kcross_vc_shift = vc_results$pval_vc_shift,
    pval_Kcross_evc_shift = vc_results$pval_evc_shift,
    
    # Variance-corrected p-values (sample-size-based)
    pval_Kcross_vc_sample = vc_results$pval_vc_sample,
    pval_Kcross_evc_sample = vc_results$pval_evc_sample,
    
    # Minus p-values (no variance correction)
    pval_Kcross_minus = pval_Kcross_minus,
    pval_auc_minus = pval_auc_minus,
    pval_median_minus = pval_median_minus,
    
    # AUC and median variance-corrected
    pval_auc_vc = auc_vc_results$pval,
    pval_median_vc = median_vc_results$pval
  ))
}

# Helper function: Compute envelope test p-value
compute_envelope_pval <- function(r, obs, sim_matrix) {
  df <- data.frame(r = r, obs = obs, t(sim_matrix))
  clean_idx <- apply(df, 1, function(x) all(is.finite(x)))
  df <- df[clean_idx, ]
  
  if (nrow(df) == 0) {
    return(NA)
  }
  
  cset <- curve_set(obs = df$obs, sim = as.matrix(df[, -c(1, 2)]), r = df$r)
  GET_result <- global_envelope_test(cset, typeone = "fwer", type = "rank",
                                     alpha = 0.05, alternative = "greater")
  return(attr(GET_result, "p"))
}

# Helper function: Select bandwidth
select_bandwidth <- function(X, method = "silverman") {
  n_sim <- ncol(X)
  d <- nrow(X)
  
  sd_X <- apply(X, 1, sd)
  sigma_X <- exp(mean(log(sd_X)))
  
  if (method == "silverman") {
    bw <- sigma_X * (n_sim * (d + 2) / 4)^(-1 / (d + 4))
  } else if (method == "scott") {
    bw <- sigma_X * n_sim^(-1 / (d + 4))
  } else {
    bw <- method  # Assume numeric bandwidth provided
  }
  
  return(bw)
}

# Helper function: Compute variance-corrected tests
compute_variance_corrected_tests <- function(K_vc_matrix, obs_K, shift_vectors, 
                                             sample_sizes, bw_shift, bw_sample, r, n_perm) {
  
  K_full <- cbind(t(K_vc_matrix), obs_K)
  K_mean <- rowMeans(K_full)
  K_centered <- sweep(K_full, 1, K_mean, "-")
  
  # Shift-based variance correction
  var_shift <- t(nadaraya_watson(X = t(shift_vectors), Y = t(K_centered^2),
                                 X_pred = t(shift_vectors), bw = bw_shift))
  var_shift[var_shift <= 1e-10] <- 1e-10
  K_std_shift <- K_centered / sqrt(var_shift)
  K_std_shift[!is.finite(K_std_shift)] <- 0
  
  CS_shift <- create_curve_set(list(r = r, obs = K_std_shift[, n_perm + 1],
                                    sim_m = K_std_shift[, 1:n_perm]))
  pval_vc_shift <- attr(rank_envelope(CS_shift, type = "erl"), "p")
  
  # Shift-based expectation-variance correction
  E_K_shift <- t(nadaraya_watson(X = t(shift_vectors), Y = t(K_full),
                                 X_pred = t(shift_vectors), bw = bw_shift))
  E2_K_shift <- t(nadaraya_watson(X = t(shift_vectors), Y = t(K_full^2),
                                  X_pred = t(shift_vectors), bw = bw_shift))
  var_K_shift <- E2_K_shift - E_K_shift^2
  var_K_shift[var_K_shift <= 1e-10] <- 1e-10
  K_std_evc_shift <- (K_full - E_K_shift) / sqrt(var_K_shift)
  K_std_evc_shift[!is.finite(K_std_evc_shift)] <- 0
  
  CS_evc_shift <- create_curve_set(list(r = r, obs = K_std_evc_shift[, n_perm + 1],
                                        sim_m = K_std_evc_shift[, 1:n_perm]))
  pval_evc_shift <- attr(rank_envelope(CS_evc_shift, type = "erl"), "p")
  
  # Sample-size-based variance correction
  var_sample <- t(nadaraya_watson(X = t(sample_sizes), Y = t(K_centered^2),
                                  X_pred = t(sample_sizes), bw = bw_sample))
  var_sample[var_sample <= 1e-10] <- 1e-10
  K_std_sample <- K_centered / sqrt(var_sample)
  K_std_sample[!is.finite(K_std_sample)] <- 0
  
  CS_sample <- create_curve_set(list(r = r, obs = K_std_sample[, n_perm + 1],
                                     sim_m = K_std_sample[, 1:n_perm]))
  pval_vc_sample <- attr(rank_envelope(CS_sample, type = "erl"), "p")
  
  # Sample-size-based expectation-variance correction
  E_K_sample <- t(nadaraya_watson(X = t(sample_sizes), Y = t(K_full),
                                  X_pred = t(sample_sizes), bw = bw_sample))
  E2_K_sample <- t(nadaraya_watson(X = t(sample_sizes), Y = t(K_full^2),
                                   X_pred = t(sample_sizes), bw = bw_sample))
  var_K_sample <- E2_K_sample - E_K_sample^2
  var_K_sample[var_K_sample <= 1e-10] <- 1e-10
  K_std_evc_sample <- (K_full - E_K_sample) / sqrt(var_K_sample)
  K_std_evc_sample[!is.finite(K_std_evc_sample)] <- 0
  
  CS_evc_sample <- create_curve_set(list(r = r, obs = K_std_evc_sample[, n_perm + 1],
                                         sim_m = K_std_evc_sample[, 1:n_perm]))
  pval_evc_sample <- attr(rank_envelope(CS_evc_sample, type = "erl"), "p")
  
  return(list(
    pval_vc_shift = pval_vc_shift,
    pval_evc_shift = pval_evc_shift,
    pval_vc_sample = pval_vc_sample,
    pval_evc_sample = pval_evc_sample
  ))
}

# Helper function: Compute scalar variance-corrected test
compute_scalar_vc_test <- function(sim_values, obs_value, shift_vectors, bw, n_perm) {
  values_full <- c(sim_values, obs_value)
  mean_value <- mean(values_full)
  centered_values <- values_full - mean_value
  
  var_values <- nadaraya_watson(X = t(shift_vectors), Y = centered_values^2,
                                X_pred = t(shift_vectors), bw = bw)
  var_values[var_values <= 1e-10] <- 1e-10
  std_values <- centered_values / sqrt(var_values)
  
  pval <- (1 + sum(std_values[1:n_perm] >= std_values[n_perm + 1])) / (1 + n_perm)
  
  return(list(pval = pval))
}