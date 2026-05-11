##all plots
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
M = 4
data <- list()
data[[1]] <- readRDS("simulation_results_irregW_1_power13.rds")
data[[2]] <- readRDS("simulation_results_irregW_3_power13.rds")
data[[3]] <- readRDS("simulation_results_irregW_5_power13.rds")



results_list <- list(
  sigma_0.01 = list(Kstar_torshift = data[[1]]$Kstar_torshift, Kcross.vc = data[[1]]$Kcross.vc.gauss, Kcor.vc = data[[1]]$Kcor.vc.gauss, NN.vc = data[[1]]$NN.vc.gauss),
  sigma_0.03 = list(Kstar_torshift = data[[2]]$Kstar_torshift, Kcross.vc = data[[2]]$Kcross.vc.gauss, Kcor.vc = data[[2]]$Kcor.vc.gauss, NN.vc = data[[2]]$NN.vc.gauss),
  sigma_0.05 = list(Kstar_torshift = data[[3]]$Kstar_torshift, Kcross.vc = data[[3]]$Kcross.vc.gauss, Kcor.vc = data[[3]]$Kcor.vc.gauss, NN.vc = data[[3]]$NN.vc.gauss)
  
)
convert_to_long <- function(results_list, sigma_values) {
  
  all_data <- map2_dfr(results_list, sigma_values, function(sigma_list, sigma_val) {
    
    map_dfr(names(sigma_list), function(method_name) {
      mat <- sigma_list[[method_name]]
      n <- nrow(mat)
      
      # Extract upper triangular indices
      upper_tri_indices <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
      
      if(nrow(upper_tri_indices) > 0) {
        tibble(
          row = upper_tri_indices[, 1],
          col = upper_tri_indices[, 2],
          value = mat[upper_tri_indices],
          method = method_name,
          sigma = sigma_val
        )
      } else {
        tibble()
      }
    })
  })
  
  return(all_data)
}
sigma_values <- c(0.01, 0.03, 0.05)
long_data <- convert_to_long(results_list, sigma_values)
long_data <- long_data %>%
  mutate(
    facet_label = paste0("(", row, ",", col, ")"),
    taxa_i = factor(row),
    taxa_j = factor(col)
  )
annotation_df <- data.frame(
  Sigma = Inf,  # or a specific x value
  `Empirical Rejection Rate` = Inf,  # or a specific y value
  taxa_j = "3",  # This should match the factor level for taxa_j: 3
  taxa_i = "1",  # This should match the factor level for taxa_i: 1
  label = "alternative"
)

# Step 4: Plot with facet_grid
p <- ggplot(long_data, aes(x = sigma, y = value, color = method, group = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(taxa_i ~ taxa_j, 
             scales = "free_y",
             labeller = label_both) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines")
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(
    x = "Sigma",
    y = "Empirical Rejection Rate",
    color = "Method",
    title = "Empirical rejection under null and alternative (1,3)"
  ) +
   geom_text(data = annotation_df,
              aes(x = Sigma, y = `Empirical.Rejection.Rate`, label = label),
              hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE)+
  scale_color_brewer(palette = "Set1")

print(p)
