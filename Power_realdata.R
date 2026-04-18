##Power, based on real data
##read a slide, calculate p values for 100 samples of a specific depth
library(RColorBrewer)
library(Polychrome)
library(corrplot)
library(ggcorrplot)
library(gridGraphics)
rm(list = ls())
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/generate_data.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/pval_v2.R", echo = TRUE)

process_df <- function(df, tile_num=NULL){
  if (is.null(tile_num)) tile <- df
  else
    tile <- df[df$tile==tile_num,]
  
  #coords <- gsub("\\[|\\]", "", tile$coord)  # Remove brackets
  coords <- gsub("\\[|\\]|\\(|\\)", "", tile$coord)
  split_coords <- strsplit(coords, ",")  # Split into x and y
  
  x <- as.numeric(sapply(split_coords, `[`, 1))  # Extract x values
  y <- as.numeric(sapply(split_coords, `[`, 2))  # Extract y values
  
  # adjust x,y coordinates
  x <- x-min(x)
  y <- y-min(y)
  
  tile.df <- data.frame(x=x,y=y,sciname=tile$sciname, tile = tile$tile)
  return(tile.df)
}

muco.list <- c("2023_02_08_hsdm_group_2_sample_06_fov_01_centroid_sciname.csv",            
               "2023_02_08_hsdm_group_2_sample_06_fov_02_centroid_sciname.csv" ,           
               "2023_02_18_hsdm_group_II_patient_6_fov_01_centroid_sciname.csv"  ,         
               "2023_10_16_hsdm_slide_IIB_fov_01_centroid_sciname.csv",                    
               "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv" ,                   
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_01_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_02_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_03_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_04_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_05_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_13_aspect_MB_fov_06_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_14_aspect_MB_fov_01_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_14_aspect_MB_fov_02_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_15_aspect_MB_fov_01_centroid_sciname.csv",
               "2024_04_19_hsdm_group_II_patient_15_aspect_MB_fov_02_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_01_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_02_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_03_centroid_sciname.csv",
               "2024_04_27_hsdm_group_II_patient_11_aspect_DL_fov_04_centroid_sciname.csv")

healthy.list <- c( "2023_02_08_hsdm_group_1_sample_06_fov_01_centroid_sciname.csv",           
                   "2023_02_08_hsdm_group_1_sample_11_fov_01_centroid_sciname.csv",           
                   "2023_02_08_hsdm_group_1_sample_12_fov_01_centroid_sciname.csv",           
                   "2023_02_18_hsdm_group_I_patient_11_fov_01_centroid_sciname.csv",          
                   "2023_02_18_hsdm_group_I_patient_11_fov_02_centroid_sciname.csv",          
                   "2023_02_18_hsdm_group_I_patient_13_fov_01_centroid_sciname.csv",          
                   "2023_02_18_hsdm_group_I_patient_6_fov_01_centroid_sciname.csv",           
                   "2023_10_16_hsdm_slide_IL_fov_01_centroid_sciname.csv",                    
                   "2023_10_16_hsdm_slide_IL_fov_02_centroid_sciname.csv",                    
                   "2023_10_16_hsdm_slide_IL_fov_03_centroid_sciname.csv",                    
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_01_centroid_sciname.csv",
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_02_centroid_sciname.csv",
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_03_centroid_sciname.csv",
                   "2024_04_24_hsdm_group_I_patient_16_aspect_MB_fov_04_centroid_sciname.csv")
path <- muco.list[5]
img <- read.csv(file=paste0("centroid_sciname_tables/",path))
img.df <- process_df(img)
ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
            marks=as.factor(img.df$sciname))
ppp_df <- data.frame(ppp_img)

split_ppp <- split.ppp(ppp_img, f = marks(ppp_img))




  load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome4/spatial_microbiome4/windows/", path, ".RData"))
X_sampled_list <- lapply(split_ppp, function(pp) {
    n <- npoints(pp)
    n_sample <- floor(p * n)
    indices <- sample(1:n, size = n_sample)
    pp[indices]
  })
  
  # Combine back into multitype ppp
  X_sampled <- do.call(superimpose, X_sampled_list)
  X_sampled$marks <- factor(X_sampled$marks[,"origMarks"], levels = levels(ppp_img$marks))
  X_sampled$window <- W
  
  all_marks <- levels(X_sampled$marks)
  rmax <- incircle(W)$r/5
  r <- seq(0, rmax, length.out = 50)
  bw <- rmax
  p = 0.1

  M = length(all_marks)
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
    results <- mclapply(1:100, function(i){
      
      X_sampled_list <- lapply(split_ppp, function(pp) {
        n <- npoints(pp)
        n_sample <- floor(p * n)
        indices <- sample(1:n, size = n_sample)
        pp[indices]
      })
      
      # Combine back into multitype ppp
      X_sampled <- do.call(superimpose, X_sampled_list)
      X_sampled$marks <- factor(X_sampled$marks[,"origMarks"], levels = levels(ppp_img$marks))
      X_sampled$window <- W
      pvals <- suppressWarnings(pval.assoc(X_sampled, base.taxa = all_marks[base.taxa], shift.taxa = all_marks[shift.taxa],r = r, n.perm = 199, bw = bw))
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
  
colnames(Kcross_torshift) <- all_marks
colnames(Kstar_torshift) <- all_marks
colnames(Kcor_torshift) <- all_marks
colnames(NN_torshift) <- all_marks

colnames(Kcross.vc.ep) <- all_marks
colnames(Kcross.vc.gauss) <- all_marks
colnames(Kcross.vc.uniform) <- all_marks
colnames(Kstar.vc.ep) <- all_marks
colnames(Kstar.vc.gauss) <- all_marks
colnames(Kstar.vc.uniform) <- all_marks
colnames(Kcor.vc.ep) <- all_marks
colnames(Kcor.vc.gauss) <- all_marks
colnames(Kcor.vc.uniform) <- all_marks
colnames(NN.vc.ep) <- all_marks
colnames(NN.vc.gauss) <- all_marks
colnames(NN.vc.uniform) <- all_marks

rownames(Kcross_torshift) <- all_marks
rownames(Kstar_torshift) <- all_marks
rownames(Kcor_torshift) <- all_marks
rownames(NN_torshift) <- all_marks
rownames(Kcross.vc.ep) <- all_marks
rownames(Kcross.vc.gauss) <- all_marks
rownames(Kcross.vc.uniform) <- all_marks
rownames(Kstar.vc.ep) <- all_marks
rownames(Kstar.vc.gauss) <- all_marks
rownames(Kstar.vc.uniform) <- all_marks
rownames(Kcor.vc.ep) <- all_marks
rownames(Kcor.vc.gauss) <- all_marks
rownames(Kcor.vc.uniform) <- all_marks
rownames(NN.vc.ep) <- all_marks
rownames(NN.vc.gauss) <- all_marks
rownames(NN.vc.uniform) <- all_marks




  saveRDS(list(Kcross_torshift = Kcross_torshift, Kstar_torshift = Kstar_torshift, Kcor_torshift = Kcor_torshift, NN_torshift = NN_torshift,
               Kcross.vc.ep = Kcross.vc.ep, Kcross.vc.gauss = Kcross.vc.gauss, Kcross.vc.uniform = Kcross.vc.uniform,
               Kstar.vc.ep = Kstar.vc.ep, Kstar.vc.gauss = Kstar.vc.gauss, Kstar.vc.uniform = Kstar.vc.uniform,
               Kcor.vc.ep = Kcor.vc.ep, Kcor.vc.gauss = Kcor.vc.gauss, Kcor.vc.uniform = Kcor.vc.uniform,
               NN.vc.ep = NN.vc.ep, NN.vc.gauss = NN.vc.gauss, NN.vc.uniform = NN.vc.uniform), file = "Power_real_p1.rds")
  beepr::beep(4)
  
  
  df <- data.frame(X_sampled)
  col_pallette <- kelly.colors(n = 18)
  col_pallette[1] <- "#e25"
  names(col_pallette) <- levels(df$marks)
  p1 <- ggplot(df, aes(x = x, y= y, color = marks, shape = marks)) + geom_point(alpha = 2) + scale_color_manual(values = col_pallette) + scale_shape_manual(values = rep(15:20,3))
 library(pheatmap)
  library(cowplot)
  library(ggplotify)
  p2 <- pheatmap(NN.vc.gauss,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("white", "red"))(100),
           display_numbers = TRUE,
           na_col = "black",# Set TRUE if you want values shown
           border_color = "grey90",
           main = "NN.vc.gauss")
p2_gg <- as.ggplot(p2[[4]])
p3 <- pheatmap(Kcross.vc.gauss,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(100),
         display_numbers = TRUE,
         na_col = "black",# Set TRUE if you want values shown
         border_color = "grey90",
         main = "Kcross.vc.gauss")
p3_gg <- as.ggplot(p3[[4]])
p4 <- pheatmap(Kcor.vc.gauss,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(100),
         display_numbers = TRUE,
         na_col = "black",# Set TRUE if you want values shown
         border_color = "grey90",
         main = "Kcor.vc.gauss")
p4_gg <- as.ggplot(p4[[4]])
library(gridExtra)
plot_grid(p1, p2_gg, p3_gg, p4_gg,
          labels = c("A)", "B)", "C)", "D)"),
          ncol = 2)
