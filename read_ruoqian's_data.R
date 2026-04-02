##read Ruoqian's data
library(spatstat)
library(ggplot2)
setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome3")

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


# load the alpha hull window for downstream analysis
path <- "2023_10_18_hsdm_slide_IIL_fov_01_centroid_sciname.csv"
load(paste0("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome3/windows/", path, ".RData"))
#list.files("/Users/asmitaroy/Library/CloudStorage/OneDrive-JohnsHopkins/Spatial_microbiome/biofilm_Ben/windows")
# sample use, construct a ppp object using window W for analysis W = c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y))
plot_muco = list()

for(i in seq_along(muco.list)){
  path <- muco.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=img.df$sciname)
  ppp_df <- data.frame(ppp_img)
  df_muco[[i]] <- ppp_df
  plot_muco[[i]] <- ggplot(img.df, aes(x = x, y = y, color = sciname)) +
    geom_point(size = 0.1) +
    coord_equal() + # Important for spatial data
    ggtitle(path)+
    theme_minimal()
  ggsave(paste0("plot_muco/",substr(muco.list[i], 1, nchar(muco.list[i])-4 ),".png"),plot_muco[[i]] )
}

plot.healthy = list()
for(i in seq_along(healthy.list)){
  path <- healthy.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=img.df$sciname)
  ppp_df <- data.frame(ppp_img)
  plot.healthy[[i]] <- ggplot(img.df, aes(x = x, y = y, color = sciname)) +
    geom_point(size = 0.1) +
    coord_equal() + # Important for spatial data
    ggtitle(path)+
    theme_minimal()
  ggsave(paste0("plot_healthy/",substr(healthy.list[i], 1, nchar(healthy.list[i])-4 ),".png"),plot.healthy[[i]] )
}




##Read all data into list
df_muco <- list()
for(i in seq_along(muco.list)){
  path <- muco.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=img.df$sciname)
  ppp_df <- data.frame(ppp_img)
  df_muco[[i]] <- ppp_df
}
df_healthy <- list()
for(i in seq_along(healthy.list)){
  path <- healthy.list[i]
  img <- read.csv(file=paste0("centroid_sciname_tables/",path))
  img.df <- process_df(img)
  ppp_img=ppp(x=img.df$x, y=img.df$y, c(min(img.df$x), max(img.df$x)), c(min(img.df$y), max(img.df$y)),
              marks=img.df$sciname)
  ppp_df <- data.frame(ppp_img)
  df_healthy[[i]] <- ppp_df
}
nnct_list <- list()
M_muco <- lapply(df_muco, function(x)length(table(x[,3])))
for(i in seq_along(df_muco)){
  nnct_list[[i]] <- fast_nnct(all_coords = df_muco[[i]][,1:2], type = df_muco[[i]][,3], M = M_muco[[i]])$NNCT
  print(i)
}

