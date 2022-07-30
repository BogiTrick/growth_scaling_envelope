remove(list=ls())
library(dplyr)
library(ggplot2)

setwd("/home/bogi/Desktop/growth_rate/data/")


#================ Geometry functions ================
capsule_vol <- function(d1, d2) {
  return((pi*(d1/2)^2)*(d2-d1/3))
}

capsule_surf <- function(d1, d2) {
  return(pi*d1*d2)
}

sphere_vol <- function(d1) {
  return(pi*(d1^3)/6)
}

sphere_surf <- function(d1) {
  return(pi*d1^2)
}

helix_volume <- function(d_cell, len_tot, pitch, d_helix) {
  actual_len <- sqrt(pitch^2 + (pi*(d_helix-d_cell))^2)*(len_tot/pitch)
  return(capsule_vol(d_cell, actual_len))
}

helix_surface <- function(d_cell, len_tot, pitch, d_helix) {
  actual_len <- sqrt(pitch^2 + (pi*(d_helix-d_cell))^2)*(len_tot/pitch)
  return(capsule_surf(d_cell, actual_len))
}

#========== Calculating cell geomtery ==========
df_size <- read.csv("./cell_traits_data/bacterial_envelopes_Sep232021 - cell_dimensions.csv",
                    skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_size <- df_size[,-c(15:ncol(df_size))]

x <- aggregate(df_size[,-c(1,2)], list(df_size[,1]), mean, na.action = na.omit) # Average all linear dimensions by species (1st column)
colnames(x) <- c("species", "D.low", "D.high", "D.mean", "L.low", "L.high", "L.mean", "P.low", "P.high", "P.mean", "Dh.low", "Dh.high", "Dh.mean")
y <- unique(df_size[,c(1,2)]) # Keep only unique entries for shape of the species...
colnames(y) <- c("species", "shape") 
df_size <- merge(y, x, by="species") #...and merge with averaged linear dimensions to get data frame containing both linear dims and shapes


#========== Processing envelope thickness data ==========
in_file_env <- "./cell_traits_data/bacterial_envelopes_Sep232021 - envelope_ultrastructure.csv"
df_env <- read.csv(in_file_env, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_env <- select(df_env, species, Total, Method)
df_env <- aggregate(df_env$Total, list(df_env$species), mean, na.action = na.omit)
colnames(df_env) <- c("species","Total")
df_env <- df_env[!grepl("\\bNitro",df_env$species),]
df_size <- merge(df_size, df_env, by = "species", all.x = TRUE)

# Since the theory states that internal volume of the cell matters as opposed to total volume
# we account for this by modifying all linear dimensions by subtracting 2x the thickness of an E. coli envelope of 30nm
#L_env_generic <- 0.03;
L_env_mollicutes <- 0.004; # From Trachtenberg et al. 2014
L_env_specific <- df_size$Total/1000;

df_size$D.low_corr <- df_size$D.low - 2*L_env_specific
df_size$D.high_corr <- df_size$D.high - 2*L_env_specific
df_size$D.mean_corr <- df_size$D.mean - 2*L_env_specific
df_size$L.low_corr <- df_size$L.low - 2*L_env_specific
df_size$L.high_corr <- df_size$L.high - 2*L_env_specific
df_size$L.mean_corr <- df_size$L.mean - 2*L_env_specific
df_size$Dh.low_corr <- df_size$Dh.low - 2*L_env_specific
df_size$Dh.high_corr <- df_size$Dh.high - 2*L_env_specific
df_size$Dh.mean_corr <- df_size$Dh.mean - 2*L_env_specific

# Subset by cell shape, and calculate volumes and surface areas of capsules and spheres
df_size_rod <- df_size[(df_size$shape=="rod"),]
df_size_sphere <- df_size[(df_size$shape=="sphere"),]
df_size_helix <- df_size[(df_size$shape=="helical"),]

# For each size feature, we calculate lower and upper bound, and mean
lo_vol_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_vol(df_size_rod$D.low_corr[x], df_size_rod$L.low_corr[x]))))
hi_vol_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_vol(df_size_rod$D.high_corr[x], df_size_rod$L.high_corr[x]))))
mean_vol_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_vol(df_size_rod$D.mean_corr[x], df_size_rod$L.mean_corr[x]))))

lo_vol_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) helix_volume(df_size_helix$D.low_corr[x], df_size_helix$L.low_corr[x], df_size_helix$P.low[x], df_size_helix$Dh.low_corr[x]))))
hi_vol_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) helix_volume(df_size_helix$D.high_corr[x], df_size_helix$L.high_corr[x], df_size_helix$P.high[x], df_size_helix$Dh.high_corr[x]))))
mean_vol_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) helix_volume(df_size_helix$D.mean_corr[x], df_size_helix$L.mean_corr[x], df_size_helix$P.mean[x], df_size_helix$Dh.mean_corr[x]))))

lo_vol_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_vol(df_size_sphere$D.low_corr[x]))))
hi_vol_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_vol(df_size_sphere$D.high_corr[x]))))
mean_vol_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_vol(df_size_sphere$D.mean_corr[x]))))

lo_surf_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_surf(df_size_rod$D.low[x], df_size_rod$L.low[x]))))
hi_surf_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_surf(df_size_rod$D.high[x], df_size_rod$L.high[x]))))
mean_surf_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_surf(df_size_rod$D.mean[x], df_size_rod$L.mean[x]))))

lo_surf_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) helix_surface(df_size_helix$D.low[x], df_size_helix$L.low[x], df_size_helix$P.low[x], df_size_helix$Dh.low[x]))))
hi_surf_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) helix_surface(df_size_helix$D.high[x], df_size_helix$L.high[x], df_size_helix$P.high[x], df_size_helix$Dh.high[x]))))
mean_surf_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) helix_surface(df_size_helix$D.mean[x], df_size_helix$L.mean[x], df_size_helix$P.mean[x], df_size_helix$Dh.mean[x]))))

lo_surf_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_surf(df_size_sphere$D.low[x]))))
hi_surf_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_surf(df_size_sphere$D.high[x]))))
mean_surf_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_surf(df_size_sphere$D.mean[x]))))

# Update rods
df_size_rod$V.low <- lo_vol_rod
df_size_rod$V.high <- hi_vol_rod
df_size_rod$V.mean <- mean_vol_rod

df_size_rod$S.low <- lo_surf_rod
df_size_rod$S.high <- hi_surf_rod
df_size_rod$S.mean <- mean_surf_rod

# Update spheres
df_size_sphere$V.low <- lo_vol_sphere
df_size_sphere$V.high <- hi_vol_sphere
df_size_sphere$V.mean <- mean_vol_sphere

df_size_sphere$S.low <- lo_surf_sphere
df_size_sphere$S.high <- hi_surf_sphere
df_size_sphere$S.mean <- mean_surf_sphere

# Update helices
df_size_helix$V.low <- lo_vol_hel
df_size_helix$V.high <- hi_vol_hel
df_size_helix$V.mean <- mean_vol_hel

df_size_helix$S.low <- lo_surf_hel
df_size_helix$S.high <- hi_surf_hel
df_size_helix$S.mean <- mean_surf_hel

# Bind together
df_size <- rbind(df_size_rod, df_size_sphere, df_size_helix)

#========== Processing growth rate data ==========
in_file_growth <- "./cell_traits_data/bacterial_envelopes_Sep232021 - growth_rates.csv"
df_growth <- read.csv(in_file_growth, header = TRUE, sep = ",", stringsAsFactors = FALSE)
in_file_genomics <- "./cell_traits_data/bacterial_envelopes_Sep232021 - genomic_properties.csv"
df_genomics <- read.csv(in_file_genomics, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_growth <- cbind(df_growth, df_genomics[,-1])
df_growth <- df_growth %>% group_by(species) %>% top_n(1, Q10.correction)
df_growth <- select(df_growth, species, Q10.correction, rRNA.genes, tRNA.genes)
colnames(df_growth) <- c("species","growth_rate","rRNA.genes","tRNA.genes")

df1 <- merge(df_size, df_growth, by="species", all.y = TRUE)
df1 <- select(df1, species, shape, V.mean, S.mean, growth_rate, Total, rRNA.genes, tRNA.genes)
colnames(df1) <- c("species", "shape", "V.mean","S.mean", "growth_rate", "total_thickness_nm", "rRNA genes", "tRNA genes")

in_file_metabolism <- "./cell_traits_data/bacterial_envelopes_Sep232021 - physiological_properties.csv"
df_metabolism <- read.csv(in_file_metabolism, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_metabolism <- select(df_metabolism, species, tca_steps)
df1 <- merge(df1, df_metabolism, all.x = TRUE)

write.csv(df1, "./main/growth_scaling.shape.envelope_corrected.csv", row.names = FALSE)

sum(duplicated(df1$species))
