#=============== USAGE =====================
# Input: 
# 1. ./cell_traits_data/bacterial_envelopes_Sep232021 - {growth, size, other properties}.csv
#
# Output:
# 1. growth_scaling.shape.genome.csv

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

helix_uncoiled_length <- function(d_cell, len_tot, pitch, d_helix) {
  return(sqrt(pitch^2 + (pi*(d_helix-d_cell))^2)*(len_tot/pitch))
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
# IMPORTANT NOTE: For S/V_{tot} and V_{tot}, set both L_env parameters to zero. This ensures that full linear dimensions are
#                 taken for both Mollicutes and the rest of bacterial species
L_env_generic <- 0.03;
L_env_mollicutes <- 0.004; # From Trachtenberg et al. 2014
#L_env_specific <- df_size$Total/1000;

mollicutes_set <- (grepl("Spiroplasma|Mycoplasma|Phytoplasma|Mesoplasma",df_size$species));

# Calculate the length of an uncoiled helical cells
# L_uncoiled -- the length of uncoiled helical cell, and will be used to calculate its internal volume
df_size$L_uncoiled.low <- NA
df_size$L_uncoiled.high <- NA
df_size$L_uncoiled.mean <- NA
df_size[df_size$shape=="helical",]$L_uncoiled.low <- with(df_size[df_size$shape=="helical",], helix_uncoiled_length(D.low, L.low, P.low, Dh.low))
df_size[df_size$shape=="helical",]$L_uncoiled.high <- with(df_size[df_size$shape=="helical",], helix_uncoiled_length(D.high, L.high, P.high, Dh.high))
df_size[df_size$shape=="helical",]$L_uncoiled.mean <- with(df_size[df_size$shape=="helical",], helix_uncoiled_length(D.mean, L.mean, P.mean, Dh.mean))

# Calculate the inner diameter, length, and uncoiled length (in the case of helical species)
# *_corr -- linear dimensions of the cytosol (i.e., linear dimensions corrected for envelope thickness)
df_size$D.low_corr <- df_size$D.low - 2*L_env_generic
df_size$D.high_corr <- df_size$D.high - 2*L_env_generic
df_size$D.mean_corr <- df_size$D.mean - 2*L_env_generic
df_size$L.low_corr <- df_size$L.low - 2*L_env_generic
df_size$L.high_corr <- df_size$L.high - 2*L_env_generic
df_size$L.mean_corr <- df_size$L.mean - 2*L_env_generic
df_size$L_uncoiled.low_corr <- df_size$L_uncoiled.low - 2*L_env_generic
df_size$L_uncoiled.high_corr <- df_size$L_uncoiled.high - 2*L_env_generic
df_size$L_uncoiled.mean_corr <- df_size$L_uncoiled.mean - 2*L_env_generic

# Repeat the same calculations for mollicutes, but use thinner envelope given that mollicutes only have a plasma membrane
df_size[mollicutes_set,]$D.low_corr <- df_size[mollicutes_set,]$D.low - 2*L_env_mollicutes
df_size[mollicutes_set,]$D.high_corr <- df_size[mollicutes_set,]$D.high - 2*L_env_mollicutes
df_size[mollicutes_set,]$D.mean_corr <- df_size[mollicutes_set,]$D.mean - 2*L_env_mollicutes
df_size[mollicutes_set,]$L.low_corr <- df_size[mollicutes_set,]$L.low - 2*L_env_mollicutes
df_size[mollicutes_set,]$L.high_corr <- df_size[mollicutes_set,]$L.high - 2*L_env_mollicutes
df_size[mollicutes_set,]$L.mean_corr <- df_size[mollicutes_set,]$L.mean - 2*L_env_mollicutes
df_size[mollicutes_set,]$L_uncoiled.low_corr <- df_size[mollicutes_set,]$L_uncoiled.low - 2*L_env_mollicutes
df_size[mollicutes_set,]$L_uncoiled.high_corr <- df_size[mollicutes_set,]$L_uncoiled.high - 2*L_env_mollicutes
df_size[mollicutes_set,]$L_uncoiled.mean_corr <- df_size[mollicutes_set,]$L_uncoiled.mean - 2*L_env_mollicutes

# Subset by cell shape, and calculate volumes and surface areas of capsules and spheres
df_size_rod <- df_size[(df_size$shape=="rod"),]
df_size_sphere <- df_size[(df_size$shape=="sphere"),]
df_size_helix <- df_size[(df_size$shape=="helical"),]

# For each size feature, we calculate lower and upper bound, and mean
lo_vol_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_vol(df_size_rod$D.low_corr[x], df_size_rod$L.low_corr[x]))))
hi_vol_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_vol(df_size_rod$D.high_corr[x], df_size_rod$L.high_corr[x]))))
mean_vol_rod <- t(as.data.frame(lapply(seq(1,nrow(df_size_rod)), function(x) capsule_vol(df_size_rod$D.mean_corr[x], df_size_rod$L.mean_corr[x]))))

# The cytoplasmic volume of a helical cell is calculated as the capsule volume, with inner diameter and uncoiled length as arguments
lo_vol_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) capsule_vol(df_size_helix$D.low_corr[x], df_size_helix$L_uncoiled.low_corr[x]))))
hi_vol_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) capsule_vol(df_size_helix$D.high_corr[x], df_size_helix$L_uncoiled.high_corr[x]))))
mean_vol_hel <- t(as.data.frame(lapply(seq(1,nrow(df_size_helix)), function(x) capsule_vol(df_size_helix$D.mean_corr[x], df_size_helix$L_uncoiled.mean_corr[x]))))

# Finally, calculate sphere volume
lo_vol_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_vol(df_size_sphere$D.low_corr[x]))))
hi_vol_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_vol(df_size_sphere$D.high_corr[x]))))
mean_vol_sphere <- t(as.data.frame(lapply(seq(1,nrow(df_size_sphere)), function(x) sphere_vol(df_size_sphere$D.mean_corr[x]))))

# Calculate total surface area; Note that we are using total linear dimensions, and not inner ones
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


#========== Scaling of surface with volume ==========
scaling_surf_vol_rod <- lm(data=df_size_rod, log10(S.mean)~log10(V.mean))
summary(scaling_surf_vol_rod)

scaling_surf_vol_sphere <- lm(data=df_size_sphere, log10(S.mean)~log10(V.mean))
summary(scaling_surf_vol_sphere)

scaling_surf_vol_helix <- lm(data=df_size_helix, log10(S.mean)~log10(V.mean))
summary(scaling_surf_vol_helix)

scaling_surf_vol_full <- lm(data=df_size, log10(S.mean)~log10(V.mean))
summary(scaling_surf_vol_full)

scaling_surf_vol_no_helix <- lm(data=df_size[df_size$shape!="helical",], log10(S.mean)~log10(V.mean))
summary(scaling_surf_vol_no_helix)


#========== Justification for fixed S/V across the cell cycle ==========
# Set the color-blind scale
colBlindScale <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7", "#DB6D00", "#EF3B2C")
colBlindScaleExtended <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                           "#0072B2", "#D55E00", "#CC79A7", "#DB6D00", "#EF3B2C",
                           "#67001F","#F7FCFD","#CB181D","#78C679","#F46D43","#A6CEE3","#FD8D3C","#A6D854","#D4B9DA","#6A51A3",
                           "#7F0000","#D9D9D9","#FFF7BC","#F0F0F0","#C7EAE5","#003C30","#F16913","#FFF7FB","#8C6BB1","#C7E9B4",
                           "#762A83","#FC9272","#AE017E","#F7F7F7","#DF65B0","#74C476")
scaleFUN <- function(x) sprintf("%.2f", x)
axes_style <- element_text(size=15, face="bold", colour = "black")

tmp_size_sphere <- df_size_sphere
tmp_size_rod <- df_size_rod
tmp_size_helix <- df_size_helix

# Spheres
birth_pi <- with(
                  df_size_sphere,
                  capsule_surf(D.mean,D.mean)/capsule_vol(D.mean_corr,D.mean_corr)
)
div_pi <- with(
              df_size_sphere,
              capsule_surf(D.mean,2*D.mean-D.mean/3)/
              capsule_vol(D.mean_corr,2*D.mean_corr-D.mean_corr/3)
)
tmp_size_sphere$delta_pi <- birth_pi-div_pi

# Rods
birth_pi <- with(
                df_size_rod,
                capsule_surf(D.mean,L.mean)/capsule_vol(D.mean_corr,L.mean_corr)
)
div_pi <- with(
              df_size_rod,
              capsule_surf(D.mean,2*L.mean-D.mean/3)/capsule_vol(D.mean_corr,2*L.mean_corr-D.mean_corr/3)
)  
tmp_size_rod$delta_pi <- birth_pi-div_pi

# Helices
birth_pi <- with(
               df_size_helix,
               capsule_surf(D.mean,L_uncoiled.mean)/capsule_vol(D.mean_corr,L_uncoiled.mean_corr)
                )
div_pi <- with(
               df_size_helix, 
               capsule_surf(D.mean,2*L_uncoiled.mean-D.mean/3)/capsule_vol(D.mean_corr,2*L_uncoiled.mean_corr-D.mean_corr/3)
              )
tmp_size_helix$delta_pi <- birth_pi-div_pi

p_dist_delta_pi <- ggplot(rbind(tmp_size_helix,tmp_size_rod,tmp_size_sphere), aes(x=log10(S.mean/V.mean), y=log10(delta_pi), col=shape)) + 
  geom_point(size=3.5, alpha=0.5) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_y_continuous(labels=scaleFUN) +
  ylab(bquote(bold('S/V difference, Log'[10]*'['*Delta*Pi*' ('*mu*''*m^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]')))
p_dist_delta_pi

p_scaling_delta_pi  <- ggplot(rbind(tmp_size_helix,tmp_size_rod,tmp_size_sphere), aes(x=log10(delta_pi), fill=shape)) + 
  geom_histogram(alpha=0.5, position="identity", color="black") +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_y_continuous(labels=scaleFUN) +
  xlab(bquote(bold('S/V difference, Log'[10]*'['*Delta*Pi*' ('*mu*''*m^-1*')]'))) + ylab(bquote(bold('Count'))) +
  xlim(-3,1)
p_scaling_delta_pi

legend_delta_pi <- ggplot(rbind(tmp_size_helix,tmp_size_rod,tmp_size_sphere), aes(x=log10(delta_pi), fill=shape)) + 
  geom_histogram(alpha=0.5, position="identity", color="black") +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_y_continuous(labels=scaleFUN) +
  xlab(bquote(bold('S/V difference, Log'[10]*'['*Delta*Pi*' ('*mu*''*m^-1*')]'))) + ylab(bquote(bold('Count'))) +
  xlim(-3,1)

ggsave(file="./main/p_dist_delta_pi.pdf", plot=p_dist_delta_pi, width = 5, height = 5)
ggsave(file="./main/p_scaling_delta_pi.pdf", plot=p_scaling_delta_pi, width = 5, height = 5)
ggsave(file="./main/legend_delta_pi.pdf", plot=legend_delta_pi, width = 5, height = 5)


#========== Correlation between volume and cell dry weight ==========
df_dens <- read.csv("./cell_traits_data/bacterial_envelopes_Sep232021 - bacterial_cell_density.csv",
                     header = TRUE, sep = ",", stringsAsFactors = FALSE)
# The slope is not significantly different from unity
model_dens <- lm(log10(dcw_ng)-offset(1*log10(v_um3))~log10(v_um3), data = df_dens)
summary(model_dens)
# Regression coefficient
model_dens <- lm(log10(dcw_ng)~log10(v_um3), data = df_dens)
summary(model_dens)

#========== Merging ribosomal abundances, quantitative proteomics, and cell size data ==========
df_abund <- read.csv("./cell_traits_data/bacterial_envelopes_Sep232021 - ribosome_abundance.csv",
                     header = TRUE, sep = ",", stringsAsFactors = FALSE)
colnames(df_abund)[4] <- "source"
df_abund$ribosome_mass_total <- (df_abund$ribosome_abundance*7336*110)/(6.022*10^23)
df_abund$ribo_frac <- df_abund$ribosome_mass_total/(
  ((10^(-9))*(0.54*10^(-3.30537))*df_abund$cell_volume..um3.^(0.92892))
)

spiroplasma_cell_mass <- 3.69*10^(-14) # grams per cell; Trachtenberg et al. 2014 Plos
df_abund$ribo_frac[grepl("\\bSpiroplasma", df_abund$species)] <- 
  df_abund$ribosome_mass_total[grepl("\\bSpiroplasma", df_abund$species)]/spiroplasma_cell_mass

df_proteomics <- read.csv("./main/quant_proteome_data.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
substr(df_proteomics$species, 1, 1) <- toupper(substr(df_proteomics$species, 1, 1))
df_proteomics <- bind_rows(df_proteomics, df_abund)
df_proteomics <- select(df_proteomics, species, lip_frac, ribo_frac, source)


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
df_proteomics <- df_proteomics[!(df_proteomics$species=="Mycoplasma pneumoniae"&df_proteomics$source=="Wang2012"),]
tmp_source <- aggregate(df_proteomics$source, list(df_proteomics$species), toString, na.rm=TRUE)
colnames(tmp_source) <- c("species","source")

df_proteomics <- aggregate(df_proteomics[,c(2:4)], list(df_proteomics$species), mean, na.rm=TRUE)
df_proteomics <- df_proteomics[,-4]
colnames(df_proteomics)[1] <- "species"
df_proteomics <- merge(df_proteomics,tmp_source,by = "species",all.y = FALSE)

df2 <- merge(df_proteomics, df_size, by="species", all.x = TRUE)
df2 <- select(df2, species, shape, V.mean, S.mean, lip_frac, ribo_frac, Total, source)
colnames(df2) <- c("species", "shape", "V.mean","S.mean","lip_frac","ribo_frac","total_thickness_nm","source")

df2 <- df2[df2$species!="Desulfovibrio vulgaris",]
df1 <- select(df1, species, shape, V.mean, S.mean, growth_rate, Total, rRNA.genes, tRNA.genes)
colnames(df1) <- c("species", "shape", "V.mean","S.mean", "growth_rate", "total_thickness_nm", "rRNA genes", "tRNA genes")

in_file_metabolism <- "./cell_traits_data/bacterial_envelopes_Sep232021 - physiological_properties.csv"
df_metabolism <- read.csv(in_file_metabolism, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_metabolism <- select(df_metabolism, species, tca_steps)
df1 <- merge(df1, df_metabolism, all.x = TRUE)

write.csv(df1, "./main/growth_scaling.shape.genome.csv", row.names = FALSE)
write.csv(df2, "./main/proteome_scaling.csv", row.names = FALSE)

nrow(df2[!is.na(df2$V.mean)|!is.na(df2$S.mean),])

sum(duplicated(df1$species))
sum(duplicated(df2$species))