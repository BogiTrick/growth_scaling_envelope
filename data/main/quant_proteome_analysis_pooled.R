#=============== USAGE =====================
# Purpose: annotated mass fractions --> total mass fractions for each sector
#    Takes data with assigned functions to each protein and associated mass fractions, and calculates total
#    mass fractions for each proteomic sector.
# Input: functionally annotated proteomes with proteomic mass fractions for each protein entry
# Output 1: proteomic_fractions_full.csv, with the following columnar structure
#    ""	"perturbation"	"strain"	"growth_medium"	"perturbant_concentration"	"growth_rates"	"phi_R"	"phi_L"	"source"
#    The data is further used to infer rate constants in the model.
# Output 2: proteomic_count.csv, contains number of proteins classified in each sector per quantitative proteomic study

remove(list=ls())
library(dplyr)

word_type <- 2
setwd("/home/bogi/Desktop/growth_rate/data/eco_proteomes/assigned_function")
colLabels <- c("perturbation", "strain", "growth_medium", "perturbant_concentration","growth_rates","phi_R","phi_L","source")


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

helix_volume <- function(d_cell, d_helix, pitch, len_tot) {
  actual_len <- sqrt(pitch^2 + (pi*(d_helix-d_cell))^2)*(len_tot/pitch)
  return(capsule_vol(d_cell, actual_len))
}

helix_surface <- function(d_cell, d_helix, pitch, len_tot) {
  actual_len <- sqrt(pitch^2 + (pi*(d_helix-d_cell))^2)*(len_tot/pitch)
  return(capsule_surf(d_cell, actual_len))
}


# Direct measurements of proteome mass fractions
#============= Schmidt2016 ======================
df_Schmidt <- read.csv("schmidt2016_data_assigned.csv",
                       skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[,-c(1:3)]

# Retain only data for proteins which are associated with target KEGG pathways
lipo_genes_Schmidt <- df_Schmidt[df_Schmidt$lipo_genes,]
fractions_lipid_Schmidt <- colSums(lipo_genes_Schmidt[,-c((ncol(df_Schmidt)-word_type+1):ncol(df_Schmidt))], na.rm = TRUE)

ribo_genes_Schmidt <- df_Schmidt[df_Schmidt$ribo_genes,]
fractions_ribo_Schmidt <- colSums(ribo_genes_Schmidt[,-c((ncol(df_Schmidt)-word_type+1):ncol(df_Schmidt))], na.rm = TRUE)

growth_rates_Schmidt <- c(
  0.58, 1.90, 1.27, 0.30, 0.42, 0.46, 0.47, 0.40, 0.50, 0.35,
  0.20, 0.12, 0.00, 0.00, 0.55, 0.66, 0.63, 0.55, 0.47, 0.26, 0.44, 0.65
)

df_fractions_Schmidt <- data.frame(matrix(nrow = length(growth_rates_Schmidt), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Schmidt) <- colLabels
df_fractions_Schmidt$perturbation <- "nutrient conditions"
df_fractions_Schmidt$strain <- "BW25113"
df_fractions_Schmidt$growth_medium <- NA
df_fractions_Schmidt$perturbant_concentration <- NA
df_fractions_Schmidt$growth_rates <- growth_rates_Schmidt
df_fractions_Schmidt$phi_R <- fractions_ribo_Schmidt
df_fractions_Schmidt$phi_L <- fractions_lipid_Schmidt
df_fractions_Schmidt$source <- "Schmidt2016"


#============= Peebo2015 ======================
df_Peebo <- read.csv("peebo2015_data_assigned.csv",
                       skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[,-c(1:3)]

# Retain only data for proteins which are associated with target KEGG pathways
lipo_genes_Peebo <- df_Peebo[df_Peebo$lipo_genes,]
fractions_lipid_Peebo <- colSums(lipo_genes_Peebo[,-c((ncol(df_Peebo)-word_type+1):ncol(df_Peebo))], na.rm = TRUE)

ribo_genes_Peebo <- df_Peebo[df_Peebo$ribo_genes,]
fractions_ribo_Peebo <- colSums(ribo_genes_Peebo[,-c((ncol(df_Peebo)-word_type+1):ncol(df_Peebo))], na.rm = TRUE)

growth_rates_Peebo <- c(
  0.21, 0.31, 0.41, 0.51, 0.22, 0.26, 0.36, 0.46, 0.51, 0.22, 0.25, 0.35,
  0.45, 0.55, 0.65, 0.74, 0.82, 0.22, 0.42, 0.53, 0.63, 0.73, 0.78
)

df_fractions_Peebo <- data.frame(matrix(nrow = length(growth_rates_Peebo), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Peebo) <- colLabels
df_fractions_Peebo$perturbation <- "nutrient conditions"
df_fractions_Peebo$strain <- "BW25113"
df_fractions_Peebo$growth_medium <- NA
df_fractions_Peebo$perturbant_concentration <- NA
df_fractions_Peebo$growth_rates <- growth_rates_Peebo
df_fractions_Peebo$phi_R <- fractions_ribo_Peebo
df_fractions_Peebo$phi_L <- fractions_lipid_Peebo
df_fractions_Peebo$source <- "Peebo2015"


#============= Erickson2017 ======================
df_Erickson <- read.csv("erickson2017_data_assigned.csv",
                     skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[,-c(1:3)]

# Retain only data for proteins which are associated with target KEGG pathways
lipo_genes_Erickson <- df_Erickson[df_Erickson$lipo_genes,]
fractions_lipid_Erickson <- colSums(lipo_genes_Erickson[,-c((ncol(df_Erickson)-word_type+1):ncol(df_Erickson))], na.rm = TRUE)

ribo_genes_Erickson <- df_Erickson[df_Erickson$ribo_genes,]
fractions_ribo_Erickson <- colSums(ribo_genes_Erickson[,-c((ncol(df_Erickson)-word_type+1):ncol(df_Erickson))], na.rm = TRUE)

growth_rates_Erickson <- c(0.98, 0.45, 0.45, 0.98)

df_fractions_Erickson <- data.frame(matrix(nrow = length(growth_rates_Erickson), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Erickson) <- colLabels
df_fractions_Erickson$perturbation <- "nutrient conditions"
df_fractions_Erickson$strain <- "NCM3722"
df_fractions_Erickson$growth_medium <- NA
df_fractions_Erickson$perturbant_concentration <- NA
df_fractions_Erickson$growth_rates <- growth_rates_Erickson
df_fractions_Erickson$phi_R <- fractions_ribo_Erickson
df_fractions_Erickson$phi_L <- fractions_lipid_Erickson
df_fractions_Erickson$source <- "Erickson2017"


#============= Valgepea2013 ======================
df_Valgepea <- read.csv("valgepea2013_data_assigned.csv",
                        skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[,-c(1:3)]

# Retain only data for proteins which are associated with target KEGG pathways
lipo_genes_Valgepea <- df_Valgepea[df_Valgepea$lipo_genes,]
fractions_lipid_Valgepea <- colSums(lipo_genes_Valgepea[,-c((ncol(df_Valgepea)-word_type+1):ncol(df_Valgepea))], na.rm = TRUE)

ribo_genes_Valgepea <- df_Valgepea[df_Valgepea$ribo_genes,]
fractions_ribo_Valgepea <- colSums(ribo_genes_Valgepea[,-c((ncol(df_Valgepea)-word_type+1):ncol(df_Valgepea))], na.rm = TRUE)

growth_rates_Valgepea <- c(0.11, 0.21, 0.31, 0.40, 0.49)

df_fractions_Valgepea <- data.frame(matrix(nrow = length(growth_rates_Valgepea), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Valgepea) <- colLabels
df_fractions_Valgepea$perturbation <- "nutrient conditions"
df_fractions_Valgepea$strain <- "MG1655"
df_fractions_Valgepea$growth_medium <- NA
df_fractions_Valgepea$perturbant_concentration <- NA
df_fractions_Valgepea$growth_rates <- growth_rates_Valgepea
df_fractions_Valgepea$phi_R <- fractions_ribo_Valgepea
df_fractions_Valgepea$phi_L <- fractions_lipid_Valgepea
df_fractions_Valgepea$source <- "Valgepea2013"


#============= Li2014 ======================
df_Li <- read.csv("li2014_data_assigned.csv",
                        skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[,-c(1:3)]

# Retain only data for proteins which are associated with target KEGG pathways
lipo_genes_Li <- df_Li[df_Li$lipo_genes,]
fractions_lipid_Li <- colSums(lipo_genes_Li[,-c((ncol(df_Li)-word_type+1):ncol(df_Li))], na.rm = TRUE)

ribo_genes_Li <- df_Li[df_Li$ribo_genes,]
fractions_ribo_Li <- colSums(ribo_genes_Li[,-c((ncol(df_Li)-word_type+1):ncol(df_Li))], na.rm = TRUE)

growth_rates_Li <- c(60*log(2)/21.5, 60*log(2)/56.3, 60*log(2)/26.5)

df_fractions_Li <- data.frame(matrix(nrow = length(growth_rates_Li), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Li) <- colLabels
df_fractions_Li$perturbation <- "nutrient conditions"
df_fractions_Li$strain <- "MG1655"
df_fractions_Li$growth_medium <- NA
df_fractions_Li$perturbant_concentration <- NA
df_fractions_Li$growth_rates <- growth_rates_Li
df_fractions_Li$phi_R <- fractions_ribo_Li
df_fractions_Li$phi_L <- fractions_lipid_Li
df_fractions_Li$source <- "Li2014"


#============= Mori2021 ======================
df_Mori <- read.csv("mori2021_data_assigned.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[-c(1,38)] # column 38 is bNum
df_growths <- read.csv("mori2021_growth.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Retain only data for proteins which are associated with target KEGG pathways
lipo_genes_Mori <- df_Mori[df_Mori$lipo_genes,]
fractions_lipid_Mori <- colSums(lipo_genes_Mori[,-c((ncol(df_Mori)-word_type+1):ncol(df_Mori))], na.rm = TRUE)

ribo_genes_Mori <- df_Mori[df_Mori$ribo_genes,]
fractions_ribo_Mori <- colSums(ribo_genes_Mori[,-c((ncol(df_Mori)-word_type+1):ncol(df_Mori))], na.rm = TRUE)

df_fractions_Mori <- data.frame(matrix(nrow = length(fractions_lipid_Mori), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Mori) <- colLabels

df_fractions_Mori$phi_R <- fractions_ribo_Mori
df_fractions_Mori$phi_L <- fractions_lipid_Mori
df_fractions_Mori$source <- "Mori2021"
df_fractions_Mori$growth_rates <- df_growths$Growth.rate..1.h.
df_fractions_Mori$strain <- df_growths$Strain
df_fractions_Mori$perturbation <- df_growths$Group
df_fractions_Mori$growth_medium <- df_growths$Growth.medium
df_fractions_Mori$perturbant_concentration <- NA

df_fractions_Mori[df_fractions_Mori$perturbation=="C-limitation",]$perturbation <- "nutrient conditions"
df_fractions_Mori[df_fractions_Mori$perturbation=="R-limitation",]$perturbation <- "chloramphenicol"
  

# Indirect measurements of proteome mass fractions
#============= Dai2018 and Si2017 ======================
rho_Dai <- 0.48 # Page 3 of the supplement in Dai2018. This is obtained from 1/sigma, where sigma=2.1 
growth_rates_Dai <- c(
  1.8,1.57,1.28,1.17,1.12,0.98,0.99,0.95,0.92,0.75,0.72,0.70,0.69,0.69,0.55,0.50,0.46,
  0.41,0.34,0.38,0.33,0.29,0.20,0.23,0.13,0.035
)
phiR <- rho_Dai*c(
  0.476,0.431,0.364,0.330,0.306,0.294,0.268,0.272,0.253,0.233,0.216,
  0.225,0.227,0.216,0.193,0.184,0.172,0.172,0.152,0.160,0.152,0.147,0.130,0.134,0.118,0.097
)

df_fractions_Dai <- data.frame(matrix(nrow = length(growth_rates_Dai), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Dai) <- colLabels
df_fractions_Dai$perturbation <- "nutrient conditions"
df_fractions_Dai$strain <- "NCM3722"
df_fractions_Dai$growth_medium <- NA
df_fractions_Dai$perturbant_concentration <- NA
df_fractions_Dai$growth_rates <- growth_rates_Dai
df_fractions_Dai$phi_R <- phiR
df_fractions_Dai$phi_L <- NA
df_fractions_Dai$source <- "Dai2018"


df_Si <- read.csv("../si2017_data.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_fractions_Si <- select(df_Si, experiment.name, type.of.perturbation, strain.type..background., growth.media,
                          inducer.or.drug.added, concentration, growth.rate..1.hours., RNA.protein)
df_fractions_Si_aggregate <- df_fractions_Si %>% group_by(type.of.perturbation, strain.type..background., 
                                                          growth.media, concentration) %>%  
  summarise(mean_growth.rate=mean(growth.rate..1.hours.), mean_RNA.protein=mean(RNA.protein), .groups='drop')
df_fractions_Si <- df_fractions_Si_aggregate[which(((df_fractions_Si_aggregate$strain.type..background.=="NCM3722")|
                                                                 (df_fractions_Si_aggregate$strain.type..background.=="MG1655"))),]

df_fractions_Si$phiR <- rho_Dai*df_fractions_Si$mean_RNA.protein
df_fractions_Si <- df_fractions_Si[,-6]
df_fractions_Si$phi_L <- NA
df_fractions_Si$source <- "Si2017"
colnames(df_fractions_Si) <- colLabels


#============= Forchhammer1971, Dennis1996, and Scott2010 ======================
growth_rates_Scott <- c(0.4,0.57,0.71,1.00,1.31,1.58)
phiR <- rho_Dai*c(0.177,0.230,0.224,0.287,0.414,0.466)
df_fractions_Scott <- data.frame(matrix(nrow = length(growth_rates_Scott), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Scott) <- colLabels
df_fractions_Scott$perturbation <- "nutrient conditions"
df_fractions_Scott$strain <- "EQ2"
df_fractions_Scott$growth_medium <- NA
df_fractions_Scott$perturbant_concentration <- NA
df_fractions_Scott$growth_rates <- growth_rates_Scott
df_fractions_Scott$phi_R <- phiR
df_fractions_Scott$phi_L <- NA
df_fractions_Scott$source <- "Scott2010"

df_Scott_Ab <- read.csv("../scott2010_data_Ab.csv",
                        skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_fractions_Scott_Ab <- data.frame(matrix(nrow = nrow(df_Scott_Ab), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Scott_Ab) <- colLabels
df_fractions_Scott_Ab$perturbation <- "chloramphenicol"
df_fractions_Scott_Ab$strain <- "EQ2"
df_fractions_Scott_Ab$growth_medium <- df_Scott_Ab$medium
df_fractions_Scott_Ab$perturbant_concentration <- NA
df_fractions_Scott_Ab$growth_rates <- df_Scott_Ab$growth_rate_mean
df_fractions_Scott_Ab$phi_R <- df_Scott_Ab$rna.protein_mean*rho_Dai
df_fractions_Scott_Ab$phi_L <- NA
df_fractions_Scott_Ab$source <- "Scott2010"

growth_rates_Forchhammer <- c(0.38,0.60,1.04,1.46,1.73)
phiR <- rho_Dai*c(0.189,0.224,0.295,0.421,0.469)
df_fractions_Forchhammer <- data.frame(matrix(nrow = length(growth_rates_Forchhammer), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Forchhammer) <- colLabels
df_fractions_Forchhammer$perturbation <- "nutrient conditions"
df_fractions_Forchhammer$strain <- "EQ2"
df_fractions_Forchhammer$growth_medium <- NA
df_fractions_Forchhammer$perturbant_concentration <- NA
df_fractions_Forchhammer$growth_rates <- growth_rates_Forchhammer
df_fractions_Forchhammer$phi_R <- phiR
df_fractions_Forchhammer$phi_L <- NA
df_fractions_Forchhammer$source <- "Forchhammer1971"

growth_rates_Bremer <- c(0.42,0.69,1.04,1.39,1.73)
phiR <- rho_Dai*c(0.20,0.255,0.331,0.391,0.471)
df_fractions_Bremer <- data.frame(matrix(nrow = length(growth_rates_Bremer), ncol = length(colLabels)), stringsAsFactors = FALSE)
colnames(df_fractions_Bremer) <- colLabels
df_fractions_Bremer$perturbation <- "nutrient conditions"
df_fractions_Bremer$strain <- "EQ2"
df_fractions_Bremer$growth_medium <- NA
df_fractions_Bremer$perturbant_concentration <- NA
df_fractions_Bremer$growth_rates <- growth_rates_Bremer
df_fractions_Bremer$phi_R <- phiR
df_fractions_Bremer$phi_L <- NA
df_fractions_Bremer$source <- "Bremer1996"


#============= Write pooled data =============
df_full_out <- rbind(df_fractions_Schmidt,df_fractions_Peebo,df_fractions_Valgepea,df_fractions_Erickson,df_fractions_Li,
      df_fractions_Si,df_fractions_Dai,df_fractions_Bremer,df_fractions_Scott,df_fractions_Forchhammer,df_fractions_Scott_Ab,df_fractions_Mori)

write.csv(df_full_out, "../../main/proteomic_fractions_full.csv")


#============= Protein number statistics ======================
study_nb <- 6
class_nb <- 2

get_protein_nb <- function(df) {
  return(nrow(df[rowSums(is.na(df[,-c((length(df)-class_nb+1):length(df))])) != ncol(df[,-c((length(df)-class_nb+1):length(df))]), ]))
}

prot_statistics <- data.frame(matrix(nrow = class_nb, ncol = study_nb), stringsAsFactors = FALSE)
colnames(prot_statistics) <- c("Schmidt2016","Peebo2015","Valgepea2013","Erickson2017","Li2014")
rownames(prot_statistics) <- c("ribo_nb","lipo_nb")
prot_statistics$Schmidt2016 <- c(get_protein_nb(ribo_genes_Schmidt), get_protein_nb(lipo_genes_Schmidt))
prot_statistics$Peebo2015 <- c(get_protein_nb(ribo_genes_Peebo), get_protein_nb(lipo_genes_Peebo))
prot_statistics$Valgepea2013 <- c(get_protein_nb(ribo_genes_Valgepea), get_protein_nb(lipo_genes_Valgepea))
prot_statistics$Erickson2017 <- c(get_protein_nb(ribo_genes_Erickson), get_protein_nb(lipo_genes_Erickson))
prot_statistics$Li2014 <- c(get_protein_nb(ribo_genes_Li), get_protein_nb(lipo_genes_Li))
prot_statistics$Mori2021 <- c(get_protein_nb(ribo_genes_Mori), get_protein_nb(lipo_genes_Mori))
write.csv(prot_statistics, "../../main/proteome_count.csv")

