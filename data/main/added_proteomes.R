#=============== USAGE =====================
# Purpose: process additional proteomes, individually collected aside from Wang2012 and Muller2020
# Input: folder containing subfolders from individual proteomic studies
# Output: processed data in the same format as given by paxDB_proteomes.R
# "species"	"lip_frac" "ribo_frac"

remove(list=ls())

library(KEGGREST)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(Peptides)

word_type <- 2
sector_classification <- function(kegg_entry) {
  # If non-existing entry is passed, set return to zero
  if(is.null(kegg_entry)) {return(rep(0,word_type+1))}
  # Set the keywords used to match BRITE hierarchy
  lipid_words <- c("\\bFatty acid biosynthesis\\b", "\\bPeptidoglycan biosynthesis\\b", "\\bLipopolysaccharide biosynthesis\\b")
  ribo_words <- c("\\bRibosomal protein")
  
  lip_flag <- sum(unlist(lapply(lipid_words, function(x) sum(grepl(x, kegg_entry[[1]]$BRITE)))))>0
  ribo_flag <- sum(unlist(lapply(ribo_words, function(x) sum(grepl(x, kegg_entry[[1]]$BRITE)))))>0
  prot_mw <- ifelse(is.null(kegg_entry[[1]]$AASEQ[[1]]),0,mw(kegg_entry[[1]]$AASEQ[[1]]))
  
  c(lip_flag, ribo_flag, prot_mw)
}


# =================== Processing Osbak2016 -- Treponema pallidum ===================
setwd("/home/bogi/Desktop/growth_rate/data/bacterial_proteomes/Osbak2016/")
df_abund <- read.csv("./treponema_pallidum_abundances.csv",
                   skip = 2, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_abund$kegg_id <- paste0("tpa:",df_abund$TP.number,sep="")

class_list <- lapply(df_abund$kegg_id, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_class <- as.data.frame(t(as.data.frame(class_list)))
# Subsetting NSAF values, given that they correspond to relative abundances
# We also use MW that authors reported as opposed to mw from sector_classigication() function
df_abund <- select(df_abund, kegg_id, Molecular.Weight..kDa., colnames(df_abund)[c(13:18)])
df_abund <- cbind(df_abund,df_class)
rownames(df_abund) <- c()
colnames(df_abund)[c(9:11)] <- c("lip_frac","ribo_frac","prot_mw")

df_mass <- data.frame(lapply(df_abund[,c(3:8)], function(x) x*df_abund$Molecular.Weight..kDa.))
df_fractions <- t(t(df_mass)/colSums(df_mass, na.rm = TRUE))

df <- cbind(df_abund$kegg_id, df_abund$Molecular.Weight..kDa., df_fractions, df_abund[c(9:11)])
write.csv(df, file="osbak2016.kegged.csv", row.names = FALSE)

# Mass fractions averaged across replicates
df_abund <- read.csv("osbak2016.kegged.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_l <- df_abund[df_abund$lip_frac==TRUE,]
df_r <- df_abund[df_abund$ribo_frac==TRUE,]

phiR <- mean(colSums(df_r[c(3:8)], na.rm = TRUE))
phiL <- mean(colSums(df_l[c(3:8)], na.rm = TRUE))
df_out <- data.frame("Treponema pallidum", phiL, phiR)
colnames(df_out) <- c("species","lip_frac", "ribo_frac")
write.csv(df_out, "../../main/osbak2016.processed.csv", row.names = FALSE)


# =================== Processing Angel2010 -- Borrelia burgdorferi ===================
setwd("/home/bogi/Desktop/growth_rate/data/bacterial_proteomes/Angel2010/")

df_abund <- read.csv("angel2010.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_abund <- df_abund[df_abund$Log=="x",-c(1:3)]
df_abund$kegg_id <- paste0("bbu:BB_",unlist(lapply(df_abund$characterized_orf, function(x) unlist(strsplit(x,"BB"))[2])),sep="")
class_list <- lapply(df_abund$kegg_id, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_class <- as.data.frame(t(as.data.frame(class_list)))
df_abund <- select(df_abund, kegg_id, total_spectral_count)
df_abund <- cbind(df_abund, df_class)
rownames(df_abund) <- c()
colnames(df_abund)[c(3:5)] <- c("lip_frac","ribo_frac","prot_mw")

# Mass calculated by normalizing spectral count with protein length
df_mass <- as.data.frame(with(df_abund, total_spectral_count*prot_mw))
# Calculate mass fractions across all growth conditions
df_mass <- t(t(df_mass)/colSums(df_mass, na.rm = TRUE))
df_out <- cbind(df_mass, df_abund[,c(3:5)])
write.csv(df_out, file="angel2010.kegged.csv", row.names = FALSE)

# Mass fractions averaged across replicates
df_abund <- read.csv("angel2010.kegged.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_l <- df_abund[df_abund$lip_frac==TRUE,]
df_r <- df_abund[df_abund$ribo_frac==TRUE,]

# Reported mass fractions are averages across all examined growth conditions
phiR <- mean(colSums(df_r[1], na.rm = TRUE))
phiL <- mean(colSums(df_l[1], na.rm = TRUE))
df_out <- data.frame("Borrelia burgdorferi", phiL, phiR)
colnames(df_out) <- c("species","lip_frac","ribo_frac")
write.csv(df_out, "../../main/angel2010.processed.csv", row.names = FALSE)


# =================== Processing Srivastava2016 -- Polynucleobacter symbioticus ===================
setwd("/home/bogi/Desktop/growth_rate/data/bacterial_proteomes/Srivastava2020/")
df_abund <- read.csv("srivastava2020.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
name_list <- lapply(paste0("up:", df_abund$X.1, sep=""), function(x) tryCatch(keggConv("genes", x)))
df_abund$kegg_id <- name_list
class_list <- lapply(df_abund$kegg_id, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_class <- as.data.frame(t(as.data.frame(class_list)))

# Convert from reported Log_2[LFQ] to LFQ, as with the data for other species
df_abund <- select(df_abund, kegg_id, X26...C_01, X26...C_02, X26...C_03)
df_abund <- cbind(df_abund, df_class)
rownames(df_abund) <- c()
colnames(df_abund)[c(5:7)] <- c("lip_frac","ribo_frac","prot_mw")
df_abund$X26...C_01 <- 2^df_abund$X26...C_01
df_abund$X26...C_02 <- 2^df_abund$X26...C_02
df_abund$X26...C_03 <- 2^df_abund$X26...C_03

df_mass <- df_abund[,c(2:4)]*df_abund$prot_mw
df_mass <- t(t(df_mass)/colSums(df_mass, na.rm = TRUE))

df_out <- cbind(df_mass, df_abund[,c(5:7)])
write.csv(df_out, file="srivastava2020.kegged.csv", row.names = FALSE)

# Mass fractions averaged across replicates
df_abund <- read.csv("srivastava2020.kegged.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_l <- df_abund[df_abund$lip_frac==TRUE,]
df_r <- df_abund[df_abund$ribo_frac==TRUE,]

phiR <- mean(colSums(df_r[1], na.rm = TRUE))
phiL <- mean(colSums(df_l[1], na.rm = TRUE))
df_out <- data.frame("Polynucleobacter asymbioticus", phiL, phiR)
colnames(df_out) <- c("species","lip_frac","ribo_frac")
write.csv(df_out, "../../main/srivastava2020.processed.csv", row.names = FALSE)


# =================== Processing Matteau2020 -- Mesoplasma florum ===================
setwd("/home/bogi/Desktop/growth_rate/data/bacterial_proteomes/Matteau2020/")
df_abund <- read.csv("matteau2020_proteomics.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[-1,]
name_frame <- data.frame(paste0("mfl:", df_abund$gene.name..RefSeq., sep = ""))
df_abund <- select(df_abund, fraction.of.total.protein.mass....)
df_abund <- cbind(name_frame, df_abund)
colnames(df_abund) <- c("kegg_id", "mass_fractions")

class_list <- lapply(df_abund$kegg_id, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_class <- as.data.frame(t(as.data.frame(class_list)))
df_abund <- cbind(df_abund, df_class)
rownames(df_abund) <- c()
colnames(df_abund)[c(3:5)] <- c("lip_frac","ribo_frac","prot_mw")
write.csv(df_abund, file="matteau2020.kegged.csv", row.names = FALSE)

# Mass fractions averaged across replicates
df_abund <- read.csv("matteau2020.kegged.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[-1]
df_abund <- df_abund[-nrow(df_abund),]

df_l <- df_abund[df_abund$lip_frac==TRUE,]
df_r <- df_abund[df_abund$ribo_frac==TRUE,]

phiR <- sum(as.numeric(df_r[[1]]), na.rm = TRUE)/100
phiL <- sum(as.numeric(df_l[[1]]), na.rm = TRUE)/100
df_out <- data.frame("Mesoplasma florum", phiL, phiR)
colnames(df_out) <- c("species","lip_frac","ribo_frac")
write.csv(df_out, "../../main/matteau2020.processed.csv", row.names = FALSE)


# =================== Processing Masson2021 -- Spiroplasma poulsonii ===================
setwd("/home/bogi/Desktop/growth_rate/data/bacterial_proteomes/Masson2021/")
df_abund <- read.csv("masson2021_proteomics.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_abund$LFQ_intensity_SPOR_01 <- df_abund$LFQ_intensity_SPOR_01 * df_abund$Mol._weight_.kDa.
df_abund$LFQ_intensity_SPOR_02 <- df_abund$LFQ_intensity_SPOR_02 * df_abund$Mol._weight_.kDa.
df_abund$LFQ_intensity_SPOR_03 <- df_abund$LFQ_intensity_SPOR_03 * df_abund$Mol._weight_.kDa.
df_abund$LFQ_intensity_SPOR_04 <- df_abund$LFQ_intensity_SPOR_04 * df_abund$Mol._weight_.kDa.
spor1_sum <- sum(df_abund$LFQ_intensity_SPOR_01)
spor2_sum <- sum(df_abund$LFQ_intensity_SPOR_02)
spor3_sum <- sum(df_abund$LFQ_intensity_SPOR_03)
spor4_sum <- sum(df_abund$LFQ_intensity_SPOR_04)

df_abund <- df_abund[grepl("ribosomal protein",df_abund$Description),c(5,6,31:34)]
ribo_frac1 <- sum(df_abund$LFQ_intensity_SPOR_01)/spor1_sum
ribo_frac2 <- sum(df_abund$LFQ_intensity_SPOR_02)/spor2_sum
ribo_frac3 <- sum(df_abund$LFQ_intensity_SPOR_03)/spor3_sum
ribo_frac4 <- sum(df_abund$LFQ_intensity_SPOR_04)/spor4_sum

df_out <- data.frame("Spiroplasma poulsonii", NA, mean(ribo_frac1, ribo_frac2, ribo_frac3, ribo_frac4))
colnames(df_out) <- c("species","lip_frac","ribo_frac")
write.csv(df_out, "../../main/masson2021.processed.csv", row.names = FALSE)

