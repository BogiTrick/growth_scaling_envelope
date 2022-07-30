#=============== USAGE =====================
# Purpose: quant. proteomics --> mass fractions
#  Takes quantitative proteomic data across different growth conditions in E. coli
#  and calculates proteomic mass fractions allocated to three sectors: 
#  envelope-producing enzymes, genome-producing enzymes, and ribosomes
#
# Input: full proteome csv files
# Output: _assigned.csv contains four additional columns, with either TRUE or FALSE for each
#         entry, depending on whether the protein has appropriate functional designation in KEGG BRITE database.
#         Protein can either be "lipo_genes", "ribo_genes".
#         All input files are transformed such that the output reports mass fractions.


remove(list=ls())
setwd("/home/bogi/Desktop/growth_rate/data/eco_proteomes/")
dir.create("assigned_function", showWarnings = FALSE)

library(KEGGREST)
library(dplyr)

# Function sector_classification takes KEGG entry and parses BRITE hierarchy for specific keywords.
# Then, it sets the flags for each function if the appropriate designation is present
word_type <- 2
sector_classification <- function(kegg_entry) {
  # Set the keywords used to match BRITE hierarchy
  # Special character '\b' is used to search for strings that begin and/or end with the exact keyword
  lipid_words <- c("\\bFatty acid biosynthesis\\b", "\\bPeptidoglycan biosynthesis\\b", "\\bLipopolysaccharide biosynthesis\\b")
  ribo_words <- c("\\bRibosomal protein")
  
  # Given multiple possible targets, grep each target in BRITE;
  # Grepl returns TRUE (or 1) for each entry in BRITE hierarchy; 
  # The first sum allows one to know whether a given keyword is present in the BRITE;
  # The second sum (after 'unlist') takes into account that there are multiple words for each function
  #     and having any one of them suffices for belonging in the sector.
  lip_flag <- sum(unlist(lapply(lipid_words, function(x) sum(grepl(x, kegg_entry[[1]]$BRITE)))))>0
  ribo_flag <- sum(unlist(lapply(ribo_words, function(x) sum(grepl(x, kegg_entry[[1]]$BRITE)))))>0
  c(lip_flag, ribo_flag)
}

colLabels <- c("gene","bNum","lipo_genes","ribo_genes")

# Direct measurements of proteome mass fractions
#============= Schmidt2016 ======================
df_Schmidt <- read.csv("schmidt2016_data.csv",
                       skip = 2, header = TRUE, sep = ",", stringsAsFactors = FALSE)
# Start and end columns are set such that data on protein mass per cell is taken
# Dataframe df_cut_ always contains all the data used for function assignment
# It contains KEGG id, gene name, and masses across conditions
start_col <- which(colnames(df_Schmidt)=="Glucose.1")
end_col <- which(colnames(df_Schmidt)=="Glucose.2")-1
df_cut_Schmidt <- cbind(df_Schmidt$Bnumber, df_Schmidt$Gene, df_Schmidt[,start_col:end_col])
colnames(df_cut_Schmidt)[1:2] <- c("bNum","gene")
df_cut_Schmidt$bNum <- paste0("eco:", df_cut_Schmidt$bNum, sep="")      # Add organism code for KEGG querying; 'eco' in this case

# class_flags_ is a list containing BRITE flags for each protein species;
class_flags_Schmidt <- lapply(df_cut_Schmidt$bNum, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_flags_Schmidt <- as.data.frame(t(as.data.frame(class_flags_Schmidt)))
rownames(df_flags_Schmidt) <- c()

# Calculate mass fractions of each reported protein; We are dividing each row with the column sum to obtain mass fractions
# Note that the first two columns with id and name are excluded. Also note the division form (i.e., df/colsum(df)[col(df)]) which
#     is equivalent to transpose(transpose(df)/colsum(df)) and returns row-wise division with denominator
# Finally, flags are appended to the dataframe
mass_fractions_Schmidt <- df_cut_Schmidt[,-c(1,2)]/colSums(df_cut_Schmidt[,-c(1,2)], na.rm = TRUE)[col(df_cut_Schmidt[,-c(1,2)])]
df_Schmidt <- cbind(df_cut_Schmidt$gene, df_cut_Schmidt$bNum, mass_fractions_Schmidt, df_flags_Schmidt)
colnames(df_Schmidt)[c(1:2,(ncol(df_Schmidt)-word_type+1):ncol(df_Schmidt))] <- colLabels

# All data is saved with _assigned extension
write.csv(df_Schmidt, "./assigned_function/schmidt2016_data_assigned.csv")


#============= Peebo2015 ======================
kDa_convert <- 1.66*10^(-21) # Used for conversion from molecs/fL to g/fL

df_Peebo <- read.csv("peebo2015_data.csv", skip = 4, header = TRUE, sep = ",", stringsAsFactors = FALSE)
# Start and end column subset molecules per fL of cell volume
start_col <- which(colnames(df_Peebo)=="X0.21")
end_col <- which(colnames(df_Peebo)=="X0.21.1")-1
# Note that we need molecular weight of each protein to convert absolute abundances to mass fractions
df_cut_Peebo <- cbind(df_Peebo$B.number.identifier, df_Peebo$Gene.name, df_Peebo$Molecular.weight..kDa., df_Peebo[,start_col:end_col])
colnames(df_cut_Peebo)[1:3] <- c("bNum","gene","mass")
df_cut_Peebo$bNum <- paste0("eco:", df_cut_Peebo$bNum, sep="")
# mass_protein = MW * abundance * constant (i.e., unit conversion to g/fL -- in retrospect not necessary given that we divide everything with the same constant to obtain mass fractions)
df_cut_Peebo <- cbind(df_cut_Peebo$bNum, df_cut_Peebo$gene, data.frame(lapply(df_cut_Peebo[,c(4:ncol(df_cut_Peebo))], function(x) x*kDa_convert*df_cut_Peebo$mass)))
colnames(df_cut_Peebo)[1:2] <- c("bNum","gene")

class_flags_Peebo <- lapply(df_cut_Peebo$bNum, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_flags_Peebo <- as.data.frame(t(as.data.frame(class_flags_Peebo)))
rownames(df_flags_Peebo) <- c()

# Calculate mass fractions of each reported protein
mass_fractions_Peebo <- df_cut_Peebo[,-c(1,2)]/colSums(df_cut_Peebo[,-c(1,2)], na.rm = TRUE)[col(df_cut_Peebo[,-c(1,2)])]
df_Peebo <- cbind(df_cut_Peebo$gene, df_cut_Peebo$bNum, mass_fractions_Peebo, df_flags_Peebo)
colnames(df_Peebo)[c(1:2,(ncol(df_Peebo)-word_type+1):ncol(df_Peebo))] <- colLabels
write.csv(df_Peebo, "./assigned_function/peebo2015_data_assigned.csv")


#============= Erickson2017 ======================
# Only gene names reported, so we need to retrieve ids from KEGG
all_genes <- keggList("eco")
all_genes <- data.frame(unname(all_genes), names(all_genes))
all_genes <- cbind(unlist(lapply(all_genes$unname.all_genes., function(x) sub(';.*$','', x))), data.frame(all_genes$names.all_genes.))
colnames(all_genes) <- c("gene","bNum")

df_Erickson <- read.csv("erickson2017_data.csv", skip = 1, header = TRUE, sep = ";", dec=",", stringsAsFactors = FALSE)
# Date for protein masses is taken for the first and the last points in both upshift and downshift experiments,
#  when protein allocation has settled into a steady-state
df_cut_Erickson <- merge(all_genes, df_Erickson[,c(1,4,11,12,19)], by = "gene")             

class_flags_Erickson <- lapply(df_cut_Erickson$bNum, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_flags_Erickson <- as.data.frame(t(as.data.frame(class_flags_Erickson)))
rownames(df_flags_Erickson) <- c()

# Calculate mass fractions of each reported protein
mass_fractions_Erickson <- df_cut_Erickson[,-c(1,2)]/colSums(df_cut_Erickson[,-c(1,2)], na.rm = TRUE)[col(df_cut_Erickson[,-c(1,2)])]
df_Erickson <- cbind(df_cut_Erickson$gene, df_cut_Erickson$bNum, mass_fractions_Erickson, df_flags_Erickson)
colnames(df_Erickson)[c(1:2,(ncol(df_Erickson)-word_type+1):ncol(df_Erickson))] <- colLabels
write.csv(df_Erickson, "./assigned_function/erickson2017_data_assigned.csv")


#============= Valgepea2013 ======================
df_Valgepea <- read.csv("valgepea2013_data.csv",
                        skip = 11, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_cut_Valgepea <- df_Valgepea[,c(1,2,19:23,82)]
df_cut_Valgepea$KEGG.ID <- paste0("eco:", df_Valgepea$KEGG.ID, sep="")
colnames(df_cut_Valgepea)[1:2] <- c("bNum","gene")

df_cut_Valgepea <- cbind(df_cut_Valgepea$bNum, df_cut_Valgepea$gene, data.frame(lapply(df_cut_Valgepea[,c(3:7)], function(x) x*kDa_convert*df_cut_Valgepea$Protein.MW)))
colnames(df_cut_Valgepea)[1:2] <- c("bNum","gene")

class_flags_Valgepea <- lapply(df_cut_Valgepea$bNum, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_flags_Valgepea <- as.data.frame(t(as.data.frame(class_flags_Valgepea)))
rownames(df_flags_Valgepea) <- c()

# Calculate mass fractions of each reported protein
mass_fractions_Valgepea <- df_cut_Valgepea[,-c(1,2)]/colSums(df_cut_Valgepea[,-c(1,2)], na.rm = TRUE)[col(df_cut_Valgepea[,-c(1,2)])]
df_Valgepea <- cbind(df_cut_Valgepea$gene, df_cut_Valgepea$bNum, mass_fractions_Valgepea, df_flags_Valgepea)
colnames(df_Valgepea)[c(1:2,(ncol(df_Valgepea)-word_type+1):ncol(df_Valgepea))] <- colLabels
write.csv(df_Valgepea, "./assigned_function/valgepea2013_data_assigned.csv")


#============= Li2014 ======================
all_genes <- keggList("eco")
all_genes <- data.frame(unname(all_genes), names(all_genes))
all_genes <- cbind(unlist(lapply(all_genes$unname.all_genes., function(x) sub(';.*$','', x))), data.frame(all_genes$names.all_genes.))
colnames(all_genes) <- c("gene","bNum")

df_Li <- read.csv("li2014_data.csv",
                  skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_Li <- cbind(df_Li[,1], data.frame(lapply(df_Li[,-1], function(x) as.numeric(x))))
colnames(df_Li)[1] <- "gene"

tmp_Peebo <- read.csv("peebo2015_data.csv",
                      skip = 4, header = TRUE, sep = ",", stringsAsFactors = FALSE)
tmp_Peebo <- select(tmp_Peebo, Molecular.weight..kDa., B.number.identifier)
tmp_Peebo$B.number.identifier <- paste0("eco:", tmp_Peebo$B.number.identifie, sep="")
colnames(tmp_Peebo) <- c("MW","bNum")

df_Li <- merge(all_genes, df_Li, by = "gene")
df_cut_Li <- merge(tmp_Peebo, df_Li, by = "bNum")
df_cut_Li <- cbind(df_cut_Li[,c(1:3)], data.frame(lapply(df_cut_Li[,-c(1:3)], function(x) x*df_cut_Li$MW)))

class_flags_Li <- lapply(df_cut_Li$bNum, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_flags_Li <- as.data.frame(t(as.data.frame(class_flags_Li)))
rownames(df_flags_Li) <- c()

# Calculate mass fractions of each reported protein
mass_fractions_Li <- df_cut_Li[,-c(1:3)]/colSums(df_cut_Li[,-c(1:3)], na.rm = TRUE)[col(df_cut_Li[,-c(1:3)])]
df_Li <- cbind(df_cut_Li$gene, df_cut_Li$bNum, mass_fractions_Li, df_flags_Li)
colnames(df_Li)[c(1:2,(ncol(df_Li)-word_type+1):ncol(df_Li))] <- colLabels
write.csv(df_Li, "./assigned_function/li2014_data_assigned.csv")


#============= Mori2021 ======================
df_Mori <- read.csv("mori2021_data.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_Mori$bNum <- paste0("eco:",df_Mori$Gene.locus)

# class_flags_ is a list containing BRITE flags for each protein species;
class_flags_Mori <- lapply(df_Mori$bNum, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
df_flags_Mori <- as.data.frame(t(as.data.frame(class_flags_Mori)))
rownames(df_flags_Mori) <- c()
colnames(df_flags_Mori) <- c("lipo_genes","ribo_genes")

df_cut_Mori <- df_Mori[,-c(1:3)]
df_cut_Mori <- cbind(df_cut_Mori, df_flags_Mori)
write.csv(df_cut_Mori, "./assigned_function/mori2021_data_assigned.csv")
