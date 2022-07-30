#=============== USAGE =====================
# Purpose: Merges the proteomic data across the bacteria into a single csv file
# Input: folder with all the proteomes produced by
#       muller_proteomes.R
#       paxDB_proteomes.R
#       added_proteomes.R
# Output: proteome_scaling.size.csv

remove(list=ls())
setwd("/home/bogi/Desktop/growth_rate/data/main/")


#========== Merging with proteomic data ==========
in_file_proteome <- "paxdb_proteomes.csv"
df_all <- read.csv(in_file_proteome, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_all$source <- "Wang2012"

in_file_proteome <- "osbak2016.processed.csv"
df_osbak <- read.csv(in_file_proteome, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_osbak$source <- "Osbak2016"

in_file_proteome <- "angel2010.processed.csv"
df_angel <- read.csv(in_file_proteome, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_angel$source <- "Angel2010"

in_file_proteome <- "srivastava2020.processed.csv"
df_srivastava <- read.csv(in_file_proteome, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_srivastava$source <- "Srivastava2020"

in_file_proteome <- "matteau2020.processed.csv"
df_matteau <- read.csv(in_file_proteome, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_matteau$source <- "Matteau2020"

in_file_proteome <- "masson2021.processed.csv"
df_masson <- read.csv(in_file_proteome, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_masson$source <- "Masson2021"

in_file_proteome <- "muller2020_processed.csv"
df_muller <- read.csv(in_file_proteome, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_muller <- df_muller[(df_muller$lip_frac!=0&df_muller$ribo_frac!=0),]
df_muller$source <- "Muller2020"
df_muller <- subset(df_muller, select = -c(tot_intensity))

df_all <- rbind(df_all, df_muller, df_osbak, df_angel, df_srivastava, df_matteau, df_masson)
write.csv(df_all, "quant_proteome_data.csv", row.names = FALSE)
