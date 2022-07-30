#=============== USAGE =====================
# Purpose: Muller et al. 2020 data --> total mass fractions for each sector across all bacterial species
# Input: folder with Muller2020 proteomes
# Output: muller2020_processed.csv
# Overall algorithm:
#     1) Convert all Uniprot ids into KEGG for each entry
#     2) Use these to assign each protein into one of the three sectors
#     3) Compute mass fractions of each protein and sum across members of each sector

remove(list=ls())

library(KEGGREST)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

setwd("/home/bogi/Desktop/growth_rate/data/bacterial_proteomes/Muller2020/")

word_type <- 2
sector_classification <- function(kegg_entry) {
  # If non-existing entry is passed, set return to zero
  if(is.null(kegg_entry)) {return(rep(0,word_type+1))}
  # Set the keywords used to match BRITE hierarchy
  lipid_words <- c("\\bFatty acid biosynthesis\\b", "\\bPeptidoglycan biosynthesis\\b", "\\bLipopolysaccharide biosynthesis\\b")
  ribo_words <- c("\\bRibosomal protein")
  
  lip_flag <- sum(unlist(lapply(lipid_words, function(x) sum(grepl(x, kegg_entry[[1]]$BRITE)))))>0
  ribo_flag <- sum(unlist(lapply(ribo_words, function(x) sum(grepl(x, kegg_entry[[1]]$BRITE)))))>0
  
  c(lip_flag, ribo_flag)
}

# Reported intensities in the paper are already mass fractions, so we are not accounting for protein length or mass
get_frac <- function(df, species, flag) {
  x <- subset(df, df$Organism==species)
  sum(x$Intensity*x[,flag], na.rm = TRUE)/sum(x$Intensity, na.rm = TRUE)
}


# ===== Filter only bacteria =====
# Subset only bacterial species from Muller2020
df_species <- read.csv("./table1_species_list_Muller.csv", header = TRUE, stringsAsFactors = FALSE)
species_list <- df_species[df_species$Superkingdom=="Bacteria",]$Spezies

df_prot <- read.csv("./table2_intensities_data_Muller.csv", header = TRUE, stringsAsFactors = FALSE)
write.csv(as.data.frame(species_list),"species_Muller.csv")

df_species <- read.csv("./species_Muller.csv", header = TRUE, stringsAsFactors = FALSE)

df <- df_prot[df_prot$Organism %in% df_species$species_list,]


# ===== Convert Uniprot to KEGG ids =====
# Each entry initially contains multiple Uniprot ids, so we convert each one into KEGG
name_list <- lapply(df$Protein_IDs, 
                    function(x) paste(
                      tryCatch(
                          keggConv("genes", paste0("up:", strsplit(x, split=';', fixed=TRUE)[[1]], sep=""))
                        ), 
                        collapse=";"
                      )
                   )
df$kegg_id <- unlist(name_list)
df <- df[,-1]
write.csv(df, "table2_intensities_data_Muller_bacteria_only_kegg.csv", row.names = FALSE)


# ===== Get BRITE hierarchies =====
df <- read.csv("table2_intensities_data_Muller_bacteria_only_kegg.csv", header = TRUE, stringsAsFactors = FALSE)
class_list <- lapply(
  df$kegg_id, function(x) 
  Reduce(`|`, lapply(unlist(strsplit(x, ";")), function(y) sector_classification(tryCatch(keggGet(y), error = function(c) NULL))))
)

df$lipo_flag <- lapply(c(1:nrow(df)), function(x) unlist(class_list[x][[1]][1]))
df$lipo_flag[sapply(df$lipo_flag, is.null)] <- NA
df$lipo_flag <- unlist(df$lipo_flag)

df$ribo_flag <- lapply(c(1:nrow(df)), function(x) unlist(class_list[x][[1]][2]))
df$ribo_flag[sapply(df$ribo_flag, is.null)] <- NA
df$ribo_flag <- unlist(df$ribo_flag)

write.csv(df, "table2_intensities_data_Muller_bacteria_kegg_flagged.csv", row.names = FALSE)


# ===== Calculate mass fractions of each sector =====
df <- read.csv("table2_intensities_data_Muller_bacteria_kegg_flagged.csv", header = TRUE)

lip_frac <- unlist(lapply(unique(df$Organism), function(x) get_frac(df, x, "lipo_flag")))
ribo_frac <- unlist(lapply(unique(df$Organism), function(x) get_frac(df, x, "ribo_flag")))
tot_intensity <- aggregate(df$Intensity, by=list(df$Organism), FUN = sum)

df_out <- as.data.frame(cbind(as.character(unique(df$Organism)), lip_frac, ribo_frac, tot_intensity))
colnames(df_out)[c(1,5)] <- c("species","tot_intensity")
df_out <- df_out[,-4]
write.csv(df_out, "../../main/muller2020_processed.csv", row.names = FALSE)
