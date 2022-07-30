#=============== USAGE =====================
# Purpose: PaxDb bacterial proteomes --> total mass fractions for each sector across all species
# Input: folder with PaxDb proteomes
# Output: paxdb_fractions.csv in the following columnar format 
# "species"	"lip_frac"	"ribo_frac"
# Overall algorithm:
#     1) Convert STRINGdb used in PaxDb into Uniprot ids
#     2) Convert Uniprot to KEGG ids, and use these to assign each protein into one of the three sectors
#     3) Compute mass fractions of each protein and sum across members of each sector

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("KEGGREST")
#BiocManager::install("STRINGdb")
#BiocManager::install("pathview")

remove(list=ls())
setwd("/home/bogi/Desktop/growth_rate/data/bacterial_proteomes/bacterial_proteomes_PaxDb/")

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
    prot_mw <- 0
    prot_mw <- mw(kegg_entry[[1]]$AASEQ[[1]])

    c(lip_flag, ribo_flag, prot_mw)
}


# =================== Mapping STRINGdb to KEGG id ===================
# PaxDb retains STRINGdb ids for each protein entry; This is part of their effort to estimate the qulity
#   of the proteomic data by inspecting stoichiometry in protein-protein interactions.
# Here, we convert STRINGdb to Uniprot by using mappings file
df_map <- read.csv("./map/full_uniprot_2_paxdb.04.2015.tsv",
                     skip = 0, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
df_string_id <- as.data.table(with(df_map, paste0(V1, ".", V3)))
df_map <- cbind(df_string_id, df_map$V2, df_map$V4, df_map$V5)
colnames(df_map) <- c("string_id", "uniprot", "identity", "bit_score")
# Multiple proteins in STRINGdb can be be mapped onto the same Uniport entry
# For each entry, we retained a single Uniprot id with the best alignment quality (i.e., highest bitscore)
df_map <- df_map[,.SD[which.max(bit_score)], by=string_id]

# With conversion map, we can now functionally annotate all PaxDb entries
all_files <- list.files(pattern = "\\.tsv$")

for(i in 1:length(all_files))
{

  df_abund <- read.csv(all_files[[i]],
                       skip = 12, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df_abund) <- c("internal_id","string_id","abundance")
  df_abund <- df_abund[,c(1:3)]
  print(all_files[[i]])
  name_list <- vector(mode = "list", length = nrow(df_abund))
  
  # Merge abundances with map by mathcing STRINGdb numbers
  df_abund <- merge(df_abund, df_map, by="string_id", all.x = TRUE)
  # Convert Uniprot number into accession and id part by separating at "|"
  df_abund <- separate(data = df_abund, col = uniprot, sep = "\\|", into = c("uniprot_ac", "uniprot_id"))
  # Prepend "up:" required for conversion into KEGG through keggConv()
  df_abund$uniprot_ac <- paste0("up:", df_abund$uniprot_ac)
  df_abund <- df_abund[!is.na(df_abund$identity),]
  
  # name_list holds Uniprot 
  name_list <- tryCatch(keggConv("genes", df_abund$uniprot_ac), error = function(e) NULL)
  # Some organisms might not be represented in KEGG, hence the condition below
  if (length(name_list)>0)
  {
    # Puts all the KEGG ids into a single column, and merges this with abundance data
    # Output is saved with .kegged.csv extension; These files are used subsequently for functional annotation
    kegg_stack <- stack(name_list)
    kegg_stack <- tryCatch(stack(name_list), error = function(e) data.frame(t(c(NA,NA))))
    colnames(kegg_stack) <- c("kegg_id", "uniprot_ac")
    
    df_abund <- merge(df_abund, kegg_stack, by = "uniprot_ac", all.x = TRUE)
    write.csv(df_abund, file=paste(all_files[[i]], ".kegged.csv", sep = ""), row.names = FALSE)
    print(paste("processed: ", all_files[[i]], sep = "")) 
  } else print(paste("failed: ", all_files[[i]], sep = ""))
  
}


# =================== KEGG functional annotation ===================
all_files <- list.files(pattern = "\\.kegged.csv$")

for(i in 1:length(all_files))
{
  df_abund <- read.csv(all_files[[i]], skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  print(all_files[[i]])
  class_list <- lapply(df_abund$kegg_id, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))
  df_class <- as.data.frame(t(as.data.frame(class_list)))
  df_abund <- cbind(df_abund,df_class)
  rownames(df_abund) <- c()
  colnames(df_abund)[c(9:11)] <- c("lip_frac","ribo_frac","prot_mw")
  
  # Convert relative abundance into mass fractions by normalizing each entry by the protein length
  df_abund$tmp <- df_abund$prot_mw*df_abund$abundance
  tot_rel_mass <- sum(df_abund$tmp, na.rm = TRUE)
  df_abund$abundance <- df_abund$tmp/tot_rel_mass
  # Remove Uniprot entries not having a corresponding KEGG engtry, and exclude temporary column
  df_abund <- df_abund[!is.na(df_abund$kegg_id),]
  df_abund <- df_abund[ , -which(names(df_abund) %in% c("tmp"))]
  
  # Output of functional annotation is saved with .processed.csv extension
  write.csv(df_abund, file=paste(all_files[[i]], ".processed.csv", sep=""), row.names = FALSE)
  print(paste("processed: ", all_files[[i]], sep = "")) 
}


# =================== Mass fraction computation ===================
# With functionally annotated proteomes, we calculate mass fractions of three sectors in our model
all_files <- list.files(pattern = "\\.processed.csv$")

fractions <- as.data.frame(matrix(data = 0, nrow = length(all_files), ncol = word_type+1))
colnames(fractions) <- c("species","lip_frac", "ribo_frac")

for(i in 1:length(all_files))
{
  df_abund <- read.csv(all_files[[i]],
                       skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  df_l <- df_abund[df_abund$lip_frac==TRUE,]
  df_r <- df_abund[df_abund$ribo_frac==TRUE,]
  
  species_name <- paste0(strsplit(gsub("\\..*$", "", all_files[[i]]), "[_]")[[1]][1]," ",strsplit(gsub("\\..*$", "", all_files[[i]]), "[_]")[[1]][2])
  fractions[i,] <- t(c(species_name, sum(df_l$abundance), sum(df_r$abundance)))
}

# Final output has the following columns of proteomic mass fractions:
# "species"	"lip_frac"	"ribo_frac"
fractions$species <- word(fractions$species, 1, 2, sep=" ")
write.csv(fractions, file="../../main/paxdb_proteomes.csv", row.names = FALSE)



