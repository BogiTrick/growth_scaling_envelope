remove(list=ls())
setwd("/home/bogi/Desktop/growth_rate/data/favate2021/")

library(KEGGREST)
library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
library(Peptides)

axes_style <- element_text(size=15, face="bold", colour = "black")

# Function sector_classification takes KEGG entry and parses BRITE hierarchy for specific keywords.
# Then, it sets the flags for each function if the appropriate designation is present
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
  if(!is.null(kegg_entry[[1]]$AASEQ[[1]])) {prot_mw <- mw(kegg_entry[[1]]$AASEQ[[1]])} else {prot_mw <- NA}
  if(!is.null(kegg_entry[[1]]$SYMBOL)) {gene_name <- kegg_entry[[1]]$SYMBOL} else {gene_name <- NA}

  c(lip_flag, ribo_flag, prot_mw, gene_name)
}

get_mass_fractions <- function(df) {
  lip_part <- df[as.logical(df$lip_genes),]
  fractions_lipid <- sum(lip_part$mass_frac, na.rm = TRUE)
  
  ribo_part <- df[as.logical(df$ribo_genes),]
  fractions_ribo <- sum(ribo_part$mass_frac, na.rm = TRUE)
  
  c(fractions_lipid, fractions_ribo)
}


#=========== Assign the proteins ===============
df_Favate <- read.csv("favate2021.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = TRUE)
df_Favate$target_id <- paste("ebr:",df_Favate$target_id,sep="")

df_cut_Favate <- df_Favate[(df_Favate$seqtype=="ribo")&(df_Favate$repl=="rep1"),]

# class_flags_ is a list containing BRITE flags for each protein species;
class_flags_Favate <- lapply(df_cut_Favate$target_id, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))

df_flags_Favate <- as.data.frame(t(stri_list2matrix(class_flags_Favate)))
colnames(df_flags_Favate) <- c("lip_genes","ribo_genes","mw","gene_name")
rownames(df_flags_Favate) <- c()

# Save the assigned data for later manipulation
df_Favate_new <- cbind(df_cut_Favate, df_flags_Favate)
write.csv(df_Favate_new, "./favate2021_data_assigned.csv")


#=========== Load full data ===============
# Load Valgepea2013 data and get the gene names from KEGG. We use mRNA abundance in this paper to convert Favate2021 RNA abundance into protein abundance
df_Valgepea <- read.csv("valgepea2013_data.csv", skip = 11, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_Valgepea$KEGG.ID <- paste0("eco:", df_Valgepea$KEGG.ID, sep="")
colnames(df_Valgepea)[1:2] <- c("bNum","gene")
class_flags_Valgepea <- lapply(df_Valgepea$bNum, function(x) sector_classification(tryCatch(keggGet(x), error = function(c) NULL)))

# Replace a few entries that did not match in KEGG (mostly hypothetical proteins) with NA vector
empty_entries_flag <- unlist(lapply(seq(1,length(class_flags_Valgepea),1), function(x) length(class_flags_Valgepea[[x]])))!=4
empty_entries_pos <- which(empty_entries_flag %in% TRUE)
for(i in 1:length(empty_entries_pos)) {class_flags_Valgepea[empty_entries_pos[i]][[1]] <- c(NA, NA, NA, NA)}
df_flags_Valgepea <- as.data.frame(t(as.data.frame(class_flags_Valgepea)))
rownames(df_flags_Valgepea) <- c()
df_Valgepea <- cbind(df_Valgepea,df_flags_Valgepea)
write.csv(df_Valgepea, "./valgepea2013_assigned_modified.csv")

df_Valgepea <- read.csv("./valgepea2013_assigned_modified.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_Favate_all <- read.csv("favate2021.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_Favate_all$target_id <- paste("ebr:",df_Favate_all$target_id,sep="")
df_Favate_assigned <- read.csv("favate2021_data_assigned.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)[,-1]
# Given that only genes from "ribo"/"rep1" have been assigned classification, we merge this with other replicates and methods.
# This assumes that the same genes picked out in the experiment above, were picked up in other experiments. The number of gene entires are identical across all experiments.
df_Favate_all_assigned <- merge(df_Favate_all, df_Favate_assigned, by=c("target_id","line"), all.x = TRUE, all.y = FALSE)
df_Favate_all_assigned <- merge(df_Favate_all_assigned, df_Valgepea, by.x="gene_name", by.y = "V4", all.x = TRUE, all.y = FALSE)

# The conversion factor used to convert mRNA abundance into protein abundance is calculated in fastest growing media in the paper, with \lambda=0.48 per hour
# This is still slower than the ancestral growth rate in LTEE, which is 0.77 per hour (See Table 1 in Vasi1994)
# In other words, 0.48 is closest to 0.77, hence our choice
df_Favate_all_assigned$conv_factor <- df_Favate_all_assigned$Protein.abun..4/df_Favate_all_assigned$mRNA.abun..4
df_Favate_all_assigned$prot_mass <- df_Favate_all_assigned$conv_factor*df_Favate_all_assigned$tpm.x*as.numeric(df_Favate_all_assigned$mw)

# The entries that are missing molecular mass are tRNA, and we remove these from analysis
test_miss_mw <- as.numeric(df_Favate_all_assigned$mw)
df_Favate_all_assigned <- df_Favate_all_assigned[!is.na(test_miss_mw),]
test_miss_name <- as.character(df_Favate_all_assigned$gene_name)
df_Favate_all_assigned<- df_Favate_all_assigned[!is.na(test_miss_name),]


# Going through each replicate and each method used, calculate the proteome mass fractions
#=========== Replicate 1, Ribo-seq ===============
df_Favate_new <- df_Favate_all_assigned[(df_Favate_all_assigned$seqtype.x=="ribo")&(df_Favate_all_assigned$repl.x=="rep1"),]
nrow(df_Favate_new)
df_Favate_Aram1 <- df_Favate_new[df_Favate_new$line=="Ara-1",]
df_Favate_Aram2 <- df_Favate_new[df_Favate_new$line=="Ara-2",]
df_Favate_Aram3 <- df_Favate_new[df_Favate_new$line=="Ara-3",]
df_Favate_Aram4 <- df_Favate_new[df_Favate_new$line=="Ara-4",]
df_Favate_Aram5 <- df_Favate_new[df_Favate_new$line=="Ara-5",]
df_Favate_Arap1 <- df_Favate_new[df_Favate_new$line=="Ara+1",]
df_Favate_Arap2 <- df_Favate_new[df_Favate_new$line=="Ara+2",]
df_Favate_Arap3 <- df_Favate_new[df_Favate_new$line=="Ara+3",]
df_Favate_Arap4 <- df_Favate_new[df_Favate_new$line=="Ara+4",]
df_Favate_Arap5 <- df_Favate_new[df_Favate_new$line=="Ara+5",]
df_Favate_REL606 <- df_Favate_new[df_Favate_new$line=="REL606",]
df_Favate_REL607 <- df_Favate_new[df_Favate_new$line=="REL607",]

df_Favate_Aram1$mass_frac <- df_Favate_Aram1$prot_mass/sum(df_Favate_Aram1$prot_mass, na.rm = TRUE)
df_Favate_Aram2$mass_frac <- df_Favate_Aram2$prot_mass/sum(df_Favate_Aram2$prot_mass, na.rm = TRUE)
df_Favate_Aram3$mass_frac <- df_Favate_Aram3$prot_mass/sum(df_Favate_Aram3$prot_mass, na.rm = TRUE)
df_Favate_Aram4$mass_frac <- df_Favate_Aram4$prot_mass/sum(df_Favate_Aram4$prot_mass, na.rm = TRUE)
df_Favate_Aram5$mass_frac <- df_Favate_Aram5$prot_mass/sum(df_Favate_Aram5$prot_mass, na.rm = TRUE)
df_Favate_Arap1$mass_frac <- df_Favate_Arap1$prot_mass/sum(df_Favate_Arap1$prot_mass, na.rm = TRUE)
df_Favate_Arap2$mass_frac <- df_Favate_Arap2$prot_mass/sum(df_Favate_Arap2$prot_mass, na.rm = TRUE)
df_Favate_Arap3$mass_frac <- df_Favate_Arap3$prot_mass/sum(df_Favate_Arap3$prot_mass, na.rm = TRUE)
df_Favate_Arap4$mass_frac <- df_Favate_Arap4$prot_mass/sum(df_Favate_Arap4$prot_mass, na.rm = TRUE)
df_Favate_Arap5$mass_frac <- df_Favate_Arap5$prot_mass/sum(df_Favate_Arap5$prot_mass, na.rm = TRUE)
df_Favate_REL606$mass_frac <- df_Favate_REL606$prot_mass/sum(df_Favate_REL606$prot_mass, na.rm = TRUE)
df_Favate_REL607$mass_frac <- df_Favate_REL607$prot_mass/sum(df_Favate_REL607$prot_mass, na.rm = TRUE)

anc_araM_lip_rep1_riboSeq <- get_mass_fractions(df_Favate_REL606)[1]
evo_araM_lip_rep1_riboSeq <- c(get_mass_fractions(df_Favate_Aram1)[1],
                  get_mass_fractions(df_Favate_Aram2)[1],
                  get_mass_fractions(df_Favate_Aram3)[1],
                  get_mass_fractions(df_Favate_Aram4)[1],
                  get_mass_fractions(df_Favate_Aram5)[1])

anc_araP_lip_rep1_riboSeq <- get_mass_fractions(df_Favate_REL607)[1]
evo_araP_lip_rep1_riboSeq <- c(get_mass_fractions(df_Favate_Arap1)[1],
                  get_mass_fractions(df_Favate_Arap2)[1],
                  get_mass_fractions(df_Favate_Arap3)[1],
                  get_mass_fractions(df_Favate_Arap4)[1],
                  get_mass_fractions(df_Favate_Arap5)[1])

anc_araM_ribo_rep1_riboSeq <- get_mass_fractions(df_Favate_REL606)[2]
evo_araM_ribo_rep1_riboSeq <- c(get_mass_fractions(df_Favate_Aram1)[2],
                   get_mass_fractions(df_Favate_Aram2)[2],
                   get_mass_fractions(df_Favate_Aram3)[2],
                   get_mass_fractions(df_Favate_Aram4)[2],
                   get_mass_fractions(df_Favate_Aram5)[2])

anc_araP_ribo_rep1_riboSeq <- get_mass_fractions(df_Favate_REL607)[2]
evo_araP_ribo_rep1_riboSeq <- c(get_mass_fractions(df_Favate_Arap1)[2],
                      get_mass_fractions(df_Favate_Arap2)[2],
                      get_mass_fractions(df_Favate_Arap3)[2],
                      get_mass_fractions(df_Favate_Arap4)[2],
                      get_mass_fractions(df_Favate_Arap5)[2])


#=========== Replicate 2, Ribo-seq ===============
df_Favate_new <- df_Favate_all_assigned[(df_Favate_all_assigned$repl.x=="rep2")&(df_Favate_all_assigned$seqtype.x=="ribo"),]
nrow(df_Favate_new)
df_Favate_Aram1 <- df_Favate_new[df_Favate_new$line=="Ara-1",]
df_Favate_Aram2 <- df_Favate_new[df_Favate_new$line=="Ara-2",]
df_Favate_Aram3 <- df_Favate_new[df_Favate_new$line=="Ara-3",]
df_Favate_Aram4 <- df_Favate_new[df_Favate_new$line=="Ara-4",]
df_Favate_Aram5 <- df_Favate_new[df_Favate_new$line=="Ara-5",]
df_Favate_Arap1 <- df_Favate_new[df_Favate_new$line=="Ara+1",]
df_Favate_Arap2 <- df_Favate_new[df_Favate_new$line=="Ara+2",]
df_Favate_Arap3 <- df_Favate_new[df_Favate_new$line=="Ara+3",]
df_Favate_Arap4 <- df_Favate_new[df_Favate_new$line=="Ara+4",]
df_Favate_Arap5 <- df_Favate_new[df_Favate_new$line=="Ara+5",]
df_Favate_REL606 <- df_Favate_new[df_Favate_new$line=="REL606",]
df_Favate_REL607 <- df_Favate_new[df_Favate_new$line=="REL607",]

df_Favate_Aram1$mass_frac <- df_Favate_Aram1$prot_mass/sum(df_Favate_Aram1$prot_mass, na.rm = TRUE)
df_Favate_Aram2$mass_frac <- df_Favate_Aram2$prot_mass/sum(df_Favate_Aram2$prot_mass, na.rm = TRUE)
df_Favate_Aram3$mass_frac <- df_Favate_Aram3$prot_mass/sum(df_Favate_Aram3$prot_mass, na.rm = TRUE)
df_Favate_Aram4$mass_frac <- df_Favate_Aram4$prot_mass/sum(df_Favate_Aram4$prot_mass, na.rm = TRUE)
df_Favate_Aram5$mass_frac <- df_Favate_Aram5$prot_mass/sum(df_Favate_Aram5$prot_mass, na.rm = TRUE)
df_Favate_Arap1$mass_frac <- df_Favate_Arap1$prot_mass/sum(df_Favate_Arap1$prot_mass, na.rm = TRUE)
df_Favate_Arap2$mass_frac <- df_Favate_Arap2$prot_mass/sum(df_Favate_Arap2$prot_mass, na.rm = TRUE)
df_Favate_Arap3$mass_frac <- df_Favate_Arap3$prot_mass/sum(df_Favate_Arap3$prot_mass, na.rm = TRUE)
df_Favate_Arap4$mass_frac <- df_Favate_Arap4$prot_mass/sum(df_Favate_Arap4$prot_mass, na.rm = TRUE)
df_Favate_Arap5$mass_frac <- df_Favate_Arap5$prot_mass/sum(df_Favate_Arap5$prot_mass, na.rm = TRUE)
df_Favate_REL606$mass_frac <- df_Favate_REL606$prot_mass/sum(df_Favate_REL606$prot_mass, na.rm = TRUE)
df_Favate_REL607$mass_frac <- df_Favate_REL607$prot_mass/sum(df_Favate_REL607$prot_mass, na.rm = TRUE)

anc_araM_lip_rep2_riboSeq <- get_mass_fractions(df_Favate_REL606)[1]
evo_araM_lip_rep2_riboSeq <- c(get_mass_fractions(df_Favate_Aram1)[1],
                  get_mass_fractions(df_Favate_Aram2)[1],
                  get_mass_fractions(df_Favate_Aram3)[1],
                  get_mass_fractions(df_Favate_Aram4)[1],
                  get_mass_fractions(df_Favate_Aram5)[1])

anc_araP_lip_rep2_riboSeq <- get_mass_fractions(df_Favate_REL607)[1]
evo_araP_lip_rep2_riboSeq <- c(get_mass_fractions(df_Favate_Arap1)[1],
                  get_mass_fractions(df_Favate_Arap2)[1],
                  get_mass_fractions(df_Favate_Arap3)[1],
                  get_mass_fractions(df_Favate_Arap4)[1],
                  get_mass_fractions(df_Favate_Arap5)[1])

anc_araM_ribo_rep2_riboSeq <- get_mass_fractions(df_Favate_REL606)[2]
evo_araM_ribo_rep2_riboSeq <- c(get_mass_fractions(df_Favate_Aram1)[2],
                   get_mass_fractions(df_Favate_Aram2)[2],
                   get_mass_fractions(df_Favate_Aram3)[2],
                   get_mass_fractions(df_Favate_Aram4)[2],
                   get_mass_fractions(df_Favate_Aram5)[2])

anc_araP_ribo_rep2_riboSeq <- get_mass_fractions(df_Favate_REL607)[2]
evo_araP_ribo_rep2_riboSeq <- c(get_mass_fractions(df_Favate_Arap1)[2],
                      get_mass_fractions(df_Favate_Arap2)[2],
                      get_mass_fractions(df_Favate_Arap3)[2],
                      get_mass_fractions(df_Favate_Arap4)[2],
                      get_mass_fractions(df_Favate_Arap5)[2])


#=========== Replicate 1, RNA-seq ===============
df_Favate_new <- df_Favate_all_assigned[(df_Favate_all_assigned$repl.x=="rep1")&(df_Favate_all_assigned$seqtype.x=="rna"),]
nrow(df_Favate_new)
df_Favate_Aram1 <- df_Favate_new[df_Favate_new$line=="Ara-1",]
df_Favate_Aram2 <- df_Favate_new[df_Favate_new$line=="Ara-2",]
df_Favate_Aram3 <- df_Favate_new[df_Favate_new$line=="Ara-3",]
df_Favate_Aram4 <- df_Favate_new[df_Favate_new$line=="Ara-4",]
df_Favate_Aram5 <- df_Favate_new[df_Favate_new$line=="Ara-5",]
df_Favate_Arap1 <- df_Favate_new[df_Favate_new$line=="Ara+1",]
df_Favate_Arap2 <- df_Favate_new[df_Favate_new$line=="Ara+2",]
df_Favate_Arap3 <- df_Favate_new[df_Favate_new$line=="Ara+3",]
df_Favate_Arap4 <- df_Favate_new[df_Favate_new$line=="Ara+4",]
df_Favate_Arap5 <- df_Favate_new[df_Favate_new$line=="Ara+5",]
df_Favate_REL606 <- df_Favate_new[df_Favate_new$line=="REL606",]
df_Favate_REL607 <- df_Favate_new[df_Favate_new$line=="REL607",]

df_Favate_Aram1$mass_frac <- df_Favate_Aram1$prot_mass/sum(df_Favate_Aram1$prot_mass, na.rm = TRUE)
df_Favate_Aram2$mass_frac <- df_Favate_Aram2$prot_mass/sum(df_Favate_Aram2$prot_mass, na.rm = TRUE)
df_Favate_Aram3$mass_frac <- df_Favate_Aram3$prot_mass/sum(df_Favate_Aram3$prot_mass, na.rm = TRUE)
df_Favate_Aram4$mass_frac <- df_Favate_Aram4$prot_mass/sum(df_Favate_Aram4$prot_mass, na.rm = TRUE)
df_Favate_Aram5$mass_frac <- df_Favate_Aram5$prot_mass/sum(df_Favate_Aram5$prot_mass, na.rm = TRUE)
df_Favate_Arap1$mass_frac <- df_Favate_Arap1$prot_mass/sum(df_Favate_Arap1$prot_mass, na.rm = TRUE)
df_Favate_Arap2$mass_frac <- df_Favate_Arap2$prot_mass/sum(df_Favate_Arap2$prot_mass, na.rm = TRUE)
df_Favate_Arap3$mass_frac <- df_Favate_Arap3$prot_mass/sum(df_Favate_Arap3$prot_mass, na.rm = TRUE)
df_Favate_Arap4$mass_frac <- df_Favate_Arap4$prot_mass/sum(df_Favate_Arap4$prot_mass, na.rm = TRUE)
df_Favate_Arap5$mass_frac <- df_Favate_Arap5$prot_mass/sum(df_Favate_Arap5$prot_mass, na.rm = TRUE)
df_Favate_REL606$mass_frac <- df_Favate_REL606$prot_mass/sum(df_Favate_REL606$prot_mass, na.rm = TRUE)
df_Favate_REL607$mass_frac <- df_Favate_REL607$prot_mass/sum(df_Favate_REL607$prot_mass, na.rm = TRUE)

anc_araM_lip_rep1_rnaSeq <- get_mass_fractions(df_Favate_REL606)[1]
evo_araM_lip_rep1_rnaSeq <- c(get_mass_fractions(df_Favate_Aram1)[1],
                  get_mass_fractions(df_Favate_Aram2)[1],
                  get_mass_fractions(df_Favate_Aram3)[1],
                  get_mass_fractions(df_Favate_Aram4)[1],
                  get_mass_fractions(df_Favate_Aram5)[1])

anc_araP_lip_rep1_rnaSeq <- get_mass_fractions(df_Favate_REL607)[1]
evo_araP_lip_rep1_rnaSeq <- c(get_mass_fractions(df_Favate_Arap1)[1],
                  get_mass_fractions(df_Favate_Arap2)[1],
                  get_mass_fractions(df_Favate_Arap3)[1],
                  get_mass_fractions(df_Favate_Arap4)[1],
                  get_mass_fractions(df_Favate_Arap5)[1])

anc_araM_ribo_rep1_rnaSeq <- get_mass_fractions(df_Favate_REL606)[2]
evo_araM_ribo_rep1_rnaSeq <- c(get_mass_fractions(df_Favate_Aram1)[2],
                   get_mass_fractions(df_Favate_Aram2)[2],
                   get_mass_fractions(df_Favate_Aram3)[2],
                   get_mass_fractions(df_Favate_Aram4)[2],
                   get_mass_fractions(df_Favate_Aram5)[2])

anc_araP_ribo_rep1_rnaSeq <- get_mass_fractions(df_Favate_REL607)[2]
evo_araP_ribo_rep1_rnaSeq <- c(get_mass_fractions(df_Favate_Arap1)[2],
                      get_mass_fractions(df_Favate_Arap2)[2],
                      get_mass_fractions(df_Favate_Arap3)[2],
                      get_mass_fractions(df_Favate_Arap4)[2],
                      get_mass_fractions(df_Favate_Arap5)[2])


#=========== Replicate 2, RNA-seq ===============
df_Favate_new <- df_Favate_all_assigned[(df_Favate_all_assigned$repl.x=="rep2")&(df_Favate_all_assigned$seqtype.x=="rna"),]
nrow(df_Favate_new)
df_Favate_Aram1 <- df_Favate_new[df_Favate_new$line=="Ara-1",]
df_Favate_Aram2 <- df_Favate_new[df_Favate_new$line=="Ara-2",]
df_Favate_Aram3 <- df_Favate_new[df_Favate_new$line=="Ara-3",]
df_Favate_Aram4 <- df_Favate_new[df_Favate_new$line=="Ara-4",]
df_Favate_Aram5 <- df_Favate_new[df_Favate_new$line=="Ara-5",]
df_Favate_Arap1 <- df_Favate_new[df_Favate_new$line=="Ara+1",]
df_Favate_Arap2 <- df_Favate_new[df_Favate_new$line=="Ara+2",]
df_Favate_Arap3 <- df_Favate_new[df_Favate_new$line=="Ara+3",]
df_Favate_Arap4 <- df_Favate_new[df_Favate_new$line=="Ara+4",]
df_Favate_Arap5 <- df_Favate_new[df_Favate_new$line=="Ara+5",]
df_Favate_REL606 <- df_Favate_new[df_Favate_new$line=="REL606",]
df_Favate_REL607 <- df_Favate_new[df_Favate_new$line=="REL607",]

df_Favate_Aram1$mass_frac <- df_Favate_Aram1$prot_mass/sum(df_Favate_Aram1$prot_mass, na.rm = TRUE)
df_Favate_Aram2$mass_frac <- df_Favate_Aram2$prot_mass/sum(df_Favate_Aram2$prot_mass, na.rm = TRUE)
df_Favate_Aram3$mass_frac <- df_Favate_Aram3$prot_mass/sum(df_Favate_Aram3$prot_mass, na.rm = TRUE)
df_Favate_Aram4$mass_frac <- df_Favate_Aram4$prot_mass/sum(df_Favate_Aram4$prot_mass, na.rm = TRUE)
df_Favate_Aram5$mass_frac <- df_Favate_Aram5$prot_mass/sum(df_Favate_Aram5$prot_mass, na.rm = TRUE)
df_Favate_Arap1$mass_frac <- df_Favate_Arap1$prot_mass/sum(df_Favate_Arap1$prot_mass, na.rm = TRUE)
df_Favate_Arap2$mass_frac <- df_Favate_Arap2$prot_mass/sum(df_Favate_Arap2$prot_mass, na.rm = TRUE)
df_Favate_Arap3$mass_frac <- df_Favate_Arap3$prot_mass/sum(df_Favate_Arap3$prot_mass, na.rm = TRUE)
df_Favate_Arap4$mass_frac <- df_Favate_Arap4$prot_mass/sum(df_Favate_Arap4$prot_mass, na.rm = TRUE)
df_Favate_Arap5$mass_frac <- df_Favate_Arap5$prot_mass/sum(df_Favate_Arap5$prot_mass, na.rm = TRUE)
df_Favate_REL606$mass_frac <- df_Favate_REL606$prot_mass/sum(df_Favate_REL606$prot_mass, na.rm = TRUE)
df_Favate_REL607$mass_frac <- df_Favate_REL607$prot_mass/sum(df_Favate_REL607$prot_mass, na.rm = TRUE)

anc_araM_lip_rep2_rnaSeq <- get_mass_fractions(df_Favate_REL606)[1]
evo_araM_lip_rep2_rnaSeq <- c(get_mass_fractions(df_Favate_Aram1)[1],
                  get_mass_fractions(df_Favate_Aram2)[1],
                  get_mass_fractions(df_Favate_Aram3)[1],
                  get_mass_fractions(df_Favate_Aram4)[1],
                  get_mass_fractions(df_Favate_Aram5)[1])

anc_araP_lip_rep2_rnaSeq <- get_mass_fractions(df_Favate_REL607)[1]
evo_araP_lip_rep2_rnaSeq <- c(get_mass_fractions(df_Favate_Arap1)[1],
                  get_mass_fractions(df_Favate_Arap2)[1],
                  get_mass_fractions(df_Favate_Arap3)[1],
                  get_mass_fractions(df_Favate_Arap4)[1],
                  get_mass_fractions(df_Favate_Arap5)[1])

anc_araM_ribo_rep2_rnaSeq <- get_mass_fractions(df_Favate_REL606)[2]
evo_araM_ribo_rep2_rnaSeq <- c(get_mass_fractions(df_Favate_Aram1)[2],
                   get_mass_fractions(df_Favate_Aram2)[2],
                   get_mass_fractions(df_Favate_Aram3)[2],
                   get_mass_fractions(df_Favate_Aram4)[2],
                   get_mass_fractions(df_Favate_Aram5)[2])

anc_araP_ribo_rep2_rnaSeq <- get_mass_fractions(df_Favate_REL607)[2]
evo_araP_ribo_rep2_rnaSeq <- c(get_mass_fractions(df_Favate_Arap1)[2],
                      get_mass_fractions(df_Favate_Arap2)[2],
                      get_mass_fractions(df_Favate_Arap3)[2],
                      get_mass_fractions(df_Favate_Arap4)[2],
                      get_mass_fractions(df_Favate_Arap5)[2])


#============ Analysis ==============
# Individual
anc_araP_lip_rep1_rnaSeq
mean(evo_araP_lip_rep1_rnaSeq)

anc_araP_ribo_rep1_rnaSeq
mean(evo_araP_ribo_rep1_rnaSeq)

anc_araM_lip_rep1_riboSeq
mean(evo_araM_lip_rep1_riboSeq)

anc_araM_ribo_rep1_riboSeq
mean(evo_araM_ribo_rep1_riboSeq)


anc_araP_lip_rep2_rnaSeq
mean(evo_araP_lip_rep2_rnaSeq)

anc_araP_ribo_rep2_rnaSeq
mean(evo_araP_ribo_rep2_rnaSeq)

anc_araM_lip_rep2_riboSeq
mean(evo_araM_lip_rep2_riboSeq)

anc_araM_ribo_rep2_riboSeq
mean(evo_araM_ribo_rep2_riboSeq)

# Pool - Ara(+) vs. Ara(-) and anc vs. evo
mean(c(anc_araP_lip_rep1_rnaSeq,
       anc_araP_lip_rep1_riboSeq,
       anc_araP_lip_rep2_rnaSeq,
       anc_araP_lip_rep2_riboSeq))
mean(c(evo_araP_lip_rep1_rnaSeq,
       evo_araP_lip_rep1_riboSeq,
       evo_araP_lip_rep2_rnaSeq,
       evo_araP_lip_rep2_riboSeq))

mean(c(anc_araP_ribo_rep1_rnaSeq,
       anc_araP_ribo_rep1_riboSeq,
       anc_araP_ribo_rep2_rnaSeq,
       anc_araP_ribo_rep2_riboSeq))
mean(c(evo_araP_ribo_rep1_rnaSeq,
       evo_araP_ribo_rep1_riboSeq,
       evo_araP_ribo_rep2_rnaSeq,
       evo_araP_ribo_rep2_riboSeq))


mean(c(anc_araM_lip_rep1_rnaSeq,
       anc_araM_lip_rep1_riboSeq,
       anc_araM_lip_rep2_rnaSeq,
       anc_araM_lip_rep2_riboSeq))
mean(c(evo_araM_lip_rep1_rnaSeq,
       evo_araM_lip_rep1_riboSeq,
       evo_araM_lip_rep2_rnaSeq,
       evo_araM_lip_rep2_riboSeq))

mean(c(anc_araM_ribo_rep1_rnaSeq,
       anc_araM_ribo_rep1_riboSeq,
       anc_araM_ribo_rep2_rnaSeq,
       anc_araM_ribo_rep2_riboSeq))
mean(c(evo_araM_ribo_rep1_rnaSeq,
       evo_araM_ribo_rep1_riboSeq,
       evo_araM_ribo_rep2_rnaSeq,
       evo_araM_ribo_rep2_riboSeq))

# Pool - anc vs. evo
mean(c(anc_araP_lip_rep1_rnaSeq,
       anc_araP_lip_rep1_riboSeq,
       anc_araP_lip_rep2_rnaSeq,
       anc_araP_lip_rep2_riboSeq,
       anc_araM_lip_rep1_rnaSeq,
       anc_araM_lip_rep1_riboSeq,
       anc_araM_lip_rep2_rnaSeq,
       anc_araM_lip_rep2_riboSeq))
mean(c(evo_araP_lip_rep1_rnaSeq,
       evo_araP_lip_rep1_riboSeq,
       evo_araP_lip_rep2_rnaSeq,
       evo_araP_lip_rep2_riboSeq,
       evo_araM_lip_rep1_rnaSeq,
       evo_araM_lip_rep1_riboSeq,
       evo_araM_lip_rep2_rnaSeq,
       evo_araM_lip_rep2_riboSeq))

mean(c(anc_araP_ribo_rep1_rnaSeq,
       anc_araP_ribo_rep1_riboSeq,
       anc_araP_ribo_rep2_rnaSeq,
       anc_araP_ribo_rep2_riboSeq,
       anc_araM_ribo_rep1_rnaSeq,
       anc_araM_ribo_rep1_riboSeq,
       anc_araM_ribo_rep2_rnaSeq,
       anc_araM_ribo_rep2_riboSeq))
mean(c(evo_araP_ribo_rep1_rnaSeq,
       evo_araP_ribo_rep1_riboSeq,
       evo_araP_ribo_rep2_rnaSeq,
       evo_araP_ribo_rep2_riboSeq,
       evo_araM_ribo_rep1_rnaSeq,
       evo_araM_ribo_rep1_riboSeq,
       evo_araM_ribo_rep2_rnaSeq,
       evo_araM_ribo_rep2_riboSeq))

# Pool - anc vs. evo
anc_lip <- c(anc_araP_lip_rep1_rnaSeq,
                  anc_araP_lip_rep1_riboSeq,
                  anc_araP_lip_rep2_rnaSeq,
                  anc_araP_lip_rep2_riboSeq,
                  anc_araM_lip_rep1_rnaSeq,
                  anc_araM_lip_rep1_riboSeq,
                  anc_araM_lip_rep2_rnaSeq,
                  anc_araM_lip_rep2_riboSeq)
evo_lip <- c(evo_araP_lip_rep1_rnaSeq,
                  evo_araP_lip_rep1_riboSeq,
                  evo_araP_lip_rep2_rnaSeq,
                  evo_araP_lip_rep2_riboSeq,
                  evo_araM_lip_rep1_rnaSeq,
                  evo_araM_lip_rep1_riboSeq,
                  evo_araM_lip_rep2_rnaSeq,
                  evo_araM_lip_rep2_riboSeq)

anc_ribo <- c(anc_araP_ribo_rep1_rnaSeq,
                   anc_araP_ribo_rep1_riboSeq,
                   anc_araP_ribo_rep2_rnaSeq,
                   anc_araP_ribo_rep2_riboSeq,
                   anc_araM_ribo_rep1_rnaSeq,
                   anc_araM_ribo_rep1_riboSeq,
                   anc_araM_ribo_rep2_rnaSeq,
                   anc_araM_ribo_rep2_riboSeq)
evo_ribo <- c(evo_araP_ribo_rep1_rnaSeq,
                   evo_araP_ribo_rep1_riboSeq,
                   evo_araP_ribo_rep2_rnaSeq,
                   evo_araP_ribo_rep2_riboSeq,
                   evo_araM_ribo_rep1_rnaSeq,
                   evo_araM_ribo_rep1_riboSeq,
                   evo_araM_ribo_rep2_rnaSeq,
                   evo_araM_ribo_rep2_riboSeq)

mean(anc_lip)
mean(evo_lip)

mean(anc_ribo)
mean(evo_ribo)

hist(anc_lip)
hist(evo_lip)

hist(anc_ribo)
hist(evo_ribo)

mean(anc_lip)
sd(anc_lip)/sqrt(length(anc_lip))

mean(evo_lip)
sd(evo_lip)/sqrt(length(evo_lip))

mean(anc_ribo)
sd(anc_ribo)/sqrt(length(anc_ribo))

mean(evo_ribo)
sd(evo_ribo)/sqrt(length(evo_ribo))
