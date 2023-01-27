# growth_scaling_envelope
# ************************************************
# UPDATE (July 30, 2022):
# 1) Files affected -- "cellular_traits_processing.R" and "cellular_traits_processing_envelope_thickness.R"
#    Bug -- When calling functions helix_volume() and helix_surface(), helix diameter and helix length were swapped when passing arguments. Lines 
#           Lines where bug is manifested are 105 to 107, and 117 to 119 in "cellular_traits_processing.R", 
#           and 82 to 84 and 94 to 96 in "cellular_traits_processing_envelope_thickness.R"
#    Status -- Fixed by swapping arguments in function definition. Statistics updated but none of our conclusions are affected.
# 2) Files affected -- "cellular_traits_processing.R"
#    Bug -- Mycoplasma pneumoniae has been accidentally omitted during data processing. We had issues assigning genes via KEGG database, 
#           which is why we excluded it from the analysis. Later we realized that the species had estimate for ribosome number from a separate study
#            that we automatically removed during data filtering.
#    Status -- "Processing growth rate data" section was re-written to an include estimate on ribosome number from Lynch et al. 2017.
#              Statistics are updated and non of the conclusions are changed.
# ************************************************
# 
# UPDATE (January 26, 2023):
# 1) Added a section of the scaling of metabolic rate in LTEE in "integrated_growth_theory.nb"
#           Plotting graphs via "master.R" requires "Respiration_scaling.csv" from Marshall et al. 2022 at
#           https://datadryad.org/stash/dataset/doi:10.5061/dryad.kd51c5b7n
# 2) Added "favate_processing.R" -- Script which computed \Phi_R and \Phi_L based on transcriptomic data from two time points in LTEE
# 3) Replaced "grant2021.csv" with "grant_wiser_full.csv" due to a mistake in computing S/V across LTEE:
#    Bug -- An error occurred while processing data from Grant et al. 2021 (i.e., Lenski's cell size measurements), and the wrong column was used as cell #           dimensions. That is, I used the length-to-width ratio as length, and length as width. This is now corrected, and the results are not
#           significantly impacted.
# 4) Modified "lennon2021.csv"
#    Bug -- Cell diameter was initially read off from Figure 4. We now took values directly reported in Table S2. Small differences
#           exist, but nothing that challenges conclusions.
# 5) Modified "master.R"
#    Bug -- Eq S.17 was erroneously transcribed from authors notes, causing misestimation of \kappa_l. This is now fixed, with the improvement of the
#           actual fit. The whole section on accounting for variation in S/V across growth conditions is elaborated on in the manuscript.
# ************************************************


# ********** Directory descriptions **********
# ./data/cell_traits_data/: Contains the data on cell dimensions and growth rates. Other data (e.g., physiological and genomic properties) also included.

# ./data/cell_size_data/: Contains cell size and growth measurements in Lenski's LTEE. Used in Figure 8.

# ./data/main/: Location for all processed data used by master.R to perform analysis and plot figures.
# ********************************************


# ********** Spreadsheets descriptions ***********
# ./data/main/proteomic_fractions_full.csv: Escherichia coli proteomic mass fractions of ribosomes and envelope-producers under different treatments.
# ./data/main/quant_proteome_data.csv: Proteomic mass fractions of ribosomes and envelope-producers under different treatments across species.
# ./data/main/growth_scaling.shape.genome.csv: Data on surface area, internal volume, growth rate, envelope thickness, and genomic features of bacteria.
# ./data/main/growth_scaling.shape.genome.envelope_corrected.csv: Same as above, but internal volume calculated using species-specific envelope thickness.
# ./data/main/proteome_scaling.csv: Contains data on proteomic mass fractions, surface area, and internal volume across bacteria.

# ./data/ltee_data/gallet2017.csv: Growth rate and cell volume in "The evolution of bacterial cell size: the internal diffusion-constraint hypothesis."
# ./data/ltee_data/grant_wiser_full.csv: Fitness and linear dimensions of LTEE lines from "Changes in Cell Size and Shape during 50,000 Generations of Experimental Evolution with Escherichia coli."
# ./data/ltee_data/lennon2021.csv: Fitness and width of Mycoplasma mycoides from "Evolution of a minimal cell."

# ./data/cell_traits_data/bacterial_envelopes_Sep232021.xlsx: Compiled data from literature on cell size, growth, genomic features, envelope thickness, etc.
# ************************************************


# ********** Script descriptions **********
#
# NOTE: Source files needed to generate "proteomic_fractions_full.csv" are omitted due to potential Copyright issues. Either download them from
#       the designated URLs or contact me for assistance. 
#
# ./data/main/quant_proteome_annotation.R
#         Description: Assigns each gene as a ribosomal protein or an envelope-producer
#                      for E. coli quantiative proteomics studies across growth conditions. It uses previosuly published studies, and outputs
#                      a modified spreadsheet for each study with three additional columns, denoting KEGG id, and the assignment to a proteome sector
#         Input:    Supplementary Table 2 (Erickson et al. 2017)[https://www.nature.com/articles/nature24299],
#                   Table S6 (Schmidt et al. 2016)[https://www.nature.com/articles/nbt.3418],
#                   Supplementary Information (Peebo et al. 2015)[https://pubs.rsc.org/en/content/articlelanding/2015/MB/C4MB00721B],
#                   Table S1 (Li et al. 2014)[https://www.sciencedirect.com/science/article/pii/S0092867414002323?via%3Dihub],
#                   Table S2 (Valgepea et al. 2013)[https://pubs.rsc.org/en/content/articlelanding/2013/MB/c3mb70119k],
#                   Dataset EV2 and Dataset EV9 (Mori et al. 2021)[https://www.embopress.org/doi/full/10.15252/msb.20209536]                
#         Output:   ./assigned_function/*_assigned.csv
#
#
# ./data/main/quant_proteome_analysis_pooled.R
#         Description: Calculates mass fractions for each proteomic class when E. coli is reared in different conditions. It outputs a joint spreadsheet
#                      with values of \Phi_R and \Phi_L for all studies used to generate Figure 4 in the manuscript.
#         Input: ./assigned_function/*_assigned.csv
#         Output: proteomic_fractions_full.csv
#
#
# ./data/main/paxDB_proteomes.R
#         Description: Assigns proteins into proteome sectors and calculates proteomic mass fractions for all species in PaxDB.
#                      Note that it requires original mapping tsv, and tsv for each species.
#         Input: ./bacterial_proteomes_PaxDb/
#         Output: paxdb_proteomes.csv
#
#
# ./data/main/muller_proteomes.R
#         Description: Processes bacterial proteomes from Muller et al. 2020 [https://www.nature.com/articles/s41586-020-2402-x]
#         Input: Supplementary Table 1 ("./table1_species_list_Muller.csv") and Supplementary Table 2 ("./table2_intensities_data_Muller.csv")
#         Output: muller2020_processed.csv
#
#
# ./data/main/added_proteomes.R
#         Description: Assigns and calculates mass fractions for additional proteomic studies not included in previouslt processed.
#         Input: Table S2 (Masson et al. 2021)[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0250524]
#                Supplementary Table S1 (Angel et al. 2010)[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0013800]
#                Supplementary Table S3 (Srivastava et al. 2020)[https://www.frontiersin.org/articles/10.3389/fmicb.2020.544785/full]
#                Dataset EV6 and EV7 (Matteau et al. 2020)[https://www.embopress.org/doi/full/10.15252/msb.20209844]
#                Table S3 (Osbak et al. 2016)[https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004988]
#         Output: *.processed.csv
#
#
# ./data/main/merge_proteomes.R
#         Description: Merges all of the above-processed proteomic data. It is used for plotting Figure 5 in the manuscript. 
#         Input: ./data/main/(with all associated directories and files)
#         Output: quant_proteome_data.csv
#
#
# ./data/main/cellular_traits_processing.R
#         Description: Integrates different cell traits and proteomic data into a common file used for downstream analysis.
#         Input  -- ./cell_traits_data/bacterial_envelopes_Sep232021 - {growth, size, other properties}.csv AND quant_proteome_data.csv
#         Output -- growth_scaling.shape.genome.csv
#
#
# ./data/main/cellular_traits_processing_envelope_thickness.R
#         Description: The same as the previous script, but uses species-specific estimates of envelope thickness.
#         Input  -- ./cell_traits_data/bacterial_envelopes_Sep232021 - {growth, size, other properties}.csv
#         Output -- growth_scaling.shape.envelope_corrected.csv
#
#
# ./data/main/favate_processing.R
#         Description: Computes \Phi_L and \Phi_R for ancestral and evolved populations
#         Input -- LTEE mRNA abundances reported in [Supplemental Table 1]([url](https://www.biorxiv.org/content/10.1101/2021.01.12.426406v1.supplementary-material)) and data from Valgepea et al. 2013 is used to convert mRNA to protein abundances
#         Outpur -- Printed proteome mass fractions
#
#
# ./data/main/master.R
#         Description: Estimates all parameters, performs regression analyses, and plots the figures in the manuscript.
# *****************************************


# ********** Notebook descriptions ***********
# integrated_growth_theory.nb: Derivation of all theoretical results, numerical and analytical.
# cellEngine.nb: Numerically integrates system of non-linear ODEs.
# cellEngine_calling.nb: Executes the cellEngine.nb and checks the analytics against numerics. Used for Figure 2.
# ********************************************


