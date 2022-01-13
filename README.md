# growth_scaling_envelope


# ********** Directory descriptions **********
# ./data/cell_traits_data/: Contains the data on cell dimensions and growth rates. Other data (e.g., physiological and genomic properties) also included.

# ./data/cell_size_data/: Contains cell size and growth measurements in Lenski's LTEE. Used in Figure 8.

# ./data/main/: Location for all processed data used by master.R to perform analysis and plot figures.
# ********************************************


# ********** Spreadsheets descriptions ***********
# ./data/main/proteomic_fractions_full.csv: Escherichia coli proteomic mass fractions of ribosomes and envelope-producers under different treatments.
# ./data/main/quant_proteome_data.csv: Proteomic mass fractions of ribosomes and envelope-producers under different treatments across species.
# ./data/main/growth_scaling.shape.genome.csv: Data on surface area, internl volume, growth rate, envelope thickness, and genomic features of bacteria.
# ./data/main/growth_scaling.shape.genome.envelope_corrected.csv: Same as above, but internal volume calculated using species-specific envelope thickness.
# ./data/main/proteome_scaling.csv: Contains data on proteomic mass fractions, surface area, and internal volume across bacteria.

# ./data/ltee_data/gallet2017.csv: Growth rate and cell volume in "The evolution of bacterial cell size: the internal diffusion-constraint hypothesis."
# ./data/ltee_data/gallet2017.csv: Fitness and linear dimensions of LTEE lines from "Changes in Cell Size and Shape during 50,000 Generations of Experimental Evolution with Escherichia coli."
# ./data/ltee_data/lennon2021.csv: Fitness and width of Mycoplasma mycoides from "Evolution of a minimal cell."

# ./data/cell_traits_data/bacterial_envelopes_Sep232021.xlsx: Compiled data from literature on cell size, growth, genomic features, envelope thickness, etc.
# ************************************************


# ********** Script descriptions **********
# ./data/main/cellular_traits_processing.R: Integrates different cell traits and proteomic data into a common file used for downstream analysis.
#         Input  -- ./cell_traits_data/bacterial_envelopes_Sep232021 - {growth, size, other properties}.csv AND quant_proteome_data.csv
#         Output -- growth_scaling.shape.genome.csv

# ./data/main/cellular_traits_processing_envelope_thickness.R: The same as the previous script, but uses species-specific estimates of envelope thickness.
#         Input  -- ./cell_traits_data/bacterial_envelopes_Sep232021 - {growth, size, other properties}.csv
#         Output -- growth_scaling.shape.envelope_corrected.csv

# ./data/main/master.R: Estimates all parameters, performs regression analyses, and plots the figures in the manuscript.
# *****************************************


# ********** Notebook descriptions ***********
# integrated_growth_theory.nb: Derivation of all theoretical results, numerical and analytical.
# cellEngine.nb: Numerically integrates system of non-linear ODEs.
# cellEngine_calling.nb: Executes the cellEngine.nb and checks the analytics against numerics. Used for Figure 2.
# ********************************************


