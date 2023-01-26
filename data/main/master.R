#=============== USAGE =====================
# Unites both data (from other R scripts) and theory (from Mathematica notebook)
# Output: All reported figures in the text

remove(list=ls())

library(dplyr)
library(ggplot2)
library(viridis)
library(scales)
library(grid)
library(gridExtra)
library(car)

setwd("/home/bogi/Desktop/growth_rate/data/main/")
# Set the color-blind scale
colBlindScale <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                   "#0072B2", "#D55E00", "#CC79A7", "#DB6D00", "#EF3B2C",
                   "#67001F","#CB181D","#78C679","#F46D43","#A6CEE3","#FD8D3C","#A6D854","#D4B9DA","#6A51A3",
                   "#7F0000","#D9D9D9","#FFF7BC","#F0F0F0","#C7EAE5","#003C30","#F16913","#FFF7FB","#8C6BB1","#C7E9B4",
                   "#762A83","#FC9272","#AE017E","#F7F7F7","#DF65B0","#74C476")
shape_scale <- c(0,16,2,5,6,7,15,1,17,18,19)
scaleFUN <- function(x) sprintf("%.2f", x)

medium_text <- c(bquote(kappa[n]*"(M63+gly)"), bquote(kappa[n]*"(M63+glc)"), bquote(kappa[n]*"(cAA+gly)"),
                 bquote(kappa[n]*"(cAA+glc)"), bquote(kappa[n]*"(RDM+gly)"), bquote(kappa[n]*"(RDM+glc)"))
plot_points <- 5000
line_thick <- 1
axes_style <- element_text(size=15, face="bold", colour = "black")
L_env <- 0.03
rho_Dai <- 0.48

#========== Growth rate and proteomic fractions equations ==========
# The growth rate under cell envelope burden
lambda_env <- function(kn, kt, kl, eps, svPar, dp, dl) {
  thetaE <- (
    eps*svPar*(kl+kn)*(kt+dp)*(kn*dl+(dl-dp+kn)*kt)
  )/(
    (kn+kt)*(kl*kt*(kn-dp)-dl*(kn+kl)*(kt+dp)*eps*svPar)
  )
  
  return(
    (kt*(kn-dp))/((kn + kt)*(1 + thetaE))
  )
}

lambda_anaero_correct <- function(kn, kt, kl, eps, svPar, dp, dl) {
  anaero_coeff <- 0.63
  a <- (kl*(dp-kn)*kt+dl*eps*(kn+kl)*(dp+kt)*svPar)*(kl*(anaero_coeff*kn+kt)+eps*(kl+anaero_coeff*kn)*(dp+kt)*svPar)
  b <- (kl*(kn+kt)+eps*(kl+kn)*(dp+kt)*svPar)*(kl*(dp-anaero_coeff*kn)*kt+dl*eps*(anaero_coeff*kn+kl)*(dp+kt)*svPar)
  return(a/b)
}

riboFrac_env <- function(kn, kt, kl, dp, dl, eps, sv) {
  phiR_opt <- (kl*(-dp + kn)*kt - dl*eps*(kl + kn)*(dp + kt)*sv)/(
    kl*kt*(kn + kt) + eps*(kl + kn)*kt*(dp + kt)*sv)
  return(
    ((dp+kt*phiR_opt))/(dp+kt)
  )
}

lipoFrac_env <- function(kn, kt, kl, dp, dl, eps, sv) {
  phiL_opt <- (eps*(dp + kt)*(dl*kn + (dl - dp + kn)*kt)*sv)/(
    kl*kt*(kn + kt) + eps*(kl + kn)*kt*(dp + kt)*sv) 
  return(
    kt*phiL_opt/(dp+kt)
  )
}

q10_corrected <- function(rate, q10_coeff, temp) {
  return(rate*q10_coeff^((20-temp)/10))
}

# d1 -- diameter; d2 -- tip-to-tip length
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


#========== Estimating slopes of E. coli proteomic responses ==========
# Load E. coli proteomes across various growth conditions
df_eco <- read.csv("proteomic_fractions_full.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# First, we focus only on scaling of proteomic composition across different media quality
# Studies by Valgepea2013 and Peebo2015 are excluded because the slopes are larger than in other papers; These are reported separately
# We also exclude two data points measured in stationary phase, given that our theory assumes expo-growing cell
# All parameters are inferred from Schmidt2016 due to superior quality relative
# Other data is plotted to make sure that results do not look widely off across studies
df_eco_nutri <- subset(df_eco, perturbation=="nutrient conditions" & growth_rates>0)
df_eco_pl_chloram <- subset(df_eco, perturbation=="chloramphenicol" & source=="Mori2021")
df_eco_nutri_all <- subset(df_eco, perturbation=="nutrient conditions" & growth_rates>0)
df_eco_nutri_Schmidt <- subset(df_eco, perturbation=="nutrient conditions" & source=="Schmidt2016" & growth_rates>0)

# olsR is used to later compute kappa_t
fit_phiR_nutri <- lm(phi_R ~ growth_rates, data=df_eco_nutri_Schmidt)
olsR <- summary(fit_phiR_nutri)

# Second, we focus on scaling of envelope-producer across different media quality
# olsL is the slope of regression, and will be used to compute kappa_l
fit_phiL_nutri <- lm(phi_L ~ growth_rates, data=df_eco_nutri_Schmidt)
olsL <- summary(fit_phiL_nutri)

# Thidly, we now focus on scaling of proteome composition with growth rate, when latter is modulated by translation inhibitor
# olsR_ab are slopes in various media, and will be used to compute kappa_n
df_eco_ab_ribo <- subset(df_eco, perturbation=="chloramphenicol" & (source=="Scott2010"))
fit_phiR_ab1 <- lm(phi_R ~ growth_rates, data=df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="M63 + Glyc",])
fit_phiR_ab2 <- lm(phi_R ~ growth_rates, data=df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="M63 + Glc",])
fit_phiR_ab3 <- lm(phi_R ~ growth_rates, data=df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="cAA + Glyc",])
fit_phiR_ab4 <- lm(phi_R ~ growth_rates, data=df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="cAA + Glc",])
fit_phiR_ab5 <- lm(phi_R ~ growth_rates, data=df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="RDM + Glyc",])
fit_phiR_ab6 <- lm(phi_R ~ growth_rates, data=df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="RDM + Glc",])
olsR_ab1 <- summary(fit_phiR_ab1)
olsR_ab2 <- summary(fit_phiR_ab2)
olsR_ab3 <- summary(fit_phiR_ab3)
olsR_ab4 <- summary(fit_phiR_ab4)
olsR_ab5 <- summary(fit_phiR_ab5)
olsR_ab6 <- summary(fit_phiR_ab6)

# Fourtly, we examine how ribosome fraction scales with growth rate, when latter is modulated by other antibiotics
# All data should roughly fall onto the same regression
df_eco_env_ribo <- subset(df_eco, (perturbation=="triclosan" | perturbation=="fosfomycin" |
                                     perturbation=="nutrient conditions") & source!="Valgepea2013" & source!="Peebo2015" & growth_rates>0)
df_eco_env_ribo$perturbation <- factor(df_eco_env_ribo$perturbation, levels = c('rifampicin', 'triclosan', 'fosfomycin', 'nutrient conditions'))


#========== Accounting for the variation in S:V across growth conditions ==========
# Load data on \Phi_R, S:V, and the growth rate from four studies
df_si <- read.csv("/home/bogi/Desktop/growth_rate/data/eco_proteomes/si2017_data.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_volkmer <- read.csv("/home/bogi/Desktop/growth_rate/data/eco_proteomes/volkmer2011_size.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_vadia <- read.csv("/home/bogi/Desktop/growth_rate/data/eco_proteomes/vadia2017_size.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_basan <- read.csv("/home/bogi/Desktop/growth_rate/data/eco_proteomes/basan2015_data.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Split Basan2015 into with and without a.b. (Chloramphenicol) treatment
# Calculate \Pi by first subtracting 2x E. coli cell envelope thickness from each linear dimension
basan_flag <- grepl("Cm",df_basan$medium)
df_basan$sv <- capsule_surf(df_basan$cell.width-2*L_env,df_basan$cell.length-2*L_env)/capsule_vol(df_basan$cell.width-2*L_env,df_basan$cell.length-2*L_env)
df_basan_Cm <- df_basan[basan_flag,]
df_basan <- df_basan[!basan_flag,]
df_vadia$sv <- capsule_surf(df_vadia$cell.width-2*L_env,df_vadia$cell.length-2*L_env)/capsule_vol(df_vadia$cell.width-2*L_env,df_vadia$cell.length-2*L_env)
df_volkmer$sv <- capsule_surf(df_volkmer$cell_width.mean-2*L_env,df_volkmer$cell_length.mean-2*L_env)/capsule_vol(df_volkmer$cell_width.mean-2*L_env,df_volkmer$cell_length.mean-2*L_env)
df_si$surface.to.volume.ratio..μm.1. <- capsule_surf(df_si$cel.width..μm.-2*L_env,df_si$cel.length..μm.-2*L_env)/
  capsule_vol(df_si$cel.width..μm.-2*L_env,df_si$cel.length..μm.-2*L_env)

# ***** 1. Process Si2017 *****
# Convert RNA/Protein ratio into \Phi_R
df_si$RNA.protein <- rho_Dai*df_si$RNA.protein

# Take only chloramphenicol-treated experiments, and average across replicates for each Cm concentration, for each medium
df_si <- df_si[df_si$type.of.perturbation=='chloramphenicol',]
df_si_processed <- aggregate(df_si, by=list(df_si$concentration,df_si$growth.media), mean, na.rm = TRUE)
media_list <- unique(df_si_processed$Group.2)
df_si_nutrients <- df_si_processed[df_si_processed$Group.1=="0 µM",]      # 0 uM is the experiment with carbon source-modulated growth rate

# \Pi~\lambda across carbon sources
lm_si_nutrients <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_nutrients)
summary(lm_si_nutrients)

lm_basan_nutrients <- lm(sv~growth.rate, data=df_basan)
summary(lm_basan_nutrients)

lm_volkmer_nutrients <- lm(sv~growth_rate.mean, data=df_volkmer)
summary(lm_volkmer_nutrients)

p_sv_growth_nut <- ggplot() + 
  geom_point(size=4.5, shape=15, 
             aes(x=df_si_nutrients$growth.rate..1.hours., 
                 y=df_si_nutrients$surface.to.volume.ratio..μm.1.,
                 col="a")) +
  geom_point(size=4.5, shape=16, 
             aes(x=df_basan$growth.rate, 
                 y=df_basan$sv,
                 col="b")) +
  geom_point(size=4.5, shape=17, 
             aes(x=df_volkmer$growth_rate.mean, 
                 y=df_volkmer$sv,
                 col="c")) +
  stat_function(fun=function(x) lm_volkmer_nutrients$coefficients[1]+lm_volkmer_nutrients$coefficients[2]*x,
                aes(col="c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) lm_basan_nutrients$coefficients[1]+lm_basan_nutrients$coefficients[2]*x,
                aes(col="b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) lm_si_nutrients$coefficients[1]+lm_si_nutrients$coefficients[2]*x,
                aes(col="a"), size=line_thick, n=plot_points) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual("Source", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3]), 
                      labels = c("Volkemer2011","Basan2015","Si2017")) +
  scale_shape_manual(name="Source", values = c(15,16,17), 
                     labels = c("Volkemer2011","Basan2015","Si2017")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('S/V, '*Pi*' ('*mu*''*m^-1*')]'))) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17),linetype=0)))
p_sv_growth_nut


# \Pi~\lambda across Cm concentration, for each carbon-source
media_list
df_si_1 <- df_si_processed[df_si_processed$Group.2==media_list[1],]
df_si_2 <- df_si_processed[df_si_processed$Group.2==media_list[2],]
df_si_3 <- df_si_processed[df_si_processed$Group.2==media_list[3],]
df_si_4 <- df_si_processed[df_si_processed$Group.2==media_list[4],]
df_si_5 <- df_si_processed[df_si_processed$Group.2==media_list[5],]
df_si_6 <- df_si_processed[df_si_processed$Group.2==media_list[6],]
df_si_7 <- df_si_processed[df_si_processed$Group.2==media_list[7],]
df_si_8 <- df_si_processed[df_si_processed$Group.2==media_list[8],]
df_si_9 <- df_si_processed[df_si_processed$Group.2==media_list[9],]

lm_si_1 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_1)
lm_si_2 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_2)
lm_si_3 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_3)
lm_si_4 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_4)
lm_si_5 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_5)
lm_si_6 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_6)
lm_si_7 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_7)
lm_si_8 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_8)
lm_si_9 <- lm(surface.to.volume.ratio..μm.1.~growth.rate..1.hours., data=df_si_9)
summary(lm_si_1)
summary(lm_si_2)
summary(lm_si_3) # ***
summary(lm_si_4) 
summary(lm_si_5) 
summary(lm_si_6) # *
summary(lm_si_7)
summary(lm_si_8)
summary(lm_si_9)

# With the exception of two media, \Pi is constant across growth conditions
mean(df_si_1$surface.to.volume.ratio..μm.1.)
mean(df_si_2$surface.to.volume.ratio..μm.1.)
mean(df_si_4$surface.to.volume.ratio..μm.1.)
mean(df_si_5$surface.to.volume.ratio..μm.1.)
mean(df_si_7$surface.to.volume.ratio..μm.1.)
mean(df_si_8$surface.to.volume.ratio..μm.1.)
mean(df_si_9$surface.to.volume.ratio..μm.1.)

p_sv_growth_ab <- ggplot() + 
  geom_point(size=4.5, 
             aes(x=df_si_1$growth.rate..1.hours., 
                 y=df_si_1$surface.to.volume.ratio..μm.1.,
                 col="a"), shape=15, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_2$growth.rate..1.hours., 
                 y=df_si_2$surface.to.volume.ratio..μm.1.,
                 col="b"), shape=16, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_3$growth.rate..1.hours., 
                 y=df_si_3$surface.to.volume.ratio..μm.1.,
                 col="c"), shape=17, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_4$growth.rate..1.hours., 
                 y=df_si_4$surface.to.volume.ratio..μm.1.,
                 col="d"), shape=18, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_5$growth.rate..1.hours., 
                 y=df_si_5$surface.to.volume.ratio..μm.1.,
                 col="e"), shape=19, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_6$growth.rate..1.hours., 
                 y=df_si_6$surface.to.volume.ratio..μm.1.,
                 col="f"), shape=0, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_7$growth.rate..1.hours., 
                 y=df_si_7$surface.to.volume.ratio..μm.1.,
                 col="g"), shape=1, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_8$growth.rate..1.hours., 
                 y=df_si_8$surface.to.volume.ratio..μm.1.,
                 col="h"), shape=2, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_9$growth.rate..1.hours., 
                 y=df_si_9$surface.to.volume.ratio..μm.1.,
                 col="i"), shape=5, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_basan_Cm$growth.rate, 
                 y=df_basan_Cm$sv,
                 col="j"), shape=6, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_vadia$growth.rate, 
                 y=df_vadia$sv,
                 col="k"), shape=7, stroke=1.5) +
  stat_function(fun=function(x) lm_si_3$coefficients[1]+lm_si_3$coefficients[2]*x,
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) lm_si_6$coefficients[1]+lm_si_6$coefficients[2]*x,
                aes(colour = "f"), size=line_thick, n=plot_points) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual("Medium", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                           "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6],
                                           "g"=colBlindScale[7], "h"=colBlindScale[8], "i"=colBlindScale[9],
                                           "j"=colBlindScale[10], "k"=colBlindScale[11]), 
                      labels = c("Glucose (Si2017)","Glucose + 12 a.a. (Si2017)","Glucose + 6 a.a. (Si2017)",
                                 "Glucose + caa (Si2017)","Glucose synth. rich (Si2017)","Glycerol (Si2017)",
                                 "Mannose (Si2017)","Sorbitol (Si2017)","TSB (Si2017)",
                                 "Glucose (Basan2015)","LB + glucose (Vadia2017)")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('S/V, '*Pi*' ('*mu*''*m^-1*')'))) +
  ylim(c(3.5,9)) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0,1,2,5,6,7),linetype=0)))
p_sv_growth_ab


#========== Inferring capacities while correcting for S:V across growth conditions ==========
# Base parameters that are collected from the literature
SAVpar <- 5 # Rough average in our dataset; More precisely, it would be 4.94 but this might change if more data is added
epsPar <- (0.30/0.55)/SAVpar
dp_est <- 0.05 # page 2896 in Nath and Koch 1970 (taken for slowly decaying component)
dl_est <- 1.00

# Phi_L and growth rates from Schmidt2016
y <- df_eco_nutri_Schmidt$phi_L
x <- df_eco_nutri_Schmidt$growth_rates

# The obtained \Pi~\lambda_N is plugged into eqn. 27, in order to account for changes in cell shape when inferring \kappa_l
# We compare results across three different datasets that perturbed \Pi by changing nutrient quality
# Non-linear regression to obtain kappaL estimate, then check for goodness-of-fit
m_Si <- nls(y~y_int+(lm_si_nutrients$coefficients[1]+lm_si_nutrients$coefficients[2]*x)*x*epsPar/kappaL)
m_Basan <- nls(y~y_int+(lm_basan_nutrients$coefficients[1]+lm_basan_nutrients$coefficients[2]*x)*x*epsPar/kappaL)
m_Volkmer <- nls(y~y_int+(lm_volkmer_nutrients$coefficients[1]+lm_volkmer_nutrients$coefficients[2]*x)*x*epsPar/kappaL)

m_summary_Si <- summary(m_Si)
m_summary_Si
kappaL_est_Si <- m_summary_Si$coefficients[2]
kappaL_int_Si <- m_summary_Si$coefficients[1]

m_summary_Basan <- summary(m_Basan)
m_summary_Basan
kappaL_est_Basan <- m_summary_Basan$coefficients[2]
kappaL_int_Basan <- m_summary_Basan$coefficients[1]

m_summary_Volkmer <- summary(m_Volkmer)
m_summary_Volkmer
kappaL_est_Volkmer <- m_summary_Volkmer$coefficients[2]
kappaL_int_Volkmer <- m_summary_Volkmer$coefficients[1]

cor(y,predict(m_Si))
cor(y,predict(m_Volkmer))
cor(y,predict(m_Basan))

# The line is a best fit, given the changes in Pi across growth conditions
df_fin_Schmidt <- df_eco_nutri[df_eco_nutri$source=="Schmidt2016",]
df_fin_Peebo <- df_eco_nutri[df_eco_nutri$source=="Peebo2015",]
df_fin_Valgepea <- df_eco_nutri[df_eco_nutri$source=="Valgepea2013",]
df_fin_Erickson <- df_eco_nutri[df_eco_nutri$source=="Erickson2017",]
df_fin_Li <- df_eco_nutri[df_eco_nutri$source=="Li2014",]
df_fin_Si <- df_eco_nutri[df_eco_nutri$source=="Si2017",]
df_fin_Dai <- df_eco_nutri[df_eco_nutri$source=="Dai2018",]
df_fin_Bremer <- df_eco_nutri[df_eco_nutri$source=="Bremer1996",]
df_fin_Scott <- df_eco_nutri[df_eco_nutri$source=="Scott2010",]
df_fin_Forchhammer <- df_eco_nutri[df_eco_nutri$source=="Forchhammer1971",]
df_fin_Mori <- df_eco_nutri[df_eco_nutri$source=="Mori2021",]

p_lipo_nut_pert <- ggplot() + 
  stat_function(fun=function(x) kappaL_int_Volkmer+(lm_volkmer_nutrients$coefficients[1]*epsPar/kappaL_est_Volkmer)*x+(lm_volkmer_nutrients$coefficients[2]*epsPar/kappaL_est_Volkmer)*x^2, aes(colour = "g"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_fin_Valgepea$growth_rates, y=df_fin_Valgepea$phi_L, colour = "d")) +
  geom_point(size=4.5, stroke=1.5, shape=19, aes(x=df_fin_Li$growth_rates, y=df_fin_Li$phi_L, colour = "e")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_fin_Peebo$growth_rates, y=df_fin_Peebo$phi_L, colour = "f")) +
  geom_point(size=4.5, stroke=1.5, shape=1, aes(x=df_fin_Schmidt$growth_rates, y=df_fin_Schmidt$phi_L, colour = "g")) +
  geom_point(size=4.5, stroke=1.5, shape=2, aes(x=df_fin_Erickson$growth_rates, y=df_fin_Erickson$phi_L, colour = "h")) +
  geom_point(size=4.5, stroke=1.5, shape=7, aes(x=df_fin_Mori$growth_rates, y=df_fin_Mori$phi_L, colour = "k")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name="Source", values = c("d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6],
                                                "g"=colBlindScale[7], "h"=colBlindScale[8], "k"=colBlindScale[11]), 
                      labels = c("Forchhammer1971","Bremer1996","Scott2010",
                                 "Valgepea2013","Li2014","Peebo2015",
                                 "Schmidt2016","Erickson2017","Si2017","Dai2018","Mori2021")) +
  scale_shape_manual(name="Medium", values = c(18,19,0,1,2,7), 
                     labels = c("Valgepea2013","Li2014","Peebo2015",
                                "Schmidt2016","Erickson2017","Mori2021")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Envelope-producer mass fraction, '*Phi[L]))) +
  guides(colour=guide_legend(override.aes=list(shape=c(18,19,0,1,2,7),linetype=0,stroke=1.5)))
p_lipo_nut_pert


p_lipo_lambda_alt <- ggplot() + 
  stat_function(fun=function(x) kappaL_int_Si+(lm_si_nutrients$coefficients[1]+lm_si_nutrients$coefficients[2]*x)*(epsPar*x/kappaL_est_Si), aes(colour = "d"), size = 1.5) +
  stat_function(fun=function(x) kappaL_int_Basan+(lm_basan_nutrients$coefficients[1]+lm_basan_nutrients$coefficients[2]*x)*(epsPar*x/kappaL_est_Basan), aes(colour = "e"), size = 1.5) +
  stat_function(fun=function(x) kappaL_int_Volkmer+(lm_volkmer_nutrients$coefficients[1]+lm_volkmer_nutrients$coefficients[2]*x)*(epsPar*x/kappaL_est_Volkmer), aes(colour = "g"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=1, aes(x=df_fin_Schmidt$growth_rates, y=df_fin_Schmidt$phi_L, colour = "g")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name=bquote("Source of "*Pi*"("*lambda[N]*")"), values = c(colBlindScale[4], colBlindScale[5], colBlindScale[7], colBlindScale[7]), 
                      labels = c("Si2017","Basan2015","Volkmer2011","Schmidt2016"),
                      guide = guide_legend(override.aes = list(linetype = c(rep("solid", 3)), shape = c(rep(NA,3))))) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Envelope-producer mass fraction, '*Phi[L])))
p_lipo_lambda_alt


# By substituting \Pi in formula for \Phi_R=f(\lambda_T), one accounts for changes in cell shape with growth rate when estimating \kappa_n
kappaL_est <- kappaL_est_Si
phiR_pertT_1 <-function(x,kN,y_int) {
  y_int+((kappaL_est+6.818*epsPar*(kappaL_est+kN))*x)/(kappaL_est*(dp_est-kN))
}

phiR_pertT_2 <-function(x,kN,y_int) {
  y_int+((kappaL_est+6.19*epsPar*(kappaL_est+kN))*x)/(kappaL_est*(dp_est-kN))
}

phiR_pertT_3 <-function(x,kN,y_int) {
  y_int+((kappaL_est+epsPar*(kappaL_est+kN)*(7.162 -0.944*x))*x)/(kappaL_est*(dp_est-kN))
}

phiR_pertT_4 <-function(x,kN,y_int) {
  y_int+((kappaL_est+4.907*epsPar*(kappaL_est+kN))*x)/(kappaL_est*(dp_est-kN))
}

phiR_pertT_5 <-function(x,kN,y_int) {
  y_int+((kappaL_est+4.838*epsPar*(kappaL_est+kN))*x)/(kappaL_est*(dp_est-kN))
}

phiR_pertT_6 <-function(x,kN,y_int) {
  y_int+(x*(kappaL_est+epsPar*(kappaL_est+kN)*(7.128 +1.862*x)))/(kappaL_est*(dp_est-kN))
}

phiR_pertT_7 <-function(x,kN,y_int) {
  y_int+((kappaL_est+8.7*epsPar*(kappaL_est+kN))*x)/(kappaL_est*(dp_est-kN))
}

phiR_pertT_8 <-function(x,kN,y_int) {
  y_int+((kappaL_est+8.22*epsPar*(kappaL_est+kN))*x)/(kappaL_est*(dp_est-kN))
}

phiR_pertT_9 <-function(x,kN,y_int) {
  y_int+((kappaL_est+4.156*epsPar*(kappaL_est+kN))*x)/(kappaL_est*(dp_est-kN))
}

# Phi_R and growth rates from Si2017; Then infer \kappa_n
y <- df_si_1$RNA.protein
x <- df_si_1$growth.rate..1.hours.
m_kN1 <- nls(y~phiR_pertT_1(x,kN,y_int))
m_kN1_summ <- summary(m_kN1)

y <- df_si_2$RNA.protein
x <- df_si_2$growth.rate..1.hours.
m_kN2 <- nls(y~phiR_pertT_2(x,kN,y_int))
m_kN2_summ <- summary(m_kN2)

y <- df_si_3$RNA.protein
x <- df_si_3$growth.rate..1.hours.
m_kN3 <- nls(y~phiR_pertT_3(x,kN,y_int))
m_kN3_summ <- summary(m_kN3)

y <- df_si_4$RNA.protein
x <- df_si_4$growth.rate..1.hours.
m_kN4 <- nls(y~phiR_pertT_4(x,kN,y_int))
m_kN4_summ <- summary(m_kN4)

y <- df_si_5$RNA.protein
x <- df_si_5$growth.rate..1.hours.
m_kN5 <- nls(y~phiR_pertT_5(x,kN,y_int))
m_kN5_summ <- summary(m_kN5)

y <- df_si_6$RNA.protein
x <- df_si_6$growth.rate..1.hours.
m_kN6 <- nls(y~phiR_pertT_6(x,kN,y_int))
m_kN6_summ <- summary(m_kN6)

y <- df_si_7$RNA.protein
x <- df_si_7$growth.rate..1.hours.
m_kN7 <- nls(y~phiR_pertT_7(x,kN,y_int))
m_kN7_summ <- summary(m_kN7)

y <- df_si_8$RNA.protein
x <- df_si_8$growth.rate..1.hours.
m_kN8 <- nls(y~phiR_pertT_8(x,kN,y_int))
m_kN8_summ <- summary(m_kN8)

y <- df_si_9$RNA.protein
x <- df_si_9$growth.rate..1.hours.
m_kN9 <- nls(y~phiR_pertT_9(x,kN,y_int))
m_kN9_summ <- summary(m_kN9)

kappaN1_est <- m_kN1_summ$coefficients[1]
kappaN2_est <- m_kN2_summ$coefficients[1]
kappaN3_est <- m_kN3_summ$coefficients[1]
kappaN4_est <- m_kN4_summ$coefficients[1]
kappaN5_est <- m_kN5_summ$coefficients[1]
kappaN6_est <- m_kN6_summ$coefficients[1]
kappaN7_est <- m_kN7_summ$coefficients[1]
kappaN8_est <- m_kN8_summ$coefficients[1]
kappaN9_est <- m_kN9_summ$coefficients[1]

y1_est <- m_kN1_summ$coefficients[2]
y2_est <- m_kN2_summ$coefficients[2]
y3_est <- m_kN3_summ$coefficients[2]
y4_est <- m_kN4_summ$coefficients[2]
y5_est <- m_kN5_summ$coefficients[2]
y6_est <- m_kN6_summ$coefficients[2]
y7_est <- m_kN7_summ$coefficients[2]
y8_est <- m_kN8_summ$coefficients[2]
y9_est <- m_kN9_summ$coefficients[2]

p_ribo_ab_pert <- ggplot() + 
  stat_function(fun=function(x) phiR_pertT_1(x,kappaN1_est,y1_est), aes(colour = "a"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_2(x,kappaN2_est,y2_est), aes(colour = "b"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_3(x,kappaN3_est,y3_est), aes(colour = "c"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_4(x,kappaN4_est,y4_est), aes(colour = "d"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_5(x,kappaN5_est,y5_est), aes(colour = "e"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_6(x,kappaN6_est,y6_est), aes(colour = "f"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_7(x,kappaN7_est,y7_est), aes(colour = "g"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_8(x,kappaN8_est,y8_est), aes(colour = "h"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_9(x,kappaN9_est,y9_est), aes(colour = "i"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=15, aes(x=df_si_1$growth.rate..1.hours., y=df_si_1$RNA.protein, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=16, aes(x=df_si_2$growth.rate..1.hours., y=df_si_2$RNA.protein, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=17, aes(x=df_si_3$growth.rate..1.hours., y=df_si_3$RNA.protein, colour = "c")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_si_4$growth.rate..1.hours., y=df_si_4$RNA.protein, colour = "d")) +
  geom_point(size=4.5, stroke=1.5, shape=19, aes(x=df_si_5$growth.rate..1.hours., y=df_si_5$RNA.protein, colour = "e")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_si_6$growth.rate..1.hours., y=df_si_6$RNA.protein, colour = "f")) +
  geom_point(size=4.5, stroke=1.5, shape=1, aes(x=df_si_7$growth.rate..1.hours., y=df_si_7$RNA.protein, colour = "g")) +
  geom_point(size=4.5, stroke=1.5, shape=2, aes(x=df_si_8$growth.rate..1.hours., y=df_si_8$RNA.protein, colour = "h")) +
  geom_point(size=4.5, stroke=1.5, shape=5, aes(x=df_si_9$growth.rate..1.hours., y=df_si_9$RNA.protein, colour = "i")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name="Medium", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                           "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6],
                                           "g"=colBlindScale[7], "h"=colBlindScale[8], "i"=colBlindScale[9]), 
                      labels = c("Glucose (Si2017)","Glucose + 12 a.a. (Si2017)","Glucose + 6 a.a. (Si2017)",
                                 "Glucose + caa (Si2017)","Glucose synth. rich (Si2017)","Glycerol (Si2017)",
                                 "Mannose (Si2017)","Sorbitol (Si2017)","TSB (Si2017)")) +
  scale_shape_manual(name="Medium", values = c(15,16,17,18,19,0,1,2,5), 
                     labels = c("Glucose (Si2017)","Glucose + 12 a.a. (Si2017)","Glucose + 6 a.a. (Si2017)",
                                "Glucose + caa (Si2017)","Glucose synth. rich (Si2017)","Glycerol (Si2017)",
                                "Mannose (Si2017)","Sorbitol (Si2017)","TSB (Si2017)")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) + ylim(c(0.05,0.40)) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0,1,2,5),linetype=0,stroke=1.5)))
p_ribo_ab_pert

# Finally, infer \kappa_t from Eq 27; Note that it is independent of \Pi so there is no need to account for cell shape
kappaT_est <- (1-dp_est*olsR$coefficients[2])/olsR$coefficients[2]


# Correct for temperature using Q20
q10_factor <- 2.5^((20-37)/10)

kN1_q10_alt <- q10_corrected(kappaN1_est, 2.5, 37)
kN2_q10_alt <- q10_corrected(kappaN2_est, 2.5, 37)
kN3_q10_alt <- q10_corrected(kappaN3_est, 2.5, 37)
kN4_q10_alt <- q10_corrected(kappaN4_est, 2.5, 37)
kN5_q10_alt <- q10_corrected(kappaN5_est, 2.5, 37)
kN6_q10_alt <- q10_corrected(kappaN6_est, 2.5, 37)
kN7_q10_alt <- q10_corrected(kappaN7_est, 2.5, 37)
kN8_q10_alt <- q10_corrected(kappaN8_est, 2.5, 37)
kN9_q10_alt <- q10_corrected(kappaN9_est, 2.5, 37)
kL_q10_alt <- q10_corrected(kappaL_est, 2.5, 37)
kT_q10 <- q10_corrected(kappaT_est, 2.5, 37)

dp_q10 <- q10_corrected(dp_est, 2.5, 37)
dl_q10 <- q10_corrected(dl_est, 2.5, 37)

kN_q10_mean_alt <- mean(c(kN1_q10_alt, kN2_q10_alt, kN3_q10_alt, kN4_q10_alt, kN5_q10_alt, kN6_q10_alt, kN7_q10_alt, kN8_q10_alt, kN9_q10_alt))

c(kN1_q10_alt, kN2_q10_alt, kN3_q10_alt, kN4_q10_alt, kN5_q10_alt, kN6_q10_alt, kN7_q10_alt, kN8_q10_alt, kN9_q10_alt)
c(kT_q10, kL_q10_alt, dp_q10, dl_q10)
c(kappaT_est, kappaL_est, epsPar)

c(kappaN1_est, kappaN2_est, kappaN3_est, kappaN4_est, kappaN5_est, kappaN6_est, kappaN7_est, kappaN8_est, kappaN9_est, kappaL_est, kappaT_est)

ggsave(file="p_phiR_ab_pert.pdf", plot=p_ribo_ab_pert, width = 5.2, height = 5.2)
ggsave(file="p_lipo_nutri.pdf", plot=p_lipo_nut_pert, width = 5.2, height = 5.2)
ggsave(file="p_sv_growth_nut.pdf", plot=p_sv_growth_nut, width = 5.2, height = 5.2)
ggsave(file="p_sv_growth_ab.pdf", plot=p_sv_growth_ab, width = 5.2, height = 5.2)
ggsave(file="p_lipo_lambda_alt.pdf", plot=p_lipo_lambda_alt, width = 5.2, height = 5.2)


#========== Plotting legends ==========
legend_si <- ggplot() + 
  stat_function(fun=function(x) phiR_pertT_1(x,kappaN1_est,y1_est), aes(colour = "a"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_2(x,kappaN2_est,y2_est), aes(colour = "b"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_3(x,kappaN3_est,y3_est), aes(colour = "c"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_4(x,kappaN4_est,y4_est), aes(colour = "d"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_5(x,kappaN5_est,y5_est), aes(colour = "e"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_6(x,kappaN6_est,y6_est), aes(colour = "f"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_7(x,kappaN7_est,y7_est), aes(colour = "g"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_8(x,kappaN8_est,y8_est), aes(colour = "h"), size = 1.5) +
  stat_function(fun=function(x) phiR_pertT_9(x,kappaN9_est,y9_est), aes(colour = "i"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=15, aes(x=df_si_1$growth.rate..1.hours., y=df_si_1$RNA.protein, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=16, aes(x=df_si_2$growth.rate..1.hours., y=df_si_2$RNA.protein, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=17, aes(x=df_si_3$growth.rate..1.hours., y=df_si_3$RNA.protein, colour = "c")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_si_4$growth.rate..1.hours., y=df_si_4$RNA.protein, colour = "d")) +
  geom_point(size=4.5, stroke=1.5, shape=19, aes(x=df_si_5$growth.rate..1.hours., y=df_si_5$RNA.protein, colour = "e")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_si_6$growth.rate..1.hours., y=df_si_6$RNA.protein, colour = "f")) +
  geom_point(size=4.5, stroke=1.5, shape=1, aes(x=df_si_7$growth.rate..1.hours., y=df_si_7$RNA.protein, colour = "g")) +
  geom_point(size=4.5, stroke=1.5, shape=2, aes(x=df_si_8$growth.rate..1.hours., y=df_si_8$RNA.protein, colour = "h")) +
  geom_point(size=4.5, stroke=1.5, shape=5, aes(x=df_si_9$growth.rate..1.hours., y=df_si_9$RNA.protein, colour = "i")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name="Medium", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                                "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6],
                                                "g"=colBlindScale[7], "h"=colBlindScale[8], "i"=colBlindScale[9]), 
                      labels = c("Glucose (Si2017)","Glucose + 12 a.a. (Si2017)","Glucose + 6 a.a. (Si2017)",
                                 "Glucose + caa (Si2017)","Glucose synth. rich (Si2017)","Glycerol (Si2017)",
                                 "Mannose (Si2017)","Sorbitol (Si2017)","TSB (Si2017)")) +
  scale_shape_manual(name="Medium", values = c(15,16,17,18,19,0,1,2,5), 
                     labels = c("Glucose (Si2017)","Glucose + 12 a.a. (Si2017)","Glucose + 6 a.a. (Si2017)",
                                "Glucose + caa (Si2017)","Glucose synth. rich (Si2017)","Glycerol (Si2017)",
                                "Mannose (Si2017)","Sorbitol (Si2017)","TSB (Si2017)")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) + ylim(c(0.05,0.40)) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0,1,2,5),linetype=0,stroke=1.5)))
legend_si <- cowplot::get_legend(legend_si)
ggsave(file="legend_si.pdf", legend_si, width = 8, height = 8)


legend_sv_source <- ggplot() + 
  geom_point(size=4.5, shape=15, 
             aes(x=df_si_nutrients$growth.rate..1.hours., 
                 y=df_si_nutrients$surface.to.volume.ratio..μm.1.,
                 col="a")) +
  geom_point(size=4.5, shape=16, 
             aes(x=df_basan$growth.rate, 
                 y=df_basan$sv,
                 col="b")) +
  geom_point(size=4.5, shape=17, 
             aes(x=df_volkmer$growth_rate.mean, 
                 y=df_volkmer$sv,
                 col="c")) +
  stat_function(fun=function(x) lm_volkmer_nutrients$coefficients[1]+lm_volkmer_nutrients$coefficients[2]*x,
                aes(col="c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) lm_basan_nutrients$coefficients[1]+lm_basan_nutrients$coefficients[2]*x,
                aes(col="b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) lm_si_nutrients$coefficients[1]+lm_si_nutrients$coefficients[2]*x,
                aes(col="a"), size=line_thick, n=plot_points) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual("Source", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3]), 
                      labels = c("Volkemer2011","Basan2015","Si2017")) +
  scale_shape_manual(name="Source", values = c(15,16,17), 
                     labels = c("Volkemer2011","Basan2015","Si2017")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('S/V, '*Pi*' ('*mu*''*m^-1*')]'))) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17),linetype=0)))
legend_sv_source <- cowplot::get_legend(legend_sv_source)
ggsave(file="legend_sv_source.pdf", legend_sv_source, width = 8, height = 8)

legend_sv_ab <- ggplot() + 
  geom_point(size=4.5, 
             aes(x=df_si_1$growth.rate..1.hours., 
                 y=df_si_1$surface.to.volume.ratio..μm.1.,
                 col="a"), shape=15, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_2$growth.rate..1.hours., 
                 y=df_si_2$surface.to.volume.ratio..μm.1.,
                 col="b"), shape=16, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_3$growth.rate..1.hours., 
                 y=df_si_3$surface.to.volume.ratio..μm.1.,
                 col="c"), shape=17, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_4$growth.rate..1.hours., 
                 y=df_si_4$surface.to.volume.ratio..μm.1.,
                 col="d"), shape=18, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_5$growth.rate..1.hours., 
                 y=df_si_5$surface.to.volume.ratio..μm.1.,
                 col="e"), shape=19, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_6$growth.rate..1.hours., 
                 y=df_si_6$surface.to.volume.ratio..μm.1.,
                 col="f"), shape=0, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_7$growth.rate..1.hours., 
                 y=df_si_7$surface.to.volume.ratio..μm.1.,
                 col="g"), shape=1, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_8$growth.rate..1.hours., 
                 y=df_si_8$surface.to.volume.ratio..μm.1.,
                 col="h"), shape=2, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_si_9$growth.rate..1.hours., 
                 y=df_si_9$surface.to.volume.ratio..μm.1.,
                 col="i"), shape=5, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_basan_Cm$growth.rate, 
                 y=df_basan_Cm$sv,
                 col="j"), shape=6, stroke=1.5) +
  geom_point(size=4.5, 
             aes(x=df_vadia$growth.rate, 
                 y=df_vadia$sv,
                 col="k"), shape=7, stroke=1.5) +
  stat_function(fun=function(x) lm_si_3$coefficients[1]+lm_si_3$coefficients[2]*x,
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) lm_si_6$coefficients[1]+lm_si_6$coefficients[2]*x,
                aes(colour = "f"), size=line_thick, n=plot_points) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual("Medium", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                           "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6],
                                           "g"=colBlindScale[7], "h"=colBlindScale[8], "i"=colBlindScale[9],
                                           "j"=colBlindScale[10], "k"=colBlindScale[11]), 
                      labels = c("Glucose (Si2017)","Glucose + 12 a.a. (Si2017)","Glucose + 6 a.a. (Si2017)",
                                 "Glucose + caa (Si2017)","Glucose synth. rich (Si2017)","Glycerol (Si2017)",
                                 "Mannose (Si2017)","Sorbitol (Si2017)","TSB (Si2017)",
                                 "Glucose (Basan2015)","LB + glucose (Vadia2017)")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('S/V, '*Pi*' ('*mu*''*m^-1*')'))) +
  ylim(c(3.5,9)) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0,1,2,5,6,7),linetype=0)))
legend_sv_ab <- cowplot::get_legend(legend_sv_ab)
ggsave(file="legend_sv_ab.pdf", legend_sv_ab, width = 8, height = 8)


#========== Inferring capacities from OLS slopes ==========
# The rational here is that \Pi is largely independent of \lambda_T, given that it only shows negative dependence in two types of media
# Therefore, one can assume that \Pi is fixed and use the value which is average for E. coli

SAVpar <- 5 # Rough average in our dataset; More precisely, it would be 4.94 but this might change if more data is added
epsPar <- (0.30/0.55)/SAVpar
dp_est <- 0.05 # page 2896 in Nath and Koch 1970 (taken for slowly decaying component)
dl_est <- 1.00

# We take \kappa_l inferred using data on \Pi~\lambda_N. This is because we use proteomic data from Schmidt2016, and they the same media as Volkmer2011
kappaL_est <- kappaL_est_Volkmer
kappaN_est <- function(x) {
  (kappaL_est*(x*dp_est-epsPar*SAVpar-1))/(x*kappaL_est+epsPar*SAVpar)
}

# Calculating standard error of the estimate by error propagation
# This is the temperature correction via Q10 method
q10_factor <- 2.5^((20-37)/10)

# Phi_R ~ lambda_N
# Get the estimate of slope and its variance
olsR
var_b <- (olsR$coefficients[2,2]*sqrt(20))^2
est_b <- olsR$coefficients[2,1]
var_kT <- (q10_factor^2/est_b^4)*var_b
kT_q10_err <- sqrt(var_kT)/sqrt(20)

# Phi_L ~ f(lambda_N)
# Get the estimate of slope and its variance
m_summary_Volkmer
var_kL <- (sqrt(20)*m_summary_Volkmer$parameters[2,2])^2
kL_q10_err <- sqrt(q10_factor^2*var_kL)/sqrt(20)

# Phi_R ~ lambda_T
# Get the estimate of slope and its variance
kN_std_q10_est <- function(x,n) {
  var_b <- (x$coefficients[2,2]*sqrt(n))^2
  est_b <- x$coefficients[2,1]
  var_kN <- (q10_factor^2*var_kL*epsPar^2*SAVpar^2*(1-est_b*dp_est+epsPar*SAVpar)^2)/((est_b*kappaL_est+epsPar*SAVpar)^4) + 
            (q10_factor^2*var_b*kappaL_est^2*(kappaL_est+epsPar*SAVpar*(dp_est+kappaL_est))^2)/((est_b*kappaL_est+epsPar*SAVpar)^4)
  return(sqrt(var_kN)/sqrt(5))
}

kappaN1_q10_err <- kN_std_q10_est(olsR_ab1, 5)
kappaN2_q10_err <- kN_std_q10_est(olsR_ab2, 5)
kappaN3_q10_err <- kN_std_q10_est(olsR_ab3, 5)
kappaN4_q10_err <- kN_std_q10_est(olsR_ab4, 5)
kappaN5_q10_err <- kN_std_q10_est(olsR_ab5, 5)
kappaN6_q10_err <- kN_std_q10_est(olsR_ab6, 5)
c(kappaN1_q10_err, kappaN2_q10_err, kappaN3_q10_err, kappaN4_q10_err, kappaN5_q10_err, kappaN6_q10_err,
  kT_q10_err, kL_q10_err)

kappaN1_est <- kappaN_est(olsR_ab1$coefficients[2])
kappaN2_est <- kappaN_est(olsR_ab2$coefficients[2])
kappaN3_est <- kappaN_est(olsR_ab3$coefficients[2])
kappaN4_est <- kappaN_est(olsR_ab4$coefficients[2])
kappaN5_est <- kappaN_est(olsR_ab5$coefficients[2])
kappaN6_est <- kappaN_est(olsR_ab6$coefficients[2])

kN1_q10 <- q10_corrected(kappaN1_est, 2.5, 37)
kN2_q10 <- q10_corrected(kappaN2_est, 2.5, 37)
kN3_q10 <- q10_corrected(kappaN3_est, 2.5, 37)
kN4_q10 <- q10_corrected(kappaN4_est, 2.5, 37)
kN5_q10 <- q10_corrected(kappaN5_est, 2.5, 37)
kN6_q10 <- q10_corrected(kappaN6_est, 2.5, 37)
kL_q10 <- q10_corrected(kappaL_est, 2.5, 37)
kT_q10 <- q10_corrected(kappaT_est, 2.5, 37)

dp_q10 <- q10_corrected(dp_est, 2.5, 37)
dl_q10 <- q10_corrected(dl_est, 2.5, 37)

kN_q10_mean <- mean(c(kN1_q10, kN2_q10, kN3_q10, kN4_q10, kN5_q10, kN6_q10))

c(kN1_q10, kN2_q10, kN3_q10, kN4_q10, kN5_q10, kN6_q10)
c(kT_q10, kL_q10, dp_q10, dl_q10)
c(kappaT_est, kappaL_est, epsPar)


#========== Plotting figures ==========
p_ribo_nutri <- ggplot() + 
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, aes(colour = "g"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=15, aes(x=df_fin_Forchhammer$growth_rates, y=df_fin_Forchhammer$phi_R, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=16, aes(x=df_fin_Bremer$growth_rates, y=df_fin_Bremer$phi_R, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=17, aes(x=df_fin_Scott$growth_rates, y=df_fin_Scott$phi_R, colour = "c")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_fin_Valgepea$growth_rates, y=df_fin_Valgepea$phi_R, colour = "d")) +
  geom_point(size=4.5, stroke=1.5, shape=19, aes(x=df_fin_Li$growth_rates, y=df_fin_Li$phi_R, colour = "e")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_fin_Peebo$growth_rates, y=df_fin_Peebo$phi_R, colour = "f")) +
  geom_point(size=4.5, stroke=1.5, shape=1, aes(x=df_fin_Schmidt$growth_rates, y=df_fin_Schmidt$phi_R, colour = "g")) +
  geom_point(size=4.5, stroke=1.5, shape=2, aes(x=df_fin_Erickson$growth_rates, y=df_fin_Erickson$phi_R, colour = "h")) +
  geom_point(size=4.5, stroke=1.5, shape=5, aes(x=df_fin_Si$growth_rates, y=df_fin_Si$phi_R, colour = "i")) +
  geom_point(size=4.5, stroke=1.5, shape=6, aes(x=df_fin_Dai$growth_rates, y=df_fin_Dai$phi_R, colour = "j")) +
  geom_point(size=4.5, stroke=1.5, shape=7, aes(x=df_fin_Mori$growth_rates, y=df_fin_Mori$phi_R, colour = "k")) +
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
  scale_colour_manual(name="Source", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                                "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6],
                                                "g"=colBlindScale[7], "h"=colBlindScale[8], "i"=colBlindScale[9],
                                                "j"=colBlindScale[10], "k"=colBlindScale[11]), 
                      labels = c("Forchhammer1971","Bremer1996","Scott2010",
                                 "Valgepea2013","Li2014","Peebo2015",
                                 "Schmidt2016","Erickson2017","Si2017","Dai2018","Mori2021")) +
  scale_shape_manual(name="Medium", values = c(15,16,17,18,19,0,1,2,5,6,7), 
                     labels = c("Forchhammer1971","Bremer1996","Scott2010",
                                "Valgepea2013","Li2014","Peebo2015",
                                "Schmidt2016","Erickson2017","Si2017","Dai2018","Mori2021")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0,1,2,5,6,7),linetype=0,stroke=1.5)))
p_ribo_nutri


df_eco_ab_ribo_1 <- df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="M63 + Glyc",]
df_eco_ab_ribo_2 <- df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="M63 + Glc",]
df_eco_ab_ribo_3 <- df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="cAA + Glyc",]
df_eco_ab_ribo_4 <- df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="cAA + Glc",]
df_eco_ab_ribo_5 <- df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="RDM + Glyc",]
df_eco_ab_ribo_6 <- df_eco_ab_ribo[df_eco_ab_ribo$growth_medium=="RDM + Glc",]
p_ribo_ab <- ggplot() + 
  stat_function(fun=function(x) olsR_ab1$coefficients[1]+olsR_ab1$coefficients[2]*x, aes(colour = "a"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab2$coefficients[1]+olsR_ab2$coefficients[2]*x, aes(colour = "b"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab3$coefficients[1]+olsR_ab3$coefficients[2]*x, aes(colour = "c"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab4$coefficients[1]+olsR_ab4$coefficients[2]*x, aes(colour = "d"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab5$coefficients[1]+olsR_ab5$coefficients[2]*x, aes(colour = "e"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab6$coefficients[1]+olsR_ab6$coefficients[2]*x, aes(colour = "f"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=15, aes(x=df_eco_ab_ribo_1$growth_rates, y=df_eco_ab_ribo_1$phi_R, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=16, aes(x=df_eco_ab_ribo_2$growth_rates, y=df_eco_ab_ribo_2$phi_R, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=17, aes(x=df_eco_ab_ribo_3$growth_rates, y=df_eco_ab_ribo_3$phi_R, colour = "c")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_eco_ab_ribo_4$growth_rates, y=df_eco_ab_ribo_4$phi_R, colour = "d")) +
  geom_point(size=4.5, stroke=1.5, shape=19, aes(x=df_eco_ab_ribo_5$growth_rates, y=df_eco_ab_ribo_5$phi_R, colour = "e")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_eco_ab_ribo_6$growth_rates, y=df_eco_ab_ribo_6$phi_R, colour = "f")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name="Source", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                                "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6]), 
                      labels = c("M63 + Glyc","M63 + Glc","cAA + Glyc",
                                 "cAA + Glc","RDM + Glyc","RDM + Glc")) +
  scale_shape_manual(name="Medium", values = c(15,16,17,18,19,0), 
                     labels = c("M63 + Glyc","M63 + Glc","cAA + Glyc",
                                "cAA + Glc","RDM + Glyc","RDM + Glc")) +
  xlim(c(0.0,2.00)) + ylim(c(0.0,0.4)) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0),linetype=0,stroke=1.5)))
p_ribo_ab


df_fin_nutrients <- df_eco_env_ribo[df_eco_env_ribo$perturbation=="nutrient conditions",]
df_fin_triclosan <- df_eco_env_ribo[df_eco_env_ribo$perturbation=="triclosan",]
df_fin_fosfomycin <- df_eco_env_ribo[df_eco_env_ribo$perturbation=="fosfomycin",]

p_ribo_env_ab <- ggplot() + 
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, aes(colour = "a"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=21, aes(x=df_fin_nutrients$growth_rates, y=df_fin_nutrients$phi_R, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_fin_triclosan$growth_rates, y=df_fin_triclosan$phi_R, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_fin_fosfomycin$growth_rates, y=df_fin_fosfomycin$phi_R, colour = "c")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name="Perturbation", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3]), 
                      labels = c("Nutrients","Triclosan","Fosfomycin")) +
  scale_shape_manual(name="Pertrubation", values = c(21,18,0), 
                     labels = c("Nutrients","Triclosan","Fosfomycin")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) +
  guides(colour=guide_legend(override.aes=list(shape=c(21,18,0),linetype=0,stroke=1.5)))
p_ribo_env_ab

ribo_legend1 <- ggplot() + 
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, aes(colour = "g"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=15, aes(x=df_fin_Forchhammer$growth_rates, y=df_fin_Forchhammer$phi_R, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=16, aes(x=df_fin_Bremer$growth_rates, y=df_fin_Bremer$phi_R, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=17, aes(x=df_fin_Scott$growth_rates, y=df_fin_Scott$phi_R, colour = "c")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_fin_Valgepea$growth_rates, y=df_fin_Valgepea$phi_R, colour = "d")) +
  geom_point(size=4.5, stroke=1.5, shape=19, aes(x=df_fin_Li$growth_rates, y=df_fin_Li$phi_R, colour = "e")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_fin_Peebo$growth_rates, y=df_fin_Peebo$phi_R, colour = "f")) +
  geom_point(size=4.5, stroke=1.5, shape=1, aes(x=df_fin_Schmidt$growth_rates, y=df_fin_Schmidt$phi_R, colour = "g")) +
  geom_point(size=4.5, stroke=1.5, shape=2, aes(x=df_fin_Erickson$growth_rates, y=df_fin_Erickson$phi_R, colour = "h")) +
  geom_point(size=4.5, stroke=1.5, shape=5, aes(x=df_fin_Si$growth_rates, y=df_fin_Si$phi_R, colour = "i")) +
  geom_point(size=4.5, stroke=1.5, shape=6, aes(x=df_fin_Dai$growth_rates, y=df_fin_Dai$phi_R, colour = "j")) +
  geom_point(size=4.5, stroke=1.5, shape=7, aes(x=df_fin_Mori$growth_rates, y=df_fin_Mori$phi_R, colour = "k")) +
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
  scale_colour_manual(name="Source", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                                "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6],
                                                "g"=colBlindScale[7], "h"=colBlindScale[8], "i"=colBlindScale[9],
                                                "j"=colBlindScale[10], "k"=colBlindScale[11]), 
                      labels = c("Forchhammer1971","Bremer1996","Scott2010",
                                 "Valgepea2013","Li2014","Peebo2015",
                                 "Schmidt2016","Erickson2017","Si2017","Dai2018","Mori2021")) +
  scale_shape_manual(name="Medium", values = c(15,16,17,18,19,0,1,2,5,6,7), 
                     labels = c("Forchhammer1971","Bremer1996","Scott2010",
                                "Valgepea2013","Li2014","Peebo2015",
                                "Schmidt2016","Erickson2017","Si2017","Dai2018","Mori2021")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0,1,2,5,6,7),linetype=0,stroke=1.5)))
legend1 <- cowplot::get_legend(ribo_legend1)
ggsave(file="legend1.pdf", legend1, width = 8, height = 8)

ribo_legend2 <- ggplot() + 
  stat_function(fun=function(x) olsR_ab1$coefficients[1]+olsR_ab1$coefficients[2]*x, aes(colour = "a"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab2$coefficients[1]+olsR_ab2$coefficients[2]*x, aes(colour = "b"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab3$coefficients[1]+olsR_ab3$coefficients[2]*x, aes(colour = "c"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab4$coefficients[1]+olsR_ab4$coefficients[2]*x, aes(colour = "d"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab5$coefficients[1]+olsR_ab5$coefficients[2]*x, aes(colour = "e"), size = 1.5) +
  stat_function(fun=function(x) olsR_ab6$coefficients[1]+olsR_ab6$coefficients[2]*x, aes(colour = "f"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=15, aes(x=df_eco_ab_ribo_1$growth_rates, y=df_eco_ab_ribo_1$phi_R, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=16, aes(x=df_eco_ab_ribo_2$growth_rates, y=df_eco_ab_ribo_2$phi_R, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=17, aes(x=df_eco_ab_ribo_3$growth_rates, y=df_eco_ab_ribo_3$phi_R, colour = "c")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_eco_ab_ribo_4$growth_rates, y=df_eco_ab_ribo_4$phi_R, colour = "d")) +
  geom_point(size=4.5, stroke=1.5, shape=19, aes(x=df_eco_ab_ribo_5$growth_rates, y=df_eco_ab_ribo_5$phi_R, colour = "e")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_eco_ab_ribo_6$growth_rates, y=df_eco_ab_ribo_6$phi_R, colour = "f")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name="Source", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                                "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6]), 
                      labels = c("M63 + Glyc","M63 + Glc","cAA + Glyc",
                                 "cAA + Glc","RDM + Glyc","RDM + Glc")) +
  scale_shape_manual(name="Medium", values = c(15,16,17,18,19,0), 
                     labels = c("M63 + Glyc","M63 + Glc","cAA + Glyc",
                                "cAA + Glc","RDM + Glyc","RDM + Glc")) +
  xlim(c(0.0,1.75)) + ylim(c(0.0,0.4)) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) +
  guides(colour=guide_legend(override.aes=list(shape=c(15,16,17,18,19,0),linetype=0,stroke=1.5)))
legend2 <- cowplot::get_legend(ribo_legend2)
ggsave(file="legend2.pdf", legend2, width = 8, height = 8)

ribo_legend3 <- ggplot() + 
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, aes(colour = "a"), size = 1.5) +
  geom_point(size=4.5, stroke=1.5, shape=21, aes(x=df_fin_nutrients$growth_rates, y=df_fin_nutrients$phi_R, colour = "a")) +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=df_fin_triclosan$growth_rates, y=df_fin_triclosan$phi_R, colour = "b")) +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=df_fin_fosfomycin$growth_rates, y=df_fin_fosfomycin$phi_R, colour = "c")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_colour_manual(name="Perturbation", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3]), 
                      labels = c("Nutrients","Triclosan","Fosfomycin")) +
  scale_shape_manual(name="Pertrubation", values = c(21,18,0), 
                     labels = c("Nutrients","Triclosan","Fosfomycin")) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R]))) +
  guides(colour=guide_legend(override.aes=list(shape=c(21,18,0),linetype=0,stroke=1.5)))
legend3 <- cowplot::get_legend(ribo_legend3)
ggsave(file="legend3.pdf", plot=legend3, width = 8, height = 8)

ggsave(file="p_ribo_nutri.pdf", plot=p_ribo_nutri, width = 5, height = 5)
ggsave(file="p_ribo_ab.pdf", plot=p_ribo_ab, width = 5, height = 5)
ggsave(file="p_ribo_env_ab.pdf", plot=p_ribo_env_ab, width = 5, height = 5)


#========== Loading cross-species data ==========
line_thick = 2.5
point_size1 = 3
point_size2 = 5.5
L_env <- 0.03
gammaShapeConst <- 2*pi
anaerobic_cutoff <- 11

df_lambda_env_corr <- read.csv("/home/bogi/Desktop/growth_rate/data/main/growth_scaling.shape.envelope_corrected.csv",
                               skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_lambda_env_corr <- df_lambda_env_corr[!is.na(df_lambda_env_corr$total_thickness_nm),]
df_lambda_env_corr[is.na(df_lambda_env_corr$tca_steps),]$tca_steps <- 12
df_lambda_env_corr$growth_rate <- df_lambda_env_corr$growth_rate/24
df_lambda_env_corr$growth_rate <- ifelse(df_lambda_env_corr$tca_steps<=anaerobic_cutoff,
                                df_lambda_env_corr$growth_rate*lambda_anaero_correct(kN6_q10,kT_q10,kL_q10_alt,epsPar,
                                                                                     df_lambda_env_corr$S.mean/df_lambda_env_corr$V.mean,dp_q10,dl_q10),
                                df_lambda_env_corr$growth_rate)

df_lambda <- read.csv("/home/bogi/Desktop/growth_rate/data/main/growth_scaling.shape.genome.csv",
                      skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_lambda[is.na(df_lambda$tca_steps),]$tca_steps <- 12
df_lambda$growth_rate <- df_lambda$growth_rate/24
df_lambda$growth_rate <- ifelse(df_lambda$tca_steps<=anaerobic_cutoff,
                                df_lambda$growth_rate*lambda_anaero_correct(kN6_q10,kT_q10,kL_q10_alt,epsPar,df_lambda$S.mean/df_lambda$V.mean,dp_q10,dl_q10),
                                df_lambda$growth_rate)

#df_lambda <- df_lambda %>% group_by(species) %>% slice(which.max(growth_rate))
#df_lambda_env_corr <- df_lambda_env_corr %>% group_by(species) %>% slice(which.max(growth_rate))
df_lambda$ID <- "Generic thickness (30 nm)"
df_lambda_env_corr$ID <- "Specific thickness"
df_lambda_full <- rbind(df_lambda, df_lambda_env_corr[!is.na(df_lambda_env_corr$total_thickness_nm),])
df_lambda_full$SV.mean <- df_lambda_full$S.mean/df_lambda_full$V.mean

# Add Grant2021 data
df_Grant <- read.csv("/home/bogi/Desktop/growth_rate/data/cell_size_data/ltee_data/grant_wiser_full.csv")
df_Grant$S.mean <- capsule_surf(df_Grant$mean.width-2*L_env, df_Grant$mean.length-2*L_env)
df_Grant$V.mean <- capsule_vol(df_Grant$mean.width-2*L_env, df_Grant$mean.length-2*L_env)
df_Grant$SV.mean <- df_Grant$S.mean/df_Grant$V.mean
df_Grant$fitnessClones <- df_Grant$fitnessClones*0.7726 # Value from Table 1 in Vasi1994
colnames(df_Grant)[3] <- "growth_rate"
df_Grant$growth_rate <- q10_corrected(df_Grant$growth_rate,2.5,37)
df_Grant_avg <- df_Grant %>% group_by(generation, ID) %>%  summarise(SV.mean=mean(SV.mean, na.rm = TRUE), growth_rate=mean(growth_rate, na.rm = TRUE), .groups='drop')

mean(capsule_vol(df_volkmer$cell_width.mean-2*L_env,df_volkmer$cell_length.mean-2*L_env))

# Add Lennon2021 data
df_Lennon <- read.csv("/home/bogi/Desktop/growth_rate/data/cell_size_data/ltee_data/lennon2021.csv")
df_Lennon$fitnessMP <- 1.5926*df_Lennon$fitnessMP # This value is from our dataset
colnames(df_Lennon)[2] <- "growth_rate"
df_Lennon$S.mean <- sphere_surf(df_Lennon$mean.width-2*0.01)
df_Lennon$V.mean <- sphere_vol(df_Lennon$mean.width-2*0.01)
df_Lennon$SV.mean <- df_Lennon$S.mean/df_Lennon$V.mean
df_Lennon$growth_rate <- q10_corrected(df_Lennon$growth_rate,2.5,37)

# Add Gallet2017 data -- Note that this is sligthly different because we are not subtracting L_env and we have S/V_{tot} not S/V_{cyt}
df_Gallet <- read.csv("/home/bogi/Desktop/growth_rate/data/cell_size_data/ltee_data/gallet2017.csv")
df_Gallet$S.mean <- gammaShapeConst*df_Gallet$V.mean^(2/3)
df_Gallet$SV.mean <- df_Gallet$S.mean/df_Gallet$V.mean
df_Gallet$growth_rate <- q10_corrected(df_Gallet$growth_rate,2.5,37)

df_Lennon_wt <- df_Lennon[df_Lennon$ID=="Mycoplasma mycoides (Moger-Reischer2021)",]
df_Lennon_jcvi <- df_Lennon[df_Lennon$ID=="JCVI Mycoplasma (Moger-Reischer2021)",]

p_scaled_growth <- ggplot() + 
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "e"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "f"), size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black", aes(x=log10(df_lambda$S.mean/df_lambda$V.mean), y=log10(df_lambda$growth_rate))) +
  geom_point(size=3.5, stroke=1.5, shape=25, aes(x=log10(df_Grant_avg$SV.mean), y=log10(df_Grant_avg$growth_rate), colour="x")) +
  geom_point(size=3.5, stroke=1.5, shape=18, aes(x=log10(df_Gallet$SV.mean), y=log10(df_Gallet$growth_rate), colour="x")) +
  geom_point(size=3.5, stroke=1.5, shape=1, aes(x=log10(df_Lennon_wt$SV.mean), y=log10(df_Lennon_wt$growth_rate), colour="x")) +
  geom_point(size=3.5, stroke=1.5, shape=0, aes(x=log10(df_Lennon_jcvi$SV.mean), y=log10(df_Lennon_jcvi$growth_rate), colour="x")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual("Nutrient capacity", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                                      "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6], "x"="red"), labels = medium_text) +
  guides(colour=guide_legend(override.aes=list(shape=c(25,18,1,0),linetype=0,stroke=1.5))) +
  ylim(c(-2.5,0.2))
p_scaled_growth

nutrient_legend <- ggplot() + 
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "e"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "f"), size=line_thick, n=plot_points) +
  geom_point() +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual("Nutrient capacity", values = c("a"=colBlindScale[1],"b"=colBlindScale[2],
                                                      "c"=colBlindScale[3],"d"=colBlindScale[4],
                                                      "e"=colBlindScale[5], "f"=colBlindScale[6]), labels = medium_text) +
  ylim(c(-2.5,0.2))
nutrient_legend <- cowplot::get_legend(nutrient_legend)


grantPlot <- ggplot() + 
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "e"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "f"), size=line_thick, n=plot_points) +
  geom_point(size=4.5, stroke=1.5, shape=25, aes(x=log10(df_Grant_avg$SV.mean), y=log10(df_Grant_avg$growth_rate)), colour = "red") +
  geom_point(size=4.5, stroke=1.5, shape=18, aes(x=log10(df_Gallet$SV.mean), y=log10(df_Gallet$growth_rate)), colour = "red") +
  geom_point(size=4.5, stroke=1.5, shape=1, aes(x=log10(df_Lennon_wt$SV.mean), y=log10(df_Lennon_wt$growth_rate)), colour = "red") +
  geom_point(size=4.5, stroke=1.5, shape=0, aes(x=log10(df_Lennon_jcvi$SV.mean), y=log10(df_Lennon_jcvi$growth_rate)), colour = "red") +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual("Source", values = c("a"=colBlindScale[1],"b"=colBlindScale[2],
                                           "c"=colBlindScale[3],"d"=colBlindScale[4],
                                           "e"=colBlindScale[5],"f"=colBlindScale[6]), labels = medium_text) +
  guides(colour=guide_legend(override.aes=list(shape=c(21,18,1,0),linetype=0,stroke=1.5))) +
  ylim(c(-1.1,0.0))
grantPlot

grant_legend <- ggplot() + 
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black", aes(x=log10(df_lambda$S.mean/df_lambda$V.mean), y=log10(df_lambda$growth_rate))) +
  geom_point(size=3.5, stroke=1.5, shape=18, aes(x=log10(df_Gallet$SV.mean), y=log10(df_Gallet$growth_rate), colour="a")) +
  geom_point(size=3.5, stroke=1.5, shape=25, aes(x=log10(df_Grant_avg$SV.mean), y=log10(df_Grant_avg$growth_rate), colour="b")) +
  geom_point(size=3.5, stroke=1.5, shape=1, aes(x=log10(df_Lennon_wt$SV.mean), y=log10(df_Lennon_wt$growth_rate), colour="c")) +
  geom_point(size=3.5, stroke=1.5, shape=0, aes(x=log10(df_Lennon_jcvi$SV.mean), y=log10(df_Lennon_jcvi$growth_rate), colour="d")) +
  geom_point(size=3.5, stroke=1.5, shape=23, aes(x=log10(df_Lennon_jcvi$SV.mean), y=log10(df_Lennon_jcvi$growth_rate), colour="e")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual(name="Source", values = c("x"="black","a"="red","b"="red","c"="red","d"="red","e"="red"),
                      labels = c("Bacterial species","Escherichia coli (Gallet2017)","Escherichia coli (Grant2021)",
                                 "M. mycoides (Moger-Reisher2021)","JCVI M. mycoides (Moger-Reisher2021)","Escherichia coli (Favate2021)")) +
  scale_shape_manual(name="Nutrient capacity", values = c(1,18,25,1,0,23), labels = c("x","a","b","c","d","e")) +
  guides(colour=guide_legend(override.aes=list(shape=c(1,18,25,1,0,23),linetype=0,stroke=1.5))) +
  ylim(c(-2.5,0.2))
grant_legend <- cowplot::get_legend(grant_legend)

df_phi <- read.csv("proteome_scaling.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
p_scaled_ribo_frac <- ggplot(df_phi, aes(x=log10(S.mean/V.mean), y=log10(ribo_frac))) + 
  stat_function(fun=function(x) log10(riboFrac_env(kN1_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), aes(colour="a"), 
                size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN2_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), aes(colour="b"), 
                size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN3_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), aes(colour="c"), 
                size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN4_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), aes(colour="d"), 
                size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN5_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), aes(colour="e"), 
                size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN6_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), aes(colour="f"), 
                size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  geom_point(size = 3.5, stroke=1.5, shape=23, alpha=0.75, aes(x=log10(df_Grant_avg[df_Grant_avg$generation==0,]$SV.mean), y=log10(0.2100)), col="red") +
  geom_point(size = 3.5, stroke=1.5, shape=23, alpha=0.75, aes(x=log10(df_Grant_avg[df_Grant_avg$generation==50,]$SV.mean), y=log10(0.2177)), col="red") +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Ribosomal mass fraction, Log'[10]*'['*Phi[R]*']'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual("Source", values = c("a"=colBlindScale[1],"b"=colBlindScale[2],
                                           "c"=colBlindScale[3],"d"=colBlindScale[4],
                                           "e"=colBlindScale[5],"f"=colBlindScale[6]), labels = medium_text) +
  ylim(c(-2.0,0.0))
p_scaled_ribo_frac


p_scaled_lipo_frac <- ggplot(df_phi, aes(x=log10(S.mean/V.mean), y=log10(lip_frac))) + 
  stat_function(fun=function(x) log10(lipoFrac_env(kN1_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                aes(colour="a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN2_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                aes(colour="b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN3_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                aes(colour="c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN4_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                aes(colour="d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN5_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                aes(colour="e"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN6_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                aes(colour="f"), size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  geom_point(size = 3.5, stroke=1.5, shape=23, alpha=0.75, aes(x=log10(df_Grant_avg[df_Grant_avg$generation==0,]$SV.mean), y=log10(0.034)), col="red") +
  geom_point(size = 3.5, stroke=1.5, shape=23, alpha=0.75, aes(x=log10(df_Grant_avg[df_Grant_avg$generation==50,]$SV.mean), y=log10(0.028)), col="red") +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Envelope-producer mass fraction, Log'[10]*'['*Phi[L]*']'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.75,-1.0)) +
  scale_colour_manual("Source", values = c("a"=colBlindScale[1],"b"=colBlindScale[2],
                                           "c"=colBlindScale[3],"d"=colBlindScale[4],
                                           "e"=colBlindScale[5],"f"=colBlindScale[6]), labels = medium_text)
p_scaled_lipo_frac

#growth_proteome_scaling <- grid.arrange(p_scaled_growth, legend_growth, p_scaled_ribo_frac, p_scaled_lipo_frac, nrow = 2)
#ggsave(file="growth_proteome_scaling.pdf", plot=growth_proteome_scaling, width = 12, height = 12)
ggsave(file="p_scaled_growth.pdf", plot=p_scaled_growth, width = 5.2, height = 5.2)
ggsave(file="p_scaled_ribo_frac.pdf", plot=p_scaled_ribo_frac, width = 5.2, height = 5.2)
ggsave(file="p_scaled_lipo_frac.pdf", plot=p_scaled_lipo_frac, width = 5.2, height = 5.2)
ggsave(file="legend_nutrients.pdf", plot=nutrient_legend, width = 7.2, height = 7.2)
ggsave(file="legend_grant.pdf", plot=grant_legend, width = 7.2, height = 7.2)
ggsave(file="grant_growth_scaling.pdf", plot=grantPlot, width = 5.2, height = 5.2)


#========== Plots of growth-scaling under alternative parameters (from Si2017) ==========
sv_array <- seq(0.5,2.0,0.01)
kappaN_est_alt <- function(x, yKappa_L) {(yKappa_L*(x*dp_est-epsPar*SAVpar-1))/(x*yKappa_L+epsPar*SAVpar)}

kappaN1_est_Basan <- kappaN_est_alt(olsR_ab1$coefficients[2], kappaL_est_Basan)
kappaN2_est_Basan <- kappaN_est_alt(olsR_ab2$coefficients[2], kappaL_est_Basan)
kappaN3_est_Basan <- kappaN_est_alt(olsR_ab3$coefficients[2], kappaL_est_Basan)
kappaN4_est_Basan <- kappaN_est_alt(olsR_ab4$coefficients[2], kappaL_est_Basan)
kappaN5_est_Basan <- kappaN_est_alt(olsR_ab5$coefficients[2], kappaL_est_Basan)
kappaN6_est_Basan <- kappaN_est_alt(olsR_ab6$coefficients[2], kappaL_est_Basan)
kN1_Basan_q10 <- q10_corrected(kappaN1_est_Basan, 2.5, 37)
kN2_Basan_q10 <- q10_corrected(kappaN2_est_Basan, 2.5, 37)
kN3_Basan_q10 <- q10_corrected(kappaN3_est_Basan, 2.5, 37)
kN4_Basan_q10 <- q10_corrected(kappaN4_est_Basan, 2.5, 37)
kN5_Basan_q10 <- q10_corrected(kappaN5_est_Basan, 2.5, 37)
kN6_Basan_q10 <- q10_corrected(kappaN6_est_Basan, 2.5, 37)
kL_Basan_q10 <- q10_corrected(kappaL_est_Basan, 2.5, 37)
kN_q10_mean_Basan <- mean(c(kN1_Basan_q10, kN2_Basan_q10, kN3_Basan_q10, kN4_Basan_q10, kN5_Basan_q10, kN6_Basan_q10))

kappaN1_est_Si <- kappaN_est_alt(olsR_ab1$coefficients[2], kappaL_est_Si)
kappaN2_est_Si <- kappaN_est_alt(olsR_ab2$coefficients[2], kappaL_est_Si)
kappaN3_est_Si <- kappaN_est_alt(olsR_ab3$coefficients[2], kappaL_est_Si)
kappaN4_est_Si <- kappaN_est_alt(olsR_ab4$coefficients[2], kappaL_est_Si)
kappaN5_est_Si <- kappaN_est_alt(olsR_ab5$coefficients[2], kappaL_est_Si)
kappaN6_est_Si <- kappaN_est_alt(olsR_ab6$coefficients[2], kappaL_est_Si)
kN1_Si_q10 <- q10_corrected(kappaN1_est_Si, 2.5, 37)
kN2_Si_q10 <- q10_corrected(kappaN2_est_Si, 2.5, 37)
kN3_Si_q10 <- q10_corrected(kappaN3_est_Si, 2.5, 37)
kN4_Si_q10 <- q10_corrected(kappaN4_est_Si, 2.5, 37)
kN5_Si_q10 <- q10_corrected(kappaN5_est_Si, 2.5, 37)
kN6_Si_q10 <- q10_corrected(kappaN6_est_Si, 2.5, 37)
kL_Si_q10 <- q10_corrected(kappaL_est_Si, 2.5, 37)
kN_q10_mean_Si <- mean(c(kN1_Si_q10, kN2_Si_q10, kN3_Si_q10, kN4_Si_q10, kN5_Si_q10, kN6_Si_q10))

kN_q10_mean
kN_q10_mean_alt
kN_q10_mean_Si
kN_q10_mean_Basan

kL_q10
kL_q10_alt
kL_Si_q10
kL_Basan_q10

p_scaled_growth_alt <- ggplot(data.frame(x = sv_array), aes(x)) + 
  stat_function(fun=function(x) log10(lambda_env(kN_q10_mean, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(col="a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN_q10_mean_Basan, kT_q10, kL_Basan_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(col="b"), size=1, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN_q10_mean_Si, kT_q10, kL_Si_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(col="c"), size=1, n=plot_points, linetype="dashed") +
  stat_function(fun=function(x) log10(lambda_env(kN_q10_mean_alt, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(col="d"), size=line_thick, n=plot_points) +
  scale_colour_manual(bquote("Source of "*Pi*"("*lambda[N]*")"), values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                                            "c"=colBlindScale[4],"d"=colBlindScale[5]), 
                      labels = c(bquote("Volkmer2011, "*Pi*"("*lambda[T]*") const."), bquote("Basan2015, "*Pi*"("*lambda[T]*") const."),
                                 bquote("Si2017, "*Pi*"("*lambda[T]*") const."), bquote("Si2017, "*Pi*"("*lambda[T]*") varies"))) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.5,0.2))
p_scaled_growth_alt

p_scaled_ribo_frac_alt <- ggplot(data.frame(x = sv_array), aes(x)) + 
  stat_function(fun=function(x) log10(riboFrac_env(kN1_q10, kT_q10, kL_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN_q10_mean_Basan, kT_q10, kL_Basan_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="b"), size=1, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN_q10_mean_Si, kT_q10, kL_Si_q10, dp_q10, dl_q10, epsPar,10^x)), 
                aes(col="c"), size=1, n=plot_points, linetype="dashed") +
  stat_function(fun=function(x) log10(riboFrac_env(kN_q10_mean_alt, kT_q10, kL_q10_alt, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="d"), size=line_thick, n=plot_points) +
  scale_colour_manual(bquote("Source of "*Pi*"("*lambda[N]*")"), values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                                            "c"=colBlindScale[4],"d"=colBlindScale[5]), 
                      labels = c(bquote("Volkmer2011, "*Pi*"("*lambda[T]*") const."), bquote("Basan2015, "*Pi*"("*lambda[T]*") const."),
                                 bquote("Si2017, "*Pi*"("*lambda[T]*") const."), bquote("Si2017, "*Pi*"("*lambda[T]*") varies"))) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Ribosomal mass fraction, Log'[10]*'['*Phi[R]*']'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.0,0.0))
p_scaled_ribo_frac_alt

p_scaled_lipo_frac_alt <- ggplot(data.frame(x = sv_array), aes(x)) + 
  stat_function(fun=function(x) log10(lipoFrac_env(kN1_q10, kT_q10, kL_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN_q10_mean_Basan,kT_q10, kL_Basan_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="b"), size=1, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN_q10_mean_Si, kT_q10, kL_Si_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="c"), size=1, n=plot_points, linetype="dashed") +
  stat_function(fun=function(x) log10(lipoFrac_env(kN_q10_mean_alt, kT_q10, kL_q10_alt, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="d"), size=line_thick, n=plot_points) +
  scale_colour_manual(bquote("Source of "*Pi*"("*lambda[N]*")"), values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                                            "c"=colBlindScale[4],"d"=colBlindScale[5]), 
                      labels = c(bquote("Volkmer2011, "*Pi*"("*lambda[T]*") const."), bquote("Basan2015, "*Pi*"("*lambda[T]*") const."),
                                 bquote("Si2017, "*Pi*"("*lambda[T]*") const."), bquote("Si2017, "*Pi*"("*lambda[T]*") varies"))) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Envelope-producer mass fraction, Log'[10]*'['*Phi[L]*']'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.75,-1.0))
p_scaled_lipo_frac_alt

p_scaled_legend_alt <- ggplot(data.frame(x = sv_array), aes(x)) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN1_q10, kT_q10, kL_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN_q10_mean_Basan,kT_q10, kL_Basan_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="b"), size=1, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN_q10_mean_Si, kT_q10, kL_Si_q10, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="c"), size=1, n=plot_points, linetype="dashed") +
  stat_function(fun=function(x) log10(lipoFrac_env(kN_q10_mean_alt, kT_q10, kL_q10_alt, dp_q10, dl_q10, epsPar, 10^x)), 
                aes(col="d"), size=line_thick, n=plot_points) +
  scale_colour_manual(bquote("Source of "*Pi*"("*lambda[N]*")"), values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                                            "c"=colBlindScale[4],"d"=colBlindScale[5]), 
                      labels = c(bquote("Volkmer2011, "*Pi*"("*lambda[T]*") const."), bquote("Basan2015, "*Pi*"("*lambda[T]*") const."),
                                 bquote("Si2017, "*Pi*"("*lambda[T]*") const."), bquote("Si2017, "*Pi*"("*lambda[T]*") varies"))) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "left",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Envelope-producer mass fraction, Log'[10]*'['*Phi[L]*']'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.75,-1.0))
p_scaled_legend_alt
p_scaled_legend_alt <- cowplot::get_legend(p_scaled_legend_alt)

ggsave(file="p_scaled_growth_alt.pdf", plot=p_scaled_growth_alt, width = 5.2, height = 5.2)
ggsave(file="p_scaled_ribo_frac_alt.pdf", plot=p_scaled_ribo_frac_alt, width = 5.2, height = 5.2)
ggsave(file="p_scaled_lipo_frac_alt.pdf", plot=p_scaled_lipo_frac_alt, width = 5.2, height = 5.2)
ggsave(file="p_scaled_legend_alt.pdf", plot=p_scaled_legend_alt, width = 5.2, height = 5.2)


#========== Regression analysis ==========
eqn_label <- function(lm_summ){
  eqn <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
                   list(a = format(unname(lm_summ$coefficients[1]), digits = 2),
                        b = format(unname(lm_summ$coefficients[2]), digits = 2),
                        r2 = format(lm_summ$r.squared, digits = 3)))
  as.character(as.expression(eqn));
}

# First, obtain slopes
df_lambda$SV.mean <- df_lambda$S.mean/df_lambda$V.mean
model_sv_growth <- lm(log10(growth_rate)~log10(S.mean/V.mean), data = df_lambda)
model_v_growth <- lm(log10(growth_rate)~log10(V.mean), data = df_lambda)
model_ribo_sv <- lm(log10(ribo_frac)~log10(S.mean/V.mean), data = df_phi)
model_ribo_v <- lm(log10(ribo_frac)~log10(V.mean), data = df_phi)
model_lipo_sv <- lm(log10(lip_frac)~log10(S.mean/V.mean), data = df_phi)
model_lipo_v <- lm(log10(lip_frac)~log10(V.mean), data = df_phi)
# Take slopes and adjusted R^2 values
summary(model_sv_growth)
summary(model_v_growth)
summary(model_ribo_sv)
summary(model_ribo_v)
summary(model_lipo_sv)
summary(model_lipo_v)

lm_summ <- summary(model_sv_growth)
p_reg_growth_sv <- ggplot() + 
  stat_function(fun=function(x) lm_summ$coefficients[1]+lm_summ$coefficients[2]*x, aes(colour = "a"), size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black", aes(x=log10(df_lambda$S.mean/df_lambda$V.mean), y=log10(df_lambda$growth_rate))) +
  annotate(geom = "text", x = 1.2, y = 0.2, label = eqn_label(lm_summ), parse = TRUE, col = "red", size = 6) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V'['tot']*', Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual(values = c("a"=colBlindScale[1]), labels = medium_text) +
  guides(colour=guide_legend(override.aes=list(shape=c(25,18,1,0),linetype=0,stroke=1.5))) +
  ylim(c(-2.5,0.2))
p_reg_growth_sv
ggsave(file="p_reg_growth_sv.pdf", plot=p_reg_growth_sv, width = 6.2, height = 6.2)

lm_summ <- summary(model_ribo_sv)
p_reg_phiR_sv <- ggplot(df_phi, aes(x=log10(S.mean/V.mean), y=log10(ribo_frac))) + 
  stat_function(fun=function(x) lm_summ$coefficients[1]+lm_summ$coefficients[2]*x, aes(colour = "a"), size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  annotate(geom = "text", x = 1.2, y = 0.00, label = eqn_label(lm_summ), parse = TRUE, col = "red", size = 6) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Ribosomal mass fraction, Log'[10]*'['*Phi[L]*']'))) + xlab(bquote(bold('S/V'['tot']*', Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual(values = c("a"=colBlindScale[1]), labels = medium_text) +
  ylim(c(-2.0,0.0))
p_reg_phiR_sv
ggsave(file="p_reg_phiR_sv.pdf", plot=p_reg_phiR_sv, width = 6.2, height = 6.2)

lm_summ <- summary(model_lipo_sv)
p_reg_phiL_sv <- ggplot(df_phi, aes(x=log10(S.mean/V.mean), y=log10(lip_frac))) + 
  stat_function(fun=function(x) lm_summ$coefficients[1]+lm_summ$coefficients[2]*x, aes(colour = "a"), size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  annotate(geom = "text", x = 1.2, y = -1, label = eqn_label(lm_summ), parse = TRUE, col = "red", size = 6) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Envelope-producer mass fraction, Log'[10]*'['*Phi[L]*']'))) + xlab(bquote(bold('S/V'['tot']*', Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.75,-1.0)) +
  scale_colour_manual(values = c("a"=colBlindScale[1]), labels = medium_text)
p_reg_phiL_sv
ggsave(file="p_reg_phiL_sv.pdf", plot=p_reg_phiL_sv, width = 6.2, height = 6.2)


lm_summ <- summary(model_v_growth)
p_reg_growth_v <- ggplot() + 
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black", aes(x=log10(df_lambda$V.mean), y=log10(df_lambda$growth_rate))) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('V'['tot']*', Log'[10]*'['*V*' ('*mu*''*m^3*')]'))) +
  scale_colour_manual(values = c("a"=colBlindScale[1]), labels = medium_text) +
  guides(colour=guide_legend(override.aes=list(shape=c(25,18,1,0),linetype=0,stroke=1.5))) +
  ylim(c(-2.5,0.2))
p_reg_growth_v
ggsave(file="p_reg_growth_v.pdf", plot=p_reg_growth_v, width = 6.2, height = 6.2)

lm_summ <- summary(model_ribo_v)
p_reg_phiR_v <- ggplot(df_phi, aes(x=log10(V.mean), y=log10(ribo_frac))) +
  stat_function(fun=function(x) lm_summ$coefficients[1]+lm_summ$coefficients[2]*x, aes(colour = "a"), size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  annotate(geom = "text", x = -0.5, y = 0, label = eqn_label(lm_summ), parse = TRUE, col = "red", size = 6) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Ribosomal mass fraction, Log'[10]*'['*Phi[R]*']'))) + xlab(bquote(bold('V'['tot']*', Log'[10]*'['*V*' ('*mu*''*m^3*')]'))) +
  scale_colour_manual(values = c("a"=colBlindScale[1]), labels = medium_text) +
  ylim(c(-2.0,0.0))
p_reg_phiR_v
ggsave(file="p_reg_phiR_v.pdf", plot=p_reg_phiR_v, width = 6.2, height = 6.2)

lm_summ <- summary(model_lipo_v)
p_reg_phiL_v <- ggplot(df_phi, aes(x=log10(V.mean), y=log10(lip_frac))) + 
  stat_function(fun=function(x) lm_summ$coefficients[1]+lm_summ$coefficients[2]*x, aes(colour = "a"), size=line_thick, n=plot_points) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  annotate(geom = "text", x = -0.95, y = -1.05, label = eqn_label(lm_summ), parse = TRUE, col = "red", size=6) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Envelope-producer mass fraction, Log'[10]*'['*Phi[L]*']'))) + xlab(bquote(bold('V'['tot']*', Log'[10]*'['*V*' ('*mu*''*m^3*')]'))) +
  ylim(c(-2.75,-1.0)) +
  scale_colour_manual(values = c("a"=colBlindScale[1]), labels = medium_text)
p_reg_phiL_v
ggsave(file="p_reg_phiL_v.pdf", plot=p_reg_phiL_v, width = 6.2, height = 6.2)

# Next, test whether the slopes are significantly different from our theoretical expectations
model_sv_growth <- lm(log10(growth_rate)+offset(1*log10(S.mean/V.mean))~log10(S.mean/V.mean), data = df_lambda)
model_v_growth <- lm(log10(growth_rate)-offset(0.33*log10(V.mean))~log10(V.mean), data = df_lambda)
model_ribo_sv <- lm(log10(ribo_frac)+offset(1*log10(S.mean/V.mean))~log10(S.mean/V.mean), data=df_phi)
model_ribo_v <- lm(log10(ribo_frac)-offset(0.33*log10(V.mean))~log10(V.mean), data=df_phi)
model_lipo_sv <- lm(log10(lip_frac)-offset(0*log10(S.mean/V.mean))~log10(S.mean/V.mean), data=df_phi)
model_lipo_v <- lm(log10(lip_frac)-offset(0*log10(V.mean))~log10(V.mean), data=df_phi)
# Take p-values
summary(model_sv_growth)
summary(model_v_growth)
summary(model_ribo_sv)
summary(model_ribo_v)
summary(model_lipo_sv)
summary(model_lipo_v)


#========== Correction for possibly confounding factors (envelope thickness and natural habitat) ==========
xx <- df_lambda_full[duplicated(df_lambda_full$species),]$species
pos <- which(df_lambda_full$species %in% xx)
xx <- df_lambda_full[pos,]
xx_raw <- xx[xx$ID=="Generic thickness (30 nm)",]
xx_transf <- xx[xx$ID=="Specific thickness",]
df <- data.frame(x=log10(xx_raw$SV.mean),
                 y=log10(xx_raw$growth_rate),
                 xend=log10(xx_transf$SV.mean),
                 yend=log10(xx_transf$growth_rate))

# Regression using the species-specific thicknesses
# First, obtain slopes
model_sv_growth <- lm(log10(growth_rate)~log10(S.mean/V.mean), data = xx[xx$ID=="Specific thickness",])
model_v_growth <- lm(log10(growth_rate)~log10(V.mean), data = xx[xx$ID=="Specific thickness",])
plot_model_specific <- model_sv_growth
# Take slopes and adjusted R^2 values
summary(model_sv_growth)
summary(model_v_growth)

# Next, test whether the slopes are significantly different from our theoretical expectations
model_sv_growth <- lm(log10(growth_rate)+offset(1*log10(S.mean/V.mean))~log10(S.mean/V.mean), data = xx[xx$ID=="Specific thickness",])
model_v_growth <- lm(log10(growth_rate)-offset(0.33*log10(V.mean))~log10(V.mean), data = xx[xx$ID=="Specific thickness",])
# Take p-values
summary(model_sv_growth)
summary(model_v_growth)

# Regression using the generic thickness of 30 nm
# First, obtain slopes
model_sv_growth <- lm(log10(growth_rate)~log10(S.mean/V.mean), data = xx[xx$ID=="Generic thickness (30 nm)",])
model_v_growth <- lm(log10(growth_rate)~log10(V.mean), data = xx[xx$ID=="Generic thickness (30 nm)",])
plot_model_generic <- model_sv_growth
# Take slopes and adjusted R^2 values
summary(model_sv_growth)
summary(model_v_growth)

# Next, test whether the slopes are significantly different from our theoretical expectations
model_sv_growth <- lm(log10(growth_rate)+offset(1*log10(S.mean/V.mean))~log10(S.mean/V.mean), data = xx[xx$ID=="Generic thickness (30 nm)",])
model_v_growth <- lm(log10(growth_rate)-offset(0.33*log10(V.mean))~log10(V.mean), data = xx[xx$ID=="Generic thickness (30 nm)",])
# Take p-values
summary(model_sv_growth)
summary(model_v_growth)


env_corr_growth <- ggplot(df_lambda_env_corr, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(species))) + 
  stat_function(fun=function(x) plot_model_generic$coefficients[1]+plot_model_generic$coefficients[2]*x, 
                colour = "black", size=line_thick, linetype="dashed") +
  stat_function(fun=function(x) plot_model_specific$coefficients[1]+plot_model_specific$coefficients[2]*x, 
                colour = "red", size=line_thick, linetype="solid") +
  geom_point(size = 4.5) +
  geom_point(data=xx[xx$ID=="Generic thickness (30 nm)",], size = 4.5, shape=1, stroke=1.5) +
  geom_segment(aes(x=df$x,y=df$y, xend=df$xend, yend=df$yend)) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.5,0.2))
env_corr_growth
ggsave(file="env_corr_growth.pdf", plot=env_corr_growth, width = 5.2, height = 5.2)


legend_growth <- ggplot(df_lambda_env_corr, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(species))) + 
  geom_point(size = 4.5) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]')))
legend_growth
legend_growth <- cowplot::get_legend(legend_growth)
ggsave(file="legend_growth_env_corr.pdf", plot=legend_growth, width = 7.2, height = 7.2)

filter_shape <- function(shape_df, shape_str) {
  if(shape_str == "sphere") {
    return(shape_df$cell_shape=="coccus")
  } else {
    return(with(shape_df, (cell_shape=="bacillus")|(cell_shape=="coccobacillus")|(cell_shape=="vibrio")))
  }
}

L_env <- 0.03 # Remember to change this, if you want the total volume
in_file <- "madin2020_growth.csv"
df_madin <- read.csv(in_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
shape_sphere <- filter_shape(df_madin, "sphere")
shape_capsule <- filter_shape(df_madin, "capsule")
df_sub <- select(df_madin, superkingdom, phylum, gram_stain, species, isolation_source, metabolism, doubling_h, doubling_h_norm, optimum_tmp, growth_tmp, cell_shape, d1_lo, d1_up, d2_lo, d2_up, data_source)
df_sub <- df_sub[complete.cases(df_sub)&(shape_capsule|shape_sphere),]
df_sub$d1_mean <- rowMeans(data.frame(cbind(df_sub$d1_lo,df_sub$d1_up)))-2*L_env
df_sub$d2_mean <- rowMeans(data.frame(cbind(df_sub$d2_lo,df_sub$d2_up)))-2*L_env
shape_sphere <- filter_shape(df_sub, "sphere")
shape_capsule <- filter_shape(df_sub, "capsule")

df_sub$V.mean <- ifelse(shape_sphere, sphere_vol(df_sub$d1_mean),
                          capsule_vol(df_sub$d1_mean,df_sub$d2_mean))
df_sub$S.mean <- ifelse(shape_sphere, sphere_surf(df_sub$d1_mean),
                           capsule_surf(df_sub$d1_mean,df_sub$d2_mean))

data_select <- df_sub[(df_sub$superkingdom == "Bacteria"),]
data_select <- data_select[!grepl("^Cyano", data_select$species),] # Remove cyanobacteria because they are not heterotrophs
data_select$growth_rate <- log(2)/data_select$doubling_h_norm

data_select$environment <- "Host-associated"
data_select[!grepl("^host", data_select$isolation_source),]$environment <- "Free-living"

#+offset(1*log10(S.mean/V.mean)) -- this part should be added as offset to test whether the slope is significantly different from the expectation
lm_madinSV <- lm(log10(growth_rate)~log10(S.mean/V.mean), data=data_select)
lm_madinV <- lm(log10(growth_rate)~log10(V.mean), data=data_select)
summary(lm_madinSV)
summary(lm_madinV)

model_sv_Madin <- lm(log10(growth_rate)+offset(1*log10(S.mean/V.mean))~log10(S.mean/V.mean), data = data_select)
model_v_Madin <- lm(log10(growth_rate)-offset(0.33*log10(V.mean))~log10(V.mean), data = data_select)
# Take p-values
summary(model_sv_Madin)
summary(model_v_Madin)

lm_madin_free <- lm(log10(growth_rate)~log10(S.mean/V.mean), data=data_select[data_select$environment!="Host-associated",])
lm_madin_host <- lm(log10(growth_rate)~log10(S.mean/V.mean), data=data_select[data_select$environment=="Host-associated",])
summary(lm_madin_free)
summary(lm_madin_host)

model_free_Madin <- lm(log10(growth_rate)+offset(1*log10(S.mean/V.mean))~log10(S.mean/V.mean), data=data_select[data_select$environment!="Host-associated",])
model_host_Madin <- lm(log10(growth_rate)+offset(1*log10(S.mean/V.mean))~log10(V.mean), data=data_select[data_select$environment=="Host-associated",])
# Take p-values
summary(model_free_Madin)
summary(model_host_Madin)

data_select_1 <- data_select[data_select$environment=="Free-living",]
data_select_2 <- data_select[data_select$environment=="Host-associated",]
madinPlot <- ggplot() + 
  stat_function(fun=function(x) lm_madin_free$coefficients[1]+lm_madin_free$coefficients[2]*x, 
                aes(colour = "a"), size=line_thick*0.75, linetype="dashed") +
  stat_function(fun=function(x) lm_madin_host$coefficients[1]+lm_madin_host$coefficients[2]*x, 
                aes(colour = "b"), size=line_thick*0.75, linetype="dashed") +
  stat_function(fun=function(x) lm_madinSV$coefficients[1]+lm_madinSV$coefficients[2]*x, 
                aes(colour = "c"), size=line_thick*0.75) +
  geom_point(size=4.5, aes(x=log10(data_select_1$S.mean/data_select_1$V.mean), y=log10(data_select_1$growth_rate), colour = "a")) +
  geom_point(size=4.5, aes(x=log10(data_select_2$S.mean/data_select_2$V.mean), y=log10(data_select_2$growth_rate), colour = "b")) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual("Environment", values = c("a"=colBlindScale[1],"b"=colBlindScale[2],"c"=colBlindScale[3]), 
                      labels = c("Free-living","Host-associated","All"))
madin_scaling <- grid.arrange(madinPlot, nrow = 1)
ggsave(file="madin_scaling.pdf", plot=madin_scaling, width = 6.5, height = 5)


#df_lambda_env_corr[df_lambda_env_corr$total_thickness_nm<10^1.25,]$species
lm_thick <- lm(log10(total_thickness_nm)~log10(V.mean), data=df_lambda_env_corr)
summary(lm_thick)

thickness_plot <- ggplot(df_lambda_env_corr, aes(x=log10(V.mean), y=log10(total_thickness_nm), color=(species))) + 
  geom_point(size = 4.5) +
  stat_function(fun=function(x) lm_thick$coefficients[1]+lm_thick$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Cell envelope thickness, Log'[10]*'['*h*' ('*mu*''*m*')]'))) + xlab(bquote(bold('Internal cell volume, Log'[10]*'['*V*' ('*mu*''*m^3*')]')))
thickness_scaling <- grid.arrange(thickness_plot, nrow = 1)
ggsave(file="thickness_scaling.pdf", plot=thickness_plot, width = 5, height = 5)


#========== Correlation between S:V and genomic properties ==========
lm_trna_sv <- lm(log10(tRNA.genes)~log10(S.mean/V.mean), data=df_lambda[!is.na(df_lambda$shape),])
lm_trna_v <- lm(log10(tRNA.genes)~log10(V.mean), data=df_lambda[!is.na(df_lambda$shape),])
lm_trna_lambda <- lm(log10(tRNA.genes)~log10(growth_rate), data=df_lambda[!is.na(df_lambda$shape),])
summary(lm_trna_sv)
summary(lm_trna_v)
summary(lm_trna_lambda)

lm_rrna_sv <- lm(log10(rRNA.genes)~log10(S.mean/V.mean), data=df_lambda[!is.na(df_lambda$shape),])
lm_rrna_v <- lm(log10(rRNA.genes)~log10(V.mean), data=df_lambda[!is.na(df_lambda$shape),])
lm_rrna_lambda <- lm(log10(rRNA.genes)~log10(growth_rate), data=df_lambda[!is.na(df_lambda$shape),])
summary(lm_rrna_sv)
summary(lm_rrna_v)
summary(lm_rrna_lambda)

p_rrna_sv <- ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(S.mean/V.mean), y=log10(rRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_rrna_sv$coefficients[1]+lm_rrna_sv$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(rRNA gene number)'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]')))

p_trna_sv <-ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(S.mean/V.mean), y=log10(tRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_trna_sv$coefficients[1]+lm_trna_sv$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(tRNA gene number)'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]')))

ggsave(file="genomic_trna_sv.pdf", plot=p_trna_sv, width = 6, height = 6)
ggsave(file="genomic_rrna_sv.pdf", plot=p_rrna_sv, width = 6, height = 6)

p_rrna_lambda <- ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(growth_rate), y=log10(rRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_rrna_lambda$coefficients[1]+lm_rrna_lambda$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(rRNA gene number)'))) + xlab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]')))

p_trna_lambda <- ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(growth_rate), y=log10(tRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_trna_lambda$coefficients[1]+lm_trna_lambda$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(tRNA gene number)'))) + xlab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]')))

ggsave(file="genomic_trna_lambda.pdf", plot=p_trna_lambda, width = 6, height = 6)
ggsave(file="genomic_rrna_lambda.pdf", plot=p_rrna_lambda, width = 6, height = 6)

p_rrna_v <- ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(V.mean), y=log10(rRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_rrna_v$coefficients[1]+lm_rrna_v$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(rRNA gene number)'))) + xlab(bquote(bold('Volume, Log'[10]*'[V ('*mu*''*m^3*')]')))

p_trna_v <-ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(V.mean), y=log10(tRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_trna_v$coefficients[1]+lm_trna_v$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(tRNA gene number)'))) + xlab(bquote(bold('Volume, Log'[10]*'[V ('*mu*''*m^3*')]')))

ggsave(file="genomic_trna_v.pdf", plot=p_trna_v, width = 6, height = 6)
ggsave(file="genomic_rrna_v.pdf", plot=p_rrna_v, width = 6, height = 6)


#========== Sensitivity analysis of protein degradation rates ==========
dp_q10_vec <- q10_corrected(
                c(log(2)/70,
                log(2)/mean(c(25,70)),
                log(2)/mean(c(5,25)),
                log(2)/mean(c(2,5)),
                0.4*log(2)/70+0.35*log(2)/mean(c(25,70))+0.1*log(2)/mean(c(5,25))+0.05*log(2)/mean(c(2,5))),
                2.5,37)
dp_q10_vec
deg_rates_text <- c(bquote("70 hours, "*d[p]*"=0.002"), 
                 bquote("25-70 hours, "*d[p]*"=0.003"),
                 bquote("5-25 hours, "*d[p]*"=0.009"),
                 bquote("2-5 hours, "*d[p]*"=0.04"),
                 bquote("Fraction weighted mean, "*d[p]*"=0.005"))

p_sensitivity_analysis_dp <- ggplot(df_lambda, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(ID))) + 
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10_vec[1], dl_q10)), 
                aes(colour = "a"), size=0.3*line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10_vec[2], dl_q10)), 
                aes(colour = "b"), size=0.3*line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10_vec[3], dl_q10)), 
                aes(colour = "c"), size=0.3*line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10_vec[4], dl_q10)), 
                aes(colour = "d"), size=0.3*line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10_vec[5], dl_q10)), 
                aes(colour = "e"), size=0.3*line_thick, n=plot_points) +
  scale_colour_manual("Protein half-life:", values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                      "c"=colBlindScale[4],"d"=colBlindScale[5],
                                                      "e"=colBlindScale[6]), labels = deg_rates_text) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.5,0.2))
p_sensitivity_analysis_dp

sensitivity_legend <- ggplot(df_lambda, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(ID))) + 
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10_vec[1], dl_q10)), 
                aes(colour = "a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10_vec[2], dl_q10)), 
                aes(colour = "b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10_vec[3], dl_q10)), 
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10_vec[4], dl_q10)), 
                aes(colour = "d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10_vec[5], dl_q10)), 
                aes(colour = "e"), size=line_thick, n=plot_points) +
  scale_colour_manual("Protein half-life:", values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                      "c"=colBlindScale[4],"d"=colBlindScale[5],
                                                      "e"=colBlindScale[6]), labels = deg_rates_text) +
  geom_point(size = 2.0, shape=1, stroke=1.5, alpha=0.75, col="black") +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "right",
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.5,0.2))
sensitivity_legend
sensitivity_legend <- cowplot::get_legend(sensitivity_legend)

ggsave(file="p_sensitivity_analysis_dp.pdf", plot=p_sensitivity_analysis_dp, width = 5.2, height = 5.2)
ggsave(file="sensitivity_legend.pdf", plot=sensitivity_legend, width = 5.2, height = 5.2)


#========== The scaling of metabolic rate in LTEE ==========
# Surface area taken to be 4.32 um^2
lipPerSurf <- (22*10^6/3)/4.32
pgnPerSurf <- (3.5*10^6)/4.32
lpsPerSurf <- 1.43*10^6/4.32
aaPerSurf <- 78*10^6/4.32
lipPerUnit <- lipPerSurf/lipPerSurf

# 4 ATP for direct costs of elongating a polypeptide (from Lynch and Marinov 2015)
# The parenthesis contains the number of NADH, NADPH (from Figure 2D in Mahmoudabadi2020). Each is energetically equivalent to 2 ATPs.
aaDirCost <- 2+2*mean(c(1,3,1,1,4,1,1,0,-2,5,1,4,8,2,3,0,3,1,1,2))
lppPerAACost <- 4+aaDirCost

pgnPerUnit <- pgnPerSurf/lipPerSurf
lpsPerUnit <- lpsPerSurf/lipPerSurf
aaPerUnit <- aaPerSurf/lipPerSurf

qlPar <- 74*3*lipPerUnit + 17*pgnPerUnit + 301*lpsPerUnit + lppPerAACost*aaPerUnit
ls <- ((27 + 4)*aaPerUnit + 223*pgnPerUnit + 232*3*lipPerUnit + 2243*lpsPerUnit)/27

mean_vol_Volkmer <- mean(capsule_vol(df_volkmer$cell_width.mean,df_volkmer$cell_length.mean))
gammaPar <- (mean_vol_Volkmer*4.4*10^6 - 30000)*325 + (30000*7336)
lp <- 325; lr <- 7336;
Qaa <- aaDirCost; Qt <- 4; Ql <- qlPar; Qdp <- 1*325; Qdl <- 0;

mr_eq <- function(kN, kT, kL, epsPar, svPar, dp, dl) {
  return(
    (8*gammaPar*pi^3*(kL*kT*ls*(dp*Qdp + kN*lp*(Qaa + Qt)) + epsPar*(kL*kN*kT*lp*(ls*Qaa + Ql) + dp*kN*kT*ls*(Qdp + lp*(Qaa + Qt)) + 
    dl*kN*ls*(dp*Qdp - kT*lp*(Qaa + Qt)) + dp*kL*kT*(-lp*Ql + ls*(Qdp + lp*Qt)) + dl*kL*(dp*ls*Qdp + kN*lp*(ls*Qaa + Qdl + Ql) + 
    kT*lp*(Qdl + Ql - ls*Qt)))*svPar + dl*(epsPar^2)*(kL + kN)*(dp + kT)*lp*Qdl*svPar^2))/(lp*ls*(svPar^3)*(kL*(kN + kT) + epsPar*(kL + kN)*(dp + kT)*svPar))
  )
}

df_mr <- read.csv("/home/bogi/Desktop/growth_rate/data/cell_traits_data/metabolic_rate/Respiration_scaling.csv")
df_mr <- df_mr[df_mr$OPTICAL.DENSITY==30,]
df_mr$LOG10RESP <- 10^df_mr$LOG10RESP
df_mr$LOG10VOLUME <- 10^df_mr$LOG10VOLUME
df_mr$LOG10RESP <- 60*(6.022*10^23*10^(-6))*((28/6)/0.512)*df_mr$LOG10RESP*(2.5^((20-37)/10))
df_mr$SV.mean <- 2*pi*(df_mr$LOG10VOLUME)^(-1/3)

p_mr_ltee <- ggplot() + 
  stat_function(fun=function(x) log10(mr_eq(kN1_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(mr_eq(kN2_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(mr_eq(kN3_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(mr_eq(kN4_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(mr_eq(kN5_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "e"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(mr_eq(kN6_q10, kT_q10, kL_q10_alt, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "f"), size=line_thick, n=plot_points) +
  geom_point(size = 4.0, shape=1, stroke=1.5, alpha=0.75, col="red", aes(x=log10(df_mr$SV.mean), y=log10(df_mr$LOG10RESP))) +
  theme(aspect.ratio = 1,
        axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = 'grey75'),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Metabolic rate, Log'[10]*'[dQ/dt (ATP cell'^-1*'h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  scale_colour_manual("Nutrient capacity", values = c("a"=colBlindScale[1], "b"=colBlindScale[2], "c"=colBlindScale[3],
                                                      "d"=colBlindScale[4], "e"=colBlindScale[5], "f"=colBlindScale[6], "x"="red"), labels = medium_text) +
  guides(colour=guide_legend(override.aes=list(shape=c(25,18,1,0),linetype=0,stroke=1.5)))
p_mr_ltee

ggsave(file="p_mr_ltee.pdf", plot=p_mr_ltee, width = 5.2, height = 5.2)
