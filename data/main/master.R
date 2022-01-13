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
colBlindScale <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7", "#DB6D00", "#EF3B2C")
colBlindScaleExtended <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7", "#DB6D00", "#EF3B2C",
                   "#67001F","#F7FCFD","#CB181D","#78C679","#F46D43","#A6CEE3","#FD8D3C","#A6D854","#D4B9DA","#6A51A3",
                   "#7F0000","#D9D9D9","#FFF7BC","#F0F0F0","#C7EAE5","#003C30","#F16913","#FFF7FB","#8C6BB1","#C7E9B4",
                   "#762A83","#FC9272","#AE017E","#F7F7F7","#DF65B0","#74C476")
scaleFUN <- function(x) sprintf("%.2f", x)


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


#========== Inferring capacities from OLS slopes ==========
SAVpar <- 5 # Rough average in our dataset; More precisely, it would be 4.94 but this might change if more data is added
epsPar <- (0.30/0.55)/SAVpar
dp_est <- 0.05 # page 2896 in Nath and Koch 1970 (taken for slowly decaying component)
dl_est <- 1.00

kappaT_est <- (1-dp_est*olsR$coefficients[2])/olsR$coefficients[2]
kappaL_est <- (epsPar*SAVpar)/olsL$coefficients[2]
kappaN_est <- function(x) {
  (kappaL_est*(x*dp_est-epsPar*SAVpar-1))/(x*kappaL_est+epsPar*SAVpar)
}

# Calculating standard error of the estimate by error propagation
q10_factor <- 2.5^((20-37)/10)
kT_std <- sqrt((((sqrt(20)*olsR$coefficients[4])^2))/(olsR$coefficients[2]^4))
kL_std <- sqrt((((epsPar*SAVpar)^2)*((sqrt(20)*olsL$coefficients[4])^2))/(olsL$coefficients[2]^4))
kT_q10_std <- sqrt((((sqrt(20)*olsR$coefficients[4])^2)*(q10_factor^2))/(olsR$coefficients[2]^4))
kL_q10_std <- sqrt((((epsPar*SAVpar)^2)*((sqrt(20)*olsL$coefficients[4])^2)*(q10_factor^2))/(olsL$coefficients[2]^4))
kT_q10_err <- kT_q10_std/sqrt(20)
kL_q10_err <- kL_q10_std/sqrt(20)

kN_std_q10_est <- function(x,n) {
  a <- ((q10_factor^2)*(epsPar*SAVpar)^2)/(epsPar*SAVpar+kappaL_est*x$coefficients[2])^4
  b <- (x$coefficients[2]^2 * dp_est^2 * kL_std^2 + kappaL_est^2 * (sqrt(n)*olsR_ab1$coefficients[4])^2 * (dp_est + kappaL_est)^2
        - 2*x$coefficients[2] * dp_est * kL_std^2 * epsPar * SAVpar + kL_std^2*(epsPar * SAVpar)^2)
  return(sqrt(a*b)/sqrt(n))
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
axes_style <- element_text(size=15, face="bold", colour = "black");

p_ribo_nutri <- ggplot(df_eco_nutri, aes(x=growth_rates, y=phi_R, col=source)) + 
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, colour = "black", size = 1.5) +
  geom_point(size=3.5) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
        scale_y_continuous(labels=scaleFUN) +
        scale_fill_manual(values=colBlindScale) +
        scale_colour_manual(values=colBlindScale) +
        xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R])))
p_ribo_nutri

p_lipo_nutri <- ggplot(df_eco_nutri, aes(x=growth_rates, y=phi_L, col=source)) + 
  stat_function(fun=function(x) olsL$coefficients[1]+olsL$coefficients[2]*x, colour = "black", size = 1.5) +
  geom_point(size=3.5) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_fill_manual(values=colBlindScale) +
  scale_colour_manual(values=colBlindScale) +
  scale_y_continuous(labels=scaleFUN) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Envelope-producer mass fraction, '*Phi[L])))
p_lipo_nutri

p_ribo_ab <- ggplot(df_eco_ab_ribo, aes(x=growth_rates, y=phi_R, col=growth_medium)) + 
  geom_smooth(method='lm', formula= y~x, se = FALSE, size = 1.5) +
  geom_point(size=3.5) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_continuous(labels=scaleFUN) +
  xlim(c(0.0,2.0)) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R])))
p_ribo_ab

p_ribo_env_ab <- ggplot(df_eco_env_ribo, aes(x=growth_rates, y=phi_R, col=perturbation)) + 
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, colour = "black", size = 1.5) +
  geom_point(size=3.5) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_fill_manual(values=colBlindScale) +
  scale_colour_manual(values=colBlindScale) +
  scale_y_continuous(labels=scaleFUN) +
  xlab(bquote(bold('Growth rate, '*lambda*' (h'^-1*')'))) + ylab(bquote(bold('Ribosomal mass fraction, '*Phi[R])))
p_ribo_env_ab

ribo_legend1 <- ggplot(df_eco_nutri_all, aes(x=growth_rates, y=phi_R, col=source)) + 
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, colour = "black", size = 1.5) +
  geom_point(size=3.5) +
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_y_continuous(labels=scaleFUN) +
  scale_fill_manual(values=colBlindScale) +
  scale_colour_manual(values=colBlindScale) +
  xlab(bquote('Growth rate, '*lambda*' (h'^-1*')')) + ylab(bquote('Ribosomal mass fraction, '*phi[R]))
legend1 <- cowplot::get_legend(ribo_legend1)
ggsave(file="legend1.pdf", legend1, width = 8, height = 8)

ribo_legend2 <- ggplot(df_eco_ab_ribo, aes(x=growth_rates, y=phi_R, col=growth_medium)) + 
  geom_smooth(method='lm', formula= y~x, se = FALSE, size = 1.5) +
  geom_point(size=3.5) +
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_continuous(labels=scaleFUN) +
  xlim(c(0.0,2.0)) +
  xlab(bquote('Growth rate, '*lambda*' (h'^-1*')')) + ylab(bquote('Ribosomal mass fraction, '*phi[R]))
legend2 <- cowplot::get_legend(ribo_legend2)
ggsave(file="legend2.pdf", legend2, width = 8, height = 8)

ribo_legend3 <- ggplot(df_eco_env_ribo, aes(x=growth_rates, y=phi_R, col=perturbation)) + 
  geom_point(size=3.5) +
  stat_function(fun=function(x) olsR$coefficients[1]+olsR$coefficients[2]*x, colour = "black", size = 1.5) +
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  scale_fill_manual(values=colBlindScale) +
  scale_colour_manual(values=colBlindScale) +
  scale_y_continuous(labels=scaleFUN) +
  xlab(bquote('Growth rate, '*lambda*' (h'^-1*')')) + ylab(bquote('Ribosomal mass fraction, '*Phi[R]))
legend3 <- cowplot::get_legend(ribo_legend3)
ggsave(file="legend3.pdf", plot=legend3, width = 8, height = 8)

#p_proteomic_response <- grid.arrange(p_ribo_nutri, p_ribo_ab, p_lipo_nutri, p_ribo_env_ab, nrow = 2)
#ggsave(file="proteome_response.pdf", plot=p_proteomic_response, width = 8, height = 8)
ggsave(file="p_ribo_nutri.pdf", plot=p_ribo_nutri, width = 5, height = 5)
ggsave(file="p_ribo_ab.pdf", plot=p_ribo_ab, width = 5, height = 5)
ggsave(file="p_lipo_nutri.pdf", plot=p_lipo_nutri, width = 5, height = 5)
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
                                df_lambda_env_corr$growth_rate*lambda_anaero_correct(kN6_q10,kT_q10,kL_q10,epsPar,
                                                                                     df_lambda_env_corr$S.mean/df_lambda_env_corr$V.mean,dp_q10,dl_q10),
                                df_lambda_env_corr$growth_rate)

df_lambda <- read.csv("/home/bogi/Desktop/growth_rate/data/main/growth_scaling.shape.genome.csv",
                      skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
df_lambda[is.na(df_lambda$tca_steps),]$tca_steps <- 12
df_lambda$growth_rate <- df_lambda$growth_rate/24
df_lambda$growth_rate <- ifelse(df_lambda$tca_steps<=anaerobic_cutoff,
                                df_lambda$growth_rate*lambda_anaero_correct(kN6_q10,kT_q10,kL_q10,epsPar,df_lambda$S.mean/df_lambda$V.mean,dp_q10,dl_q10),
                                df_lambda$growth_rate)

#df_lambda <- df_lambda %>% group_by(species) %>% slice(which.max(growth_rate))
#df_lambda_env_corr <- df_lambda_env_corr %>% group_by(species) %>% slice(which.max(growth_rate))
df_lambda$ID <- "Generic thickness (30 nm)"
df_lambda_env_corr$ID <- "Specific thickness"
df_lambda_full <- rbind(df_lambda, df_lambda_env_corr[!is.na(df_lambda_env_corr$total_thickness_nm),])
df_lambda_full$SV.mean <- df_lambda_full$S.mean/df_lambda_full$V.mean

# Add Grant2021 data
df_Grant <- read.csv("/home/bogi/Desktop/growth_rate/data/cell_size_data/ltee_data/grant2021.csv")
df_Grant$S.mean <- capsule_surf(df_Grant$mean.width-2*L_env, df_Grant$mean.length-2*L_env)
df_Grant$V.mean <- capsule_vol(df_Grant$mean.width-2*L_env, df_Grant$mean.length-2*L_env)
df_Grant$SV.mean <- df_Grant$S.mean/df_Grant$V.mean
df_Grant$fitnessMP <- df_Grant$fitnessMP*0.7726 # Value from Table 1 in Vasi1994
colnames(df_Grant)[2] <- "growth_rate"
df_Grant$growth_rate <- q10_corrected(df_Grant$growth_rate,2.5,37)

# Add Lennon2021 data
df_Lennon <- read.csv("/home/bogi/Desktop/growth_rate/data/cell_size_data/ltee_data/lennon2021.csv")
df_Lennon$fitnessMP <- 1.5926*df_Lennon$fitnessMP # This value is from our dataset
colnames(df_Lennon)[2] <- "growth_rate"
df_Lennon$S.mean <- sphere_surf(df_Lennon$mean.width-2*0.01)
df_Lennon$V.mean <- sphere_vol(df_Lennon$mean.width-2*0.01)
df_Lennon$SV.mean <- df_Lennon$S.mean/df_Lennon$V.mean
df_Lennon$growth_rate <- q10_corrected(df_Lennon$growth_rate,2.5,37)

# Add Gallet2017 data -- Note that this is sligthly different because we are not subtracting L_env and we have S/V_{tot} not S/v_{cyt}
df_Gallet <- read.csv("/home/bogi/Desktop/growth_rate/data/cell_size_data/ltee_data/gallet2017.csv")
df_Gallet$S.mean <- gammaShapeConst*df_Gallet$V.mean^(2/3)
df_Gallet$SV.mean <- df_Gallet$S.mean/df_Gallet$V.mean
df_Gallet$growth_rate <- q10_corrected(df_Gallet$growth_rate,2.5,37)

medium_text <- c(bquote(kappa[n]*"(M63+gly)"), bquote(kappa[n]*"(M63+glc)"), bquote(kappa[n]*"(cAA+gly)"),
                 bquote(kappa[n]*"(cAA+glc)"), bquote(kappa[n]*"(RDM+gly)"), bquote(kappa[n]*"(RDM+glc)"))
plot_points <- 5000

p_scaled_growth <- ggplot(df_lambda, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(ID))) + 
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "e"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "f"), size=line_thick, n=plot_points) +
  scale_colour_manual("Nutrient capacity", values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                      "c"=colBlindScale[4],"d"=colBlindScale[5],
                                                      "e"=colBlindScale[6],"f"=colBlindScale[7]), labels = medium_text) +
  geom_point(col="black", size = point_size1) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.5,0.2))
p_scaled_growth

nutrient_legend <- ggplot(df_lambda, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(ID))) + 
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "a"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "b"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "c"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "d"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "e"), size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                aes(colour = "f"), size=line_thick, n=plot_points) +
  scale_colour_manual("Nutrient capacity", values = c("a"=colBlindScale[2],"b"=colBlindScale[3],
                                                      "c"=colBlindScale[4],"d"=colBlindScale[5],
                                                      "e"=colBlindScale[6],"f"=colBlindScale[7]), labels = medium_text) +
  geom_point(col="black", size = point_size1) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]')))+
  ylim(c(-2.5,0.2))
nutrient_legend <- cowplot::get_legend(nutrient_legend)

xx <- dplyr::bind_rows(df_Grant, df_Gallet, df_Lennon)
yy <- aggregate(xx, by=list(xx$generation, xx$ID), mean, na.rm = TRUE)
colnames(yy)[2] <- "species" 
grantPlot <- ggplot(yy, aes(x=log10(SV.mean), y=log10(growth_rate), col=species)) + 
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[2], size=0.75*line_thick) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[3], size=0.75*line_thick) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[4], size=0.75*line_thick) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[5], size=0.75*line_thick) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[6], size=0.75*line_thick) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[7], size=0.75*line_thick) +
  geom_point(size = 4) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-1,-0.1))
grantPlot

df_phi <- read.csv("proteome_scaling.csv", skip = 0, header = TRUE, sep = ",", stringsAsFactors = FALSE)
p_scaled_ribo_frac <- ggplot(df_phi, aes(x=log10(S.mean/V.mean), y=log10(ribo_frac))) + 
  stat_function(fun=function(x) log10(riboFrac_env(kN1_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[2], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN2_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[3], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN3_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[4], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN4_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[5], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN5_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[6], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(riboFrac_env(kN6_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[7], size=line_thick, n=plot_points) +
  geom_point(color = 'black', size = point_size1) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Ribosomal mass fraction, Log'[10]*'['*Phi[R]*']'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.0,0.0))
p_scaled_ribo_frac

p_scaled_lipo_frac <- ggplot(df_phi, aes(x=log10(S.mean/V.mean), y=log10(lip_frac))) + 
  stat_function(fun=function(x) log10(lipoFrac_env(kN1_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[2], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN2_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[3], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN3_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[4], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN4_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[5], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN5_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[6], size=line_thick, n=plot_points) +
  stat_function(fun=function(x) log10(lipoFrac_env(kN6_q10,kT_q10,kL_q10,dp_q10,dl_q10,epsPar,10^x)), 
                colour = colBlindScale[7], size=line_thick, n=plot_points) +
  geom_point(color = 'black', size = point_size1) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Envelope-producer mass fraction, Log'[10]*'['*Phi[L]*']'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.75,-1.0))
p_scaled_lipo_frac

#growth_proteome_scaling <- grid.arrange(p_scaled_growth, legend_growth, p_scaled_ribo_frac, p_scaled_lipo_frac, nrow = 2)
#ggsave(file="growth_proteome_scaling.pdf", plot=growth_proteome_scaling, width = 12, height = 12)
ggsave(file="p_scaled_growth.pdf", plot=p_scaled_growth, width = 5.2, height = 5.2)
ggsave(file="p_scaled_ribo_frac.pdf", plot=p_scaled_ribo_frac, width = 5.2, height = 5.2)
ggsave(file="p_scaled_lipo_frac.pdf", plot=p_scaled_lipo_frac, width = 5.2, height = 5.2)
ggsave(file="legend_nutrients.pdf", plot=nutrient_legend, width = 7.2, height = 7.2)
ggsave(file="grant_growth_scaling.pdf", plot=grantPlot, width = 9.25, height = 6)


#========== Regression analysis ==========
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
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[2], size=1.5, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[3], size=1.5, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[4], size=1.5, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[5], size=1.5, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[6], size=1.5, n=plot_points) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[7], size=1.5, n=plot_points) +
  stat_function(fun=function(x) plot_model_generic$coefficients[1]+plot_model_generic$coefficients[2]*x, 
                colour = "black", size=line_thick*0.50, linetype="dashed") +
  stat_function(fun=function(x) plot_model_specific$coefficients[1]+plot_model_specific$coefficients[2]*x, 
                colour = "red", size=line_thick*0.50, linetype="solid") +
  geom_point(size = 3) +
  geom_point(data=xx[xx$ID=="Generic thickness (30 nm)",], size = 3, shape=1) +
  geom_segment(aes(x=df$x,y=df$y, xend=df$xend, yend=df$yend)) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]'))) +
  ylim(c(-2.5,0.2))
env_corr_growth
ggsave(file="env_corr_growth.pdf", plot=env_corr_growth, width = 5.2, height = 5.2)


legend_growth <- ggplot(df_lambda_env_corr, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(species))) + 
  stat_function(fun=function(x) log10(lambda_env(kN1_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[2], size=1.5) +
  stat_function(fun=function(x) log10(lambda_env(kN2_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[3], size=1.5) +
  stat_function(fun=function(x) log10(lambda_env(kN3_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[4], size=1.5) +
  stat_function(fun=function(x) log10(lambda_env(kN4_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[5], size=1.5) +
  stat_function(fun=function(x) log10(lambda_env(kN5_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[6], size=1.5) +
  stat_function(fun=function(x) log10(lambda_env(kN6_q10, kT_q10, kL_q10, epsPar, 10^x, dp_q10, dl_q10)), 
                colour = colBlindScale[7], size=1.5) +
  geom_point(size = point_size2) +
  theme(axis.text.y   = axes_style,
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
ggsave(file="legend_growth.pdf", plot=legend_growth, width = 7.2, height = 7.2)

filter_shape <- function(shape_df, shape_str) {
  if(shape_str == "sphere") {
    return(shape_df$cell_shape=="coccus")
  } else {
    return(with(shape_df, (cell_shape=="bacillus")|(cell_shape=="coccobacillus")|(cell_shape=="vibrio")))
  }
}

L_env <- 0.03 # Remember to change this, if you want total volume
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

madinPlot <- ggplot(data_select, aes(x=log10(S.mean/V.mean), y=log10(growth_rate),col=environment)) + 
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_madin_free$coefficients[1]+lm_madin_free$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75, linetype="dashed") +
  stat_function(fun=function(x) lm_madin_host$coefficients[1]+lm_madin_host$coefficients[2]*x, 
                colour = colBlindScale[2], size=line_thick*0.75, linetype="dashed") +
  stat_function(fun=function(x) lm_madinSV$coefficients[1]+lm_madinSV$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]')))

madin_scaling <- grid.arrange(madinPlot, nrow = 1)
ggsave(file="madin_scaling.pdf", plot=madin_scaling, width = 6.5, height = 5)


#df_lambda_env_corr[df_lambda_env_corr$total_thickness_nm<10^1.25,]$species
lm_thick <- lm(log10(total_thickness_nm)~log10(V.mean), data=df_lambda_env_corr)
summary(lm_thick)

thickness_plot <- ggplot(df_lambda_env_corr, aes(x=log10(V.mean), y=log10(total_thickness_nm), color=(species))) + 
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_thick$coefficients[1]+lm_thick$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
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
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(rRNA gene number)'))) + xlab(bquote(bold('S/V, Log'[10]*'['*Pi*' ('*mu*''*m^-1*')]')))

p_trna_sv <-ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(S.mean/V.mean), y=log10(tRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_trna_sv$coefficients[1]+lm_trna_sv$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
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
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(rRNA gene number)'))) + xlab(bquote(bold('Growth rate, Log'[10]*'['*lambda*' (h'^-1*')]')))

p_trna_lambda <- ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(growth_rate), y=log10(tRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_trna_lambda$coefficients[1]+lm_trna_lambda$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
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
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin=unit(c(1,1,1.5,1.2),"cm"),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab(bquote(bold('Log'[10]*'(rRNA gene number)'))) + xlab(bquote(bold('Volume, Log'[10]*'[V ('*mu*''*m^3*')]')))

p_trna_v <-ggplot(df_lambda[!is.na(df_lambda$shape),], aes(x=log10(V.mean), y=log10(tRNA.genes), col=shape)) +
  geom_point(size = point_size1) +
  stat_function(fun=function(x) lm_trna_v$coefficients[1]+lm_trna_v$coefficients[2]*x, 
                colour = colBlindScale[1], size=line_thick*0.75) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        axis.line = element_line(colour = "black"),
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
medium_text <- c(bquote("70 hours, "*d[p]*"=0.002"), 
                 bquote("25-70 hours, "*d[p]*"=0.003"),
                 bquote("5-25 hours, "*d[p]*"=0.009"),
                 bquote("2-5 hours, "*d[p]*"=0.04"),
                 bquote("Fraction weighted mean, "*d[p]*"=0.005"))
line_thick <- 1
p_sensitivity_analysis_dp <- ggplot(df_lambda, aes(x=log10(S.mean/V.mean), y=log10(growth_rate), color=(ID))) + 
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
                                                      "e"=colBlindScale[6]), labels = medium_text) +
  geom_point(col="black", size = point_size1) +
  theme(axis.text.y   = axes_style,
        axis.text.x   = axes_style,
        axis.title.y  = axes_style,
        axis.title.x  = axes_style,
        legend.position = "none",
        axis.line = element_line(colour = "black"),
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
                                                      "e"=colBlindScale[6]), labels = medium_text) +
  geom_point(col="black", size = point_size1) +
  theme(axis.text.y   = axes_style,
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
