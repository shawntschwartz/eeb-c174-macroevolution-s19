#
# Shawn Schwartz, 2019
# EEB C174 UCLA Spring 2019
# Lab 4 HW - Evolution of Carnivoran Dentition and Body Size
#

#clean up workspace
rm(list=ls())

#includes
library(phytools)
library(geiger)
library(phylolm)

cwd <- "~/Developer/EEB-C174-Labs/Lab4"
resources_path <- "resources/"
output_path <- "output/"

setwd(cwd)

#### Read in Tree and Trait Data ####
carnivoran_tree <- read.tree(paste0(resources_path,"carnivoran_tree.tre"))
carnivoran_trait <- read.csv(paste0(resources_path,"carnivoran_trait.csv"), header = TRUE, sep = ",")
carnivoran_SE <- read.csv(paste0(resources_path,"carnivoran_SE.csv"), header = TRUE, sep = ",")

pdf(paste0(output_path,"1_carnivoran_trait.pdf"))
plot(carnivoran_tree, cex = 0.25)
  axisPhylo()
  add.scale.bar()
dev.off()

#### Fit Models ####
#c1: fit models to cross sectional shape of the upper canine
C1_data <- carnivoran_trait$C1
names(C1_data) <- carnivoran_trait$Taxa
C1_SE <- carnivoran_SE$C1
names(C1_SE) <- carnivoran_SE$Taxa

carniv_fit_c1.bm <- fitContinuous(phy = carnivoran_tree, dat = C1_data, SE = C1_SE, model = c("BM"))
carniv_fit_c1.ou <- fitContinuous(phy = carnivoran_tree, dat = C1_data, SE = C1_SE, model = c("OU"))
carniv_fit_c1.eb <- fitContinuous(phy = carnivoran_tree, dat = C1_data, SE = C1_SE, model = c("EB"))
carniv_fit_c1.white <- fitContinuous(phy = carnivoran_tree, dat = C1_data, SE = C1_SE, model = c("white"))

#RLGA: fit models to the relative upper grinding area
RLGA_data <- carnivoran_trait$RLGA
names(RLGA_data) <- carnivoran_trait$Taxa
RLGA_SE <- carnivoran_SE$RLGA
names(RLGA_SE) <- carnivoran_SE$Taxa

carniv_fit_RLGA.bm <- fitContinuous(phy = carnivoran_tree, dat = RLGA_data, SE = RLGA_SE, model = c("BM"))
carniv_fit_RLGA.ou <- fitContinuous(phy = carnivoran_tree, dat = RLGA_data, SE = RLGA_SE, model = c("OU"))
carniv_fit_RLGA.eb <- fitContinuous(phy = carnivoran_tree, dat = RLGA_data, SE = RLGA_SE, model = c("EB"))
carniv_fit_RLGA.white <- fitContinuous(phy = carnivoran_tree, dat = RLGA_data, SE = RLGA_SE, model = c("white"))

#ln-cbrt-bodymass: fit models to the log cuberoot of body mass
lncbrt_bodymass_data <- carnivoran_trait$log_cuberoot_mass
names(lncbrt_bodymass_data) <- carnivoran_trait$Taxa
lncbrt_bodymass_SE <- carnivoran_SE$logMass
names(lncbrt_bodymass_SE) <- carnivoran_SE$Taxa

carniv_fit_lncbrt_bodymass.bm <- fitContinuous(phy = carnivoran_tree, dat = lncbrt_bodymass_data, SE = lncbrt_bodymass_SE, model = c("BM"))
carniv_fit_lncbrt_bodymass.ou <- fitContinuous(phy = carnivoran_tree, dat = lncbrt_bodymass_data, SE = lncbrt_bodymass_SE, model = c("OU"))
carniv_fit_lncbrt_bodymass.eb <- fitContinuous(phy = carnivoran_tree, dat = lncbrt_bodymass_data, SE = lncbrt_bodymass_SE, model = c("EB"))
carniv_fit_lncbrt_bodymass.white <- fitContinuous(phy = carnivoran_tree, dat = lncbrt_bodymass_data, SE = lncbrt_bodymass_SE, model = c("white"))

#### Compare Results Using AIC ####
c1.data <- c(carniv_fit_c1.bm$opt$aic, carniv_fit_c1.ou$opt$aic, carniv_fit_c1.eb$opt$aic, carniv_fit_c1.white$opt$aic)
rlga.data <- c(carniv_fit_RLGA.bm$opt$aic, carniv_fit_RLGA.ou$opt$aic, carniv_fit_RLGA.eb$opt$aic, carniv_fit_RLGA.white$opt$aic)
lncbrt.data <- c(carniv_fit_lncbrt_bodymass.bm$opt$aic, carniv_fit_lncbrt_bodymass.ou$opt$aic, carniv_fit_lncbrt_bodymass.eb$opt$aic, carniv_fit_lncbrt_bodymass.white$opt$aic)
AICresults <- rbind(c1.data, rlga.data, lncbrt.data)
colnames(AICresults) <- c("Brownian Motion", "Ornstein Uhlenbeck", "Early Burst", "White noise")
AICresults

write.csv(AICresults, file = paste0(output_path,"AICresults.csv"))

#transform AIC values into Akaike weights
##which can be interpreted as conditional probabilities for each model
AICw_c1_vals <- aicw(c1.data)
AICw_rlga_vals <- aicw(rlga.data)
AICw_lncbrt_vals <- aicw(lncbrt.data)

write.csv(AICw_c1_vals, file = paste0(output_path,"AICw_c1_vals.csv"))
write.csv(AICw_rlga_vals, file = paste0(output_path,"AICw_rlga_vals.csv"))
write.csv(AICw_lncbrt_vals, file = paste0(output_path,"AICw_lncbrt_vals.csv"))

pdf(paste0(output_path,"c1_AICw_barplot.pdf"))
  barplot(AICw_c1_vals[,3], names.arg = c("Brownian Motion", "Ornstein Uhlenbeck", "Early Burst", "White noise"), ylim = c(0,1), main = "Cross Sectional Shape of Upper Canine", cex.names = 0.7)
  abline(h = 1, lty = "dashed")
dev.off()

pdf(paste0(output_path,"rlga_AICw_barplot.pdf"))
  barplot(AICw_rlga_vals[,3], names.arg = c("Brownian Motion", "Ornstein Uhlenbeck", "Early Burst", "White noise"), ylim = c(0,1), main = "Relative Upper Grinding Area", cex.names = 0.7)
  abline(h = 1, lty = "dashed")
dev.off()

pdf(paste0(output_path,"lncbrt_AICw_barplot.pdf"))
  barplot(AICw_lncbrt_vals[,3], names.arg = c("Brownian Motion", "Ornstein Uhlenbeck", "Early Burst", "White noise"), ylim = c(0,1), main = "Log Cuberoot of Body Mass", cex.names = 0.7)
  abline(h = 1, lty = "dashed")
dev.off()

#### Disparity Through Time Plots ####
pdf(paste0(output_path,"c1_dtt_plot.pdf"))
  c1.dtt <- dtt(phy = carnivoran_tree, data = C1_data, nsim = 100, index = c("avg.sq"), plot = TRUE)
  title("Cross Sectional Shape of Upper Canine")
dev.off()

pdf(paste0(output_path,"rlga_dtt_plot.pdf"))
  rlga.dtt <- dtt(phy = carnivoran_tree, data = RLGA_data, nsim = 100, index = c("avg.sq"), plot = TRUE)
  title("Relative Upper Grinding Area")
dev.off()

pdf(paste0(output_path,"lncbrt_dtt_plot.pdf"))
  lncbrt.dtt <- dtt(phy = carnivoran_tree, data = lncbrt_bodymass_data, nsim = 100, index = c("avg.sq"), plot = TRUE)
  title("Log Cuberoot of Body Mass")
dev.off()


