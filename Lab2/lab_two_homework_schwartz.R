#
# Shawn Schwartz, 2019
# EEB C174 UCLA Spring 2019
# Lab 2 HW - Bat Testes
#

#clean up workspace
rm(list=ls())

#includes
library(phytools)
library(dplyr)

#directories
wd_path <- "~/Developer/EEB-C174-Labs/Lab2"
setwd(wd_path)
output_path <- "output/"
resources_path <- "resources/"

#### Import Data ####
bat_tree <- read.tree(paste0(resources_path, "bat_tree.tre"))
bat_traits <- read.csv(paste0(resources_path, "batTraits.csv"))

#### Clean Up Datasets ####
rownames(bat_traits) <- bat_traits$Species
cbind(rownames(bat_traits), bat_tree$tip.label) #check order
bat_traits <- bat_traits[match(bat_tree$tip.label, rownames(bat_traits)), ]
cbind(rownames(bat_traits), bat_tree$tip.label) #confirm matching order

#### Naive Analysis of Trait Data ####
trait_bodymass <- bat_traits$Body.mass..g.
trait_bodymass

trait_groupsize <- bat_traits$Group.size
trait_groupsize

trait_testesmass <- bat_traits$Testes.mass..g.
trait_testesmass

##compute relative traits
###testes mass adjusted for body mass
relative_testes <- trait_testesmass/trait_bodymass
relative_testes

plot(trait_bodymass, relative_testes)

##transform trait data with natural log
ln_trait_bodymass <- log(trait_bodymass)
ln_trait_bodymass

ln_trait_testesmass <- log(relative_testes)
ln_trait_testesmass

ln_trait_groupsize <- log(trait_groupsize)
ln_trait_groupsize

plot(ln_trait_bodymass, ln_trait_testesmass)
plot(ln_trait_groupsize, ln_trait_testesmass)

##fit a linear model to the ln(data)
fit <- lm(ln_trait_testesmass ~ ln_trait_bodymass)
plot(ln_trait_bodymass, ln_trait_testesmass, xlab = "ln(Body Mass)", ylab = "ln(Relative Testes Mass)")
abline(fit)

summary(fit)
lm_fit_summary <- capture.output(print(summary(fit)))
writeLines(lm_fit_summary, con = file(paste(output_path, "lm_summary_fit.txt")))

group_fit <- lm(ln_trait_testesmass ~ ln_trait_groupsize)
summary(group_fit)
lm_group_fit_summary <- capture.output(print(summary(group_fit)))
writeLines(lm_group_fit_summary, con = file(paste(output_path, "lm_summary_group_fit.txt")))

#### Analysis with PICs ####
PIC_ln_trait_bodymass <- pic(ln_trait_bodymass, bat_tree)
PIC_ln_trait_testesmass <- pic(ln_trait_testesmass, bat_tree)
PIC_fit <- lm(PIC_ln_trait_testesmass ~ PIC_ln_trait_bodymass - 1) #-1 to remove intercept
plot(PIC_ln_trait_bodymass, PIC_ln_trait_testesmass, xlab = "Contrasts in ln(Body Mass)", ylab = "Contrasts in ln(Relative Testes Mass)")
abline(PIC_fit)

summary(PIC_fit)
lm_PIC_fit_summary <- capture.output(print(summary(PIC_fit)))
writeLines(lm_PIC_fit_summary, con = file(paste(output_path, "lm_summary_PIC_fit.txt")))

PIC_ln_trait_groupsize <- pic(ln_trait_groupsize, bat_tree)
PIC_group_fit <- lm(PIC_ln_trait_testesmass ~ PIC_ln_trait_groupsize - 1)
summary(PIC_group_fit)
lm_PIC_group_fit_summary <- capture.output(print(summary(PIC_group_fit)))
writeLines(lm_PIC_group_fit_summary, con = file(paste(output_path, "lm_summary_PIC_group_fit.txt")))

#### Generate 2D Phylomorphospace ####
ln_traits <- cbind(ln_trait_bodymass, ln_trait_testesmass)
rownames(ln_traits) <- bat_traits$Species
pdf(paste0(output_path,"bat_testes_phylomorphospace.pdf"))
  phylomorphospace(tree = bat_tree, ln_traits, xlab = "ln(Body Mass)", ylab = "ln(Relative Testes Mass)")
dev.off()

ln_traits_group <- cbind(ln_trait_groupsize, ln_trait_testesmass)
rownames(ln_traits_group) <- bat_traits$Species
pdf(paste0(output_path,"bat_testes_phylomorphospace_groupsize.pdf"))
  phylomorphospace(tree = bat_tree, ln_traits_group, xlab = "ln(Group Size)", ylab = "ln(Relative Testes Mass)")
dev.off()

#### Prepare Plots for Report ####
pdf(paste0(output_path,"bat_fit_analyses_GROUP_combined.pdf"))
  par(mfrow = c(1,2))
  #Plotting Naive Analysis
  plot(ln_trait_groupsize, ln_trait_testesmass, xlab = "ln(Group Size)", ylab = "ln(Relative Testes Mass)", main = "Naive Analysis")
  abline(fit)
  
  #Plotting Analysis using PICs
  plot(PIC_ln_trait_groupsize, PIC_ln_trait_testesmass, xlab = "PICs for ln(Group Size)", ylab = "PICs for ln(Relative Testes Mass)", main = "Analysis using PICs")
  abline(PIC_fit)
dev.off()

pdf(paste0(output_path,"bat_fit_analyses_GROUP_separate.pdf"))
  #Plotting Naive Analysis
  plot(ln_trait_groupsize, ln_trait_testesmass, xlab = "ln(Group Size)", ylab = "ln(Relative Testes Mass)", main = "Naive Analysis")
  abline(fit)
  
  #Plotting Analysis using PICs
  plot(PIC_ln_trait_groupsize, PIC_ln_trait_testesmass, xlab = "PICs for ln(Group Size)", ylab = "PICs for ln(Relative Testes Mass)", main = "Analysis using PICs")
  abline(PIC_fit)
dev.off()

pdf(paste0(output_path,"bat_fit_analyses_separate.pdf"))
  plot(ln_trait_bodymass, ln_trait_testesmass, xlab = "ln(Body Mass)", ylab = "ln(Relative Testes Mass)", main = "Naive Analysis")
  abline(fit)

  #Plotting Analysis using PICs
  plot(PIC_ln_trait_bodymass, PIC_ln_trait_testesmass, xlab = "PICs for ln(Body Mass)", ylab = "PICs for ln(Relative Testes Mass)", main = "Analysis using PICs")
  abline(PIC_fit)
dev.off()

pdf(paste0(output_path,"bat_fit_analyses_combined.pdf"))
  par(mfrow = c(1,2))
  #Plotting Naive Analysis
  plot(ln_trait_bodymass, ln_trait_testesmass, xlab = "ln(Body Mass)", ylab = "ln(Relative Testes Mass)", main = "Naive Analysis")
  abline(fit)

  #Plotting Analysis using PICs
  plot(PIC_ln_trait_bodymass, PIC_ln_trait_testesmass, xlab = "PICs for ln(Body Mass)", ylab = "PICs for ln(Relative Testes Mass)", main = "Analysis using PICs")
  abline(PIC_fit)
dev.off()

#### Get the Bat Tree ####
ladder_bat_tree <- ladderize(bat_tree)
pdf(paste0(output_path,"bat_tree_visual.pdf"))
  plot(ladder_bat_tree, cex = 0.6)
dev.off()