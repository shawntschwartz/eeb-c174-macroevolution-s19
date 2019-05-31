#
# EEB C174 Week 9, Lab 9
# Comparative Biology and Macroevolution
# Biogeographic Distribution
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(ape)

cwd <- "~/Developer/EEB-C174-Labs/Lab9/"
setwd(cwd)

#read in data
sample_tree <- read.tree("sample_tree.tre")
sample_tree <- ladderize(sample_tree)
plot(sample_tree)
axisPhylo()

#check tree assumptions
is.ultrametric(sample_tree)
is.binary(sample_tree)
min(sample_tree$edge.length) > 0

sample_data <- read.csv("sample_data.csv", header = TRUE)
write.table(sample_data, "sample_data.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#### plot the biogeographic range ####
row.names(sample_data) <- sample_data$X
sample_data <- sample_data[sample_tree$tip.label,]

plot(sample_tree, cex = 0.001)
tiplabels(text = sample_tree$tip.label, cex = 0.5, offset = -0.2)
tiplabels(pie = sample_data$A, piecol = c("black", "white"), cex = 0.2, offset = 0)
tiplabels(pie = sample_data$B, piecol = c("black", "white"), cex = 0.2, offset = 0.1) 
tiplabels(pie = sample_data$C, piecol = c("black", "white"), cex = 0.2, offset = 0.2) 
tiplabels(pie = sample_data$D, piecol = c("black", "white"), cex = 0.2, offset = 0.3) 
axisPhylo()