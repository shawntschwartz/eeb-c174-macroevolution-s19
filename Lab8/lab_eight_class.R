#
# EEB C174 Week 8, Lab 8
# Comparative Biology and Macroevolution
# BAMM and FiSSE
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(BAMMtools)
library(coda)

cwd <- "~/Developer/EEB-C174-Labs/Lab8/BAMM_Run"
labeight_dir <- "~/Developer/EEB-C174-Labs/Lab8"
setwd(cwd)

resources_path <- "resources/"
output_path <- "output/"
source_path <- "source/"

tree <- read.tree("primatetree.tre")
tree <- ladderize(tree)
plot(x = tree, cex = 0.2, no.margin = T)
axisPhylo()

#check bamm assumptions
is.ultrametric(tree)
is.binary(tree)
min(tree$edge.length) > 0

#generate bamm control file
# estimate priors and create control files
priors <- setBAMMpriors(phy = tree, traits = "primates_logmass.txt", outfile = NULL)
generateControlFile(file = "controlfile.txt", type = "trait", params = list(
  treefile = "primatetree.tre",
  traitfile = "primates_logmass.txt",
  seed = sample(1:1000000, 1),
  overwrite = "0",
  expectedNumberOfShifts = "30",
  betaInitPrior = as.numeric(priors["betaInitPrior"]),
  betaShiftPrior = as.numeric(priors["betaShiftPrior"]),
  numberOfGenerations = "90000000",
  mcmcWriteFreq = "30000",
  eventDataWriteFreq = "30000",
  printFreq = "30000",
  acceptanceResetFreq = "30000",
  outName = "classwork",
  numberOfChains = "8"
))

#### import BAMM results ####
edata <- getEventData(tree, eventdata = "classwork_event_data.txt", burnin = 0.1, type = "trait")
summary(edata)
#number of analyzed xxxx posterior samples are the number of possible rate shift combinations that were tested
#a lot might be the same and that's okay

#### check quality of BAMM results ####
mcmcout <- read.csv("classwork_mcmc_out.txt")
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout)) # Discard the first 10% of samples as burnin 
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
#there is still some possible of auto correlation between subsequent steps because they are following each other

#### analyze BAMM results with BAMMtools

#analysis of rate shifts
plotPrior("classwork_mcmc_out.txt", expectedNumberOfShifts = 30)

#analysis of rates
##phylorate plot displays model-averaged rates of trait evolution (i.e., Brownian rate parameter) on branches of your tree
s <- plot.bammdata(edata, labels = T, font = 3, cex = 0.1, logcolor = T)
title(main = "Mean phenotypic rate", sub = "time before present")
addBAMMlegend(s, location = "left", nTicks = 1)
axisPhylo()

css <- credibleShiftSet(edata, expectedNumberOfShifts = 30, threshold = 5, set.limit = 0.95)
css$number.distinct
summary(css)

sss <- plot.credibleshiftset(css, border = F, logcolor = T)

plot.new()
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts = 30)
ss <- plot.bammdata(best, labels = T, font = 3, cex = 0.1, logcolor = T)
title(main = "Best shift configuration", sub = "time before present")
addBAMMlegend(ss, location = "left", nTicks = 1)
addBAMMshifts(best, cex = 3, pch = 1)
axisPhylo()

#This should correspond to the plot in the credible shift set (CSS) with the highest frequency, 
#but will look slightly different because it is just a single sample from the posterior whereas 
#CSS is averaged across samples in the posterior.

##IMPORTANT NOTES FOR FINAL PROJECT: ####
####Justification for sampling interval (3000 => less likely for autocorrelation vs 1000 having a higher chance of being correlated with one another)
#marginal probabilities tree (scaled)

#Put best shift configuration graph with circles on lightening talk
#if they are going to cooler colors, we see that the rate is increasing through time in this plot, but there isn't a strong enough signal to denote a rate shift


par(font = 1)
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, cex = 0.2, no.margin = T)
title(sub = "Marginal shift probability")
add.scale.bar(x = 0.5, y = 0.5, font = 1)
#the longer the branch, the higher the likelihood that there was a higher increase in body size (already log-transformed...same for the data for the homework)

#### rate through time analysis ####
###THIS PLOT WILL BE VERY USEFUL FOR INTERPRETING THE DATA###
#questions: have changes in the environment through time impacted the diversification of my group
#example: Eocene (it got really hot in terms of the climate, and a lot of carnivores/animals got really small)
#we see that monkeys don't really care and they chug along through time and there isn't a big change in trait evolution
#so maybe this tells us that monkey trait evolution isn't really effected by the heat through time, but in fact might be decreasing their speciation rate with the global cooling event that has been occuring
#Thus, you can correlate these with known geological events that were occuring through time (questions that can be asked with rate through time plots)
#asteroid killed everything over 1 KG (KT extinction)

#THIS IS FOR CONTINUOUS TRAITS
plotRateThroughTime(edata,
                    ratetype = "auto",
                    avgCol = "black",
                    intervalCol = "gray80",
                    intervals = c(0.05, 0.95),
                    opacity = 1)

plotRateThroughTime(edata,
                    ratetype = "auto",
                    avgCol = "black",
                    intervalCol = "gray80",
                    intervals = c(0.05, 0.95),
                    opacity = 1)

#FiSSE is for DISCRETE TRAITS

#FiSSE
library(ape)
library(phangorn)
library(diversitree)

setwd(labeight_dir)
source(paste0(source_path,"traitDependent_functions.R")) # functions for FiSSE method

example_tree <- read.tree(paste0(source_path,"example_tree.tre"))
example_tree <- ladderize(example_tree)
xx <- read.csv(paste0(source_path,"example_trait.csv"), header = F)
traits <- xx$V2
names(traits) <- xx$V1
traits
traits <- traits[example_tree$tip.label]

#plot tree with trait labels on tips
colvec <- rep("white", length(traits))
colvec[traits == 1] <- "black"
plot.phylo(example_tree, show.tip.label = F, no.margin = T)
plot.phylo(example_tree, show.tip.label = F, no.margin = T, type = "fan")
tiplabels(pch = 21, bg = colvec, cex = 0.6)

res <- FISSE.binary(example_tree, traits)
res

#we want to take the other direction of the p-value test
pval_2tailed <- min(res$pval, 1 - res$pval) * 2
pval_2tailed
#this will be a little different when we do it because we are doing a bunch of simulations and throwing out the bad ones
##to figure out what this specific relationship is, we use res

#trait1 -> black dots has a higher speciation rate than trait 0 (open circle trait)
