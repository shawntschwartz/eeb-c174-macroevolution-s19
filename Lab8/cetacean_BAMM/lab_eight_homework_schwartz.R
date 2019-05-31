#
# Shawn Schwartz, 2019
# EEB C174 UCLA Spring 2019
# Lab 8 HW - BAMM and FiSSE
#

rm(list=ls())

#imports
library(BAMMtools)
library(coda)

cwd <- "~/Developer/EEB-C174-Labs/Lab8/cetacean_BAMM"
setwd(cwd)

tree <- read.tree("whaleTree.tre")
tree <- ladderize(tree)
plot(x = tree, cex = 0.2, no.margin = T)
axisPhylo()

#check bamm assumptions
is.ultrametric(tree)
is.binary(tree)
min(tree$edge.length) > 0


#generate bamm control file
# estimate priors and create control files
priors <- setBAMMpriors(phy = tree, traits = "whaleSize.txt", outfile = NULL)
generateControlFile(file = "cetacean_bamm_controlfile.txt", type = "trait", params = list(
  treefile = "whaleTree.tre",
  traitfile = "whaleSize.txt",
  seed = sample(1:1000000, 1),
  overwrite = "0",
  expectedNumberOfShifts = "2",
  betaInitPrior = as.numeric(priors["betaInitPrior"]),
  betaShiftPrior = as.numeric(priors["betaShiftPrior"]),
  numberOfGenerations = "120000000",
  mcmcWriteFreq = "30000",
  eventDataWriteFreq = "30000",
  printFreq = "15000",
  acceptanceResetFreq = "30000",
  outName = "cetacean_bamm_homework",
  numberOfChains = "8"
))

edata <- getEventData(tree, eventdata = "2_expected_shifts/cetacean_bamm_homework_event_data.txt", burnin = 0.1, type = "trait")
summary(edata)

#### Assess MCMC convergence ####
mcmcout <- read.csv("2_expected_shifts/cetacean_bamm_homework_mcmc_out.txt")
plot(mcmcout$logLik ~ mcmcout$generation)

#### Test for convergence of the MCMC chains with the coda package for R ####
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#### Analysis of rate shifts ####
plotPrior("2_expected_shifts/cetacean_bamm_homework_mcmc_out.txt", expectedNumberOfShifts = 2)
s <- plot.bammdata(edata, labels = T, font = 3, cex = 0.1, logcolor = T)
title(main = "Mean phenotypic rate", sub = "time before present")
addBAMMlegend(s, location = "left", nTicks = 1)
axisPhylo()

css <- credibleShiftSet(edata, expectedNumberOfShifts = 30, threshold = 5, set.limit = 0.95)
css$number.distinct
summary(css)

sss <- plot.credibleshiftset(css, border = F, logcolor = T)

plot.new()
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts = 2)

ss <- plot.bammdata(best, labels = T, font = 3, cex = 0.1, logcolor = T) 
title(main = "Best shift configuration", sub = "time before present") 
addBAMMlegend(ss, location = "left", nTicks = 1)
addBAMMshifts(best, cex = 3, pch = 1)
axisPhylo()

#### scaled branch lengths ####
par(font = 1)
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, cex = 0.45, no.margin = T)
title(sub = "Marginal shift probability")
add.scale.bar(x = 0.5, y = 0.5, font = 1)


#### rate through time analysis ####
par(font = 1)
plotRateThroughTime(edata,
                    ratetype = "auto",
                    avgCol = "black", intervalCol = "gray80", intervals = c(0.05, 0.95), opacity = 1)

#### FiSSE ####
library(phangorn)
library(diversitree)
library(ape)

source("traitDependent_functions.R")

diet <- read.csv(file = "whaleDiet.csv", stringsAsFactors = F)
generalist <- diet$generalist
names(generalist) <- diet$species
generalist <- generalist[tree$tip.label]
colvec <- rep("white", length(generalist))
colvec[generalist == 1] <- "black"
plot.phylo(tree, show.tip.label = F, cex = 0.4, no.margin = T)
tiplabels(pch = 21, bg = colvec, cex = 0.8)
add.scale.bar(x = 0.5, y = 0.5, font = 1)
title(sub = "time before present (my)")

res <- FISSE.binary(tree, generalist)
res

# two-tailed pvalue is obtained as
pval_2tailed <- min(res$pval, 1 - res$pval) * 2
pval_2tailed

#two planel plot for lab report
par(mfrow = c(1,2))
par(mar=c(5, 1.0, 2, 0.0))
#plotPrior("2_expected_shifts/cetacean_bamm_homework_mcmc_out.txt", expectedNumberOfShifts = 2)
s <- plot.bammdata(edata, labels = T, font = 3, cex = 0.5, logcolor = T)
title(main = "Mean phenotypic rate", sub = "time before present (my)")
addBAMMlegend(s, location = "left", nTicks = 1)
axisPhylo()

par(font = 1)
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, cex = 0.5, no.margin = F)
title(main = "Marginal shift probability")
add.scale.bar(x = 0.75, y = 0.5, font = 1)


library(geiger)
library(apTreeshape)
library(phytools)

ltt.plot(tree, log = "y")
lines(c(0, 1), c(log(2), log(100)), lty = "dashed", lwd = 2, col = "red")
gammaStat(tree)
