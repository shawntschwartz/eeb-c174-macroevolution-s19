#
# Shawn Schwartz, 2019
# EEB C174 UCLA Spring 2019
# Lab 7 HW - BAMM
#

rm(list=ls())

#imports
library(BAMMtools)
library(coda)

cwd <- "~/Developer/EEB-C174-Labs/Lab7/"
setwd(cwd)

resources_path <- "resources/"
output_path <- "output/"

tree <- read.tree(paste0(resources_path,"whaleTree.tre"))
plot(tree, cex = 0.7)
axisPhylo()

#### Verify BAMM assumptions ####

# Tree Must Be Ultrametric
is.ultrametric(tree)

# Tree Must Be Binary (all of the nodes need to lead to 2 branches, and no more than 2)
is.binary.tree(tree)

# Check to make sure all branch lengths are greater than 0
min(tree$edge.length)

#### Create BAMM Control File ####
#estimate priors and create control files
priors <- setBAMMpriors(tree, outfile = NULL)
priors

generateControlFile(file = paste0(resources_path,"whaleBAMMcontrolfile.txt"), params = list(
  treefile = "whaleTree.tre",
  globalSamplingFraction = "0.8", #This is 1 for complete sampling (for homework, put in a fraction like 0.8)
  seed = sample(1:5000000, 1),
  overwrite = "0",
  expectedNumberOfShifts = "1",
  lambdaInitPrior = as.numeric(priors["lambdaInitPrior"]),
  lambdaShiftPrior = as.numeric(priors["lambdaShiftPrior"]),
  muInitPrior = as.numeric(priors["muInitPrior"]),
  numberOfGenerations = "5000000",
  mcmcWriteFreq = "1000",
  eventDataWriteFreq = "1000",
  printFreq = "1000",
  acceptanceResetFreq = "1000",
  outName = "whaleBAMMhomeworkoutput",
  ##numberOfChains = "2", #set to the number of CPUs
  numberOfChains = "4", #set to the number of CPUs
  deltaT = "0.01"
))


labrid_tree <- read.tree(paste0(resources_path,"Labridae.tre"))
#### Verify BAMM assumptions ####

# Tree Must Be Ultrametric
is.ultrametric(labrid_tree)

library(phytools)
labrid_tree <- force.ultrametric(labrid_tree, method=c("nnls","extend"))
is.ultrametric(labrid_tree)

# Tree Must Be Binary (all of the nodes need to lead to 2 branches, and no more than 2)
is.binary.tree(labrid_tree)

# Check to make sure all branch lengths are greater than 0
min(labrid_tree$edge.length)

#### Create BAMM Control File ####
#estimate priors and create control files
labrid_priors <- setBAMMpriors(labrid_tree, outfile = NULL)
labrid_priors

generateControlFile(file = paste0(resources_path,"LABRIDBAMMcontrolfile.txt"), params = list(
  treefile = "Labridae.tre",
  globalSamplingFraction = "0.8", #This is 1 for complete sampling (for homework, put in a fraction like 0.8)
  seed = sample(1:5000000, 1),
  overwrite = "0",
  expectedNumberOfShifts = "1",
  lambdaInitPrior = as.numeric(priors["lambdaInitPrior"]),
  lambdaShiftPrior = as.numeric(priors["lambdaShiftPrior"]),
  muInitPrior = as.numeric(priors["muInitPrior"]),
  numberOfGenerations = "5000000",
  mcmcWriteFreq = "1000",
  eventDataWriteFreq = "1000",
  printFreq = "1000",
  acceptanceResetFreq = "1000",
  outName = "LABRIDBAMMhomeworkoutput",
  ##numberOfChains = "2", #set to the number of CPUs
  numberOfChains = "4", #set to the number of CPUs
  deltaT = "0.01"
))

plot(labrid_tree, cex = 0.3)

whale_data <- getEventData(tree, eventdata = paste0(resources_path,"whaleBAMMhomeworkoutput_event_data.txt"), burnin = 0.1)
summary(whale_data)

#### check quality of BAMM results ####
mcmcout <- read.csv(paste0(resources_path,"whaleBAMMhomeworkoutput_mcmc_out.txt"))
plot(mcmcout$logLik ~ mcmcout$generation)

# test for convergence of the MCMC chains #
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

computeBayesFactors(paste0(resources_path,"whaleBAMMhomeworkoutput_mcmc_out.txt"), expectedNumberOfShifts = 1, burnin = 0.1)
plotPrior(paste0(resources_path,"whaleBAMMhomeworkoutput_mcmc_out.txt"), expectedNumberOfShifts = 1)

s <- plot.bammdata(whale_data, spex = "s", labels = T, font = 3, cex = 0.5)
title(main = "Mean speciation rate", sub = "time before present")
addBAMMlegend(s, location = "topleft", nTicks = 1)
axisPhylo()

css <- credibleShiftSet(whale_data, expectedNumberOfShifts = 1, threshold = 5, set.limit = 0.95)
css$number.distinct
summary(css)
sss <- plot.credibleshiftset(css, border = F)

plot.new()
best <- getBestShiftConfiguration(whale_data, expectedNumberOfShifts = 1)
par(mfrow = c(1,1))
par(mar=c(5, 4.5, 2, 0.5))
ss <- plot.bammdata(best, labels = T, font = 3, cex = 0.55)
title(main = "Best shift configuration", sub = "time before present")
addBAMMlegend(ss, location = "topleft", nTicks = 1)
addBAMMshifts(best, cex = 3, pch = 1)
axisPhylo()

par(font = 1)
marg_probs <- marginalShiftProbsTree(whale_data)
plot.phylo(marg_probs, cex = 0.55)
title(main = "Marginal shift probability")
add.scale.bar(x = 0.5, y = 0.5, font = 1)


#### clade-specific evolutionary rates ####
global_rates <- getCladeRates(whale_data)
mean(global_rates$lambda)
quantile(global_rates$lambda, c(0.05, 0.95))
rateshift_MRCA <- getMRCA(tree, tip = c("Orcinus_orca_AF084061", "Stenella_longirostris_AF084103"))
clade_rates <- getCladeRates(whale_data, node = rateshift_MRCA)
mean(clade_rates$lambda)
quantile(clade_rates$lambda, c(0.05, 0.95))

smaller_clade_MRCA <- getMRCA(tree, tip = c("Delphinapterus_leucas_DLU72037","Phocoenoides_dalli_PDU09679"))
smaller_clade_rates <- getCladeRates(whale_data, node = smaller_clade_MRCA)
mean(smaller_clade_rates$lambda)
quantile(smaller_clade_rates$lambda, c(0.05, 0.95))

lowest_clade_rates_MRCA <- getMRCA(tree, tip = c("Caperea_marginata_X75586","Balaenoptera_acutorostrata"))
lowest_clade_rates_MRCA <- getMRCA(tree, tip = c("Ziphius_cavirostris_AF304075_","Mesoplodon_europaeus_X92537"))
lowest_clade_rates <- getCladeRates(whale_data, node = lowest_clade_rates_MRCA)
mean(lowest_clade_rates$lambda)
quantile(lowest_clade_rates$lambda, c(0.05, 0.95))

#### Rate-through-time analysis ####
par(font = 1)
plotRateThroughTime(whale_data,
                    ratetype = "speciation",
                    avgCol = "black",
                    intervalCol = "gray80",
                    intervals = c(0.05, 0.95),
                    opacity = 1)



#two planel plot for lab report
par(mfrow = c(1,2))
par(mar=c(5, 1.0, 2, 0.0))
s <- plot.bammdata(whale_data, spex = "s", labels = T, font = 3, cex = 0.5)
title(main = "Mean speciation rate", sub = "time before present (my)")
addBAMMlegend(s, location = "topleft", nTicks = 1)
axisPhylo()

par(font = 1)
marg_probs <- marginalShiftProbsTree(whale_data)
plot.phylo(marg_probs, cex = 0.55)
title(main = "Marginal shift probability")
add.scale.bar(x = 1.5, y = 0.0, font = 1)

