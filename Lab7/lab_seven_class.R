#
# EEB C174 Week 7, Lab 7
# Comparative Biology and Macroevolution
# BAMM
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(BAMMtools)
library(coda)

cwd <- "~/Developer/EEB-C174-Labs/Lab7/"
setwd(cwd)

resources_path <- "resources/"
output_path <- "output/"

#NOTES:
#this tree looks like there is some different process going on here: 
##large division/rate shift potentially going on at 25MY
tree <- read.tree(paste0(resources_path,"tree.tre"))
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

generateControlFile(file = paste0(resources_path,"controlfile.txt"), params = list(
  treefile = "tree.tre",
  globalSamplingFraction = "1", #This is 1 for complete sampling (for homework, put in a fraction like 0.8)
  seed = sample(1:1000000, 1),
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
  outName = "classwork",
  ##numberOfChains = "2", #set to the number of CPUs
  numberOfChains = "8", #set to the number of CPUs
  deltaT = "0.01"
))

#### Import BAMM results into R ####
edata <- getEventData(tree, eventdata = paste0(resources_path,"classwork_event_data.txt"), burnin = 0.1)
#burnin -> we are throwing out the first 10 percent of the data b/c of the initial 
##getting into the topology with the initial chain runs
summary(edata)

#### Check quality of BAMM results ####
##Assess MCMC convergence
mcmcout <- read.csv(paste0(resources_path,"classwork_mcmc_out.txt"))
plot(mcmcout$logLik ~ mcmcout$generation)
#note: if you get a trace plot that doesn't look like a fuzzy catterpillar, and instead is jumping around a lot,
##you will need to increase the generation time so that it has more runs to get to the higher probabilities (higher on the topography)
###(i.e., blind robot climbing the mountain example)

#b/c these runs are all running in space, that means that they are going to be correlated one another
##b/c teleportation between chains isn't possible
###thus, "effectivesize" of > 200 is good to assess if the MCMC has run long enough and account for the inherent relatedness of the chains (paths)

#in the methods: "I gave it this control file and ran it for x number of generations and then use bamtools to confirm that the algorithm ran for long enought and i looked at the trace plot and there wasn't a lot of deviation from the expected value, and the effective size 
##was greater than 200, therefore I can use BAMM to analyze the rate shifts of the tree

##Test for convergence of the MCMC chains with the coda package
burnstart <- floor(0.1 * nrow(mcmcout)) # discard the first 10% of samples as burnin
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

#### Analyze BAMM results with BAMMtools ####
##Analysis of rate shifts
computeBayesFactors(paste0(resources_path,"classwork_mcmc_out.txt"), expectedNumberOfShifts = 1, burnin = 0.1)
#these give the likelihoods as well as thee plot of thee prior and posterior probabiltiese
##to be able to say if there has been evidence of there being 1 rate shift
### i.e., 0:1, the difference between the prior mode of 0 rate shifts and the posterior mode of 1 rate shifts is additional evidence for a model with a single rate shift
plotPrior(paste0(resources_path,"classwork_mcmc_out.txt"), expectedNumberOfShifts = 1)


##Analysis of rates
#phylorate plot: displays mean, model-averaged diversification rates 
#(i.e., speciation "s", extinction "e", or net diversification "netdiv") on branches of your tree with colors
s <- plot.bammdata(edata, spex = "s", labels = T, font = 3, cex = 0.5)
title(main = "Mean speciation rate", sub = "time before present")
addBAMMlegend(s, location = "topleft", nTicks = 1)
axisPhylo()
#this figure doesn't specifically tell us about the number of rate shifts, regardless of it being a pretty figure
##don't just look at this plot and say that there is a rate shift
###SO we will generate figures that have the actual number of rate shifts on our trees

#figure out the number of distinct rate shift configurations
css <- credibleShiftSet(edata, expectedNumberOfShifts = 1, threshold = 5, set.limit = 0.95)
css$number.distinct #this is number of distinct shift configurations in the data
summary(css)

sss <- plot.credibleshiftset(css, border = F)

#we just viewed a summary of the distinct rate shift configurations with the highest probabilities
#now lets find and plot the single best shift configuration (i.e, the one we see most often in the posterior distribution)
plot.new()
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts = 1)
ss <- plot.bammdata(best, labels = T, font = 3, cex = 0.7)
title(main = "Best shift configuration", sub = "time before present")
addBAMMlegend(ss, location = "topleft", nTicks = 1)
addBAMMshifts(best, cex = 3, pch = 1)
axisPhylo()
#you could say (this rate shift occurs in 86% of the posterior distribution (so there is something about certainty), meaning that there is 4% chance there is a rate shift somewhere else)
#whereas MEDUSA wouldnt tell us about alternative options of rate shifts instead of just saying there is a good chance there is a rate shift here
#MEDUSA (we can throw a rate shift down here, but we might not be entirely sure that there is a rate shift right now)
#BAMM tells us something about certainty when we are talking about our clades


#This should correspond to the polt in the credible shift set (CSS) with the highest frequency
#but will look slightly different because it is just a single sample from the posterior whereas CSS
#is averaged across samples in the posterior

#next we will plot a tree with branch lengths scaled by the probabilty that they contain a raate shift
#the longer the branches, the greater probability that there was a rate shift somewhere along that branch
###THUS,...
#plot marginal probability as a function of branch lengths (think about it as the probability of the rate shift occuring (long branch = high chance that a rate shift occured, whereas a low bvranch would = a low chance that a rate shift occured))
##^^^USE THIS INFO FOR A FIGURE LABEL IN THE FINAL PAPERr
par(font = 1)
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, cex = 0.7)
title(main = "Marginal shift probability")
add.scale.bar(x = 0.5, y = 0.5, font = 1)
#the other one is still slightly long (but not too long) because there was a 14% chance whereas the 80% chance from Bayes is much longer

## Clade-specific evolutionary rates
#we can compute clade-specific marginal distributions of rates with getCladeRates()
#below we estimate an overall speciation rate and 90% credible interval, then rates and intervals for the T and t clades separately
global_rates <- getCladeRates(edata)
mean(global_rates$lambda)

#the speciation rate estimate is 0.19 new species per million years
quantile(global_rates$lambda, c(0.05, 0.95))

#there is 90% probability that the speciation rate of this clade is between 0.11 and 0.29
T_MRCA <- getMRCA(tree, tip = c("T1", "T10"))
T_rates <- getCladeRates(edata, node = T_MRCA)
mean(T_rates$lambda)

quantile(T_rates$lambda, c(0.05, 0.95))

t_MRCA <- getMRCA(tree, tip = c("t1", "t20"))
t_rates <- getCladeRates(edata, node = t_MRCA)
mean(t_rates$lambda)

quantile(t_rates$lambda, c(0.05, 0.95))

##Rate-through-time analysis
#we can plot speciation (default) of extinction rats through timee with plotRateThroughTime()
#This plot can display dynamics in speciation or extinction rates
#The black line is the mean speciation rate, and the grey area is the 90% credible interval for the speciation rate
par(font = 1)
plotRateThroughTime(edata, 
                    ratetype = "speciation",
                    avgCol = "black",
                    intervalCol = "gray80",
                    intervals = c(0.05, 0.95),
                    opacity = 1)
