#
# EEB C174 Week 5, Lab 5
# Comparative Biology and Macroevolution
# Introduction to Brownian Motion
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(phytools)
library(geiger)
library(ape)

cwd <- "~/Developer/EEB-C174-Labs/Lab5"
resources_path <- "resources/"
output_path <- "output/"

setwd(cwd)

#### read in tree and data ####
simpleTree <- read.tree(paste0(resources_path,"simpleTree.tre"))
simpleTree <- ladderize(simpleTree)
plot(simpleTree)

simpleData <- read.csv(paste0(resources_path,"simpleData.csv"), stringsAsFactors = F)
simpleData

#cbind binds two vectors as columns
cbind(simpleData$Taxon, simpleTree$tip.label) #check trait and tip order

###Question: What is wrong with the order of the data?
rownames(simpleData) <- simpleData$Taxon #set rownames of traits
simpleData <- simpleData[match(simpleTree$tip.label, rownames(simpleData)),] #match traits and tips
cbind(simpleData$Taxon, simpleTree$tip.label) #verify trait and tips are in same order

State <- simpleData$State #create vector of trait states
names(State) <- simpleData$Taxon #setnames for trait states
State <- as.factor(State) #convert discrete trait into factor
State

##Note: We can see the strings are now factors with 2 Levels: Non-Venomous and Venomous

#### Fit models using fitDiscrete()
fitER <- fitDiscrete(phy = simpleTree, dat = State, model = "ER") #equal transition rates
fitSYM <- fitDiscrete(phy = simpleTree, dat = State, model = "SYM") #symmetric transition rates
fitARD <- fitDiscrete(phy = simpleTree, dat = State, model = "ARD") #all rates different

##Note: #sym: transitions between categories will not be the same/not allowed to vary, but they are symmetric in that back-and-forth rates are the same b/w a pair of traits, but not all individual traits necessarily have the same rate

fitER #take a look inside the returned list

###Question: How many parameters are in the equal-rates (ER) model?
Venom <- c(fitER$opt$aic, fitSYM$opt$aic, fitARD$opt$aic)
names(Venom) <- c("Equal Rates", "Symmetric Rates", "All Rates Different")
Venom

###Question: Why are the AIC Scores identical for Equal Rates and Symmetric Rates?
fitER
fitSYM

#### Marginal ancestral state reconstruction ####
# marginal ancestral state estimation for each internal node of the tree using maximum likelihood
marginal_ER_fit <- rerootingMethod(tree = simpleTree, x = State, model = "ER")
marginal_ER_fit

#### Plot estimated marginal ancestral states on the tree ####
plot(simpleTree, show.tip.label = F)
tiplabels(simpleTree$tip.label, adj = -0.5, frame = "none")
nodelabels(node = as.numeric(rownames(marginal_ER_fit$marginal.anc)), pie = marginal_ER_fit$marginal.anc, piecol = c("black", "red"), cex = 0.6)
tiplabels(pie = to.matrix(State, sort(unique(State))), piecol = c("black", "red"), cex = 0.3)


#### Stochastic Character Mapping ####
# Use AIC-selected model for stochastic character mapping. 
# Simulate and plot 100 character histories
mtrees <- make.simmap(tree = simpleTree, x = State, model = "ER", nsim = 100)

# get colors for the states
cols <- setNames(object = palette()[1:length(unique(State))], nm = sort(unique(State)))
cols
cols <- setNames(object = c("black", "red"), nm = unique(State)) # This does the same thing
cols

par(mfrow = c(10,10)) #set plot window to 10 rows by 10 columns
null <- sapply(X = mtrees, FUN = plotSimmap, colors = cols, lwd = 1, ftype = "off") #plot

# With an aggregate of stochastic character maps, we can estimate...
## 1) the number of changes of each type,
## 2) the proportion of time spent in each state,
## 3) and the posterior probabilities that each internal node is in each state
pd <- describe.simmap(tree = mtrees, plot = FALSE)
pd #this has important information about the ancestral state

# Plot posterior probability of states on nodes with legend
par(mfrow = c(1,1)) #reset graphing parameters
plot(pd)
add.simmap.legend(colors = cols, prompt = F, x = 0, y = 2, fsize = 0.8)

# Plot posterior probability of states on branches
densityMap(mtrees)



