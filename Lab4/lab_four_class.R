#
# EEB C174 Week 4, Lab 4
# Comparative Biology and Macroevolution
# Introduction to Brownian Motion
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(phytools)
library(geiger)
library(phylolm)

cwd <- "~/Developer/EEB-C174-Labs/Lab4"
setwd(cwd)

#import tree
tree <- read.tree("resources/lab4_part1_tree.tre")
plot(tree, show.tip.label = F)

#### Simulate Data Under Several Different Models ####
# Simulate traits under a Brownian Motion Model
data1 <- rTrait(phy = tree, model = "BM", parameters = list(ancestral.state = 0, sigma2 = 0.1))

# Simulate traits under an Ornstein-Uhlenbeck Model
data2 <- rTrait(phy = tree, model = "OU", parameters = list(ancestral.state = 0, sigma2 = 0.1, optimal.value = 10, alpha = 0.001))

# Simulate traits under an Early Burst Model
data3 <- rTrait(phy = tree, model = "EB", parameters = list(ancestral.state = 0, sigma2 = 0.1, rate = -0.0001))

# Simulate traits under a White Noise Model
data4 <- rnorm(n = length(tree$tip.label), mean = 0, sd = exp(1))
names(data4) <- tree$tip.label #Need to give names to the data4 vector, since rTrait() already does that for us

# TO visualize the data
cbind(data1, data2, data3, data4)

#### Model Fitting ####
# Fit models to data1 (BM)
data1.bm <- fitContinuous(phy = tree, dat = data1, model = "BM") 
data1.ou <- fitContinuous(phy = tree, dat = data1, model = "OU") 
data1.eb <- fitContinuous(phy = tree, dat = data1, model = "EB") 
data1.wn <- fitContinuous(phy = tree, dat = data1, model = "white")

# Fit models to data2 (OU)
data2.bm <- fitContinuous(phy = tree, dat = data2, model = "BM") 
data2.ou <- fitContinuous(phy = tree, dat = data2, model = "OU") 
data2.eb <- fitContinuous(phy = tree, dat = data2, model = "EB") 
data2.wn <- fitContinuous(phy = tree, dat = data2, model = "white")

# Fit models to data3 (EB)
data3.bm <- fitContinuous(phy = tree, dat = data3, model = "BM") 
data3.ou <- fitContinuous(phy = tree, dat = data3, model = "OU") 
data3.eb <- fitContinuous(phy = tree, dat = data3, model = "EB") 
data3.wn <- fitContinuous(phy = tree, dat = data3, model = "white")

# Fit models to data4 (WN)
data4.bm <- fitContinuous(phy = tree, dat = data4, model = "BM") 
data4.ou <- fitContinuous(phy = tree, dat = data4, model = "OU") 
data4.eb <- fitContinuous(phy = tree, dat = data4, model = "EB") 
data4.wn <- fitContinuous(phy = tree, dat = data4, model = "white")

#### Question: How many params are in the 'BM' model? ####
#Answer: 2 (sigma square (rate), and z0 (root state))

#### Compare Results Using AIC ####
# AIC = 2k-2ln(L)
#Make a matrix of all the AIC scores
BM.Data <- c(data1.bm$opt$aic, data1.ou$opt$aic, data1.eb$opt$aic, data1.wn$opt$aic)
OU.Data <- c(data2.bm$opt$aic, data2.ou$opt$aic, data2.eb$opt$aic, data2.wn$opt$aic)
EB.Data <- c(data3.bm$opt$aic, data3.ou$opt$aic, data3.eb$opt$aic, data3.wn$opt$aic)
WN.Data <- c(data4.bm$opt$aic, data4.ou$opt$aic, data4.eb$opt$aic, data4.wn$opt$aic)
AICresults <- rbind(BM.Data, OU.Data, EB.Data, WN.Data)
colnames(AICresults) <- c("Brownian Motion", "Ornstein-Uhlenbeck", "Early Burst", "White Nosie")
rownames(AICresults) <- c("From BM", "From OU", "From EB", "From WN")
AICresults
#NOTE: According to this, the Brownian Motion (BM) model fits the data the best (however, it is not as strong as it could be: thus there is some support for the model but not fully)...
#The lowest AIC Score is the best fit. The score difference of at least 2-3 is what is needed to have a difference between the models.
#You would expect that the generating model would have the beste support (but isn't always the case; matching columns and rows for the models in the matrix that we constructed)
#In our matrix, we have a type 1 error (false positive) for BM for the Early Burst (EB) model.
#B/c we didn't give strong parameter differences for the models, we don't get much difference in the signals for the current model states

#### Introducing Disparity Through Time ####
#There are multiple ways to calculate disparity.
#The most popular approach is to use the average pairwise Euclidean distances between species.
#This measure relates to the variance and estimates the dispersion of the points in the space (Harmon et al. 2003)

#nsmi = number of null simulations that you want (the cloud that you see: they simulate these null distributions with the idea being that your dtt falling different from the null cloud (to see how different the disparity of your clade is compared to some null))
ddt1 <- dtt(phy = tree, data = data1, nsim = 100, index = c("avg.sq"), plot = TRUE)
title("Brownian Motion")

ddt2 <- dtt(phy = tree, data = data2, nsim = 100, index = c("avg.sq"), plot = TRUE)
title("Ornstein Uhlenbeck")

ddt3 <- dtt(phy = tree, data = data3, nsim = 100, index = c("avg.sq"), plot = TRUE)
title("Early Burst")

#### Inferring Mode of Trait Evolution ####
#Now try to infer the node of trait evolution on an unknown dataset
lab4_tree <- read.tree("resources/lab4_tree.tre")
lab4_data <- read.csv("resources/lab4_trait_data.csv", stringsAsFactors = F)
head(lab4_data)

lab4_tree

plot(lab4_tree, show.tip.label = F)

#Question: What is interesting about this tree?
#Answer: Something interesting about this tree is that it has some extinct fossils (tips that don't go all the way to the end)
#Thus,
# Remove the extinct tips
extant_only_lab4_tree <- drop.extinct(lab4_tree)
plot(extant_only_lab4_tree, show.tip.label = F)

# Pull out the trait value column and name it
trait_values <- lab4_data$trait_value
names(trait_values) <- lab4_data$X
head(trait_values)

# Extract the trait data for only the living species
extant_only_trait_values <- trait_values[intersect(names(trait_values), extant_only_lab4_tree$tip.label)]
length(extant_only_trait_values) #Should be 100 to match the number of extant species within the tree after pruning the extinct species from the tree

# Fit models to extant only tree
extant.bm <- fitContinuous(phy = extant_only_lab4_tree, dat = extant_only_trait_values, model = "BM")
extant.ou <- fitContinuous(phy = extant_only_lab4_tree, dat = extant_only_trait_values, model = "OU")
extant.eb <- fitContinuous(phy = extant_only_lab4_tree, dat = extant_only_trait_values, model = "EB")
extant.wn <- fitContinuous(phy = extant_only_lab4_tree, dat = extant_only_trait_values, model = "white")

# Fit models to total tree
total.bm <- fitContinuous(phy = lab4_tree, dat = trait_values, model = "BM")
total.ou <- fitContinuous(phy = lab4_tree, dat = trait_values, model = "OU")
total.eb <- fitContinuous(phy = lab4_tree, dat = trait_values, model = "EB")
total.wn <- fitContinuous(phy = lab4_tree, dat = trait_values, model = "white")

Extant.Data <- c(extant.bm$opt$aic, extant.ou$opt$aic, extant.eb$opt$aic, extant.wn$opt$aic)
Total.Data <- c(total.bm$opt$aic, total.ou$opt$aic, total.eb$opt$aic, total.wn$opt$aic)
AIC_results <- rbind(Extant.Data, Total.Data)
colnames(AIC_results) <- c("Brownian Motion", "Orstein Uhlenbeck", "Early Burst", "White Noise")
AIC_results

# Transform AIC values into Akaike weights, which can be interpreted as conditional probabilities for each model
#NOTE: a way of comparing AIC scores to each other (they have to add up to one:-> making them conditional probabilities)
Extant_aicw <- aicw(Extant.Data)
Total_aicw <- aicw(Total.Data)
par(mfrow = c(2,1))
barplot(Extant_aicw[,3], names.arg = c("Brownian Motion", "Ornstein Uhlenbeck", "Early Burst", "White Noise"), ylim = c(0,1), main = "Extant Only", cex.names = 0.7)
abline(h = 1, lty = "dashed") #makes a dashed line at height 1 (y = 1)

barplot(Total_aicw[,3], names.arg = c("Brownian Motion", "Ornstein Uhlenbeck", "Early Burst", "White Noise"), ylim = c(0,1), main = "Total", cex.names = 0.7)
abline(h = 1, lty = "dashed") #makes a dashed line at height 1 (y = 1)

#Question: How does our incorporation of fossil data influence of inference of trait evolution?
#Answer: When we look at only the living species, we don't have that good of support for any of the models, but it isn't really clear that there is good support for BM,
#however, when we put all of thee fossils in the tree (including all the tips), we see that there is an overwhelming support for the Early Burst (EB) model (which is what we expect since these data were simulated with a high early burst (-r) paramater in the model);
#but when we remove the fossil data, we don't really good that much support for the distinct mode of trait evolution like we can see with the fossil data (which gives rise to EB as the main mode of trait evolution)
