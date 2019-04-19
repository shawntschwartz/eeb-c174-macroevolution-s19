#
# EEB C174 Week 2, Lab 2
# Comparative Biology and Macroevolution
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(phytools)

#paths
wd_path <- "~/Developer/EEB-C174-Labs/Lab2"
setwd(wd_path)
resources_path <- "resources/"
output_path <- "output/"

#### Import Data ####
tree <- read.tree(paste0(resources_path,"lab2_tree_new.tre"))
plot(tree)

traits <- read.csv(paste0(resources_path,"lab2_data1.csv"), stringsAsFactors = F)
traits

#check if the tip.label in the phylo object have the same row order as traits$Taxon
cbind(traits$Taxon, tree$tip.label)
##reorder
rownames(traits) <- traits$Taxon
#NOTE: "match(tree$tip.label, rownames(traits))" gives us the indices of the matching locations
#NOTE: then, by putting it within the traits bracket, we are able to actually mutate the traits dataframe to be in ordered
#and then restoring it into traits with (traits <- ...)
traits <- traits[match(tree$tip.label, rownames(traits)), ]
traits
#NOTE: check again to make sure that everything worked:
cbind(traits$Taxon, tree$tip.label)
##NOTE: all should be good now

#### Naive Analysis of Trait Data ####
#extract trait data from traits data.frame
trait1 <- traits$trait1
trait1

trait2 <- traits$trait2
trait2

#compute relative traits
##(e.g., brain size adjusted for body size)
###this is done easily by dividing columns
relative_trait1 <- trait1/trait2
relative_trait1

#plot trait data to visualize the data
plot(trait1, trait2)

#log transformations look at relative differences
##(e.g., 1kg increase for mouse vs. 1kg increase for elephant)
##mouse 1kg increase would be an order of magnitude increase compared to that for an elephant

#transform the trait data with natural log. Plot the ln(data)
ln_trait1 <- log(trait1)
ln_trait1

ln_trait2 <- log(trait2)
ln_trait2

plot(ln_trait1, ln_trait2)

#fit a linear model to the ln(data)
fit <- lm(ln_trait2 ~ ln_trait1)
plot(ln_trait1, ln_trait2, xlab = "ln(Trait1)", ylab = "ln(Trait2)")
abline(fit)

#### Question: What is the relationship between trait 1 and trait 2? ####
# There is a positive relationship between trait 1 and trait 2.
summary(fit)

#### Question: What can we conclude from the naive analysis about the relationship between trait 1 and trait 2? #### 
# There is a significant positive relationship (p = 0.003, Adj. R^2 = 0.29)
# Approx. 30% of the variance in trait 2 can be accounted by the variance in trait 1
# Slope, call: "fit" and take the model (beta) -> intercept
fit

#### Analysis with PICs ####
#from library(ape)
#"Phylogenetic Independent Contrasts"
#takes in a vector and a phylogeny
#NOTE: This is the non-naive phylogenetic analysis#
length(ln_trait1) #N => 25
PIC_ln_trait1 <- pic(ln_trait1, tree) #N-1 => 24 independent contrasts
PIC_ln_trait2 <- pic(ln_trait2, tree)
PIC_fit <- lm(PIC_ln_trait2 ~ PIC_ln_trait1 - 1) # "-1 " removes intercept term
plot(PIC_ln_trait1, PIC_ln_trait2,
     xlab = "Contrasts in ln(Trait1)",
     ylab = "Contrasts in ln(Trait2)")
abline(PIC_fit)

#### Question: What can we conclude about the relationship between the PICs for trait 1 and trait 2? ####
summary(PIC_fit)
PIC_fit
#There is now a negative relationship:

#### Question: Why do we fit the linear model without an intercept term? ####
##why did we remove the intercept term (-1)?:
###because it forces the line to go through (0,0)
#these magnitudes of points give us insight into how much evolution has taken place (indepdent of time)
#thus, at a value of 0 (i.e., 0 for the contrast -> such as when Taxa A and B have the same trait value at the tips: i.e., no evolution has occured)
#thus, when there is no evolution, you have to have evolution in another trait (which is why we remove the intercept point
#b/c we care about if there is a relationship between the indepedent differences between these trait values, not the actual trait values)
#(p = 0.03, Adj. R^2 = 0.14)
#independent contrasts are really powerful, but they are just relationshis and are meaningless outside the context of that independent contrasts (since we are talking about the relative differences in trait values, not the actual trait values)
#thus, you couldn't just combine the PICs from two different studies to make one combined trait value relationship, since they weren't calculated together because of their context dependence on the other points

#### Prepare Plots for a Report ####
par(mfrow = c(1, 2)) #Allows you to plot two figures side-by-side
##NOTE: reset by clicking the broom in the plot window

#Plotting naive analysis
plot(ln_trait1, ln_trait2,
     xlab = "ln(Trait1)",ylab = "ln(Trait2)",
     main = "Naive analysis")
abline(fit)

#Plotting Analysis using PICs
plot(PIC_ln_trait1, PIC_ln_trait2,
     #xlab = "Contrasts in ln(Trait1)", ylab = "Contrasts in ln(Trait2)",
     xlab = "PICs for ln(Trait1)", ylab = "PICs for ln(Trait2)",
     main = "Analysis using PICs")
abline(PIC_fit)

#### Understanding the Intuition ####
#we see that when we do the phylogenetic independent contrasts that there is a switch in the slopes
#thus, we see that this analysis is powerful enought to change relationships and still find significant relationships
#the negative relationship b/w the green group is captured by the phylogenetic contrasts (by their close evolutionary history)

#### Some Fancy Stuff You Can Do ####
ln_traits = cbind(ln_trait1, ln_trait2)
rownames(ln_traits) <- traits$Taxon
#superimposes phylogenetic tree ontop of naive trait data analysis
phylomorphospace(tree = tree, ln_traits, xlab = "ln(trait1)", ylab = "ln(trait2)")
#example method section: "i read this data into r and conducted and plotted linear relationships between the trait data using the package ape"
#give results: "give the p-value and adj. r^2
#discussion: "what you found and why do you think it is real?"
#citations: I used the bat traits data from this paper (cite) and the bat trees data from this paper (cite)
#you can find additional papers to cite
#note: there are more bats in the bat tree than there are in the bat traits
#there are more bats in the bat trees than in the bat traits, so you need to prune the tree dataset to match for the hw
