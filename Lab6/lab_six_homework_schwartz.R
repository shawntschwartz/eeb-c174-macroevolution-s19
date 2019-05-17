#
# Shawn Schwartz, 2019
# EEB C174 UCLA Spring 2019
# Lab 6 HW - Diversity Dynamics of Labrids
#

rm(list=ls())

#imports
library(phytools)
library(geiger)
library(ape)
library(apTreeshape)

cwd <- "~/Developer/EEB-C174-Labs/Lab6"
resources_path <- "resources/"
output_path <- "output/"

setwd(cwd)

labrid_tree <- read.tree(paste0(resources_path,"Labridae.tre"))
plot(labrid_tree, type = "fan", show.tip.label = F)
plot(labrid_tree, cex = 0.3)

#### Diversification Methods ####
#Turn phylo into tree shape
labrid_treeshape <- as.treeshape(labrid_tree)

## run colless test ##
labrid_colless_test <- colless.test(tree = labrid_treeshape, model = "yule", alternative = "less", n.mc = 1000)
labrid_colless_test

## gamma stat ##
makeTransparent <- function(someColor,alpha=10)
{
  newColor <- col2rgb(someColor)
  apply(newColor,2,function(curcoldata)
  {
    rgb(red=curcoldata[1],
        green=curcoldata[2],
        blue=curcoldata[3],
        alpha=alpha,
        maxColorValue=255)
  })
}
ltt.plot(labrid_tree, log = "y")
plotTree(labrid_tree, color = makeTransparent("blue", alpha = 50), ftype = "off", add = TRUE, mar = par()$mar)

# sample 50% of the tree
half_num <- (length(labrid_tree$tip.label))/2
half_num <- half_num + .5
half_num
incomplete_labrid_tree <- drop.random(phy = labrid_tree, n = half_num)
ltt.plot(incomplete_labrid_tree, log = "y")
plotTree(incomplete_labrid_tree, color = makeTransparent("blue", alpha = 50), ftype = "off", add = TRUE, mar = par()$mar)

labrid_tree_gammastat <- gammaStat(labrid_tree)
incomplete_labrid_tree_gammastat <- gammaStat(incomplete_labrid_tree)
labrid_tree_gammastat
incomplete_labrid_tree_gammastat

## Get the age of the clade by finding the branching time at the root.
age <- branching.times(labrid_tree)[1]

## Manually specify total richness
richness <- 600
labrid_birth <- (log(richness) - log(2)) / age
labrid_birth

num_simulations <- 1000
gamma_null <- numeric(num_simulations)
# gamma_null will hold the simulated gamma values
# for the trees that have been pruned down
for (i in 1:num_simulations)
{
  sim_labrid_tree <- sim.bdtree(b = labrid_birth, d = 0, stop = "taxa", n = richness)
  prune <- drop.random(sim_labrid_tree, 136) #here we drop 136 species randomly from the tree (i.e., ~40%)
  gamma_null[i] <- gammaStat(prune)
}

hist(gamma_null, xlim = c(-3.5, 3.5))
arrows(x0 = labrid_tree_gammastat, y0 = 100, x1 = labrid_tree_gammastat, y1 = 0, col = "red", lwd = 2, xlab = "Null Gammas", main = "Incomplete Sampling")

mean(gamma_null)
sd(gamma_null)


#which of the null values is smaller (more negative) than the data?
smallerNull <- gamma_null <= labrid_tree_gammastat
smallerNull
count <- sum(smallerNull)
count
mccr_pval <- (count + 1)/(num_simulations+1)
mccr_pval

## MEDUSA: Investigating Rate Shifts ##
#Question: Have different Labridae traits evolved at different rates?
run1 <- medusa(phy = labrid_tree)
run1

shift_nodes <- run1$model$split.at 
shift_nodes
plot(run1, show.tip.label = FALSE)
extracted_clade_1 <- extract.clade(phy = labrid_tree, node = shift_nodes[1])
extracted_clade_2 <- extract.clade(phy = labrid_tree, node = shift_nodes[2])
extracted_clade_3 <- extract.clade(phy = labrid_tree, node = shift_nodes[3])

plot(extracted_clade_1, cex = 0.2)
  axisPhylo()

plot(extracted_clade_2, cex = 0.6)
  axisPhylo()

plot(extracted_clade_3, cex = 0.3)
  axisPhylo()
