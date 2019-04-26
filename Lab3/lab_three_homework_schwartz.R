#
# Shawn Schwartz, 2019
# EEB C174 UCLA Spring 2019
# Lab 3 HW - Simulating Brownian Motion
#

#clean up workspace
rm(list=ls())

#includes
library(phytools)
library(dplyr)
library(car)
library(geiger)

#directories
wd_path <- "~/Developer/EEB-C174-Labs/Lab3"
setwd(wd_path)
output_path <- "output/"
resources_path <- "resources/"

#### 1. Read in Time-Calibrated Tree ####
tree <- read.tree(paste0(resources_path,"Labridae.tre"))
tree <- ladderize(tree)
pdf(paste0(output_path,"1_tc_tree_labridae.pdf"))
  plot(tree, cex = 0.25)
  axisPhylo()
  add.scale.bar()
dev.off()

#### 2. Simulate Brownian Motion ####
sig2 <- 0.01 # sigma^2
n_steps <- 1:100 # number of steps (time)

## simulate Brownian evolution on a tree with fastBM
#x <- fastBM(tree, sig2 = sig2, internal = TRUE)
x <- fastBM(tree, sig2 = sig2, a = 0)
## visualize Brownian evolution on a tree
pdf(paste0(output_path,"2_tc_fast_bm_a0.pdf"))
  phenogram(tree, x, spread.labels = TRUE, spread.cost = c(1, 0), fsize = 0.3)
dev.off()

x <- fastBM(tree, sig2 = sig2, internal = TRUE)
pdf(paste0(output_path,"2_tc_fast_bm_internal.pdf"))
  phenogram(tree, x, spread.labels = TRUE, spread.cost = c(1, 0), fsize = 0.3)
dev.off()

#### 3. Visualize Tree with Simulated Trait Values ####
pdf(paste0(output_path,"3_fast_bm_viz.pdf"))
  plotTree.wBars(tree, x, tip.labels = TRUE, fsize = 0.3)
dev.off()

#### 4. BMlk() ####
#### Fitting BM Models ####
BMlk <- function(C, inv.C, sigmasq, root.state, data) {
  N <- length(data); # the number of tips
  EX <- rep(root.state, N) # creates a vector of the expected trait value - which under V <- C * sigmasq; # multiply the entries in C by the BM rate
  V <- C * sigmasq;
  inv.V <- inv.C * sigmasq ^-1; # do the same for the inverted matrix using the inverse of the rate
  lnlNum<- -0.5*(data - EX) %*% inv.V %*% (data - EX) 
  lnlDen<- log(sqrt((2*pi)^N*det(V))) 
  L<-lnlNum-lnlDen
  return(L);
}

c <- vcvPhylo(tree, anc.nodes = FALSE)
inv.c <- solve(c) #inverse of phylogenetic variance-covariance matrix

data <- fastBM(tree, sig2 = sig2, a = 0)

estimate_1 <- BMlk(c, inv.c, sigmasq = sig2, root.state = 0, data = data) #true rate
estimate_2 <- BMlk(c, inv.c, sigmasq = 0.02, root.state = 0, data = data)
estimate_3 <- BMlk(c, inv.c, sigmasq = 0.03, root.state = 0, data = data)
estimate_4 <- BMlk(c, inv.c, sigmasq = 0.03, root.state = 0, data = data)
estimate_5 <- BMlk(c, inv.c, sigmasq = 0.05, root.state = 0, data = data)

estimates <- c(estimate_1, estimate_2, estimate_3, estimate_4, estimate_5)
estimates

#### 5. Visualize Likelihood Surface ###
vals <- numeric(100)
sigmas <- (1:100)/1000

for(i in 1:100)
{
  vals[i] <- BMlk(c, inv.c, sigmasq = sigmas[i], root.state = 0, data = data) 
}

vals[!is.finite(vals)] <- 0

pdf(paste0(output_path,"5_likelihood_surface.pdf"))
  plot(sigmas, vals, xlim = c(0.011, 0.079), type = "l", ylab = "ln(L)")
dev.off()

mod_estimates <- vals[vals < max(vals)]
best_estimate <- which.max(mod_estimates)
best_estimate
vals_estimate <- mod_estimates[best_estimate]
vals_estimate
best_estimate_sigma <- sigmas[best_estimate]
best_estimate_sigma

# check to see what the best estimate is (actual value compared to surface plot)
bm_fit <- fitContinuous(tree, data, model = "BM")
bm_fit
ML_z0 <- bm_fit$opt$z0
ML_sig2 <- bm_fit$opt$sigsq
BMlk(c, inv.c, sigmasq = ML_sig2, root.state = ML_z0, data = data)
