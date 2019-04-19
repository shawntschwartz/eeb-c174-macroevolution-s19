#
# EEB C174 Week 3, Lab 3
# Comparative Biology and Macroevolution
# Introduction to Brownian Motion
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(ape)
library(phytools)
library(car)

#### Simulating Browninan Motion ####
# Setting up our parameters
sig2 <- 0.01 #sigma^2
n_steps <- 1:100 #number of steps (time)

# Fill a vector with a bunch of small steps
## single run
vec <- numeric(100) #going to fill this with each time step
for (i in n_steps)
{
  small_step <- rnorm(n = 1, sd = sqrt(sig2))
  vec[i + 1] <- vec[i] + small_step
}
plot(n_steps, vec[1:100], type = "l", ylim = c(-2, 2), ylab = "x", xlab = "t")
# each different trajectory gives us a relaization (a run of brownian motion --> for 100 timesteps)
## which is why it is different every time

## multiple runs
# Plot an empty plot
plot(NULL, xlim = c(0, 100), ylim = c(-2, 2), ylab = "x", xlab = "t")
# Put the entire BM for loop inside another for loop
for (j in 1:100)
{
  vec <- numeric(100)
  for (i in n_steps)
  {
    small_step <- rnorm(n = 1, sd = sqrt(sig2))
    vec[i + 1] <- vec[i] + small_step
  }
  lines(n_steps, vec[1:100])
}

# QUESTION: What would happen if you made the rate parameter smaller?
## The variance in the last time step would be less (i.e., less variance in the trait value) -->
## a more constrained variance.
# Plot an empty plot
plot(NULL, xlim = c(0, 100), ylim = c(-2, 2), ylab = "x", xlab = "t")
sig2 <- 0.001
# Put the entire BM for loop inside another for loop
for (j in 1:99)
{
  vec <- numeric(100)
  for (i in n_steps)
  {
    small_step <- rnorm(n = 1, sd = sqrt(sig2))
    vec[i + 1] <- vec[i] + small_step
  }
  lines(n_steps, vec[1:100])
}

#### Simulating BM on a Tree ####
# Here we want to simulate a tree that is 100 myrs old and has 5 taxa
t <- 100 # total time
n <- 5 # total taxa
b <- (log(n) - log(2))/t
tree <- pbtree(b = b, n = n, t = t, type = "discrete", tip.label = LETTERS[5:1])
plot(tree)
axisPhylo()
edgelabels(round(tree$edge.length, 4), cex = 1)

## From Liam Revell
## simulate evolution along each edge
X <- lapply(tree$edge.length, function(x) c(0, cumsum(rnorm(n = x, sd = sqrt(sig2)))))

## reorder the edges of the tree for pre-order traversal
cw <- reorder(tree)
ll <- tree$edge.length + 1
for (i in 1:nrow(cw$edge)) {
  pp <- which(cw$edge[, 2] == cw$edge[i, 1]) 
  if (length(pp) > 0)
    X[[i]] <- X[[i]] + X[[pp]][ll[pp]] else X[[i]] <- X[[i]] + X[[1]][1]
}

## get the starting and ending points of each edge for plotting
H <- nodeHeights(tree)

## plot the simulation
plot(H[1, 1], X[[1]][1], ylim = range(X), xlim = range(H), xlab = "time", ylab = "phenotype")
for (i in 1:length(X)) lines(H[i, 1]:H[i, 2], X[[i]])

## add tip labels if desired
yy <- sapply(1:length(tree$tip.label), function(x, y) which(x == y), y = tree$edge[,2])
yy <- sapply(yy, function(x, y) y[[x]][length(y[[x]])], y = X)
text(x = max(H), y = yy, tree$tip.label, cex = 1.5)


##NOTE: once a node appears on the tree is when the BM motion line breaks off
##This simulation tracks the displacement of the phenotype through every point in time from the root to the tips.
##However, if we don;'t care about this history of displacement, we can easily simulate BM evolution on any tree following the steps given in lecture:
##For each branch:
###Get the length of the branch
###Draw from a normal distribution with mean at 0 and variance = sigma^2 * branch length . This is the amount of trait evolution that has occured on the branch.

#### fastBM ####
## simulate Brownina evolution on a tree with fastBM
x <- fastBM(tree, sig2 = sig2 , internal = TRUE)
## visualize Brownian evolution on a tree
phenogram(tree, x, spread.labels = TRUE, spread.cost = c(1, 0))

##QUESTION: Is this okay to do?
#Less computationally expensive, but doesn't allow us to see the trajectory along each node splitting event but end value (same outcome, shortcut method)

#### Asside on the CLT ####
#What if at each step, the traits were evolving under different distributions?
#Here we will simulate the traits evolving under a uniform distribution
hist(runif(10000), main = "Uniform Distribution", xlab = "x")

# Plot an empty plot
plot(NULL, xlim = c(0, 100), ylim = c(-10, 10), ylab = "x", xlab = "t")
final_vec <- numeric(100)

# Plot the entire BM for loop inside another for loop
for (j in 1:100)
{
  vec <- numeric(100)
  for (i in n_steps)
  {
    # Have each step be drawn from a uniform distribution instead of a normal distribution
    small_step <- runif(n = 1, min = -1, max = 1)
    vec[i + 1] <- vec[i] + small_step
  }
  lines(n_steps, vec[1:100])
  final_vec[j] <- vec[100] # Save last value of each realization
}
#This distribution looks normal:
hist(final_vec, xlab = "x", main = "Final Trait Distribution from Uniform")

#Exponential Distribution:
##everything will be increasing because we are adding a whole bunch of small positive numbers togteher
## Skewed right because we are adding a bunch of small positive numbers together
hist(rexp(1000, rate = 5), main = "Exponential Distribution", xlab = "x")

# Plot an empty plot
plot(NULL, xlim = c(0, 100), ylim = c(0, 30), ylab = "x",xlab = "t") 
final_vec <- numeric(100)

# Put the entire BM for loop inside another for loop
for (j in 1:100) 
{
  vec <- numeric(100)
  for (i in n_steps)
  {
    # Have each step be drawn from a uniform distribution instead of a normal distribution
    small_step <- rexp(1, rate = 5)
    vec[i+1] <- vec[i] + small_step
  }
  lines(n_steps, vec[1:100])
  final_vec[j] <- vec[100] # Save last value of each realization
}

hist(final_vec, breaks = 7, main = "Final Trait Distribution from Exponential", xlab = "x")

#### Trait Covariance and the Phylogenetic Variance - Covariance Matrix ####
# Simulate a tree 10 myrs old with 5 taxa
t <- 10 #total time
n <- 5 #total taxa
b <- (log(n) - log(2))/t
tree <- pbtree(b = b, n = n, t = t, type = "discrete", tip.label = LETTERS[5:1])
plot(tree)  
axisPhylo()  
# A and B would be expcted to covary more b/c they share a common ancestor.
# D and E would covary the most given Mark's version because they have the most recent common ancestor (in terms of evolutionary time)

## Simulate BM on the tree
X <- fastBM(tree, nsim = 500)
X <- t(X)
colnames(X) <- tree$tip.label
scatterplotMatrix(X, col = "darkblue")


# For a phylogenetic tree with n taxa:
## the phylogenetic variance - covariance matrix will be an n by n matrix
## Each row and column of the matrix correspons to one of the n taxa in the matrix
## The diagonals of the matrix are the total distances from that tip to the root of the tree (will all be the same). --> variance (they vary together)
## The off diagonals are the total branch lengths shared by a pair of taxa. (the off-diagonals will scale to how much they covary together)
plot(tree)
edgelabels(round(tree$edge.length, 4), cex = 1)
vcv.phylo(tree, anc.nodes = FALSE)
#the row and columns values (off-diagonal) are equivalent to the length of the shared branch
#0 means that they don't share a branch at all

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

# 1. simulate a 5 taxon tree
# 2. simulate BM on that tree with a known sigma^2 of 0.5
# 3. compare the likelihood of generating sigma^2 to a smaller sigma^2 of 0.005
t <- 100 #total time
n <- 5 #total time
b <- (log(n) - log(2))/t
tree <- pbtree(b = b, n = n, t = t)
# Simulate BM on the tree with a sig2 of 0.5 and z0 of -
data <- fastBM(tree, sig2 = 0.5, a = 0)
data

#calculate the phylogenetic vcm
c <- vcvPhylo(tree, anc.nodes = FALSE)

# Get inverse of phylogenetic variance-covariance matrix
inv.c <- solve(c)

estimate_1 <- BMlk(c, inv.c, sigmasq = 0.5, root.state = 0, data = data) # the true value for sigsquared
estimate_1

estimate_2 <- BMlk(c, inv.c, sigmasq = 0.005, root.state = 0, data = data)
estimate_2

##NOTE: ensure that both estimates are negative numbers (estimates < 0)

#### Question: Which parameter value fits better? ####
#estimate2 is further away from 0, so it has a lower likelihood value --> (i.e., it is less likely b/c of the lower log likelihood value)

# we want better likelihood values, but when you take the log of a likelihood value, it becomes negative
# b/c likelihoods are a value between 0 and 1

#thus, if you take e^estimates, then you get a bigger number that aren't negative values
#thus, we want to be closer to 0 with log-likelihood values

##Thus, the parameter with 0.5 is a better likelihood estimate parameter
#b/c this is a pretty good model is what it tells us compared to the one with the sigma that is really different (i.e., the sigmasquared == 0.005)
exp(estimate_1)
exp(estimate_2)

#we can also try other values of sigmasq and calculate the log likelhiood
sigma_0.3 <- BMlk(c, inv.c, sigmasq = 0.3, root.state = 0, data = data)
sigma_0.4 <- BMlk(c, inv.c, sigmasq = 0.4, root.state = 0, data = data)
sigma_0.5 <- BMlk(c, inv.c, sigmasq = 0.5, root.state = 0, data = data)
sigma_0.6 <- BMlk(c, inv.c, sigmasq = 0.6, root.state = 0, data = data)
sigma_0.7 <- BMlk(c, inv.c, sigmasq = 0.7, root.state = 0, data = data)

sigma_0.3 #this value has the highest log likelihood
sigma_0.4
sigma_0.5
sigma_0.6
sigma_0.7


#### Which value has the highest log likelihood? ####
##now if you want to know what the best likelihood value (i.e., the best parameters for our data), that is where the maximum likelihood comes in
###so we will visualize this with a max likelihood surface to see what the best value is
vals <- numeric(100)
sigmas <- (1:100)/100
for(i in 1:100)
{
  vals[i] <- BMlk(c, inv.c, sigmasq = sigmas[i], root.state = 0, data = data)
}
plot(sigmas, vals, type = "l", ylab = "ln(L)")

##soom in around the generating sigmasq value
plot(sigmas, vals, xlim = c(0.2, 0.8), type = "l", ylab = "ln(L)")
#the likelihood surface is quite flat around the generating value!

#### BONUS ####
#Is the log likelihood value we calculated actually the maximum log likelihood value?
library(geiger)
bm_fit <- fitContinuous(tree, data, model = "BM")
bm_fit
ML_z0 <- bm_fit$opt$z0
ML_sig2 <- bm_fit$opt$sigsq
BMlk(c, inv.c, sigmasq = ML_sig2, root.state = ML_z0, data = data)

