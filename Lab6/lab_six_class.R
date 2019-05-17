#
# EEB C174 Week 6, Lab 6
# Comparative Biology and Macroevolution
# Diversification Methods
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())

#imports
library(phytools)
library(geiger)
library(ape)
library(apTreeshape)

cwd <- "~/Developer/EEB-C174-Labs/Lab6/"
setwd(cwd)

resources_path <- "resources/"
output_path <- "output/"

#calculate Colless' index
my_tree <- sim.bdtree(b = 2, d = 0.5, stop = "taxa", n = 100, extinct = FALSE)
my_tree <- drop.extinct(my_tree)

plot(my_tree, show.tip.label = F)

#Turn phylo into treeshape
my_treeshape <- as.treeshape(my_tree)

#This calculates the sum of the number of left - right
#(the top part of the equation)
my_colless <- colless(my_treeshape, norm = NULL)
my_colless

#this calculates the total score (bottom part of eqn)
ntips <- length(my_tree$tip.label)
total_score <- ((ntips-1)*(ntips-2))/2
total_score

my_colless/total_score

#NOTE:p-value (whether my tree is more balanced (less unbalanced) than predicted by the null model)
##out of your x number of simulations, the proportion given were at least just as balanced or more balanced 
###if you do method = "Greater" then that would say whether they were just as unbalanced (pecetinate) or more unbalanced
colless.test(tree = my_treeshape, model = "yule", alternative = "less", n.mc = 1000)

#Exploring using the gamma statistic
tree <- pbtree(n = 100)
ltt.plot(tree, log = "y") #log transform the y-axis
lines(c(0,1), c(log(2), log(100)), lty = "dashed", lwd = 2, col = "red")

#now with multiple realizations
trees <- pbtree(n = 100, nsim = 10, scale = 1)
ltt_plots <- ltt(trees)
lines(c(0,1), c(log(2), log(100)), lty = "dashed", lwd = 2, col = "red")

#If we give ltt multiple trees
#Then applies the ltt function to all of the trees
#It then returns a list of "ltt" objects

ltt_plots[[1]]

#now lets pull out the gamma-stat for ltt_plots[[1]]
ltt_plots[[1]]$gamma
n <- length(ltt_plots)

#create an empty vector
gamma_vector <- numeric(n)

#fill the vector with a for loop
for (i in 1:n)
{
  gamma_vector[i] <- ltt_plots[[i]]$gamma
}

#histogram of all the gamma values
hist(gamma_vector)

#NOTE: we would expect a distribution centered around 0 for the gamma statistic

#### Plotting tree and LTT on the same plot ####
#our make transparent function
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

ltt.plot(tree, log = "y")
plotTree(tree,
         color=makeTransparent("blue",alpha=50),
         ftype="off",
         add=TRUE,
         mar=par()$mar)

#lets look at what happens when you incompletely sample a tree
#sample 50% of the tree
incomplete_tree <- drop.random(phy = tree, n = 50)

ltt.plot(incomplete_tree, log = "y")
plotTree(incomplete_tree,
         color=makeTransparent("blue",alpha=50),
         ftype="off",
         add=TRUE,
         mar=par()$mar)

##What would happen to thee gamma-stat if we incompletely sampled a tree?
#Another function to calculate the gamma statistic
gammaStat(tree)
gammaStat(incomplete_tree)

#### Diversification in Homalopsid Snakes ####
snake_tree <- read.tree(paste0(resources_path,"homalops.tre"))
snake_tree

snake_gamma <- gammaStat(snake_tree)
snake_gamma

#Get the age of the clade by finding the branching time at the root
age <- branching.times(snake_tree)[1]

#Manually specify total richness
richness <- 34 #species in this group
snakebirth = (log(richness) - log(2))/age
snakebirth #speciation rate estimate for all species in the tree

num_simulations <- 1000 #number of simulations
gamma_null <- numeric(num_simulations)
#gamma_null will hold the simulated gamma values
#for the trees that have been pruned down
for (i in 1:num_simulations)
{
  sim_tree <- sim.bdtree(b = snakebirth, d = 0, stop = "taxa", n = 34)
  prune <- drop.random(sim_tree, 13)
  #here we drop 13 species randomly from the tree
  gamma_null[i] <- gammaStat(prune)
  #this storesr the gamma values from the pruneed trees
}

hist(gamma_null, xlim = c(-3.5, 3.5))
arrows(x0 = snake_gamma, y0 = 100,
       x1 = snake_gamma, y1 = 0,
       col = "red", lwd = 2, xlab = "null gammas",
       main = "Incomplete Sampling")

mean(gamma_null) #find the mean gamma of nulls

#Which of the null values are smaller (more negative) than the data?
smallerNull <- gamma_null <= snake_gamma
smallerNull

#How many TRUEs are there?
count <- sum(smallerNull)
count

#We can now calculate a p-value for our empical gamma-stat using the simulations
mccr_pval <- (count+1)/(num_simulations+1)
mccr_pval
#the probability of us observing a gamma-stat of "snake_gamma" is the "p-value"

#CAN WE REJECT THE HYPOTHESIS OF CONSTANT RATE: WE CAN AND WE CAN SAY THAT THESE SNAKES ARE DECREASING IN RATE THROUGH TIME

#### Magallon and Sanderson method ####
#QUESTION: are hyenas especially species poor?
carnivora <- read.tree(paste0(resources_path,"full_carnivoran_tree.tre"))
carnivora

n <- length(carnivora$tip.label)
t <- max(nodeHeights(carnivora)) # total age of the tree

r_stem_age <- log(n)/t
r_crown_age <- (log(n)-log(2))/t

r_stem_age

#these are the actual functions built in to the package to do it for you
bd.ms(time = t, n = n, crown = FALSE, epsilon = 0)
r_crown_age

bd.ms(time = t, n = n, crown = TRUE, epsilon = 0)

net_r_e0 <- r_crown_age
net_r_e9 <- bd.ms(time = t, n = n, crown = TRUE, epsilon = 0.9)

#now we need to find the stem and crown age for the Hyenas
crown_node <- getMRCA(phy = carnivora, tip = c("Proteles_cristata", "Crocuta_crocuta"))
stem_node <- getMRCA(phy = carnivora, tip = c("Proteles_cristata", "Salanoia_concolor"))

#nodeheight computes the height above the root for a node
#Thereforer, to get the age, you need to subtract it from the total height
crown_node_age <- t - nodeheight(tree = carnivora, node = crown_node)
stem_node_age <- t - nodeheight(tree = carnivora, node = stem_node)

#prob of getting a clade as big as n (exceptionally diverse) (or 1-for especially depaupaurate)
crown.p(time = crown_node_age, n = 4, r = net_r_e0, epsilon = 0)

#prob of getting a clade as small as n (exceptionally depauparate)
1 - crown.p(time = crown_node_age, n = 4, r = net_r_e0, epsilon = 0)

#prob of gtting a clade as big as n
stem.p(time = stem_node_age, n = 4, r = net_r_e9, epsilon = 0.9)

#prob of gteting a clade as small as n
1 - stem.p(time = stem_node_age, n = 4, r = net_r_e9, epsilon = 0.9)

#Are hyenas especially species poor? --> NO

#we can try visualizing the Magallon and Sanderson method
crown_bounds <- function(max_age = 100, r = 0.1, epsilon = 0, CI = 0.95)
{
  times <- 1:max_age
  lower_bound <- numeric(max_age)
  upper_bound <- numeric(max_age)
  for (i in 1:max_age)
  {
    lower_bound[i] <- crown.limits(time = times[i],
                                   r = r,
                                   epsilon = epsilon,
                                   CI = CI)[1]
    upper_bound[i] <- crown.limits(time = times[i],
                                   r = r,
                                   epsilon = epsilon,
                                   CI = CI)[2]
  }
  crown_bounds <- data.frame(lower_bound, upper_bound)
  return(crown_bounds)
}

stem_bounds <- function(max_age = 100, r = 0.1, epsilon = 0, CI = 0.95)
{
  times <- 1:max_age
  lower_bound <- numeric(max_age)
  upper_bound <- numeric(max_age)
  for (i in 1:max_age)
  {
    lower_bound[i] <- stem.limits(time = times[i],
                                  r = r,
                                  epsilon = epsilon,
                                  CI = CI)[1]
    upper_bound[i] <- stem.limits(time = times[i],
                                  r = r,
                                  epsilon = epsilon,
                                  CI = CI)[2]
  }
  stem_bounds <- data.frame(lower_bound, upper_bound)
  return(stem_bounds)
}

crown_no_extinction <- crown_bounds(max_age = t, r = net_r_e0, epsilon = 0)
crown_high_extinction <- crown_bounds(max_age = t, r = net_r_e9, epsilon = 0.9)

stem_no_extinction <- stem_bounds(max_age = t, r = net_r_e0, epsilon = 0) 
stem_high_extinction <- stem_bounds(max_age = t, r = net_r_e9, epsilon = 0.9)

# one row, two columns
par(mfrow = c(1,2))

plot(x = crown_node_age, y = 4,
     xlim = c(0, 10), ylim = c(0, 12), main = "Crown Age for Hyenas", xlab = "MYR", ylab = "Ntaxa")
lines(1:t, crown_no_extinction$lower_bound, col = "black") 
lines(1:t, crown_no_extinction$upper_bound, col = "black") 
lines(1:t, crown_high_extinction$lower_bound, col = "red") 
lines(1:t, crown_high_extinction$upper_bound_bound, col = "red")
legend(x = -0.5, y = 12.5, legend=c("epsilon = 0", "epsilon = 0.9"),
       col=c("black", "red"), lty=1, cex=0.8, bty = "n")

plot(x = stem_node_age, y = 4,
     xlim = c(0, 28), ylim = c(0, 40), main = "Stem Age for Hyenas", xlab = "MYR", ylab = "Ntaxa")
lines(1:t, stem_no_extinction$lower_bound, col = "black")
lines(1:t, stem_no_extinction$upper_bound, col = "black")
lines(1:t, stem_high_extinction$lower_bound, col = "red")
lines(1:t, stem_high_extinction$upper_bound_bound, col = "red")
legend(x = -0.5, y = 40.5, legend=c("epsilon = 0", "epsilon = 0.9"),
       col=c("black", "red"), lty=1, cex=0.8, bty = "n")


#### MEDUSA ####
## Modeling evolutionary diversification using stepwise AIC ##
#the AIC threshold is the first number given to thrown down a rate shift

run1 <- medusa(phy = carnivora)
run1

# r is the per-lineeage net diversification rate
# epsilon is the relative extinctino rates
shift_nodes <- run1$model$split.at # get nodes with rate shifts
shift_nodes


# plot medusa rate shifts
plot(run1, show.tip.label = F)

low_div_tree <- extract.clade(phy = carnivora, node = shift_nodes[2])
plot(low_div_tree, cex = 0.8)
axisPhylo()

#NOTE:There appreas to be a rate shift in the family Procyonidae


#### Using richness data for incomplete trees ####
#simulate a 10 taxon tree
test_tree <- pbtree(b = 1, d = 0, nsim = 1, n = 10)
plot(test_tree)

##Where would we expect a rate shift?
sort(test_tree$tip.label)

taxon_richness = data.frame(taxon = sort(test_tree$tip.label), n.taxa = c(50, rep(1,9)))
taxon_richness

#we can now give the taxon_richness as an argument to richness in the medusa function
run2 <- medusa(phy = test_tree, richness = taxon_richness)

#plot medusa rate shifts
plot(run2, show.tip.label = T)





