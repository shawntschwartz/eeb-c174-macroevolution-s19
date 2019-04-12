#
# EEB C174 Week 1, Lab 1
# Comparative Biology and Macroevolution
# By Shawn T. Schwartz
# UCLA Spring Quarter 2019
#

rm(list=ls())
# call necessary libs
library(ape)
library(ggtree)
library(ggplot2)
library(ggimage)

setwd("~/Developer/EEB-C174-Labs/Lab1")

load("resources/Carnivoran_Tree")
Carnivoran_Tree
str(Carnivoran_Tree)

## Using ape package
Carnivoran_Tree <- ladderize(Carnivoran_Tree)
plot(Carnivoran_Tree)
plot(Carnivoran_Tree, cex = 0.2)
plot(Carnivoran_Tree, show.tip.label = F)
add.scale.bar()

plot(Carnivoran_Tree, show.tip.label = F, type = "unrooted")
plot(Carnivoran_Tree, show.tip.label = F, type = "fan")

par(mar = c(1,1,1,1))
zoom(Carnivoran_Tree, grep("Canidae", Carnivoran_Tree$Family), cex = 0.6)
zoom(Carnivoran_Tree, grep("Vulpes", Carnivoran_Tree$tip.label))

MRCA <- getMRCA(Carnivoran_Tree, c("Nyctereutes_procyonoides", "Pseudalopex_griseus"))
MRCA

Canidae_Tree <- extract.clade(Carnivoran_Tree, node = MRCA)
plot(Canidae_Tree, cex = 0.7)

plot(Canidae_Tree, cex = 0.7)
nodelabels()

plot(Canidae_Tree, cex = 0.7)
axisPhylo()

write.tree(Canidae_Tree, "Canidae_Tree")
Canidae_Data <- Canidae_Tree$tip.label
write.csv(Canidae_Data, file = "Canidae_Data.csv")


## using ggtree package
CT <- ggtree(Canidae_Tree) +
  geom_tiplab() +
  geom_treescale(x = 20, y = 1, offset = 1) +
  geom_text2(aes(subset = !isTip, label = node), hjust = -.3, size = 3)
print(CT)

CladeA <- getMRCA(Canidae_Tree, tip = c("Canis_adustus", "Pseudalopex_vetulus")) 
CladeA

CladeB <- getMRCA(Canidae_Tree, tip = c("Vulpes_vulpes","Vulpes_velox"))
CladeB

### hilight parts of the tree for each Clade
CT <- ggtree(Canidae_Tree) +
  geom_hilight(CladeA, fill = "darkgreen", alpha = 0.6) +
  geom_hilight(CladeB, fill = "steelblue", alpha = 0.6)
print(CT)

### annotate clades of the tree
CT <- ggtree(Canidae_Tree) +
  geom_cladelabel(node = CladeA, label = "Canini", align = T, angle = 270, hjust = 'center') +
  geom_cladelabel(node = CladeB, label = "Vulpini", align = T, angle = 270, hjust = 'center')
print(CT)

### add image to the tree
CT <- ggtree(Canidae_Tree) +
  geom_image(image = "resources/canidae.png",
  size = Inf, alpha = .5,
  color = 'firebrick') +
  ggtitle("Canidae Phylogeny")
print(CT)
