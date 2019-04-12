#
# Shawn Schwartz, 2019
# EEB C174 UCLA Spring 2019
# Lab 1 HW
#

#clean up workspace
rm(list=ls())

#includes
library(phytools)
library(ggtree)
library(ggimage)
library(tidyr)

#define dir paths
resources_path <- "~/Developer/EEB-C174-Labs/Lab1/resources/"
output_path <- "output/"

#set working directory
setwd("~/Developer/EEB-C174-Labs/Lab1")

#read in fish tree (Rabosky et al. 2019)
fish_megatree <- read.tree(paste0(resources_path,"fish.tre"))

#prep large tree for visualization
fish_megatree <- ladderize(fish_megatree)

#get the rainbow wrasses (Genus => Coris)
##extracted clade from the megaphylogeny
Coris_MRCA <- getMRCA(fish_megatree, c("Coris_julis", "Coris_pictoides"))
extracted_clade_tree <- extract.clade(fish_megatree, node = Coris_MRCA)
extracted_clade_tree_data <- extracted_clade_tree$tip.label
write.tree(extracted_clade_tree, paste0(output_path,"Labrid_Tree"))
write.csv(extracted_clade_tree_data, file = paste0(output_path,"Labridae_Data.csv"))

pdf(paste0(output_path,"fig_1_extracted_clade.pdf"))
  plot(extracted_clade_tree, cex = 0.3) #labrid tree from fish megaphylogeny
dev.off()

##tree with tip labels and a scale bar
pdf(paste0(output_path,"fig_2_extracted_tiplabs_axis.pdf"))
  plot(extracted_clade_tree, cex = 0.4)
  nodelabels(cex = 0.3)
  axisPhylo()
dev.off()

##annotate at least 2 clades on the tree with some biologically relevant information

##large labrid tree with phylopic
labrid_tree <- ggtree(extracted_clade_tree) +
  geom_image(image = paste0(resources_path,"labrid.png"), size = Inf, alpha = 0.5, color = 'steelblue') +
  ggtitle("Labridae Phylogeny")

pdf(paste0(output_path,"fig_4_labrid_tree.pdf"))
  print(labrid_tree)
dev.off()



##zoomed tree to Coris
pdf(paste0(output_path,"fig_5extra_zoomed_tree.pdf"))
  zoom(fish_megatree, grep("Coris", fish_megatree$tip.label), cex = 0.6)
dev.off()

#test new zoom
pdf(paste0(output_path,"fig_newzoomextra_zoomed_tree.pdf"))
zoom(fish_megatree, extracted_clade_tree$tip.label, cex = 0.3)
dev.off()

zoom(fish_megatree, grep("Scarus", fish_megatree$tip.label), cex = 0.6)
scarus_MRCA <- getMRCA(fish_megatree, c("Scarus_longipinnis","Scarus_compressus"))
extracted_parrotfish_tree <- extract.clade(fish_megatree, node = scarus_MRCA)
plot(extracted_parrotfish_tree, cex = 0.6)

other_tree <- read.tree(paste0(output_path,"Labridae.tre"))
other_tree <- ladderize(other_tree)
plot(other_tree)

#other_tree_MRCA <- getMRCA(other_tree, c)
#extracted_other_tree <- extract.clade(other_tree, other_tree_MRCA)
scarus_other_MRCA <- getMRCA(fish_megatree, c("Halichoeres_leucurus","Lachnolaimus_maximus"))
new_extracted <- extract.clade(fish_megatree, node = scarus_other_MRCA)
plot(new_extracted, cex = 0.5)
#par(mar = c(1,1,1,1))

#par(mar=c(5.1,4.1,4.1,2.1))


#labrid_tree <- ggtree(extracted_clade_tree) +
 # geom_hilight(Coris_MRCA, fill = "steelblue", alpha = 0.5)
#print(labrid_tree)
