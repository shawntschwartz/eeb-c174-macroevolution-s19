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

#define dir paths
resources_path <- "~/Developer/EEB-C174-Labs/Lab1/resources/"
output_path <- "output/"

#set working directory
setwd("~/Developer/EEB-C174-Labs/Lab1")

#read in fish tree (Rabosky et al. 2019)
fish_megatree <- read.tree(paste0(resources_path,"fish.tre"))

#prep large tree for visualization
fish_megatree <- ladderize(fish_megatree)

#extracted Labridae from the megaphylogeny
labrid_MRCA <- getMRCA(fish_megatree, c("Halichoeres_leucurus","Lachnolaimus_maximus"))
labrid_extracted <- extract.clade(fish_megatree, node = labrid_MRCA)
labrid_extracted_data <- labrid_extracted$tip.label

##save tree and tip labels to files
write.tree(labrid_extracted, paste0(output_path,"Labrid_Tree"))
write.csv(labrid_extracted_data, file = paste0(output_path,"Labridae_Data.csv"))

pdf(paste0(output_path,"fig_1_extracted_clade.pdf"))
  plot(labrid_extracted, cex = 0.1) #labrid tree from fish megaphylogeny
dev.off()

##tree with tip labels and a scale bar
pdf(paste0(output_path,"fig_2_extracted_tiplabsno_axis_scale.pdf"))
  plot(labrid_extracted, cex = 0.1)
  #nodelabels(cex = 0.2)
  #axisPhylo()
  add.scale.bar()
dev.off()

##annotate at least 2 clades on the tree with some biologically relevant information
parrot_fish_clade_scrape_feed <- getMRCA(labrid_extracted, c("Scarus_longipinnis","Scarus_compressus"))
parrot_fish_clade_browse_feed <- getMRCA(labrid_extracted, c("Sparisoma_cretense", "Sparisoma_axillare"))
clade_label_labrid_tree <- ggtree(labrid_extracted) +
  geom_cladelabel(node = parrot_fish_clade_scrape_feed, label = "Scrapers", align = T, angle = 270, hjust = 'center', fontsize = 3) +
  geom_cladelabel(node = parrot_fish_clade_browse_feed, label = "Browsers", align = T, angle = 270, hjust = 'center', fontsize = 3) +
  geom_hilight(parrot_fish_clade_scrape_feed, fill = 'firebrick', alpha = 0.5) +
  geom_hilight(parrot_fish_clade_browse_feed, fill = 'darkgreen', alpha = 0.5) +
  ggtitle("Labridae Phylogeny")

pdf(paste0(output_path,"fig_3_extracted_clade_labels.pdf"))
  print(clade_label_labrid_tree)
dev.off()

##large labrid tree with phylopic
labrid_tree <- ggtree(labrid_extracted) +
  geom_image(image = paste0(resources_path,"labrid.png"), size = Inf, alpha = 0.5, color = 'steelblue') +
  ggtitle("Labridae Phylogeny")

pdf(paste0(output_path,"fig_4_labrid_tree.pdf"))
  print(labrid_tree)
dev.off()

pdf(paste0(output_path,"fig_zoomed_tree.pdf"))
  zoom(fish_megatree, labrid_extracted$tip.label, cex = 0.1)
dev.off()
