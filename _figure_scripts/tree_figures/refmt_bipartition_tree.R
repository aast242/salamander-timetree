library(ggpubr)
library(ape)
library(dplyr)
library(purrr)
library(ggtree)
library(tidyverse)
library(deeptime)
library(treeio)
library(tidyverse)
library(tidytree)
library(phylobase)
library(dispRity)
library(phytools)

curr_tree = read.newick("N:/research/wiens/addressing_reviews/making_bipartition/RAxML_bipartitions.adding_bootstraps.tre")
curr_tree = root.phylo(curr_tree, "Latimeria_chalumnae")
curr_tree = ladderize(curr_tree, right = FALSE)
write.nexus(curr_tree, file="N:/research/wiens/addressing_reviews/making_bipartition/bipartition_tree.nex")
?ladderize
