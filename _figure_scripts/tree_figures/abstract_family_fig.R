library(ggpubr)
library(ape)
library(dplyr)
library(purrr)
library(patchwork)
library(ggtree)
library(tidyverse)
library(deeptime)
library(treeio)
library(tidyverse)
library(tidytree)
library(phylobase)
library(dispRity)
library(phytools)

#### RADIAL FAMILY TREE ####

# Desert Sunset
fill1 = rgb(112/255, 77/255, 158/255)
fill2 = rgb(242/255, 136/255, 145/255)
fill3 = rgb(247/255, 208/255, 135/255)

# Bi flag
fill1 = rgb(214/255, 2/255, 112/255)
fill2 = rgb(155/255, 79/255, 150/255)
fill3 = rgb(0/255, 56/255, 168/255)

# Hemidactylium circles
fill1 = rgb(183/255, 151/255, 98/255)
fill2 = rgb(0/255, 132/255, 111/255)
fill3 = rgb(79/255, 198/255, 1/255)
fill4 = rgb(255/255, 170/255, 146/255)
fill5 = rgb(163/255, 0/255, 89/255)
fill6 = rgb(44/255, 112/255, 0/255)
fill7 = rgb(255/255, 176/255, 72/255)
fill8 = rgb(236/255, 216/255, 21/255)
fill9 = rgb(255/255, 52/255, 255/255)
fill10 = rgb(10/255, 166/255, 216/255)

extend_val = 3
bar_sz = 2
txt_sz = 10
st_num = 2
timetree = read.beast(file="boostrap_nexus.nex")
tt = ggtree(timetree, size = 1.75, ladderize = FALSE, layout="fan") +
  theme_transparent() + 
  # Amniota: 4 taxa + Latimeria
  # geom_cladelab(node=756, label="Amniota", align=TRUE, angle="auto") +
  
  # Gymnophiona: 15 taxa
  #geom_cladelab(node=760, label="Gymnophiona", align=TRUE, angle="auto",
  #              barsize=bar_sz, fontsize=txt_sz) +
  # Anura: 11 taxa
  #geom_cladelab(node=775, label="Anura", align=TRUE, angle="auto",
  #              barsize=bar_sz, fontsize=txt_sz) + 
  
  #geom_cladelab(node=787, label="Cryptobranchidae", barcolor=fill1,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill1, fontsize=txt_sz) +  
  geom_cladelab(node=787, label="", barcolor=fill1,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill1, fontsize=txt_sz) +
  geom_highlight(node=787, fill=fill1, extend=extend_val) +
  
  
  #geom_cladelab(node=789, label="Hynobiidae", barcolor=fill2,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill2, fontsize=txt_sz) +  
  geom_cladelab(node=789, label="", barcolor=fill2,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill2, fontsize=txt_sz) +  
  geom_highlight(node=789, fill=fill2, extend=extend_val) + 
  
  
  #geom_cladelab(node=879, label="Sirenidae", barcolor=fill3,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill3, fontsize=txt_sz) +  
  geom_cladelab(node=879, label="", barcolor=fill3,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill3, fontsize=txt_sz) +  
  geom_highlight(node=879, fill=fill3, extend=extend_val) + 
  
  
  #geom_cladelab(node=887, label="Dicamptodontidae", barcolor=fill4,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill4, fontsize=txt_sz) +  
  geom_cladelab(node=887, label="", barcolor=fill4,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill4, fontsize=txt_sz) +  
  geom_highlight(node=887, fill=fill4, extend=extend_val) + 
  
  
  #geom_cladelab(node=890, label="Ambystomatidae", barcolor=fill5,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill5, fontsize=txt_sz) +  
  geom_cladelab(node=890, label="", barcolor=fill5,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill5, fontsize=txt_sz) +  
  geom_highlight(node=890, fill=fill5, extend=extend_val) + 
  
  
  #geom_cladelab(node=922, label="Salamandridae", barcolor=fill6,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill6, fontsize=txt_sz) +  
  geom_cladelab(node=922, label="", barcolor=fill6,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill6, fontsize=txt_sz) +  
  geom_highlight(node=922, fill=fill6, extend=extend_val) + 
  
  
  #geom_cladelab(node=1044, label="Proteidae", barcolor=fill7,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill7, fontsize=txt_sz) +  
  geom_cladelab(node=1044, label="", barcolor=fill7,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill7, fontsize=txt_sz) +  
  geom_highlight(node=1044, fill=fill7, extend=extend_val) + 
  
  
  #geom_cladelab(node=1050, label="Rhyacotritonidae", barcolor=fill8,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill8, fontsize=txt_sz) +  
  geom_cladelab(node=1050, label="", barcolor=fill8,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill8, fontsize=txt_sz) +  
  geom_highlight(node=1050, fill=fill8, extend=extend_val) + 
  
  
  #geom_cladelab(node=1054, label="Amphiumidae", barcolor=fill9,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill9, fontsize=txt_sz) +  
  geom_cladelab(node=1054, label="", barcolor=fill9,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill9, fontsize=txt_sz) +  
  geom_highlight(node=1054, fill=fill9, extend=extend_val) +   
  
  
  #geom_cladelab(node=1056, label="Plethodontidae", barcolor=fill10,
  #              barsize=bar_sz, align=TRUE, angle="auto",
  #              fontface="bold", textcolor=fill10, fontsize=txt_sz) +  
  geom_cladelab(node=1056, label="", barcolor=fill10,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill10, fontsize=txt_sz) +  
  geom_highlight(node=1056, fill=fill10, extend=extend_val) + 
  
  
  xlim(NA, 500)

# at 753 tips! all accounted for
tt
ggsave("abstract_circle.png", plot=tt, device="png", limitsize = FALSE, 
       width=32, height=32, dpi=600, bg="transparent", scale=1)
