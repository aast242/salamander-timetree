library(ape)
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

#### Main text figure ####
timetree = read.beast(file="figtree_reroot.nex")

fam_spp = c("Latimeria_chalumnae", "Homo_sapiens",
            "Rhinatrema_bivittatum", "Ascaphus_montanus",
            "Andrias_davidianus",
            "Hynobius_retardatus", "Siren_intermedia", "Taricha_torosa",
            "Ambystoma_mexicanum", "Dicamptodon_copei", "Proteus_anguinus",
            "Rhyacotriton_olympicus", "Amphiuma_means", "Plethodon_cinereus")

# only include species that are in the species file
timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, fam_spp))
timetree@phylo = ladderize(timetree@phylo)

# rename things to families
family_names = c(
  "Plethodontidae", "Amphiumidae", "Rhyacotritonidae", "Proteidae",
  "Salamandridae","Ambystomatidae", "Dicamptodontidae",  "Sirenidae",
  "Hynobiidae", "Cryptobranchidae",  
  "Anura",
  "Gymnophiona",
  "Amniota", "Coelacanthiformes"
)
timetree@phylo$tip.label = rev(family_names)

outgroup_names = c(
  "Chrysemys picta", "Gallus gallus", "Anolis carolinensis",
  "Homo sapiens", "Latimeria chalumnae"
)

# ggtree method
tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 40, size=1) +
  geom_tiplab(aes(label=label),
              fontface=1, angle=0, offset = 0.5) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.25) +
  #geom_text2(aes(label=node), hjust=-.3) +
  
  geom_text2(aes(label="Tetrapoda", subset=node==16, fontface = "bold"),
             nudge_x=-23, nudge_y=-0.25, color="black") +
 
  #geom_text2(aes(label="Batrachia", subset=node==18),
  #           nudge_x=-55, nudge_y=0.1, color="grey40") +
  geom_text2(aes(label="Batrachia", subset=node==18, fontface = "bold"),
             nudge_x=-21, nudge_y=-0.25, color="black") +
  geom_text2(aes(label="Caudata", subset=node==19, fontface = "bold"),
             nudge_x=-20, nudge_y=-0.25, color="black") +  
  
  geom_text2(aes(label="Cryptobranchoidea", subset=node==20, fontface = "bold"),
             nudge_x=-40, nudge_y=0.25, color="black") +  
  geom_text2(aes(label="Salamandroidea", subset=node==22, fontface = "bold"),
             nudge_x=-34, nudge_y=-0.25, color="black") +
  
  
  coord_geo(
    xlim = c(-450, 70), ylim = c(0, Ntip(timetree) + 0.5),
    neg = TRUE, dat = list("periods"), abbrv = FALSE,
    center_end_labels = TRUE, pos = "bottom", size = 3,
    skip=c("Ordovician","Neogene", "Quaternary")
  ) +
  # Exact breaks
  #breaks=c(-440, -419.2, -358.9, -298.9, -251.90, -201.3,-145.0, -66, -23.03, -2.58, 0)
  scale_x_continuous(breaks=c(-450, -425, -400, -375, -350, -325,
                              -300, -275, -250, -225, -200, -175,
                              -150, -125, -100, -75, -50, -25, 0),
                     #labels=c("450", "", "400", "", "350",
                     #        "", "300", "", "250", "", "200",
                     #        "","150","","100","","50","","0"),
                     labels=c("", "", "400", "", "",
                              "", "300", "", "", "", "200",
                              "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(periods$max_age),
                                         labels = rev(round(periods$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))


tt = revts(tt) %>% ggtree::rotate(20) %>% ggtree::rotate(27) %>% ggtree::rotate(24)
tt$layers = c(geom_rect(aes(xmax = -443.8, 
                            xmin = -450, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -358.9, 
                            xmin = -419.2, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -251.902, 
                            xmin = -298.9, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -145.0, 
                            xmin = -201.4, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -23.030, 
                            xmin = -66.0, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.580, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt
ggsave("fig1_family_tree.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=10)


#### Supplementary figure w/ 95% conf ####
timetree = read.beast(file="supplementary_files/S5_95conf_ML_tree.nex")

fam_spp = c("Latimeria_chalumnae", "Homo_sapiens",
            "Rhinatrema_bivittatum", "Ascaphus_montanus",
            "Andrias_davidianus",
            "Hynobius_retardatus", "Siren_intermedia", "Taricha_torosa",
            "Ambystoma_mexicanum", "Dicamptodon_copei", "Proteus_anguinus",
            "Rhyacotriton_olympicus", "Amphiuma_means", "Plethodon_cinereus")

# only include species that are in the species file
timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, fam_spp))
timetree@phylo = ladderize(timetree@phylo)

# rename things to families
family_names = c(
  "Plethodontidae", "Amphiumidae","Rhyacotritonidae", "Proteidae",
  "Salamandridae", "Ambystomatidae", "Dicamptodontidae", 
  "Sirenidae",
  "Hynobiidae", "Cryptobranchidae",  
  "Anura",
  "Gymnophiona",
  "Amniota", "Coelacanthiformes"
)
timetree@phylo$tip.label = rev(family_names)

outgroup_names = c(
  "Chrysemys picta", "Gallus gallus", "Anolis carolinensis",
  "Homo sapiens", "Latimeria chalumnae"
)

# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 1, ladderize = FALSE) ) %>%
  ggtree::rotate(20) %>% ggtree::rotate(27) %>% ggtree::rotate(24)

minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))

tt = tree_img +
  geom_rootedge(rootedge = 40, size=1) +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  theme_tree2() +
  geom_tiplab(aes(label=label),
              fontface=1, angle=0, offset = 0.5) +
  #geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.25) +
  #geom_text2(aes(label=node), hjust=-.3) +
  
  geom_text2(aes(label="Tetrapoda", subset=node==16, fontface = "bold"),
             nudge_x=-23, nudge_y=-0.25, color="black") +
  
  #geom_text2(aes(label="Batrachia", subset=node==18),
  #           nudge_x=-55, nudge_y=0.1, color="grey40") +
  geom_text2(aes(label="Batrachia", subset=node==18, fontface = "bold"),
             nudge_x=-21, nudge_y=-0.25, color="black") +
  geom_text2(aes(label="Caudata", subset=node==19, fontface = "bold"),
             nudge_x=-20, nudge_y=-0.25, color="black") +  
  
  geom_text2(aes(label="Cryptobranchoidea", subset=node==20, fontface = "bold"),
             nudge_x=-40, nudge_y=0.25, color="black") +  
  geom_text2(aes(label="Salamandroidea", subset=node==22, fontface = "bold"),
             nudge_x=-34, nudge_y=-0.25, color="black") +
  
  
  coord_geo(
    xlim = c(-450, 70), ylim = c(0, Ntip(timetree) + 0.5),
    neg = TRUE, dat = list("periods"), abbrv = FALSE,
    center_end_labels = TRUE, pos = "bottom", size = 3,
    skip=c("Ordovician","Neogene", "Quaternary")
  ) +
  # Exact breaks
  #breaks=c(-440, -419.2, -358.9, -298.9, -251.90, -201.3,-145.0, -66, -23.03, -2.58, 0)
  scale_x_continuous(breaks=c(-450, -425, -400, -375, -350, -325,
                              -300, -275, -250, -225, -200, -175,
                              -150, -125, -100, -75, -50, -25, 0),
                     #labels=c("450", "", "400", "", "350",
                     #        "", "300", "", "250", "", "200",
                     #        "","150","","100","","50","","0"),
                     labels=c("", "", "400", "", "",
                              "", "300", "", "", "", "200",
                              "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(periods$max_age),
                                         labels = rev(round(periods$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))

tt$layers = c(geom_rect(aes(xmax = -443.8, 
                            xmin = -450, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -358.9, 
                            xmin = -419.2, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -251.902, 
                            xmin = -298.9, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -145.0, 
                            xmin = -201.4, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -23.030, 
                            xmin = -66.0, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.580, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt
ggsave("S1_95conf_family_tree.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=10 )
