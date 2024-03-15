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

?geom_cladelab
#### RADIAL FAMILY TREE ####
?rgb


fill1 = rgb(112/255, 77/255, 158/255)
fill2 = rgb(242/255, 136/255, 145/255)
fill3 = rgb(247/255, 208/255, 135/255)


fill1 = rgb(214/255, 2/255, 112/255)
fill2 = rgb(155/255, 79/255, 150/255)
fill3 = rgb(0/255, 56/255, 168/255)


#fill1 = rgb(250/255, 135/255, 117/255)
#fill2 = rgb(157/255, 2/255, 215/255)
#fill3 = rgb(0/255, 0/255, 255/255)

extend_val = 3
bar_sz = 2
txt_sz = 10
st_num = 2
timetree = read.beast(file="boostrap_nexus.nex")
tt = ggtree(timetree, size = 2, ladderize = FALSE, layout="fan") +
  theme_transparent() + 
  # Amniota: 4 taxa + Latimeria
  # geom_cladelab(node=756, label="Amniota", align=TRUE, angle="auto") +
  
  # Gymnophiona: 15 taxa
  geom_cladelab(node=760, label="Gymnophiona", align=TRUE, angle="auto",
                barsize=bar_sz, fontsize=txt_sz) +
  # Anura: 11 taxa
  geom_cladelab(node=775, label="Anura", align=TRUE, angle="auto",
                barsize=bar_sz, fontsize=txt_sz) + 
  
  geom_cladelab(node=787, label="Cryptobranchidae", barcolor=fill1,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill1, fontsize=txt_sz) +
  geom_highlight(node=787, fill=fill1, extend=extend_val) +
  
  
  geom_cladelab(node=789, label="Hynobiidae", barcolor=fill2,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill2, fontsize=txt_sz) +  
  geom_highlight(node=789, fill=fill2, extend=extend_val) + 
  
  
  geom_cladelab(node=879, label="Sirenidae", barcolor=fill3,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill3, fontsize=txt_sz) +  
  geom_highlight(node=879, fill=fill3, extend=extend_val) + 
  
  
  geom_cladelab(node=887, label="Dicamptodontidae", barcolor=fill1,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill1, fontsize=txt_sz) +  
  geom_highlight(node=887, fill=fill1, extend=extend_val) + 
  
  
  geom_cladelab(node=890, label="Ambystomatidae", barcolor=fill2,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill2, fontsize=txt_sz) +  
  geom_highlight(node=890, fill=fill2, extend=extend_val) + 
  
  
  geom_cladelab(node=922, label="Salamandridae", barcolor=fill3,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill3, fontsize=txt_sz) +  
  geom_highlight(node=922, fill=fill3, extend=extend_val) + 
  
  
  geom_cladelab(node=1044, label="Proteidae", barcolor=fill1,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill1, fontsize=txt_sz) +  
  geom_highlight(node=1044, fill=fill1, extend=extend_val) + 
  
  
  geom_cladelab(node=1050, label="Rhyacotritonidae", barcolor=fill2,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill2, fontsize=txt_sz) +  
  geom_highlight(node=1050, fill=fill2, extend=extend_val) + 
  
  
  geom_cladelab(node=1054, label="Amphiumidae", barcolor=fill3,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill3, fontsize=txt_sz) +  
  geom_highlight(node=1054, fill=fill3, extend=extend_val) +   
  
  
  geom_cladelab(node=1056, label="Plethodontidae", barcolor=fill1,
                barsize=bar_sz, align=TRUE, angle="auto",
                fontface="bold", textcolor=fill1, fontsize=txt_sz) +  
  geom_highlight(node=1056, fill=fill1, extend=extend_val) + 
  
  
  xlim(NA, 500)

# at 753 tips! all accounted for
tt
ggsave("family_circle.png", plot=tt, device="png", limitsize = FALSE, 
       width=32, height=32, dpi=600, bg="transparent", scale=1)

?ggsave

#### DEFINE SUBTREES ####
cryptobranchoidea_spp = c(
  "Hynobius", "Pseudohynobius", "Liua",
  "Batrachuperus", "Salamandrella",
  "Pachyhynobius", "Paradactylodon",
  "Ranodon", "Onychodactylus", "Andrias",
  "Cryptobranchus",
  "Siren", "Pseudobranchus",
  "Ambystoma_mexicanum"
)

dicamp_amby_spp = c(
  "Dicamptodon", "Ambystoma", "Unisexual"
)
salamandridae = c(
  "Paramesotriton", "Pachytriton", "Laotriton",
  "Cynops", "Triturus", "Calotriton",
  "Neurergus", "Ommatotriton", "Lissotriton",
  "Ichthyosaura", "Euproctus", "Taricha", "Notophthalmus",
  "Tylotriton", "Echinotriton", "Pleurodeles", "Lyciasalamandra",
  "Salamandra", "Chioglossa", "Mertensiella", "Salamandrina"
)
prot_rhyach_amph = c(
  "Siren", "Pseudobranchus",
  "Ambystoma_mexicanum",
  "Proteus", "Necturus", "Rhyacotriton", "Amphiuma"
)

plethodontinae = c(
  "Aneides", "Desmognathus", "Ensatina",
  "Hydromantes", "Speleomantes", "Karsenia",
  "Phaeognathus", "Plethodon"
)

hemidactylinae_1 = c(
  "Stereochilus", "Pseudotriton", "Gyrinophilus",
  "Urspelerpes", "Eurycea", "Hemidactylium", "Batrachoseps",
  "Oedipina_taylori"
)
hemidactylinae_2 = c(
  "Cryptotriton", "Dendrotriton", "Nyctanolis",
  "Nototriton", "Bradytriton", "Oedipina",
  "Thorius_rex"
)
hemidactylinae_3 = c(
  "Chiropterotriton", "Thorius", "Parvimolge", "Aquiloeurycea",
  "Isthmura", "Pseudoeurycea", "Ixalotriton",
  "Bolitoglossa_cuna"
)
hemidactylinae_4 = c("Bolitoglossa")


#### CRYPTOBRANCHOIDEA TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = cryptobranchoidea_spp
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=1) +
  geom_tiplab(fontface=3, angle=0, offset = 3) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>80), hjust=-.3) + 
  coord_geo(
    xlim = c(-170, 85), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-120, -100.5, -66.0, -56.0, -33.9, -23.03,
                              -5.333, -2.58, 0),
                     labels=c("120", "100","66", "56", "33.9", "23",
                              "5.3", "2.6", ""),
                     position = "top") +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt
#ggsave("cryptobranchoidea.pdf", plot=tt, device="pdf", limitsize = FALSE, 
#       width=10, height=10 )


#### SIREN + PROTEIDAE + RHYAC + AMPHIUMA TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = prot_rhyach_amph
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)

# rename Ambystoma
timetree@phylo$tip.label[[7]] = "Figures 3 & 4"

tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=1) +
  geom_tiplab(fontface=3, angle=0, offset = 7) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>80), hjust=-.3) + 
  coord_geo(
    xlim = c(-180, 85), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-95, -66.0, -56.0, -33.9, -23.03,
                              -5.333, -2.58, 0),
                     labels=c("95", "66", "56", "33.9", "23",
                              "5.3", "2.6", ""),
                     position = "top") +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
#ggsave("salamandridae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
#       width=10, height=14)
tt

#### DICAMPTODON AMBY TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = dicamp_amby_spp
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=1) +
  geom_tiplab(fontface=3, angle=0, offset = 3) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>80), hjust=-.3) + 
  coord_geo(
    xlim = c(-95, 85), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-95, -66.0, -56.0, -33.9, -23.03,
                              -5.333, -2.58, 0),
                     labels=c("95", "66", "56", "33.9", "23",
                              "5.3", "2.6", ""),
                     position = "top") +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt
#ggsave("ambystoma_dicamptodon.pdf", plot=tt, device="pdf", limitsize = FALSE, 
#       width=10, height=10 )



#### SALAMANDRIDAE TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = salamandridae
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=1) +
  geom_tiplab(fontface=3, angle=0, offset = 7) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>80), hjust=-.3) + 
  coord_geo(
    xlim = c(-95, 85), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-95, -66.0, -56.0, -33.9, -23.03,
                              -5.333, -2.58, 0),
                     labels=c("95", "66", "56", "33.9", "23",
                              "5.3", "2.6", ""),
                     position = "top") +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt
#ggsave("salamandridae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
#       width=10, height=14)


#### PLETHODONTINAE TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = plethodontinae
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=1) +
  geom_tiplab(fontface=3, angle=0, offset = 7) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>80), hjust=-.3) + 
  coord_geo(
    xlim = c(-70, 85), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-95, -66.0, -56.0, -33.9, -23.03,
                              -5.333, -2.58, 0),
                     labels=c("95", "66", "56", "33.9", "23",
                              "5.3", "2.6", ""),
                     position = "top") +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt
#ggsave("salamandridae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
#       width=10, height=14)



#### HEMIDACTYLINAE 1 TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_1
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=1) +
  geom_tiplab(fontface=3, angle=0, offset = 7) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>80), hjust=-.3) + 
  coord_geo(
    xlim = c(-80, 85), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-95, -66.0, -56.0, -33.9, -23.03,
                              -5.333, -2.58, 0),
                     labels=c("95", "66", "56", "33.9", "23",
                              "5.3", "2.6", ""),
                     position = "top") +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt
#ggsave("salamandridae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
#       width=10, height=14)
