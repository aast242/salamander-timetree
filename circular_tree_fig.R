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
#### RADIAL STARTING TREE ####

fill1 = "burlywood4"
fill2 = "blue4"
extend_val = 3
bar_sz = 2
txt_sz = 5.5
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

  # at 31 tips
  # Cryptobranchoidea + Sirenidae: 99 tips
  geom_cladelab(node=81,
                label="Fig. 2: Cryptobranchoidea + Sirenidae",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill1, fontsize=txt_sz) + 
  
  geom_cladelab(node=786, label="", barcolor=fill1, barsize=bar_sz) +
  geom_cladelab(node=879, label="", barcolor=fill1, barsize=bar_sz) +  
  geom_highlight(node=786, fill=fill1, extend=extend_val) +
  geom_highlight(node=879, fill=fill1, extend=extend_val) + 
  
  # at 130 tips
  # Ambystomatidae + Dicamptodontidae: 37 tips
  geom_cladelab(node=149,
                label="Fig. 3: Amby. + Dicamp.",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill2, fontsize=txt_sz) + 
  geom_cladelab(node=886, label="",
                barcolor=fill2, barsize=bar_sz) +
  geom_highlight(node=886, fill=fill2, extend=extend_val) +
  
  # at 167 tips
  # Salamandridae: 122 tips
  geom_cladelab(node=228,
                label="Salamandridae: Fig. 4",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill1, fontsize=txt_sz) + 
  geom_cladelab(node=922, label="", barcolor=fill1, barsize=bar_sz) +
  geom_highlight(node=922, fill=fill1, extend=extend_val) +
  
  # at 289 tips
  # Proteidae, Rhyacotritonidae, Amphiumidae, Plethodontinae: 130 tips
  geom_cladelab(node=354,
                label="Pr. + R. + Amp. + Pleth.: Fig. 5",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill2, fontsize=txt_sz) +
  geom_cladelab(node=1044, label="", barcolor=fill2, barsize=bar_sz) +
  geom_cladelab(node=1050, label="", barcolor=fill2, barsize=bar_sz) +
  geom_cladelab(node=1054, label="", barcolor=fill2, barsize=bar_sz) +
  geom_cladelab(node=1057, label="", barcolor=fill2, barsize=bar_sz) +
  
  geom_highlight(node=1044, fill=fill2, extend=extend_val) +
  geom_highlight(node=1050, fill=fill2, extend=extend_val) +
  geom_highlight(node=1054, fill=fill2, extend=extend_val) +
  geom_highlight(node=1057, fill=fill2, extend=extend_val) +
  
  # at 419 tips
  # Spelerpini, Hemidactylium, Batrachoseps: 62 tips
  # Hemidactylinae I
  geom_cladelab(node=450,
                label="Hemidactylinae I: Fig. 6",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill1, fontsize=txt_sz) +
  geom_cladelab(node=1174, label="", barcolor=fill1, barsize=bar_sz) +
  geom_cladelab(node=459, label="", barcolor=fill1,
                extend=0.3, barsize=bar_sz) +
  geom_cladelab(node=1214, label="", barcolor=fill1, barsize=bar_sz) +
  
  geom_highlight(node=1174, fill=fill1, extend=extend_val) +
  geom_highlight(node=459, fill=fill1, extend=extend_val) +
  geom_highlight(node=1214, fill=fill1, extend=extend_val) +
  
  # at 481 tips
  # Hemidactylinae II: 64 tips
  geom_cladelab(node=513,
                label="Hemidactylinae II: Fig. 7",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill2, fontsize=txt_sz) +
  geom_cladelab(node=1236, label="", barcolor=fill2, barsize=bar_sz) +
  geom_highlight(node=1236, fill=fill2, extend=extend_val) +
  
  #at 545 tips
  # Hemidactylinae III: 93 tips
  geom_cladelab(node=592,
                label="Fig. 8: Hemidactylinae III",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill1, fontsize=txt_sz) +
  geom_cladelab(node=1344, label="", barcolor=fill1, barsize=bar_sz) +
  geom_cladelab(node=1300, label="", barcolor=fill1, barsize=bar_sz) +
  geom_highlight(node=1344, fill=fill1, extend=extend_val) +
  geom_highlight(node=1300, fill=fill1, extend=extend_val) +
  
  # at 638 tips
  # Hemidactylinae IV: 115 tips
  geom_cladelab(node=696,
                label="Fig. 9: Hemidactylinae IV",
                align=TRUE, angle="auto", barsize=0, barcolor=NA,
                fontface="bold", textcolor=fill2, fontsize=txt_sz) +
  geom_cladelab(node=1392, label="", barcolor=fill2, barsize=bar_sz) +
  geom_highlight(node=1392, fill=fill2, extend=extend_val) +
  xlim(NA, 500)
  
  # at 753 tips! all accounted for
tt
ggsave("big_circle.png", plot=tt, device="png", limitsize = FALSE, 
       width=20, height=20, dpi=600,bg="transparent")

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
  "Tylototriton", "Echinotriton", "Pleurodeles", "Lyciasalamandra",
  "Salamandra", "Chioglossa", "Mertensiella", "Salamandrina"
)
prot_rhyach_amph = c(
  "Proteus", "Necturus", "Rhyacotriton", "Amphiuma",
  "Aneides", "Desmognathus", "Ensatina",
  "Hydromantes", "Speleomantes", "Karsenia",
  "Phaeognathus", "Plethodon",
  "Hemidactylium_scutatum"
)

hemidactylinae_1 = c(
  "Stereochilus", "Pseudotriton", "Gyrinophilus",
  "Urspelerpes", "Eurycea", "Hemidactylium", "Batrachoseps",
  "Oedipina_taylori"
)
hemidactylinae_2 = c(
  "Cryptotriton", "Dendrotriton", "Nyctanolis",
  "Nototriton", "Bradytriton", "Oedipina",
  "Thorius_munificus"
)
hemidactylinae_3 = c(
  "Chiropterotriton", "Thorius", "Parvimolge", "Aquiloeurycea",
  "Isthmura", "Pseudoeurycea", "Ixalotriton",
  "Bolitoglossa_cuna"
)
hemidactylinae_4 = c("Bolitoglossa")


#### CRYPTOBRANCHOIDEA + SIRENIDAE TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = cryptobranchoidea_spp
epoch_labels = as.character(round(epochs$max_age, 1))
epoch_labels[[2]] = ""
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
# shorten the ambystoma branch
timetree@phylo$edge.length[[198]] = 20
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=0.75) +
  # annotate the path to other trees when the subset is not monophyletic
  geom_tiplab(aes(label=label, subset=(label!="Ambystoma mexicanum")),
              fontface=3, angle=0, offset = 0.5) +
  geom_tippoint(aes(subset=grepl("Ambystoma", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Figs. 3-9", subset=(label=="Ambystoma mexicanum")),
              fontface=1, offset = 1) +
  
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  coord_geo(
    xlim = c(-170, 50), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(0.75, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-450, -425, -400, -375, -350, -325,
                              -300, -275, -250, -225, -200, -175,
                              -150, -125, -100, -75, -50, -25, 0),
                     labels=c("450", "", "400", "", "350",
                             "", "300", "", "250", "", "200",
                             "","150","","100","","50","","0"),
                     #labels=c("", "", "400", "", "",
                    #          "", "300", "", "", "", "200",
                    #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(epoch_labels))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt$layers = c(geom_rect(aes(xmax = -161.5, 
                            xmin = -174.7, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -100.5, 
                            xmin = -145, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -56.0, 
                            xmin = -66.0, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -23.03, 
                            xmin = -33.9, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -2.58, 
                            xmin = -5.333, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 

tt
ggsave("cryptobranchoidea.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=11, height=14)
print()

#### DICAMPTODON AMBY TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = dicamp_amby_spp
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_tiplab(aes(label=label, subset=(label!="Unisexual ambystoma"))
              ,fontface=3, angle=0, offset = 0.5) +
  geom_tiplab(aes(label="Unisexual Ambystoma",
                  subset=(label=="Unisexual ambystoma"))
              ,fontface=1, angle=0, offset = 0.5) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  coord_geo(
    xlim = c(-87, 25), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt$layers = c(geom_rect(aes(xmax = -66, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -5.333, 
                            xmin = -23.03, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.58, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt
ggsave("ambystoma_dicamptodon.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=10 )
print()

#### SALAMANDRIDAE TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = salamandridae
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_tiplab(fontface=3, angle=0, offset = 0.5, size=3) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  coord_geo(
    xlim = c(-82, 20), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt$layers = c(geom_rect(aes(xmax = -66, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -5.333, 
                            xmin = -23.03, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.58, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt
ggsave("salamandridae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=14)
print()


#### PROTEIDAE + RHYAC + AMPHIUMA + PLETHODONTINAE TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = prot_rhyach_amph
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
# shorten the Hemidactylium branch
timetree@phylo$edge.length[[237]] = 10
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tt = ggtree(timetree, size = 1, ladderize = FALSE) +
  theme_tree2() +
  geom_tiplab(aes(label=label, subset=(label!="Hemidactylium scutatum")),
              fontface=3, angle=0, offset = 0.5, size=3) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  geom_tippoint(aes(subset=grepl("Hemidactylium", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Figs. 6-9", subset=(label=="Hemidactylium scutatum")),
              fontface=1, offset = 1, size=3) +
  #geom_text2(aes(label=node), hjust=-.3, size=2) + 
  coord_geo(
    xlim = c(-140, 32), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-140,-130,-120,-110,
                              -100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("140", "130", "120", "110",
                              "100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(epoch_labels))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt$layers = c(geom_rect(aes(xmax = -100.5, 
                            xmin = -145, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -56.0, 
                            xmin = -66.0, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -23.03, 
                            xmin = -33.9, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -2.58, 
                            xmin = -5.333, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)  
tt
ggsave("prot_rhyaco_amph_plethodontinae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=15)
print()

#### HEMIDACTYLINAE 1 TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_1
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[123] = 15
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  geom_tippoint(aes(subset=grepl("Oedipina", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Figs. 7-9", subset=(label=="Oedipina taylori")),
              fontface=1, offset = 1, size=3) +
  geom_tiplab(aes(label=label, subset=(label!="Oedipina taylori")),
              fontface=3, offset = 1, size=3) +
  coord_geo(
    xlim = c(-82, 20), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt) %>% rotate(104) %>% rotate(105)
tt$layers = c(geom_rect(aes(xmax = -66, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -5.333, 
                            xmin = -23.03, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.58, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt
ggsave("hemi_1.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print()

#### HEMIDACTYLINAE 2 TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_2
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[128] = 15
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  geom_tippoint(aes(subset=grepl("Thorius", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Figs. 8-9", subset=(label=="Thorius munificus")),
              fontface=1, offset = 0.5, size=3) +
  geom_tiplab(aes(label=label, subset=(label!="Thorius munificus")),
              fontface=3, offset = 0.25, size=3) +
  coord_geo(
    xlim = c(-80, 16), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt) %>% rotate(66) %>% rotate(67)
tt$layers = c(geom_rect(aes(xmax = -66, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -5.333, 
                            xmin = -23.03, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.58, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt
ggsave("hemi_2.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print()

#### HEMIDACTYLINAE 3 TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_3
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[99] = 15
#timetree@phylo$edge.length
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  geom_tippoint(aes(subset=grepl("Bolitoglossa", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Fig. 9", subset=(label=="Bolitoglossa cuna")),
              fontface=1, offset = 1, size=3) +
  geom_tiplab(aes(label=label, subset=(label!="Bolitoglossa cuna")),
              fontface=3, offset = 0.5, size=3) +
  coord_geo(
    xlim = c(-77, 20), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt) %>% rotate(139) %>% rotate(140)
tt$layers = c(geom_rect(aes(xmax = -66, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -5.333, 
                            xmin = -23.03, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.58, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt
ggsave("hemi_3.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print()

#### HEMIDACTYLINAE 4 TREE ####
timetree = read.beast(file="boostrap_nexus.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_4
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  geom_tiplab(aes(label=label), fontface=3, offset = 0.5, size=3) +
  coord_geo(
    xlim = c(-52, 14), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -5.333, 
                            xmin = -23.03, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt$layers = c(geom_rect(aes(xmax = 0, 
                            xmin = -2.58, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers)
tt
ggsave("hemi_4.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=15)
