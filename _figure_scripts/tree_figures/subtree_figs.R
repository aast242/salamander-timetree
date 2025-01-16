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

setwd("N:/research/wiens/new_submission_docs/supplementary_files/")

extend_val = 3
bar_sz = 2
txt_sz = 5.5
st_num = 2
tt_fp = "N:/research/wiens/new_submission_docs/supplementary_files/supplementary_files/S3_Dated_ML_tree.nex"
timetree = read.beast(file="N:/research/wiens/addressing_reviews/making_bipartition/figtree_reroot.nex")

#### EPOCH OFFSET num 1 ####
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

#### EPOCH OFFSET num 2 ####

tt$layers = c(geom_rect(aes(xmax = -174.7, 
                            xmin = -201.4, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -145, 
                            xmin = -161.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -66.0, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56.0, 
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
timetree = read.beast(file=tt_fp)
taxa_names = timetree@phylo$tip.label
target_tree = cryptobranchoidea_spp
epoch_labels = as.character(round(epochs$max_age, 1))
epoch_labels[[2]] = ""
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
# shorten the ambystoma branch
timetree@phylo$edge.length[[222]] = 13.13
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 10, size=0.75) +
  # annotate the path to other trees when the subset is not monophyletic
  geom_tiplab(aes(label=label, subset=(label!="Ambystoma mexicanum")),
              fontface=3, angle=0, offset = 0.5, size=3.33) +
  geom_tippoint(aes(subset=grepl("Ambystoma", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Figs. 3-9", subset=(label=="Ambystoma mexicanum")),
              fontface=1, offset = 1) +
  
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  
  # label families
  geom_text2(aes(label="Cryptobranchidae", subset=node==115, fontface = "bold"),
             nudge_x=-12, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Hynobiidae", subset=node==120, fontface = "bold"),
             nudge_x=-8, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Sirenidae", subset=node==219, fontface = "bold"),
             nudge_x=-7, nudge_y=0.75, color="black", size=3) + 
  
  # label sub-families
  geom_text2(aes(label="Onychodactylinae", subset=node==121),
             nudge_x=-11.5, nudge_y=-0.75, color="grey35", size=3) + 
  geom_text2(aes(label="Hynobiinae", subset=node==132),
             nudge_x=-8, nudge_y=-0.75, color="grey35", size=3) + 
  coord_geo(
    xlim = c(-180, 50), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(0.75, "lines")),
    lab = list(FALSE, TRUE), skip=c("Jurassic", "Quaternary"),
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-180, -170, -160, -150, -140, -130, -120,
                              -110, -100, -90, -80, -70, -60,
                              -50, -40, -30, -20, -10, 0),
                     labels=c("","", "", "150", "", "",
                             "", "", "100", "", "", "",
                             "","50","","","","","0"),
                     #labels=c("", "", "400", "", "",
                    #          "", "300", "", "", "", "200",
                    #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(epoch_labels))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt = revts(tt)
tt$layers = c(geom_rect(aes(xmax = -174.7, 
                            xmin = -201.4, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -145, 
                            xmin = -161.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -66.0, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56.0, 
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
ggsave("fig2_cryptobranchoidea.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=12, height=15)
print("")

#### SALAMANDRIDAE TREE ####
timetree = read.beast(file=tt_fp)
taxa_names = timetree@phylo$tip.label
target_tree = salamandridae
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_tiplab(fontface=3, angle=0, offset = 0.5, size=2.9) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  # label sub-families
  geom_text2(aes(label="Salamandrininae", subset=node==138),
             nudge_x=-71, nudge_y=-0.75, color="grey35", size=2.75) + 
  geom_text2(aes(label="Salamandrinae", subset=node==140),
             nudge_x=-4.25, nudge_y=-0.75, color="grey35", size=2.75) +   
  geom_text2(aes(label="Pleurodelinae", subset=node==154),
             nudge_x=-4, nudge_y=-0.75, color="grey35", size=2.75) + 
  coord_geo(
    xlim = c(-80, 20), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
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
ggsave("fig3_salamandridae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=16)
print("")


#### DICAMPTODON AMBY TREE ####
timetree = read.beast(file=tt_fp)
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
  geom_tiplab(aes(label="paste('Unisexual ', italic('Ambystoma'))",
                  subset=(label=="Unisexual ambystoma"))
              ,fontface=1, angle=0, offset = 0.5, parse=T) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  # label families
  geom_text2(aes(label="Dicamptodontidae", subset=node==33, fontface = "bold"),
             nudge_x=-8.5, nudge_y=0.33, color="black", size=3) + 
  geom_text2(aes(label="Ambystomatidae", subset=node==36, fontface = "bold"),
             nudge_x=-8, nudge_y=0.33, color="black", size=3) + 
  coord_geo(
    xlim = c(-91, 25), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
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
ggsave("fig4_ambystoma_dicamptodon.pdf", plot=tt, device="pdf",
       limitsize = FALSE, width=10, height=10 )
print("")


#### PROTEIDAE + RHYAC + AMPHIUMA + PLETHODONTINAE TREE ####
timetree = read.beast(file=tt_fp)
taxa_names = timetree@phylo$tip.label
target_tree = prot_rhyach_amph
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
# shorten the Hemidactylium branch
timetree@phylo$edge.length[[247]] = 5.504
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
  #geom_text2(aes(label=node), hjust=-.3) + 
  
  # label families
  geom_text2(aes(label="Proteidae", subset=node==140, fontface = "bold"),
             nudge_x=-39, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Rhyacotritonidae", subset=node==148, fontface = "bold"),
             nudge_x=-100, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Amphiumidae", subset=node==152, fontface = "bold"),
             nudge_x=-84.5, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Plethodontidae", subset=node==154, fontface = "bold"),
             nudge_x=-19.25, nudge_y=0.75, color="black", size=3) + 
  
  # label sub-families
  geom_text2(aes(label="Plethodontinae", subset=node==155),
             nudge_x=-19.5, nudge_y=-0.75, color="grey35", size=3) +   
  geom_text2(aes(label="Hemidactyliinae", subset=node==138),
             nudge_x=3, nudge_y=-1.15, color="grey35", size=3) +  
  coord_geo(
    xlim = c(-150, 32), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-150, -140,-130,-120,-110,
                              -100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("150","140", "130", "120", "110",
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
tt$layers = c(geom_rect(aes(xmax = -145, 
                            xmin = -161.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -66.0, 
                            xmin = -100.5, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -33.9, 
                            xmin = -56.0, 
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
ggsave("fig5_prot_rhyaco_amph_plethodontinae.pdf", plot=tt, device="pdf",
       limitsize = FALSE, width=11, height=16)
print("")

#### HEMIDACTYLINAE 1 TREE ####
timetree = read.beast(file=tt_fp)
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_1
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[127] = 3.289
tt = ggtree(timetree, size = 0.75, ladderize = FALSE) +
  theme_tree2() +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  
  # label sub-genera
  geom_text2(aes(label="Plethopsis", subset=node==109),
             nudge_x=-3.9, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") +   
  geom_text2(aes(label="Batrachoseps", subset=node==111),
             nudge_x=-5, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") +
  
  geom_tippoint(aes(subset=grepl("Oedipina", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Figs. 7-9", subset=(label=="Oedipina taylori")),
              fontface=1, offset = 1, size=3) +
  geom_tiplab(aes(label=label, subset=(label!="Oedipina taylori")),
              fontface=3, offset = 1, size=3) +
  coord_geo(
    xlim = c(-90, 20), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
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
tt = revts(tt) %>% ggtree::rotate(107) %>% ggtree::rotate(108)
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
ggsave("fig6_hemi_1.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print("")

#### HEMIDACTYLINAE 2 TREE ####
timetree = read.beast(file=tt_fp)
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_2
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[138] = 1.423
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
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
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
tt = revts(tt) %>% ggtree::rotate(71) %>% ggtree::rotate(72)
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
ggsave("fig7_hemi_2.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print("")

#### HEMIDACTYLINAE 3 TREE ####
timetree = read.beast(file=tt_fp)
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_3
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[103] = 5.322
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
    xlim = c(-80, 20), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
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
tt = revts(tt) %>% ggtree::rotate(151) %>% ggtree::rotate(152)
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
ggsave("fig8_hemi_3.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print("")

#### HEMIDACTYLINAE 4 TREE ####
timetree = read.beast(file=tt_fp)
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
  
  # label subgenera
  geom_text2(aes(label="Oaxakia", subset=node==172),
             nudge_x=-2, nudge_y=-0.75, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="Magnadigita", subset=node==199),
             nudge_x=-2.5, nudge_y=-0.75, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="Pachymandra", subset=node==198),
             nudge_x=-2.8, nudge_y=-0.75, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="paste(italic('Mayamandra'), ' + ', italic('Nanotriton'))",
                 subset=node==179),
             nudge_x=-4, nudge_y=3, color="grey35", size=2.5, parse = T) + 
  geom_text2(aes(label="Bolitoglossa", subset=node==186),
             nudge_x=-2.5, nudge_y=-0.75, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="Eladinea", subset=node==119),
             nudge_x=-1.9, nudge_y=-0.75, color="grey35", size=2.5,
             fontface="italic") + 
  
  # Add line to connect mayamandra + nanotriton
  #geom_segment(aes(x = -46, y = 69.7, xend = -42.9, yend = 67.7),
  geom_segment(aes(x = -54, y = 57.9, xend = -51, yend = 55.7),
               color = "grey35", size=0.33,
               lineend="round", linejoin="round",
               arrow = arrow(length = unit(0.1,"cm"))) +
  
  geom_tiplab(aes(label=label), fontface=3, offset = 0.5, size=3) +
  coord_geo(
    xlim = c(-64, 16.5), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(0, 2.5), abbrv = FALSE, neg = TRUE,
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
ggsave("fig9_hemi_4.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=15)

