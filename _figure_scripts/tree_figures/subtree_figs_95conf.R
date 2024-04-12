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
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = cryptobranchoidea_spp
epoch_labels = as.character(round(epochs$max_age, 1))
epoch_labels[[2]] = ""

keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
# shorten the ambystoma branch
timetree@phylo$edge.length[[198]] = 11.9564
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE))

# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))

tt = tree_img +
  theme_tree2() +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  #geom_text2(aes(label=height_0.95_HPD, subset=(node>length(keep_spp))), hjust=0, size=2) + 
  # annotate the path to other trees when the subset is not monophyletic
  geom_tiplab(aes(label=label, subset=(label!="Ambystoma mexicanum")),
              fontface=3, angle=0, offset = 0.5) +
  geom_tippoint(aes(subset=grepl("Ambystoma", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Supplementary Figs. 3-9", subset=(label=="Ambystoma mexicanum")),
              fontface=1, offset = 1) +
  #geom_text2(aes(label=node), hjust=-.3) + 
  
  # label families
  geom_text2(aes(label="Cryptobranchidae", subset=node==103, fontface = "bold"),
             nudge_x=-14.25, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Hynobiidae", subset=node==105, fontface = "bold"),
             nudge_x=-9, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Sirenidae", subset=node==195, fontface = "bold"),
             nudge_x=-8, nudge_y=0.75, color="black", size=3) + 
  
  # label sub-families
  geom_text2(aes(label="Onychodactylinae", subset=node==184),
             nudge_x=-12.5, nudge_y=-0.75, color="grey35", size=3) + 
  geom_text2(aes(label="Hynobiinae", subset=node==106),
             nudge_x=-8.5, nudge_y=-0.75, color="grey35", size=3) +
  
  coord_geo(
    xlim = c(-200, 55), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(0.75, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-200, -190, -180,
                              -170, -160, -150, -140, -130,
                              -120, -110, -100, -90, -80, -70,
                              -60, -50, -40, -30, -20, -10, 0),
                     labels=c("200","","", "", "", "150", "",
                             "", "", "", "100", "", "",
                             "","","50","","","","","0"),
                     #labels=c("", "", "400", "", "",
                    #          "", "300", "", "", "", "200",
                    #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(epoch_labels))) +
  theme(axis.line = element_blank())
tt$layers = c(geom_rect(aes(xmax = -174.7, 
                            xmin = -201.4, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -145.0, 
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

ggsave("S2_95conf_cryptobranchoidea.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=11, height=14)
print("")

#### SALAMANDRIDAE TREE ####
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = salamandridae
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)

# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE))
minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))

tt = tree_img + 
  theme_tree2() +
  geom_tiplab(fontface=3, angle=0, offset = 0.5, size=3) +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) + 
  #geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  # label sub-families
  geom_text2(aes(label="Salamandrininae", subset=node==124),
             nudge_x=-75, nudge_y=-0.75, color="grey35", size=2.5) + 
  geom_text2(aes(label="Salamandrinae", subset=node==126),
             nudge_x=-4.35, nudge_y=-0.75, color="grey35", size=2.5) +   
  geom_text2(aes(label="Pleurodelinae", subset=node==140),
             nudge_x=-5, nudge_y=-0.75, color="grey35", size=2.5) + 
  coord_geo(
    xlim = c(-100.5, 23), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-110,-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("110","100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
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
ggsave("S3_95conf_salamandridae.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=14)
print("")



#### DICAMPTODON AMBY TREE ####
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = dicamp_amby_spp
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE))
minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))
tt = tree_img + 
  theme_tree2() +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  geom_tiplab(aes(label=label, subset=(label!="Unisexual ambystoma"))
              ,fontface=3, angle=0, offset = 0.5) +
  geom_tiplab(aes(label="paste('Unisexual ', italic('Ambystoma'))",
                  subset=(label=="Unisexual ambystoma"))
              ,fontface=1, angle=0, offset = 0.5, parse=T) +
  #geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  # label families
  geom_text2(aes(label="Dicamptodontidae", subset=node==39, fontface = "bold"),
             nudge_x=-8.5, nudge_y=0.40, color="black", size=3) + 
  geom_text2(aes(label="Ambystomatidae", subset=node==42, fontface = "bold"),
             nudge_x=-8, nudge_y=0.40, color="black", size=3) + 
  coord_geo(
    xlim = c(-110, 28), ylim = c(0, Ntip(timetree) + 0.5),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE), skip=c("Quaternary" ,"Cretaceous"),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-110, -100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("110", "100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
tt$layers = c(geom_rect(aes(xmax = -100.5, 
                            xmin = -145, 
                            ymin = 0, 
                            ymax = Inf), alpha = 0.6, fill="grey"), tt$layers) 
tt$layers = c(geom_rect(aes(xmax = -56, 
                            xmin = -66, 
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
ggsave("S4_95conf_ambystoma_dicamptodon.pdf", plot=tt, device="pdf",
       limitsize = FALSE, width=10, height=10)
print("")

#### PROTEIDAE + RHYAC + AMPHIUMA + PLETHODONTINAE TREE ####
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = prot_rhyach_amph
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)
# shorten the Hemidactylium branch
timetree@phylo$edge.length[[237]] = 3.5117
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE))
minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))
tt = tree_img + 
  theme_tree2() +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  geom_tiplab(aes(label=label, subset=(label!="Hemidactylium scutatum")),
              fontface=3, angle=0, offset = 0.5, size=3) +
  geom_tippoint(aes(subset=grepl("Hemidactylium", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Supplementary Figs. 6-9",
                  subset=(label=="Hemidactylium scutatum")),
              fontface=1, offset = 1, size=3) +
  #geom_text2(aes(label=node), hjust=-.3) + 
  
  # label families
  geom_text2(aes(label="Proteidae", subset=node==257, fontface = "bold"),
             nudge_x=-38, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Rhyacotritonidae", subset=node==134, fontface = "bold"),
             nudge_x=-91.5, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Amphiumidae", subset=node==138, fontface = "bold"),
             nudge_x=-78, nudge_y=0.75, color="black", size=3) + 
  geom_text2(aes(label="Plethodontidae", subset=node==140, fontface = "bold"),
             nudge_x=-18, nudge_y=0.75, color="black", size=3) + 
  
  # label sub-families
  geom_text2(aes(label="Plethodontinae", subset=node==141),
             nudge_x=-15, nudge_y=-0.75, color="grey35", size=3) +   
  geom_text2(aes(label="Hemidactyliinae", subset=node==125),
             nudge_x=6.5, nudge_y=-1.15, color="grey35", size=3) +  
  coord_geo(
    xlim = c(-170, 35), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE), skip=c("Jurassic", "Quaternary"),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-170, -160, -150,
                              -140,-130,-120,-110,
                              -100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("170", "160", "150",
                              "140", "130", "120", "110",
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
ggsave("S5_95conf_prot_rhyaco_amph_plethodontinae.pdf", plot=tt, device="pdf",
       limitsize = FALSE, width=10, height=15)
print("")

#### HEMIDACTYLINAE 1 TREE ####
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_1
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[123] = 1.8996

# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE)) %>%
  ggtree::rotate(66) %>% ggtree::rotate(67)
  
minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))

tt = tree_img +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  theme_tree2() +
  #geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) +
  
  # label sub-genera
  geom_text2(aes(label="Plethopsis", subset=node==68),
             nudge_x=-4.9, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") +   
  geom_text2(aes(label="Batrachoseps", subset=node==70),
             nudge_x=-6, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") +
  
  geom_tippoint(aes(subset=grepl("Oedipina", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Supplementary Figs. 7-9",
                  subset=(label=="Oedipina taylori")),
              fontface=1, offset = 1, size=3) +
  geom_tiplab(aes(label=label, subset=(label!="Oedipina taylori")),
              fontface=3, offset = 1, size=3) +
  coord_geo(
    xlim = c(-120, 25), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE), skip=c("Cretaceous", "Quaternary"),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-120, -110, -100,-90,
                              -80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("120", "110", "100","90",
                              "80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
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
ggsave("S6_95conf_hemi_1.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print("")

#### HEMIDACTYLINAE 2 TREE ####
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_2
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[128] = 0.97


# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE)) %>%
  ggtree::rotate(66) %>% ggtree::rotate(67)

minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))

tt = tree_img +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  theme_tree2() +
  #geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  geom_tippoint(aes(subset=grepl("Thorius", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Supplementary Figs. 8-9", subset=(label=="Thorius munificus")),
              fontface=1, offset = 0.5, size=3) +
  geom_tiplab(aes(label=label, subset=(label!="Thorius munificus")),
              fontface=3, offset = 0.25, size=3) +
  coord_geo(
    xlim = c(-115, 23), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE), skip=c("Cretaceous", "Quaternary"),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-120, -110, -100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("120", "110", "100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
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
ggsave("S7_95conf_hemi_2.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print("")

#### HEMIDACTYLINAE 3 TREE ####
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_3
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)
timetree@phylo$edge.length[99] = 10.0449

# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE)) %>% 
  ggtree::rotate(96) %>% ggtree::rotate(97)

minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))

tt = tree_img +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  theme_tree2() +
  #geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  geom_tippoint(aes(subset=grepl("Bolitoglossa", label)),
                size=3, fill="black", color="black", shape=21) +
  geom_tiplab(aes(label="Supplementary Fig. 9",
                  subset=(label=="Bolitoglossa cuna")),
              fontface=1, offset = 1, size=3) +
  geom_tiplab(aes(label=label, subset=(label!="Bolitoglossa cuna")),
              fontface=3, offset = 0.5, size=3) +
  coord_geo(
    xlim = c(-115, 26), ylim = c(0, Ntip(timetree) + 1),
    pos = as.list(rep("bottom", 2)),
    dat = list("epochs", "periods"),
    height = list(unit(0.5, "lines"), unit(1, "lines")),
    lab = list(FALSE, TRUE), skip=c("Cretaceous", "Quaternary"),
    rot = list(0, 0), size = list(FALSE, 2.5), abbrv = FALSE, neg = TRUE,
  ) +
  scale_x_continuous(breaks=c(-110,-100,-90,-80,-70,-60,
                              -50,-40,-30,-20,-10,0),
                     labels=c("110","100","90","80","70","60",
                              "50","40","30","20","10","0"),
                     #labels=c("", "", "400", "", "",
                     #          "", "300", "", "", "", "200",
                     #          "","","","100","","","","0"),
                     position = "top", 
                     sec.axis = sec_axis(trans = ~.,
                                         breaks = -rev(epochs$max_age),
                                         labels = rev(round(epochs$max_age, 1)))) +
  theme(axis.line = element_blank(), axis.text.x = element_text(angle = -0))
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
ggsave("S8_95conf_hemi_3.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=12)
print("")

#### HEMIDACTYLINAE 4 TREE ####
timetree = read.beast(file="fullTree_95conf_dated.nex")
taxa_names = timetree@phylo$tip.label
target_tree = hemidactylinae_4
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
timetree@phylo = ladderize(timetree@phylo)  
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)


# get the correct bars. the geom_range() function centers bars at the node,
# which is not how this tree should be displayed
# solution by brunoasm in the ggtree github
# https://github.com/YuLab-SMU/ggtree/issues/306
tree_img = revts(ggtree(timetree, size = 0.75, ladderize = FALSE))

minmax = t(matrix(unlist(tree_img$data$height_0.95_HPD),nrow=2))
bar_df = as.data.frame(minmax) %>%
  rename(min = 1,
         max = 2) %>%
  mutate_all(~-.x) %>%
  bind_cols(select(tree_img$data, y))

tt = tree_img +
  geom_segment(aes(x=min, y=y, xend=max, yend=y), 
               data=bar_df, 
               color='red3',
               alpha = 0.6,
               size = 1.5) +
  theme_tree2() +
  #geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.3, size=2) + 
  #geom_text2(aes(label=node), hjust=-.3) + 
  
  # label subgenera
  geom_text2(aes(label="Oaxakia", subset=node==225),
             nudge_x=-3, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="Magnadigita", subset=node==119),
             nudge_x=-3.5, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="Pachymandra", subset=node==154),
             nudge_x=-3.8, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="paste(italic('Mayamandra'), ' + ', italic('Nanotriton'))",
                 subset=node==207),
             nudge_x=1, nudge_y=3, color="grey35", size=2.5, parse = T) + 
  geom_text2(aes(label="Bolitoglossa", subset=node==214),
             nudge_x=-4.5, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") + 
  geom_text2(aes(label="Eladinea", subset=node==157),
             nudge_x=-2.9, nudge_y=-0.5, color="grey35", size=2.5,
             fontface="italic") + 
  
  # Add line to connect mayamandra + nanotriton
  #geom_segment(aes(x = -46, y = 69.7, xend = -42.9, yend = 67.7),
  geom_segment(aes(x = -45, y = 69.7, xend = -42.9, yend = 67.7),
               color = "grey35", size=0.33,
               lineend="round", linejoin="round",
               arrow = arrow(length = unit(0.1,"cm"))) +
  
  geom_tiplab(aes(label=label), fontface=3, offset = 0.5, size=3) +
  coord_geo(
    xlim = c(-85, 22), ylim = c(0, Ntip(timetree) + 1),
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
ggsave("S9_95conf_hemi_4.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=15)
