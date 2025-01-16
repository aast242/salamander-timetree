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

setwd("N:/research/wiens/post_review_SuperCRUNCH/05_nuc_only_analyses/06_raxml/03_bipartition/")

extend_val = 3
bar_sz = 2
txt_sz = 5.5
st_num = 2
tt_fp = "nuc_only_bipartition_figtree_rr_bsrename.nex"
timetree = read.beast(file=tt_fp)
fam_spp = c(
  #"Latimeria_chalumnae", "Homo_sapiens",
  #"Rhinatrema_bivittatum", "Ascaphus_montanus",
            "Cryptobranchus_alleganiensis",
            "Andrias_davidianus",
            "Onychodactylus_fischeri",
            "Paradactylodon_persicus",
            "Ranodon_sibiricus",
            "Pachyhynobius_shangchengensis",
            "Hynobius_nigrescens",
            "Salamandrella_keyserlingii",
            "Liua_shihi",
            "Pseudohynobius_flavomaculatus",
            "Batrachuperus_yenyuanensis",
            "Pseudobranchus_axanthus",
            "Siren_intermedia",
            "Ambystoma_mexicanum",
            "Dicamptodon_copei",
            "Salamandrina_terdigitata",
            "Mertensiella_caucasica",
            "Chioglossa_lusitanica",
            "Lyciasalamandra_atifi",
            "Salamandra_salamandra",
            "Pleurodeles_waltl",
            "Tylototriton_wenxianensis",
            "Echinotriton_chinhaiensis",
            "Tylototriton_kweichowensis",
            "Notophthalmus_viridescens",
            "Taricha_granulosa",
            "Euproctus_platycephalus",
            "Laotriton_laoensis",
            "Cynops_pyrrhogaster",
            "Cynops_cyanurus",
            "Paramesotriton_hongkongensis",
            "Pachytriton_feii",
            "Calotriton_asper",
            "Triturus_marmoratus",
            "Lissotriton_vulgaris",
            "Neurergus_crocatus",
            "Ichthyosaura_alpestris",
            "Ommatotriton_ophryticus",
            "Necturus_maculosus",
            "Proteus_anguinus",
            "Rhyacotriton_olympicus",
            "Amphiuma_tridactylum",
            "Plethodon_jordani",
            "Karsenia_koreana",
            "Speleomantes_italicus",
            "Hydromantes_shastae",
            "Desmognathus_wrighti",
            "Phaeognathus_hubrichti",
            "Aneides_flavipunctatus",
            "Ensatina_eschscholtzii",
            "Eurycea_lucifuga",
            "Urspelerpes_brucei",
            "Stereochilus_marginatus",
            "Gyrinophilus_porphyriticus",
            "Pseudotriton_ruber",
            "Hemidactylium_scutatum",
            "Batrachoseps_nigriventris",
            "Cryptotriton_nasalis",
            "Dendrotriton_bromeliacius",
            "Bradytriton_silus",
            "Oedipina_nica",
            "Nyctanolis_pernix",
            "Nototriton_abscondens",
            "Thorius_troglodytes",
            "Chiropterotriton_arboreus",
            "Aquiloeurycea_cephalica",
            "Isthmura_sierraoccidentalis",
            "Bolitoglossa_riletti",
            "Ixalotriton_niger",
            "Parvimolge_townsendi",
            "Pseudoeurycea_rex",
            "Pseudoeurycea_leprosa")

# only include species that are in the species file
timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, fam_spp))
timetree@phylo = ladderize(timetree@phylo)

#### TREE ####
#timetree = read.beast(file=tt_fp)
#taxa_names = timetree@phylo$tip.label
#target_tree = cryptobranchoidea_spp
#epoch_labels = as.character(round(epochs$max_age, 1))
#epoch_labels[[2]] = ""
#keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

#timetree = drop.tip(timetree, setdiff(timetree@phylo$tip.label, keep_spp))
#timetree@phylo = ladderize(timetree@phylo)
timetree@phylo$tip.label = sub("_", " ", timetree@phylo$tip.label)

tt = ggtree(timetree, size = 0.5, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 0.01, size=0.5) +
  geom_tiplab(aes(label=label),
              fontface=3, angle=0, offset = 0.003, size=3) +
  geom_text2(aes(label=bootstrap, subset=bootstrap>0), hjust=-.35, size=2.5) + 
  geom_text2(aes(label="Cryptobranchidae", subset=node==75),
             nudge_x=-0.017, nudge_y=-0.4, color="grey35", size=2.5) +  
  geom_text2(aes(label="Hynobiidae", subset=node==76),
             nudge_x=-0.011, nudge_y=-0.45, color="grey35", size=2.5) +
  geom_text2(aes(label="Sirenidae", subset=node==85),
             nudge_x=-0.011, nudge_y=-0.4, color="grey35", size=2.5) +
  geom_text2(aes(label="Dicamptodontidae", subset=node==15),
             nudge_x=-0.016, nudge_y=-0.4, color="grey35", size=2.5) +  
  geom_text2(aes(label="Ambystomatidae", subset=node==14),
             nudge_x=-0.015, nudge_y=-0.4, color="grey35", size=2.5) +
  geom_text2(aes(label="Salamandridae", subset=node==89),
             nudge_x=-0.013, nudge_y=-0.4, color="grey35", size=2.5) +  
  geom_text2(aes(label="Proteidae", subset=node==112),
             nudge_x=-0.010, nudge_y=-0.4, color="grey35", size=2.5) +
  geom_text2(aes(label="Rhyacotritonidae", subset=node==41),
             nudge_x=-0.015, nudge_y=-0.4, color="grey35", size=2.5) +  
  geom_text2(aes(label="Amphiumidae", subset=node==42),
             nudge_x=-0.012, nudge_y=-0.4, color="grey35", size=2.5) +
  geom_text2(aes(label="Plethodontidae", subset=node==115),
             nudge_x=-0.015, nudge_y=-0.4, color="grey35", size=2.5) +
  
  #geom_cladelab(node=75, label="Cryptobranchidae", offset = 0.08) + 
  #geom_cladelab(node=76, label="Hynobiidae", offset = 0.11) + 
  #geom_cladelab(node=85, label="Sirenidae", offset = 0.08) + 
  #geom_cladelab(node=15, label="Dicamptodontidae", offset = 0.08) + 
  #geom_cladelab(node=14, label="Ambystomatidae", offset = 0.08) + 
  #geom_cladelab(node=89, label="Salamandridae", offset = 0.095) + 
  #geom_cladelab(node=112, label="Proteidae", offset = 0.08) + 
  #geom_cladelab(node=41, label="Rhyacotritonidae", offset = 0.08) + 
  #geom_cladelab(node=42, label="Amphiumidae", offset = 0.08) + 
  #geom_cladelab(node=115, label="Plethodontidae", offset = 0.07) + 
  geom_treescale() + 
  theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  xlim(c(-0.01,0.32))
tt
ggsave("FigX_nuconly_genera_tree_named.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=11)
