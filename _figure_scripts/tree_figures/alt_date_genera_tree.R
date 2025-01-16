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
tt_fp = "S3_Dated_ML_tree.nex"
fam_spp = c("Latimeria_chalumnae", "Homo_sapiens",
            "Rhinatrema_bivittatum", "Ascaphus_montanus",
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


#### TREE ####
timetree = read.nexus(file=tt_fp)
taxa_names = timetree$tip.label
target_tree = fam_spp
epoch_labels = as.character(round(epochs$max_age, 1))
epoch_labels[[2]] = ""
keep_spp = taxa_names[grep(paste(target_tree, collapse="|"), taxa_names)]

timetree = drop.tip(timetree, setdiff(timetree$tip.label, keep_spp))
timetree = ladderize(timetree)
timetree$tip.label = sub("_", " ", timetree$tip.label)

tt = ggtree(timetree, size = 0.5, ladderize = FALSE) +
  theme_tree2() +
  geom_rootedge(rootedge = 40, size=0.5) +
  geom_tiplab(aes(label=label),
              fontface=3, angle=0, offset = 0.001) + 
  coord_geo(
    xlim = c(-450, 145), ylim = c(0, Ntip(timetree) + 0.5),
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

tt = revts(tt)
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
ggsave("SX_main_date_genera_tree.pdf", plot=tt, device="pdf", limitsize = FALSE, 
       width=10, height=14)
