library(ggtree)
library(tidyverse)
library(treeio)

tree <- read.tree("Desktop/mag_paper/phylogenetic_placement/gtdbtk.bac120.classify.tree")

limno_subset <- tree_subset(tree, "daphnia_limnohabitans_sp1", levels_back = 10)


limno_subset$tip.label[limno_subset$tip.label=="daphnia_limnohabitans_sp1"] <- "Daphnia Limnohabitans MAG 1"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001412535.1"] <- "Limnohabitans sp. 63ED37-2 (GCF_001412535.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_000293865.1"] <- "Limnohabitans sp. Rim28 (GCF_000293865.1)"
limno_subset$tip.label[limno_subset$tip.label=="daphnia_limnohabitans_sp2"] <- "Daphnia Limnohabitans MAG 2"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001269345.1"] <- "Limnohabitans sp. 2KL-27 (GCF_001269345.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001517545.1"] <- "Limnohabitans sp. Rim47 (GCF_001517545.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_002778285.1"] <- "Limnohabitans sp. 15K (GCF_002778285.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_002778275.1"] <- "Limnohabitans sp. G3-2 (GCF_002778275.1)"
limno_subset$tip.label[limno_subset$tip.label=="GB_GCA_002256145.1"] <- "Comamonadaceae bacterium PBBC1 (GCA_002256145.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001269405.1"] <- "Limnohabitans sp. 2KL-3 (GCF_001269405.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001269385.1"] <- "Limnohabitans sp. DM1 (GCF_001269385.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001270065.1"] <- "Limnohabitans planktonicus II-D5 (GCF_001270065.1)"
limno_subset$tip.label[limno_subset$tip.label=="GB_GCA_002255935.1"] <- "Comamonadaceae bacterium PBBC2 (GCA_002255935.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001412575.1"] <- "Limnohabitans sp. 103DPR2 (GCF_001412575.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_001269365.1"] <- "Limnohabitans sp. Rim11 (GCF_001269365.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_002778325.1"] <- "Limnohabitans sp. B9-3 (GCF_002778325.1)"
limno_subset$tip.label[limno_subset$tip.label=="RS_GCF_002778315.1"] <- "Limnohabitans sp. JirII-31 (GCF_002778315.1)"
limno_subset$tip.label[limno_subset$tip.label=="GB_GCA_002282445.1"] <- "Burkholderiales bacterium 39-55-53 (GCA_002282445.1)"
# limno_subset$tip.label


pedobacter_subset <- tree_subset(tree, "daphnia_pedobacter", levels_back = 5)

pedobacter_subset$tip.label[pedobacter_subset$tip.label=="daphnia_pedobacter"] <- "Daphnia Pedobacter MAG"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="RS_GCF_900103545.1"] <- "Pedobacter ruber (GCF_900103545.1)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="GB_GCA_002257025.1"] <- "Sphingobacteriales bacterium 16-39-50 (GCA_002257025.1)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="RS_GCF_900168015.1"] <- "Pedobacter luteus (GCF_900168015.1)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="RS_GCF_000422945.1"] <- "Pedobacter oryzae DSM 19973 (GCF_000422945.1)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="UBA10854"] <- "Sphingobacteriaceae bacterium UBA10854 (UBA10854)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="GB_GCA_002332615.1"] <- "Sphingobacteriaceae bacterium UBA2059 (GCA_002332615.1)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="RS_GCF_000403135.1"] <- "Arcticibacter svalbardensis MN12-7 (GCF_000403135.1)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="GB_GCA_002454425.1"] <- "Sphingobacteriaceae bacterium UBA6709 (GCA_002454425.1)"
pedobacter_subset$tip.label[pedobacter_subset$tip.label=="RS_GCF_001730525.1"] <- "Arcticibacter eurypsychrophilus (GCF_001730525.1)"
#pedobacter_subset$tip.label


polaromonas_subset <- tree_subset(tree, "daphnia_polaromonas", levels_back = 3)
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="daphnia_polaromonas"] <- "Daphnia Polaromonas MAG"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_001770935.1"] <- "Burkholderiales bacterium RIFCSPLOWO2_12_FULL_61_40 (GCA_001770935.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_002163715.1"] <- "Curvibacter sp. AEP1-3 (GCF_002163715.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002390825.1"] <- "Rhodoferax sp. UBA4409 (GCA_002390825.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_001795465.1"] <- "Curvibacter sp. RIFCSPHIGHO2_12_FULL_63_18 (GCA_001795465.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_002251855.1"] <- "Rhodoferax sp. TH121 (GCF_002251855.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_001955715.1"] <- "Rhodoferax saidenbachensis (GCF_001955715.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002256115.1"] <- "Burkholderiales bacterium PBB3 (GCA_002256115.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002422455.1"] <- "Rhodoferax sp. UBA6134 (GCA_002422455.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_000013605.1"] <- "Rhodoferax ferrireducens T118 (GCF_000013605.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002413825.1"] <- "Rhodoferax sp. UBA5149 (GCA_002413825.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002381045.1"] <- "Rhodoferax sp. UBA4127 (GCA_002381045.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002789135.1"] <- "Comamonadaceae bacterium CG_4_9_14_0_8_um_filter_57_21 (GCA_002789135.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002842265.1"] <- "Betaproteobacteria bacterium HGW-Betaproteobacteria-18 (GCA_002842265.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_001771065.1"] <- "Burkholderiales bacterium RIFOXYD12_FULL_59_19 (GCA_001771065.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002789915.1"] <- "Comamonadaceae bacterium CG_4_9_14_3_um_filter_60_33 (GCA_002789915.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_002017865.1"] <- "Rhodoferax fermentans (GCF_002017865.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002083775.1"] <- "Rhodoferax ferrireducens (GCA_002083775.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_001938565.1"] <- "Rhodoferax antarcticus ANT.BR (GCF_001938565.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002786915.1"] <- "Comamonadaceae bacterium CG12_big_fil_rev_8_21_14_0_65_59_15 (GCA_002786915.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002781645.1"] <- "Comamonadaceae bacterium CG17_big_fil_post_rev_8_21_14_2_50_60_13 (GCA_002781645.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_002001015.1"] <- "Polaromonas sp. A23 (GCF_002001015.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002281715.1"] <- "Polaromonas sp. 16-63-31 (GCA_002281715.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_002379095.1"] <- "Polaromonas sp. AER18D-145 (GCF_002379095.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_002379085.1"] <- "Polaromonas sp. AET17H-212 (GCF_002379085.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002256925.1"] <- "Polaromonas sp. 28-63-22 (GCA_002256925.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_900116715.1"] <- "Polaromonas sp. YR568 (GCF_900116715.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_900103405.1"] <- "Polaromonas sp. JS666 (GCF_900103405.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_000282655.1"] <- "Polaromonas sp. CF318 (GCF_000282655.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_000013865.1"] <- "Polaromonas sp. JS666 (GCF_000013865.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_001598235.1"] <- "Polaromonas jejuensis NBRC 106434 (GCF_001598235.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002413395.1"] <- "Polaromonas sp. UBA5171 (GCA_002413395.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_002381615.1"] <- "Polaromonas sp. UBA4122 (GCA_002381615.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="GB_GCA_001770785.1"] <- "Burkholderiales bacterium RIFCSPHIGHO2_12_FULL_61_11 (GCA_001770785.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_000688115.1"] <- "Polaromonas sp. EUR3 1.2.1 (GCF_000688115.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_900112285.1"] <- "Polaromonas sp. OV174 (GCF_900112285.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_000015505.1"] <- "Polaromonas naphthalenivorans CJ2 (GCF_000015505.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_000751355.1"] <- "Polaromonas sp. CG9_12 (GCF_000751355.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_000709345.1"] <- "Polaromonas glacialis (GCF_000709345.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_001422405.1"] <- "Pseudorhodoferax sp. Leaf267 (GCF_001422405.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_001422445.1"] <- "Pseudorhodoferax sp. Leaf274 (GCF_001422445.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_001422365.1"] <- "Pseudorhodoferax sp. Leaf265 (GCF_001422365.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_001955695.1"] <- "Rhodoferax koreense (GCF_001955695.1)"
polaromonas_subset$tip.label[polaromonas_subset$tip.label=="RS_GCF_900104385.1"] <- "Albidiferax sp. OV413 (GCF_900104385.1)"
#polaromonas_subset$tip.label

limno_tree <- ggtree(limno_subset) + geom_tiplab(size = 3, color = "blue") + xlim_tree(.3)+ geom_text2(aes(label=label, subset=!isTip), hjust=-.3,size = 2) +
  geom_point2(aes(subset=!isTip), color="black", size=0) + geom_treescale(x = 0)
ggsave(filename = "limno_tree.pdf" , plot = limno_tree, device = "pdf", width = 15, height = 30 , units = "in" , limitsize = FALSE)

pedo_tree <- ggtree(pedobacter_subset) + geom_tiplab(size = 3, color = "blue") + xlim_tree(.4) + geom_text2(aes(label=label, subset=!isTip), hjust=-.3,size = 2) +
  geom_point2(aes(subset=!isTip), color="black", size=0) + geom_treescale(x = 0)
ggsave(filename = "pedo_tree.pdf" , plot = pedo_tree, device = "pdf", width = 15, height = 30 , units = "in" , limitsize = FALSE)

polaro_tree <- ggtree(polaromonas_subset) + geom_tiplab(size = 3, color = "blue") + xlim_tree(.4) + geom_text2(aes(label=label, subset=!isTip), hjust=-.3,size = 2) +
  geom_point2(aes(subset=!isTip), color="black", size=0) + geom_treescale(x = 0)
ggsave(filename = "polaro_tree.pdf" , plot = polaro_tree, device = "pdf", width = 15, height = 30 , units = "in" , limitsize = FALSE)

tree2 <- tree
tree2$tip.label[tree2$tip.label=="daphnia_limnohabitans_sp1"] <- "Limnohabitans MAG 1"
tree2$tip.label[tree2$tip.label=="daphnia_limnohabitans_sp2"] <- "Limnohabitans MAG 2"
tree2$tip.label[tree2$tip.label=="daphnia_pedobacter"] <- "Pedobacter MAG"
tree2$tip.label[tree2$tip.label=="daphnia_polaromonas"] <- "Polaromonas MAG"
tree2$tip.label[tree2$tip.label=="daphnia_unknown_burkholderiaceae"] <- "Unknown Burkholderiaceae MAG"
tree2$tip.label[tree2$tip.label=="daphnia_unknown_chitinphagales"] <- "Unknown Chitinphagales MAG"
tree2$tip.label[tree2$tip.label=="daphnia_emticicia"] <- "Emticicia MAG"

limited_tree2 <- ggtree(tree2, layout = "circular") + geom_tiplab2(aes(
  subset=(grepl('MAG',label,fixed=TRUE)==TRUE), angle=angle), size = 1)
ggsave(filename = "full_tree.pdf" , plot = limited_tree2, device = "pdf", width = 15, height = 30 , units = "in" , limitsize = FALSE)
