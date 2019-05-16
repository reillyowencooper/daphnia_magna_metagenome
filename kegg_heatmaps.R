## Creating heatmaps of pathways identified using KEGG output.

library(ggplot2)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(forcats)

kegg <- read.csv("/Users/work/Desktop/mag_paper/for_heatmap_long_format.csv", header = F)
colnames(kegg) <- c("Pathway", "MAG", "num_genes", "proportion_total", "Subtype", "Type")

amino_acid_metabolism <- ggplot(subset(kegg, Subtype == "Amino acid metabolism"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 24))

secondary_metabolism <- ggplot(subset(kegg, Subtype == "Biosynthesis of other secondary metabolities"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

carbohydrate_metabolism <- ggplot(subset(kegg, Subtype == "Carbohydrate metabolism"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 24))

energy_metabolism <- ggplot(subset(kegg, Subtype == "Energy metabolism"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

glycan_metabolism <- ggplot(subset(kegg, Subtype == "Glycan biosynthesis and metabolism"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

lipid_metabolism <- ggplot(subset(kegg, Subtype == "Lipid metabolism"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

vitamin_metabolism <- ggplot(subset(kegg, Subtype == "Metabolism of cofactors and vitamins"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

terpenoid_metabolism <- ggplot(subset(kegg, Subtype == "Metabolism of terpenoids and polyketides"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

nucleotide_metabolism <- ggplot(subset(kegg, Subtype == "Nucleotide metabolism"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

xenobiotic_metabolism <- ggplot(subset(kegg, Subtype == "Xenobiotics biodegradation and metabolism"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

env_info <- ggplot(subset(kegg, Type == "Environmental information processing"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

gene_info <- ggplot(subset(kegg, Type == "Genetic information processing"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

cell_process <- ggplot(subset(kegg, Type == "Cellular processes"), aes(MAG, fct_rev(Pathway), fill = num_genes)) + 
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

org_processing <- ggplot(subset(kegg, Type == "Organismal systems"), aes(MAG, fct_rev(Pathway), fill = num_genes)) +
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

resistance <- ggplot(subset(kegg, Type == "Human diseases"), aes(MAG, fct_rev(Pathway), fill = num_genes)) +
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1))

a1 <- amino_acid_metabolism + theme(axis.title.y = element_blank()) 
ggsave("aa_met_pdf_for_fig.pdf", a1, units="in", width=12, height=12, dpi=300)
secondary_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Biosynthesis of other secondary metabolities")
ggsave("secondary_genes.png")
a2 <- carbohydrate_metabolism + theme(axis.title.y = element_blank()) 
ggsave("carb_met_pdf_for_fig.pdf", a2, units="in", width=12, height=12, dpi=300)
energy_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Energy metabolism")
ggsave("energy_genes.png")
glycan_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Glycan biosynthesis and metabolism")
ggsave("glycan_genes.png")
lipid_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Lipid metabolism")
ggsave("lipid_genes.png")
vitamin_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Metabolism of cofactors and vitamins")
ggsave("vitamin_genes.png")
terpenoid_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Metabolism of terpenoids and polyketides")
ggsave("terp_genes.png")
nucleotide_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Nucleotide metabolism")
ggsave("nuc_genes.png")
xenobiotic_metabolism + theme(axis.title.y = element_blank()) + ggtitle("Xenobiotics biodegradation and metabolism")
ggsave("xeno_genes.png")
env_info + theme(axis.title.y = element_blank()) + ggtitle("Environmental information processing")
ggsave("envinfo_genes.png")
gene_info + theme(axis.title.y = element_blank()) + ggtitle("Genetic information processing")
ggsave("geneinfo_genes.png")
cell_process + theme(axis.title.y = element_blank()) + ggtitle("Cellular processes")
ggsave("cellprocess_genes.png")
org_processing + theme(axis.title.y = element_blank()) + ggtitle("Organismal processes")
ggsave("orgprocess_genes.png")
resistance + theme(axis.title.y = element_blank()) + ggtitle("Human diseases and resistance")
ggsave("resistance_genes.png")

trimmed_kegg <-kegg[!(kegg$Pathway=="Biosynthesis of antibiotics"),]
trimmed_kegg <- trimmed_kegg[!(trimmed_kegg$Pathway=="Biosynthesis of amino acids"),]
trimmed_kegg <- na.omit(trimmed_kegg)
all <- ggplot(trimmed_kegg, aes(MAG, fct_rev(Pathway), fill = num_genes)) +
  geom_tile() + 
  scale_fill_gradient(name = "", low = "#ffffff",high = "#2171b5") + theme_classic() + theme(axis.title.x=element_blank(),
                                                                                             axis.ticks.x=element_blank(),
                                                                                             axis.text.x = element_text(angle = 45, hjust = 1)) 
all + theme(axis.title.y = element_blank(), axis.text.y = element_text(size=3), axis.ticks.y=element_blank())