## This is another way of looking at unique genes in each of the five MAGs. Output generated can be found in Supplementary Table 3.

library(tidyverse)

drop.cols <- c("locus_tag", "ftype", "length_bp")
lim1 <- read.csv("Desktop/mag_paper/gene_func_analysis/limno1_prokka.tsv", sep = "\t") %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  separate(gene, into = c("gene", "copy"), sep="_(?=[^_]+$)") %>%
  mutate(MAG = rep("Limnohabitans sp. 1")) %>%
  subset(copy < 2 | is.na(copy)) %>%
  select(-one_of(drop.cols))
lim2 <- read.csv("Desktop/mag_paper/gene_func_analysis/limno2_prokka.tsv", sep = "\t") %>% distinct(gene, .keep_all = TRUE) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  separate(gene, into = c("gene", "copy"), sep="_(?=[^_]+$)") %>%
  mutate(MAG = rep("Limnohabitans sp. 2")) %>%
  subset(copy < 2 | is.na(copy)) %>%
  select(-one_of(drop.cols))
burk <- read.csv("Desktop/mag_paper/gene_func_analysis/burkholderia_prokka.tsv", sep = "\t") %>% distinct(gene, .keep_all = TRUE) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  separate(gene, into = c("gene", "copy"), sep="_(?=[^_]+$)") %>%
  mutate(MAG = rep("Unknown Burkholderiaceae")) %>%
  subset(copy < 2 | is.na(copy)) %>%
  select(-one_of(drop.cols))
pedo <- read.csv("Desktop/mag_paper/gene_func_analysis/pedobacter_prokka.tsv", sep = "\t") %>% distinct(gene, .keep_all = TRUE) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  separate(gene, into = c("gene", "copy"), sep="_(?=[^_]+$)") %>%
  mutate(MAG = rep("Pedobacter sp.")) %>%
  subset(copy < 2 | is.na(copy)) %>%
  select(-one_of(drop.cols))
pola <- read.csv("Desktop/mag_paper/gene_func_analysis/polaromonas_prokka.tsv", sep = "\t") %>% distinct(gene, .keep_all = TRUE) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  separate(gene, into = c("gene", "copy"), sep="_(?=[^_]+$)") %>%
  mutate(MAG = rep("Polaromonas sp.")) %>%
  subset(copy < 2 | is.na(copy)) %>%
  select(-one_of(drop.cols))

everything <- rbind(lim1, lim2, burk, pedo, pola) %>% select(-one_of("copy")) %>% subset(gene != "")
everything_unique <- distinct(everything, gene, .keep_all = T)
write.csv(everything_unique, "unique_genes.csv")
