## Exploring genes in each MAG annotated with Prokka. This cleans up the Prokka output and adds KO numbers, then compares
## among the MAGs to find which are shared and which are unique to single MAGs.

library(tidyverse)
library(KEGGREST)

setwd("/Users/reillycooper/Desktop/Experiments/Metagenome/binned_species/gene_func_analysis")

burk <- read_tsv("burkholderia_prokka.tsv")
limno1 <- read_tsv("limno1_prokka.tsv")
limno2 <- read_tsv("limno2_prokka.tsv")
pedobac <- read_tsv("pedobacter_prokka.tsv")
polaro <- read_tsv("polaromonas_prokka.tsv")

burk_clean <- subset(burk, ftype == "CDS") %>% 
    subset(product != "hypothetical protein") %>%
    subset(EC_number != "NA")
burk_clean <- burk_clean[,c(1,3,4,5,7)]

limno1_clean <- subset(limno1, ftype == "CDS") %>% 
    subset(product != "hypothetical protein") %>%
    subset(EC_number != "NA")
limno1_clean <- limno1_clean[,c(1,3,4,5,7)]

limno2_clean <- subset(limno2, ftype == "CDS") %>% 
    subset(product != "hypothetical protein") %>%
    subset(EC_number != "NA")
limno2_clean <- limno2_clean[,c(1,3,4,5,7)]

pedobac_clean <- subset(pedobac, ftype == "CDS") %>% 
    subset(product != "hypothetical protein") %>%
    subset(EC_number != "NA")
pedobac_clean <- pedobac_clean[,c(1,3,4,5,7)]

polaro_clean <- subset(polaro, ftype == "CDS") %>% 
    subset(product != "hypothetical protein") %>%
    subset(EC_number != "NA")
polaro_clean <- polaro_clean[,c(1,3,4,5,7)]

write.table(burk_clean, 
            file = "burkholderia_EC_func_only.tsv", 
            quote = F, sep = '\t', 
            col.names = T, row.names = F)

write.table(limno1_clean, 
            file = "limno1_EC_func_only.tsv", 
            quote = F, sep = '\t', 
            col.names = T, row.names = F)

write.table(limno2_clean, 
            file = "limno2_EC_func_only.tsv", 
            quote = F, sep = '\t', 
            col.names = T, row.names = F)

write.table(pedobac_clean, 
            file = "pedobacter_EC_func_only.tsv", 
            quote = F, sep = '\t', 
            col.names = T, row.names = F)

write.table(polaro_clean, 
            file = "polaromonas_EC_func_only.tsv", 
            quote = F, sep = '\t', 
            col.names = T, row.names = F)

burk_kos <- read.table("../ghostkoala_output/burkholderiaceae_func_cat.txt", fill = T, sep = "\t", row.names = NULL, col.names = F) %>% subset(FALSE. != "")
limno1_kos <- read.table("../ghostkoala_output/limno1_func_cat.txt", fill = T, sep = "\t", row.names = NULL, col.names = F) %>% subset(FALSE. != "")
limno2_kos <- read.table("../ghostkoala_output/limno2_func_cat.txt", fill = T, sep = "\t", row.names = NULL, col.names = F) %>% subset(FALSE. != "")
pedobac_kos <- read.table("../ghostkoala_output/pedobacter_func_cat.txt", fill = T, sep = "\t", row.names = NULL, col.names = F) %>% subset(FALSE. != "")
polaro_kos <- read.table("../ghostkoala_output/polaromonas_func_cat.txt", fill = T, sep = "\t", row.names = NULL, header = F) %>% subset(V2 != "")

burk_kos$row.names <- NULL
limno1_kos$row.names <- NULL
limno2_kos$row.names <- NULL
pedobac_kos$row.names <- NULL
polaro_kos$V1 <- NULL

burk_kos$color <- rep('red')
limno1_kos$color <- rep('red')
limno2_kos$color <- rep('red')
pedobac_kos$color <- rep('red')
polaro_kos$color <- rep('red')

write.table(burk_kos, file = "burkholderia_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(limno1_kos, file = "limno1_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(limno2_kos, file = "limno2_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(pedobac_kos, file = "pedobacter_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(polaro_kos, file = "polaromonas_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')

# Getting unique K numbers in each organism
burk_list <- pull(burk_kos, "FALSE.") %>% as.character()
limno1_list <- pull(limno1_kos, "FALSE.") %>% as.character()
limno2_list <- pull(limno2_kos, "FALSE.") %>% as.character()
pedobac_list <- pull(pedobac_kos, "FALSE.") %>% as.character()
polaro_list <- pull(polaro_kos, V2) %>% as.character()

burk_unique <- burk_list[!(burk_list %in% limno1_list)]
    burk_unique <- burk_unique[!(burk_unique %in% limno2_list)]
    burk_unique <- burk_unique[!(burk_unique %in% pedobac_list)]
    burk_unique <- burk_unique[!(burk_unique %in% polaro_list)]    

limno1_unique <- limno1_list[!(limno1_list %in% burk_list)]
    limno1_unique <- limno1_unique[!(limno1_unique %in% limno2_list)]
    limno1_unique <- limno1_unique[!(limno1_unique %in% pedobac_list)]
    limno1_unique <- limno1_unique[!(limno1_unique %in% polaro_list)]     
    
limno2_unique <- limno2_list[!(limno2_list %in% burk_list)]
    limno2_unique <- limno2_unique[!(limno2_unique %in% limno1_list)]
    limno2_unique <- limno2_unique[!(limno2_unique %in% pedobac_list)]
    limno2_unique <- limno2_unique[!(limno2_unique %in% polaro_list)]    
    
pedobac_unique <- pedobac_list[!(pedobac_list %in% burk_list)]
    pedobac_unique <- pedobac_unique[!(pedobac_unique %in% limno1_list)]
    pedobac_unique <- pedobac_unique[!(pedobac_unique %in% limno2_list)]
    pedobac_unique <- pedobac_unique[!(pedobac_unique %in% polaro_list)]   
    
polaro_unique <- polaro_list[!(polaro_list %in% burk_list)]
    polaro_unique <- polaro_unique[!(polaro_unique %in% limno1_list)]
    polaro_unique <- polaro_unique[!(polaro_unique %in% limno2_list)]
    polaro_unique <- polaro_unique[!(polaro_unique %in% pedobac_list)]  
    
unique_burk <- as.data.frame(burk_unique) %>% mutate(color = 'red')
unique_limno1 <- as.data.frame(limno1_unique) %>% mutate(color = 'red')
unique_limno2 <- as.data.frame(limno2_unique) %>% mutate(color = 'red')
unique_pedobac <- as.data.frame(pedobac_unique) %>% mutate(color = 'red')
unique_polaro <- as.data.frame(polaro_unique) %>% mutate(color = 'red')

write.table(unique_burk, file = "burkholderia_unique_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(unique_limno1, file = "limno1_unique_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(unique_limno2, file = "limno2_unique_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(unique_pedobac, file = "pedobacter_unique_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(unique_polaro, file = "polaromonas_unique_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')

# Creating list of K numbers in each organism but with unique KOs colored red and shared colored blue
colnames(burk_kos) <- c("k_num", "color")
colnames(limno1_kos) <- c("k_num", "color")
colnames(limno2_kos) <- c("k_num", "color")
colnames(pedobac_kos) <- c("k_num", "color")
colnames(polaro_kos) <- c("k_num", "color")
colnames(unique_burk) <- c("k_num", "color")
colnames(unique_limno1) <- c("k_num", "color")
colnames(unique_limno2) <- c("k_num", "color")
colnames(unique_pedobac) <- c("k_num", "color")
colnames(unique_polaro) <- c("k_num", "color")

burk_all_list <- pull(burk_kos, "k_num") %>% as.character()
burk_same <- burk_all_list[!(burk_all_list %in% burk_unique)] %>% as.data.frame() %>% mutate(color = 'blue')
colnames(burk_same) <- c("k_num", "color")
burk_all <- rbind(burk_same, unique_burk)

limno1_all_list <- pull(limno1_kos, "k_num") %>% as.character()
limno1_same <- limno1_all_list[!(limno1_all_list %in% limno1_unique)] %>% as.data.frame() %>% mutate(color = 'blue')
colnames(limno1_same) <- c("k_num", "color")
limno1_all <- rbind(limno1_same, unique_limno1)

limno2_all_list <- pull(limno2_kos, "k_num") %>% as.character()
limno2_same <- limno2_all_list[!(limno2_all_list %in% limno2_unique)] %>% as.data.frame() %>% mutate(color = 'blue')
colnames(limno2_same) <- c("k_num", "color")
limno2_all <- rbind(limno2_same, unique_limno2)

pedobac_all_list <- pull(pedobac_kos, "k_num") %>% as.character()
pedobac_same <- pedobac_all_list[!(pedobac_all_list %in% pedobac_unique)] %>% as.data.frame() %>% mutate(color = 'blue')
colnames(pedobac_same) <- c("k_num", "color")
pedobac_all <- rbind(pedobac_same, unique_pedobac)

polaro_all_list <- pull(polaro_kos, "k_num") %>% as.character()
polaro_same <- polaro_all_list[!(polaro_all_list %in% polaro_unique)] %>% as.data.frame() %>% mutate(color = 'blue')
colnames(polaro_same) <- c("k_num", "color")
polaro_all <- rbind(polaro_same, unique_polaro)

# red = unique, blue = shared
write.table(burk_all, file = "burkholderia_all_kos_colored.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(limno1_all, file = "limno1_all_kos_colored.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(limno2_all, file = "limno2_all_kos_colored.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(pedobac_all, file = "pedobacter_all_kos_colored.txt", col.names = F, row.names = F, quote = F, sep = '\t')
write.table(polaro_all, file = "polaromonas_all_kos_colored.txt", col.names = F, row.names = F, quote = F, sep = '\t')

shared_in_all <- Reduce(intersect, list(burk_list, limno1_list, limno2_list, pedobac_list, polaro_list)) %>% as.data.frame() %>% mutate(color = 'green')
colnames(shared_in_all) <- c("k_num", "color")
write.table(shared_in_all, file = "shared_in_all_kos.txt", col.names = F, row.names = F, quote = F, sep = '\t')


burk_module <- burk_all
burk_module$color <- NULL
burk_module$gene <- 1:nrow(burk_module)
burk_module <- rev(burk_module)
write.table(burk_module, file = "burkholderia_module.txt", col.names = F, row.names = F, quote = F, sep = '\t')

limno1_module <- limno1_all
limno1_module$color <- NULL
limno1_module$gene <- 1:nrow(limno1_module)
limno1_module <- rev(limno1_module)
write.table(limno1_module, file = "limno1_module.txt", col.names = F, row.names = F, quote = F, sep = '\t')

limno2_module <- limno2_all
limno2_module$color <- NULL
limno2_module$gene <- 1:nrow(limno2_module)
limno2_module <- rev(limno2_module)
write.table(limno2_module, file = "limno2_module.txt", col.names = F, row.names = F, quote = F, sep = '\t')

pedobac_module <- pedobac_all
pedobac_module$color <- NULL
pedobac_module$gene <- 1:nrow(pedobac_module)
pedobac_module <- rev(pedobac_module)
write.table(pedobac_module, file = "pedobacter_module.txt", col.names = F, row.names = F, quote = F, sep = '\t')

polaro_module <- polaro_all
polaro_module$color <- NULL
polaro_module$gene <- 1:nrow(polaro_module)
polaro_module <- rev(polaro_module)
write.table(polaro_module, file = "polaromonas_module.txt", col.names = F, row.names = F, quote = F, sep = '\t')