## This script generates the 16S profiles of healthy adult Daphnia magna, the food source Chlamydomonas reinhardtii,
## and the COMBO medium Daphnia are cultured in. If you don't want to go through the read denoising steps with dada2,
## you can just download the 'rc.rds' file and run from below the #### VISUALIZATION #### line.

library(phyloseq)
library(tidyverse)
library(dada2)
library(vegan)
library(plyr)
library(tidyverse)
library(DECIPHER)
library(phangorn)

path <- "/Users/work/Desktop/mag_paper/16s_comparison/16S_files/"
fnFs <- sort(list.files(path, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(17, 21), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

taxa <- assignTaxonomy(seqtab.nochim, "/Users/work/Desktop/RefSeq-RDP16S_v2_May2018.fa.gz", multithread=TRUE)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = F)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

samples.out <- data.frame(Sample = rownames(seqtab.nochim))
sample.info <- read.csv('/Users/work/Desktop/mag_paper/16s_comparison/16S_files/metadata.csv')
rownames(sample.info)<-sample.info$Sample

rc <- phyloseq(tax_table(taxa),
               sample_data(sample.info),
               otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               phy_tree(fitGTR$tree))

saveRDS(rc, "/Users/work/Desktop/mag_paper/16S_comparison/rc.rds")

#### VISUALIZATION ####

rc <- readRDS("/Users/work/Desktop/mag_paper/16s_comparison/rc.rds")

rc_metagenome <- subset_samples(rc, type == "Chlamydomonas reinhardtii" | 
                                  type == "Kit (negative control)" | type == "Culture media")
rc_metagenome= subset_taxa(rc_metagenome, (Phylum != "Cyanobacteria/Chloroplast") | is.na(Class))
rc_metagenome_phylum = tax_glom(rc_metagenome, "Phylum")
a <- plot_bar(rc_metagenome_phylum, fill = "Phylum") + theme_classic() + ylab("Number of reads") + theme(legend.position = "none", text = element_text(size = 24)) 
a
ggsave("/Users/work/Desktop/mag_paper/figs/16s_profiles/control_16s.pdf", a, device = "pdf", units="in", width=16, height=12, dpi=300)
# Image then edited in Adobe Illustrator to provide correct sample names.


rc.metagenome <- subset_samples(rc, type == "Chlamydomonas reinhardtii" | 
                                type == "Culture media" |
                                  type == "Adult Daphnia magna")
rc.metagenome.merged = merge_samples(rc.metagenome, "type")
rc.metagenome.merged.transformed = transform_sample_counts(rc.metagenome.merged, function(x) x/sum(x))
rc.metagenome.merged.transformed.glommed = tax_glom(rc.metagenome.merged.transformed, "Phylum")
b <- plot_bar(rc.metagenome.merged.transformed.glommed, fill = "Phylum") + theme_classic() + xlab("Sample Type") + ylab("Relative Abundance") + 
  theme(legend.position = "none", text = element_text(size = 24)) 
b
ggsave("/Users/work/Desktop/mag_paper/figs/16s_profiles/daphnia_vs_positivectrls_relabund.pdf", b, device = "pdf", units="in", width=16, height=12, dpi=300)
# Image edited in Adobe Illustrator to 
c.x <- plot_bar(rc_metagenome_phylum, fill = "Phylum") + theme_classic() + ylab("Number of reads") + theme(legend.position = "bottom", text = element_text(size = 36), legend.spacing.x = unit(1, 'cm'))
c.1 <- ggpubr::get_legend(c.x)
c <- ggpubr::as_ggplot(c.1)
c
ggsave("/Users/work/Desktop/mag_paper/figs/16s_profiles/16s_legend.pdf", c, device = "pdf", units="in", width=30, height=12, dpi=300)
ord = ordinate(rc.metagenome, "PCoA", "unifrac", weighted = T)
d <- plot_ordination(rc.metagenome, ord, type = "samples", color = "type") + stat_ellipse(linetype=2) + theme_classic() + geom_point(size = 4)
d
ggsave("/Users/work/Desktop/mag_paper/figs/16s_profiles/daphnia_vs_positivectrls_ordination.pdf", c, device = "pdf", dpi = 300)


