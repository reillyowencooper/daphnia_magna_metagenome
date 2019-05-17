# Step 1: Trimmed, filtered reads were assembled into metagenomes using metaSPAdes.
# 7 metagenomes were created: a master coassembly of all four samples, a coassembly of the adult samples, a coassembly of the juvenile samples, 
# and assemblies of all four samples individually.
# After assembly, contig files were renamed according to their assembly and moved out of assembly folder for read mapping.
cd $WORK/metagenomes
module load spades

spades.py --meta -1 s1bac_r1.fastq -2 s1bac_r2.fastq -o s1_assembly/
spades.py --meta -1 s2bac_r1.fastq -2 s2bac_r2.fastq -o s2_assembly/
spades.py --meta -1 s3bac_r1.fastq -2 s3bac_r2.fastq -o s3_assembly/
spades.py --meta -1 s4bac_r1.fastq -2 s4bac_r2.fastq -o s4_assembly/

cat s1bac_r1.fastq s2bac_r1.fastq > adultf.fastq
cat s1bac_r2.fastq s2bac_r2.fastq > adultr.fastq

spades.py --meta -1 adultf.fastq -2 adultr.fastq -o adult_assembly/

cat s3bac_r1.fastq s3bac_r1.fastq > juvf.fastq
cat s4bac_r2.fastq s4bac_r2.fastq > juvr.fastq

spades.py --meta -1 juvf.fastq -2 juvr.fastq -o juv_assembly/

cat adultf.fastq juvf.fastq > masterf.fastq
cat adultr.fastq juvr.fastq > masterr.fastq

spades.py --meta -1 masterf.fastq -2 masterr.fastq -o master_assembly/

cp s1_assembly/contigs.fasta s1_assembly.fasta
cp s2_assembly/contigs.fasta s2_assembly.fasta
cp s3_assembly/contigs.fasta s3_assembly.fasta
cp s4_assembly/contigs.fasta s4_assembly.fasta
cp adult_assembly/contigs.fasta adult_assembly.fasta
cp juv_assembly/contigs.fasta juv_assembly.fasta
cp master_assembly/contigs.fasta master_assembly.fasta

# Step 2: Reads from samples were mapped to specific assemblies. Juvenile reads were mapped to the juvenile coassembly, adult to the adult coassembly, and all to the master coassembly.

module load bwa samtools

bwa index master_assembly.fasta
bwa index adult_assembly.fasta
bwa index juv_assembly.fasta

bwa mem master_assembly.fasta s1bac_r1.fastq s1bac_r2.fastq > adult1_master.sam
bwa mem master_assembly.fasta s2bac_r1.fastq s2bac_r2.fastq > adult2_master.sam
bwa mem master_assembly.fasta s3bac_r1.fastq s3bac_r2.fastq > juv1_master.sam
bwa mem master_assembly.fasta s4bac_r1.fastq s4bac_r2.fastq > juv2_master.sam
bwa mem adult_assembly.fasta s1bac_r1.fastq s1bac_r2.fastq > adult1_adult.sam
bwa mem adult_assembly.fasta s2bac_r1.fastq s2bac_r2.fastq > adult2_adult.sam
bwa mem juv_assembly.fasta s3bac_r1.fastq s3bac_r2.fastq > juv1_juv.sam
bwa mem juv_assembly.fasta s4bac_r1.fastq s4bac_r2.fastq > juv2_juv.sam

samtools view -F 4 -bS adult1_master.sam > adult1_master.bam
samtools view -F 4 -bS adult2_master.sam > adult2_master.bam
samtools view -F 4 -bS juv1_master.sam > juv1_master.bam
samtools view -F 4 -bS juv2_master.sam > juv2_master.bam
samtools view -F 4 -bS adult1_adult.sam > adult1_adult.bam
samtools view -F 4 -bS adult2_adult.sam > adult2_adult.bam
samtools view -F 4 -bS juv1_juv.sam > juv1_juv.bam
samtools view -F 4 -bS juv2_juv.sam > juv2_juv.bam

samtools sort adult1_master.bam > adult1_master_sorted.bam
samtools sort adult2_master.bam > adult2_master_sorted.bam
samtools sort juv1_master.bam > juv1_master_sorted.bam
samtools sort juv2_master.bam > juv2_master_sorted.bam
samtools sort adult1_adult.bam > adult1_adult_sorted.bam
samtools sort adult2_adult.bam > adult2_adult_sorted.bam
samtools sort juv1_juv.bam > juv1_juv_sorted.bam
samtools sort juv2_juv.bam > juv2_juv_sorted.bam

samtools index adult1_master_sorted.bam
samtools index adult2_master_sorted.bam
samtools index juv1_master_sorted.bam
samtools index juv2_master_sorted.bam
samtools index adult1_adult_sorted.bam
samtools index adult2_adult_sorted.bam
samtools index juv1_juv_sorted.bam
samtools index juv2_juv_sorted.bam

# Step 3: Move all sorted .bam files and .fasta files to a folder, then download to own computer for work in Anvi'o.

mkdir for_anvio
cp *sorted.bam for_anvio/
cp *.fasta for_anvio/

scp -r rcooper@crane.unl.edu:/work/cresslerlab/rcooper/metagenomes/for_anvio/* for_anvio/
# Will be prompted for password
# Will take forEVER to download

# Step 4: Anvi'o pre-processing. Here, each assembly contig file gets trimmed for contigs > 2500bp and other identification steps occur. Also make backups of contig databases!
cd Dekstop/metagenomes/for_anvio

anvi-script-reformat-fasta master_assembly.fasta -l 2500 -o master.fasta
anvi-script-reformat-fasta adult_assembly.fasta -l 2500 -o ad.fasta
anvi-script-reformat-fasta juv_assembly.fasta -l 2500 -o juv.fasta

anvi-gen-contigs-database -f master.fasta -o master.db -n 'Master assembly contigs'
anvi-gen-contigs-database -f ad.fasta -o ad.db -n 'Adult assembly contigs'
anvi-gen-contigs-database -f juv.fasta -o juv.db -n 'Juvenile assembly contigs'

anvi-run-hmms -c master.db
anvi-run-hmms -c ad.db
anvi-run-hmms -c juv.db

anvi-get-sequences-for-gene-calls -c master.db -o master_gene_calls.fa
anvi-get-sequences-for-gene-calls -c ad.db -o ad_gene_calls.fa
anvi-get-sequences-for-gene-calls -c juv.db -o juv_gene_calls.fa

cp master.db master_backup.db
cp ad.db ad_backup.db
cp juv.db juv_backup.db

# Here, sorted read .bams are mapped onto their respective binned assemblies.

scp -r rcooper@crane.unl.edu:/work/cresslerlab/rcooper/metagenomes/for_anvio/*.names .

anvi-profile -i adult1_master_sorted.bam -c master.db
anvi-profile -i adult2_master_sorted.bam -c master.db
anvi-profile -i juv1_master_sorted.bam -c master.db
anvi-profile -i juv2_master_sorted.bam -c master.db

anvi-merge adult1_master_sorted.bam-ANVIO_PROFILE/PROFILE.db adult2_master_sorted.bam-ANVIO_PROFILE/PROFILE.db juv1_master_sorted.bam-ANVIO_PROFILE/PROFILE.db juv2_master_sorted.bam-ANVIO_PROFILE/PROFILE.db -o master_merged -c master.db

anvi-profile -i adult1_adult_sorted.bam -c ad.db
anvi-profile -i adult2_adult_sorted.bam -c ad.db

anvi-merge adult1_adult_sorted.bam-ANVIO_PROFILE/PROFILE.db adult2_adult_sorted.bam-ANVIO_PROFILE/PROFILE.db -o ad_merged -c ad.db

anvi-profile -i juv1_juv_sorted.bam -c juv.db
anvi-profile -i juv2_juv_sorted.bam -c juv.db

anvi-merge juv1_juv_sorted.bam-ANVIO_PROFILE/PROFILE.db juv2_juv_sorted.bam-ANVIO_PROFILE/PROFILE.db -o juv_merged -c juv.db

# Step 7: Anvi'o CONCOCT binning and refinement. Refining bins is a manual process, so it can't really be shown in code. My guidelines for refining bins are:
# Check CONCOCT-generated bins. If completeness > 90% and redundancy < 10%, that's a good bin. Leave it.
# If completeness > 90% and redundancy > 10%, refine bin. Probably will have to split bin, which will lower completeness (that's ok).
# If completeness < 90% and redundancy < 10%, try to merge with other closely related bins. If merging increases completeness by > 10%, keep merge. 
# Don't even try if completeness is <20%.

anvi-interactive -p master_merged/PROFILE.db -c master.db -C CONCOCT
# Example for anvi-refine -- do this for bins you want to refine
anvi-refine -p master_merged/PROFILE.db -c master.db -C CONCOCT -b Bin_6
# Example for anvi-merge-bins -- do this for bins that you think should go together
anvi-merge-bins -p master_merged/PROFILE.db -c master.db -C CONCOCT -b Bin_1,Bin_2 -B Bin_A

anvi-interactive -p ad_merged/PROFILE.db -c ad.db -C CONCOCT
anvi-interactive -p juv_merged/PROFILE.db -c juv.db -C CONCOCT

# Step 8: Once you are satisfied with bins, summarize your results. This will generate a folder with a static .html summary of your bins as well as .fasta files of each bin.
# Use these .fasta files for GTDB-Tk.

anvi-summarize -p master_merged/PROFILE.db -c master.db -C CONCOCT --report-aa-seqs-for-gene-calls -o master_bin_summary/
anvi-summarize -p ad_merged/PROFILE.db -c ad.db -C CONCOCT --report-aa-seqs-for-gene-calls -o ad_bin_summary/
anvi-summarize -p juv_merged/PROFILE.db -c juv.db -C CONCOCT --report-aa-seqs-for-gene-calls -o juv_bin_summary/

# Rename all species bins according to assembly they came from -- i.e. species 1 from master assembly should be 'sp1_master.fa'
# Move all renamed .fa bins up one level (into assembly-specific 'bin_by_bin' folder)

mkdir for_prokka
cp master_assembly.fasta for_prokka/master_assembly.fa
cp adult_assembly.fasta for_prokka/adult_assembly.fa
cp juv_assembly.fasta for_prokka/juv_assembly.fa
cp master_bin_summary/bin_by_bin/*.fa for_prokka/
cp ad_bin_summary/bin_by_bin/*.fa for_prokka/
cp juv_bin_summary/bin_by_bin/*.fa for_prokka/

# Step 8.5: Assign taxonomy to binned genomes using GTDB-Tk. This should be done on individual species.
# Needs to be done on HCC.

scp -r for_prokka/*.fa rcooper@crane.unl.edu:/work/cresslerlab/rcooper/metagenomes/
mkdir bin_taxonomy
cp sp* bin_taxonomy/

module load gtdbtk/0.1
gtdbtk classify_wf --genome_dir named_mag_identities --extension fasta --out_dir master_mag_gtdbtk

# ...Continue for number of species binned

# Step 9: Identify potential genes using Prokka for the five identified MAGs.
# Load files onto HCC


module load prokka
prokka --compliant --centre X --outdir sp1_master_prokka/ --prefix sp1_master sp1_master.fa
# ...Continue for number of species binned

# Prokka output will contain date of job start. Take all the .faa (amino acid) files and make sure they have clear names.
# Step 10. These files will be input to GhostKOALA for K Number assignment (https://www.kegg.jp/ghostkoala/).
# Unfortunately, each job has to be run sequentially.

# Files from GhostKOALA will just be named 'user_ko.txt'. Rename to appropriate assembly or species -- i.e. 'master_ko.txt', 'sp1_ko.txt'

# Step 11. In R, clean up the annotated files so only genes with assigned K numbers occur.
# extract_kos.R

sp1 <- read.csv("../sp1_ko.txt", sep = "\t", header = FALSE)
sp1$V2[sp1$V2==""] <- NA
sp1 <- subset(sp1, V2 != "NA")
sp1$V1 <- seq(1:1428)
write.table(sp1, "../sp1_ko_trimmed.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Repeat for each assembly/species. I know this could be a for loop, but I was lazy.

# Step 12: Cleaned files (named with _trimmed) should be submitted to the KEGG Module Mapping tool (https://www.genome.jp/kegg/tool/map_module.html).
# This will show complete or nearly-complete metabolic modules in pathways.

