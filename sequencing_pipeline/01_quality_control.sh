# Demultiplexed raw reads came in four sets of forward/reverse paired reads:
# s1f.fastq/s1r.fastq ... s4f.fastq/s4r.fastq

# If any step starts with 'module load', then it was performed on UNL HCC.

# Step 0: Look at read quality using FastQC
fastqc *.fastq

# Step 1: Trim using Trimmomatic with the Nextera primers. This will trim the raw reads according to quality score and will deposit trimmed reads into paired and unpaired files.

module load trimmomatic/0.36

trimmomatic PE s1f.fastq s1r.fastq s1fp.fastq s1fu.fastq s1rp.fastq s1ru.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
trimmomatic PE s2f.fastq s2r.fastq s2fp.fastq s2fu.fastq s2rp.fastq s2ru.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
trimmomatic PE s3f.fastq s3r.fastq s3fp.fastq s3fu.fastq s3rp.fastq s3ru.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
trimmomatic PE s4f.fastq s4r.fastq s4fp.fastq s4fu.fastq s4rp.fastq s4ru.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 2: Map paired, trimmed reads to the Daphnia magna genome, then remove all reads that do map using BWA. 
# This will leave behind bacteria/archaea/viruses/eukaryotes (but mostly bacteria).
# Daphnia magna genome was downloaded from: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA298946

module load bwa
module load samtools
module load bedtools

bwa index magna.fasta

bwa mem magna.fasta s1fp.fastq s1rp.fastq > s1.sam
bwa mem magna.fasta s2fp.fastq s2rp.fastq > s2.sam
bwa mem magna.fasta s3fp.fastq s3rp.fastq > s3.sam
bwa mem magna.fasta s4fp.fastq s4rp.fastq > s4.sam

samtools view -b s1.sam > s1.bam
samtools view -b s2.sam > s2.bam
samtools view -b s3.sam > s3.bam
samtools view -b s4.sam > s4.bam

samtools sort s1.bam > s1s.bam
samtools sort s2.bam > s2s.bam
samtools sort s3.bam > s3s.bam
samtools sort s4.bam > s4s.bam

samtools index s1s.bam
samtools index s2s.bam
samtools index s3s.bam
samtools index s4s.bam

samtools view -b -f 4 s1s.bam > s1_not_mag.bam
samtools view -b -f 4 s2s.bam > s2_not_mag.bam
samtools view -b -f 4 s3s.bam > s3_not_mag.bam
samtools view -b -f 4 s4s.bam > s4_not_mag.bam
samtools view -b -F 4 s1s.bam > s1_mag.bam
samtools view -b -F 4 s2s.bam > s2_mag.bam
samtools view -b -F 4 s3s.bam > s3_mag.bam
samtools view -b -F 4 s4s.bam > s4_mag.bam

samtools sort -n s1_not_mag.bam -o s1_not_mag_sorted.bam
samtools sort -n s2_not_mag.bam -o s2_not_mag_sorted.bam
samtools sort -n s3_not_mag.bam -o s3_not_mag_sorted.bam
samtools sort -n s4_not_mag.bam -o s4_not_mag_sorted.bam

bedtools bamtofastq -i s1_not_mag_sorted.bam -fq s1bac_r1.fastq -fq2 s1bac_r2.fastq
bedtools bamtofastq -i s2_not_mag_sorted.bam -fq s2bac_r1.fastq -fq2 s2bac_r2.fastq
bedtools bamtofastq -i s3_not_mag_sorted.bam -fq s3bac_r1.fastq -fq2 s3bac_r2.fastq
bedtools bamtofastq -i s4_not_mag_sorted.bam -fq s4bac_r1.fastq -fq2 s4bac_r2.fastq

# Reads not mapping to the Daphnia magna genome are titled *bac_r*.fastq.
# These reads were then run through two pipelines. One assembled the filtered bacterial reads into metagenomes using metaSPAdes for contig binning. That can be found in 02_binning.sh. 
# The second pipeline uses the filtered bacterial reads for identifying taxa using Kaiju, Kraken, and MetaPhlAn2. That can be found in 03_taxa.sh.
