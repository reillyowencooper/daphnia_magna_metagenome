# Filtered, trimmed reads can be identified using three programs: Kraken, Kaiju, and MetaPhlAn2. Each has pros and cons: 
# Kraken taxes a long time, Kaiju has a high error rate, and MetaPhlAn2 doesn't represent all taxa very well. All three will be used here, with
# Kraken as the main result and Kaiju/MetaPhlAn2 as cross-references.

# All programs must be run on the HCC.

# Step 1. Run Kaiju. Kaiju will identify exact matches and, in greedy mode, anything with up to 3 mismatches.

module load compiler/gcc/4.9
module load kaiju/1.5

kaiju -t $NODES -f $KAIJU_DB -i s1bac_r1.fastq -j s1bac_r2.fastq -a greedy -e 3 -o adult1_kaiju.out 
kaiju -t $NODES -f $KAIJU_DB -i s2bac_r1.fastq -j s2bac_r2.fastq -a greedy -e 3 -o adult2_kaiju.out
kaiju -t $NODES -f $KAIJU_DB -i s3bac_r1.fastq -j s3bac_r2.fastq -a greedy -e 3 -o juv1_kaiju.out
kaiju -t $NODES -f $KAIJU_DB -i s4bac_r1.fastq -j s4bac_r2.fastq -a greedy -e 3 -o juv2_kaiju.out

kaijuReport -t $NODES -n $NAMES -p -r species -i adult1_kaiju.out -o adult1_kaiju.report
kaijuReport -t $NODES -n $NAMES -p -r species -i adult2_kaiju.out -o adult2_kaiju.report
kaijuReport -t $NODES -n $NAMES -p -r species -i juv1_kaiju.out -o juv1_kaiju.report
kaijuReport -t $NODES -n $NAMES -p -r species -i juv2_kaiju.out -o juv2_kaiju.report

# Step 2. Run Kraken

module load kraken

kraken --preload --db $KRACKEN_DB

kraken --fastq-input --threads 8 --db $KRACKEN_DB --paired s1bac_r1.fastq s1bac_r2.fastq --output adult1_kraken.out 
kraken --fastq-input --threads 8 --db $KRACKEN_DB --paired s2bac_r1.fastq s2bac_r2.fastq --output adult2_kraken.out 
kraken --fastq-input --threads 8 --db $KRACKEN_DB --paired s3bac_r1.fastq s3bac_r2.fastq --output juv1_kraken.out 
kraken --fastq-input --threads 8 --db $KRACKEN_DB --paired s4bac_r1.fastq s4bac_r2.fastq --output juv2_kraken.out 

kraken-report --db $KRACKEN_DB adult1_kraken.out > adult1_kraken.report
kraken-report --db $KRACKEN_DB adult2_kraken.out > adult2_kraken.report
kraken-report --db $KRACKEN_DB juv1_kraken.out > juv1_kraken.report
kraken-report --db $KRACKEN_DB juv2_kraken.out > juv2_kraken.report

# Step 3. Run MetaPhlAn2. MetaPhlAn2 only runs on one file at a time, so forward and reverse need to be concatenated afterwards. The subsequent file will have two numbers,
# which correspond to the % reads classified in forward and reverse. Average the two to get total % classified at that taxonomic rank.

module load metaphlan/2.6

metaphlan2.py s1bac_r1.fastq --input_type fastq > ad1f.txt
metaphlan2.py s1bac_r2.fastq --input_type fastq > ad1r.txt
metaphlan2.py s2bac_r1.fastq --input_type fastq > s2f.txt
metaphlan2.py s2bac_r2.fastq --input_type fastq > s2r.txt
metaphlan2.py s3bac_r1.fastq --input_type fastq > s3f.txt
metaphlan2.py s3bac_r2.fastq --input_type fastq > s3r.txt
metaphlan2.py s4bac_r1.fastq --input_type fastq > s4f.txt
metaphlan2.py s4bac_r2.fastq --input_type fastq > s4r.txt