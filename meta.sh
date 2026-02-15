#File preparation
#metagenmic flies stored in /seq
#results stored in /Result
mkdir seq
mkdir Result

#Count sequence length and base number
seqkit stat seq/*.fastq

#Step1:Quality control
#Trimmomatic
java -jar Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 \
-threads 8 seq/sample_1.fastq seq/sample_2.fastq \
Result/sample_1.fq seq/sample_1.unpaired.fastq \
Result/sample_2.fq Result/sample_2.unpaired.fastq \
LEADING:20 TRAILING:20 \
SLIDINGWINDOW:4:20 MINLEN:50 \
ILLUMINACLIP:Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10

#fastqc/multiqc
fastqc result/*.fq -t 4
multiqc -d seq/ -o Result/

#Step2:Romove of the host sequences
#Bowtie2
bowtie2-build --threads 20 host_genomic.fna host_genome
bowtie2 -p 20 -x host_genome -1 Result/sample_1.fq -2 Result/sample_2.fq -S Result/samplemeta.sam
samtools view -bS --threads 20 Result/samplemeta.sam > Result/samplemeta.bam
samtools view -b -f 12 -F 256 Result/samplemeta.bam > Result/samplemeta_unmapped.bam
samtools sort Result/samplemeta_unmapped.bam -o Result/samplemeta_unmapped_sorted.bam --threads 20
bedtools bamtofastq -i Result/samplemeta_unmapped_sorted.bam -fq Result/sample_host_removed_1.fq -fq2 Result/sample_host_removed_2.fq

#Step3:Classification 
#Kraken & Bracken
mkdir Microbe
kraken2 --db kraken_database --threads 56 --report Microbe/TEST.report --output Microbe/TEST.output  --paired Result/sample_host_removed_1.fq Result/sample_host_removed_2.fq
#Phylum level
bracken -d kraken_database -i Microbe/TEST.report -t 10 -l P -o Microbe/TEST.P.bracken
#Genus level
bracken -d kraken_database -i Microbe/TEST.report -t 10 -l G -o Microbe/TEST.G.bracken

#Step4:Assembly
#Megahit
time megahit -t 50 
-1 `ls Result/sample_host_removed_1.fq|tr '\n' ','|sed 's/,$//'` 
-2 `ls Result/sample_host_removed_2.fq|tr '\n' ','|sed 's/,$//'` 
--min-contig-len 300 
-o Result/Megahit
#Assess assembly quality
ln Result/Megahit/final.contigs.fa Result/Megahit/
quast.py Result/Megahit/final.contigs.fa -o Result/Megahit/

#Step5:KEGG annotation
#Prokka
mkdir Result/prokka
time prokka Result/Megahit/final.contigs.fa --outdir Result/prokka \
--prefix mg --metagenome --kingdom Bacteria \
--force --cpus 50
#eggNOG
emapper.py \
  -i Result/prokka/protein.fa \
  -o Result/eggnog_annotation \
  --itype proteins \
  -m diamond \
  --cpu 20 \
  --data_dir eggnog #database

#Step6:Binning
mkdir Result/MAG
time bowtie2-build -f Result/Megahit/final.contigs.fa Result/MAG/sample_final --threads 16
time bowtie2 -1 Result/sample_host_removed_1.fq  -2 Result/sample_host_removed_2.fq  -p 16 -x Result/MAG/sample_final -S Result/MAG/sample_final.sam
time samtools view -@ 16 -b -S Result/MAG/sample_final.sam -o Result/MAG/sample_final.bam
time samtools sort -@ 16 -l 9 Result/MAG/sample_final.bam -o Result/MAG/sample_final.sorted.bam
#MetaBAT
time jgi_summarize_bam_contig_depths --outputDepth Result/MAG/sample_final.depth.txt Result/MAG/sample_final.sorted.bam
time metabat2 -m 1500 -t 16 -i Result/Megahit/final.contigs.fa -a Result/MAG/sample_final.depth.txt -o Result/MAG/sample_binning/sample -v
#Checkm
mkdir Result/MAG/checkm
checkm lineage_wf -t 2 -x fa Result/MAG/sample_binning/sample/ Result/MAG/checkm

#Step7:KEGG annotation
#Prokka
mkdir Result/MAG/prokka
for filename in Result/MAG/sample_binning/sample/*.fa
do
  base=$(basename $filename .fa)
  echo $base
  prokka Result/MAG/sample_binning/sample/${base}.fa --outdir Result/MAG/prokka/${base}prokka prefix goldbacteria --kingdom Bacteria
done
#KOfam
mkdir Result/MAG/ko_tmp
mkdir Result/MAG/KEGG
exec_annotation \
-f detail-tsv \
-E 1e-3 \
--profile KEGG/profiles \
--ko-list KEGG/ko_list \
--cpu 6 \
--tmp-dir Result/MAG/ko_tmp \
-o Result/MAG/KEGG/ko_out.txt \
Result/MAG/prokka/prokka/MAG.faa









