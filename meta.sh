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

#Step3:Assembly












