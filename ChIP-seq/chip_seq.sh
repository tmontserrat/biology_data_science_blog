#!/bin/bash

SECONDS=0

# Change to the working directory
cd ..


### QUALITY CONTROL ###

# Quality control for the FNR reads
fastqc data/SRR576933.fastq -o quality/

# Quality control for the control reads
fastqc data/SRR576938.fastq -o quality/

# Trimming the FNR reads
java -jar ${HOME}/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-threads 4 \
data/SRR576933.fastq data/SRR576933_trimmed.fastq \
ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
TRAILING:20 \
MINLEN:20 \
-phred33 

# Trimming the control reads
java -jar ${HOME}/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-threads 4 \
data/SRR576938.fastq data/SRR576938_trimmed.fastq \
ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 \
TRAILING:20 \
MINLEN:20 \
-phred33 

# Quality control of the trimmed reads
fastqc data/SRR576933_trimmed.fastq -o quality/
fastqc data/SRR576938_trimmed.fastq -o quality/


### MAPPING READS ###

# Create index 
bowtie2-build reference/e_coli_genome.fa reference/e_coli_genome

# Mapping FNR
bowtie2 -p 5 -q \
-x reference/e_coli_genome \
-U data/SRR576933_trimmed.fastq \
-S alignment/SRR576933_unsorted.sam

# Mapping the control
bowtie2 -p 5 -q \
-x reference/e_coli_genome \
-U data/SRR576938_trimmed.fastq \
-S alignment/SRR576938_unsorted.sam

# From sam to bam
samtools view -b alignment/SRR576933_unsorted.sam | samtools sort -o alignment/SRR576933_sorted.bam
samtools view -b alignment/SRR576938_unsorted.sam | samtools sort -o alignment/SRR576938_sorted.bam

# Index bam file
samtools index -b alignment/SRR576933_sorted.bam
samtools index -b alignment/SRR576938_sorted.bam

### ChIP-seq quality control
plotFingerprint -p 5 -b alignment/SRR576933_sorted.bam alignment/SRR576938_sorted.bam \
-plot quality/fingerprint.png --ignoreDuplicates


### PEAK CALLING ###

# Genome size
gsize=$"4641652"
echo "Genome size: ${gsize}"

# Peak calling
macs3 callpeak -t alignment/SRR576933_sorted.bam \
-c alignment/SRR576938_sorted.bam \
--keep-dup 1 \
--g $gsize \
-n  macs/macs \
-f BAM \
-B \
-q 0.05 \
--nomodel \
--extsize 211

# Print number of peaks found
number_peaks=$(wc -l macs/macs_peaks.narrowPeak)
echo "Number of peaks: ${number_peaks}"

# Remove duplicates from FNR data to scale it
samtools rmdup -s alignment/SRR576933_sorted.bam alignment/SRR576933_sorted_nodup.bam

# Index the bam file
samtools index alignment/SRR576933_sorted_nodup.bam alignment/SRR576933_sorted_nodup.bai

# Remove duplicates from input data to scale it
samtools rmdup -s alignment/SRR576938_sorted.bam alignment/SRR576938_sorted_nodup.bam

# Index the bam file
samtools index alignment/SRR576938_sorted_nodup.bam alignment/SRR576938_sorted_nodup.bai

# Convert bam files to scaled bedgraph files
bamCoverage --bam alignment/SRR576933_sorted_nodup.bam \
--outFileName alignment/SRR576933_sorted_nodup.bedgraph \
--outFileFormat bedgraph --normalizeUsing RPGC \
--effectiveGenomeSize $gsize

bamCoverage --bam alignment/SRR576938_sorted_nodup.bam \
--outFileName alignment/SRR576938_sorted_nodup.bedgraph \
--outFileFormat bedgraph --normalizeUsing RPGC \
--effectiveGenomeSize $gsize


### MOTIF ANALYSIS ###

# Retrieve the peak sequences corresponding to the peak coordinates
bedtools getfasta -fi reference/e_coli_genome.fa \
-bed macs/macs_peaks.narrowPeak \
-fo macs/macs_peaks.fa


duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
