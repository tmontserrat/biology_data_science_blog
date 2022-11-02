#!/bin/bash

SECONDS=0

# Change to the working directory
cd ..

# Quality control
fastqc data/exomeSeq_chr22.fastq -o quality/

echo "Quality control on raw data completed."

# Trimming of the reads
java -jar ${HOME}/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-threads 4 \
data/exomeSeq_chr22.fastq data/exomeSeq_chr22_trimmed.fastq \
TRAILING:15 \
MINLEN:50 \
-phred33 

echo "Trimmomatic has finished running."

# Quality control post-trimmed reads
fastqc data/exomeSeq_chr22_trimmed.fastq -o quality/

echo "Quality control on trimmed data finished."

# Change to the reference directory
# cd reference

# wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz' \
# -O chr22.fa.gz

# gunzip chr22.fa.gz

# Create the index of the reference genome
# bwa index -p chr22 chr22.fa

# Go back to the working directory
# cd ..

echo "Aligning against the reference genome."

# Align the reads against the chromosome 22
bwa mem -t 4 \
reference/chr22 \
data/exomeSeq_chr22_trimmed.fastq \
| samtools sort -o alignment/exomeSeq_chr22.bam 

echo "Reads aligned. Building the index of the bam file."

# Create the index of the bam file
samtools index -b alignment/exomeSeq_chr22.bam

echo "Marking duplicates."

# Mark duplicate reads
java -jar ${HOME}/picard/picard.jar MarkDuplicates \
INPUT=alignment/exomeSeq_chr22.bam \
OUTPUT=alignment/exomeSeq_chr22_marked.bam \
METRICS_FILE=metrics.txt \
ASSUME_SORTED=true

echo "Building the index of the bam file."

# Build the index of the not duplicate reads
samtools index -b alignment/exomeSeq_chr22_marked.bam

# Bam file statistics
samtools stats alignment/exomeSeq_chr22_marked.bam > quality/bam_file_stats.txt

echo "Variant calling."

# Variant calling
freebayes -f reference/chr22.fa \
--min-coverage 7 alignment/exomeSeq_chr22_marked.bam > results/variants.vcf

echo "Filtering and annotating variants."

# Filter variants with a probability to be false positives higher than 1%
java -jar ${HOME}/vcffilter/vcffilter-assembly-0.2.jar \
-I results/variants.vcf \
-o results/variants_filtered.vcf \
--minQualScore 20

# Filter those variants not belonging to the gene APOBEC3H
cat results/variants_filtered.vcf | \
java -jar ${HOME}/snpEff/SnpSift.jar filter \
"(CHROM = 'chr22') & (POS > 39493228) & (POS < 39500073)" \
> results/variants_apobec3h.vcf

# Predict the effects of the variants
java -jar ${HOME}/snpEff/snpEff.jar -v hg19 \
-stats results/snpEff_apobec3h.html \
results/variants_apobec3h.vcf > results/variants_apobec3h_eff.vcf

# Annotate the variants cross-referencing with a database
java -jar ${HOME}/snpEff/SnpSift.jar annotate \
data/filtered_dbsnp_chr22_hg19.vcf \
results/variants_apobec3h_eff.vcf > results/variants_apobec3h_ann.vcf

# Filter those variants present in the database
grep -v "#" results/variants_apobec3h_ann.vcf | cut -f 1-3 > results/apobec3h_variants_in_database.txt

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
