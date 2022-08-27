#!/bin/bash

SECONDS=0

# Canvi al directori de treball
cd ..

# Control de qualitat
fastqc data/exomeSeq_chr22.fastq -o quality/

echo "Quality control on raw data completed."

# Retall dels reads per la cua
java -jar ${HOME}/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-threads 4 \
data/exomeSeq_chr22.fastq data/exomeSeq_chr22_trimmed.fastq \
TRAILING:15 \
MINLEN:50 \
-phred33 
# Especifica la codificació del 'base quality score'

echo "Trimmomatic has finished running."

# # Segon control de qualitat per comprobar les millores
fastqc data/exomeSeq_chr22_trimmed.fastq -o quality/

echo "Quality control on trimmed data finished."

# Canvi al directori del genoma de referència
# cd reference

# wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz' \
# -O chr22.fa.gz

# gunzip chr22.fa.gz

# Creació de l'índex del genoma de referència
# bwa index -p chr22 chr22.fa

# Canvi al directori de treball
# cd ..

echo "Aligning against the reference genome."

# Alineament dels reads amb el cromosoma 22
bwa mem -t 4 \
reference/chr22 \
data/exomeSeq_chr22_trimmed.fastq \
| samtools sort -o alignment/exomeSeq_chr22.bam 

echo "Reads aligned. Building the index of the bam file."

# Construcció de l'índex de l'alineament
samtools index -b alignment/exomeSeq_chr22.bam

echo "Marking duplicates."

# Marcat de duplicats
java -jar ${HOME}/picard/picard.jar MarkDuplicates \
INPUT=alignment/exomeSeq_chr22.bam \
OUTPUT=alignment/exomeSeq_chr22_marked.bam \
METRICS_FILE=metrics.txt \
ASSUME_SORTED=true

echo "Building the index of the bam file."

# Construcció de l'índex de l'alineament de les seqüències no duplicades
samtools index -b alignment/exomeSeq_chr22_marked.bam

# Estadístiques el fitxer bam
samtools stats alignment/exomeSeq_chr22_marked.bam > quality/bam_file_stats.txt

echo "Variant calling."

# Descobriment de variants
freebayes -f reference/chr22.fa \
--min-coverage 7 alignment/exomeSeq_chr22_marked.bam > results/variants.vcf

echo "Filtering and annotating variants."

# Filtració de les variants amb una probabilitat de ser errònia superior a l'1%
java -jar ${HOME}/vcffilter/vcffilter-assembly-0.2.jar \
-I results/variants.vcf \
-o results/variants_filtered.vcf \
--minQualScore 20

# Filtració de les variants que no pertanyen al gen APOBEC3H
cat results/variants_filtered.vcf | \
java -jar ${HOME}/snpEff/SnpSift.jar filter \
"(CHROM = 'chr22') & (POS > 39493228) & (POS < 39500073)" \
> results/variants_apobec3h.vcf

# Predicció dels efectes de les variants
java -jar ${HOME}/snpEff/snpEff.jar -v hg19 \
-stats results/snpEff_apobec3h.html \
results/variants_apobec3h.vcf > results/variants_apobec3h_eff.vcf

# Annotació de les variants d'acord amb la base de dades
java -jar ${HOME}/snpEff/SnpSift.jar annotate \
data/filtered_dbsnp_chr22_hg19.vcf \
results/variants_apobec3h_eff.vcf > results/variants_apobec3h_ann.vcf

# Variants trobades presents a la base de dades
grep -v "#" results/variants_apobec3h_ann.vcf | cut -f 1-3 > results/apobec3h_variants_in_database.txt

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."