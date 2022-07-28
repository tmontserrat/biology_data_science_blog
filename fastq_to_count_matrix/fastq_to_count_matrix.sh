#!/bin/bash

SECONDS=0

# Canvi al directory de treball
cd ..

# Pas 1: control de qualitat
fastqc data/demo.fastq -o results/

echo "Quality control on raw data completed."

# Retall dels reads per la cua
java -jar trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-threads 4 \
data/demo.fastq data/demo_trimmed.fastq TRAILING:10 \
MINLEN:50 \
-phred33 # Especifica la codificació del 'base quality score'

echo "Trimmomatic has finished running."

# # Segon control de qualitat per comprobar les millores
fastqc data/demo_trimmed.fastq -o results/

echo "Quality control on trimmed data finished."

# Pas 2: alineament de les seqüències curtes amb HISAT2

echo "Downloading genome index."

# Descàrrega de l'índex del genoma humà
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz \
# -P hisat2/
# tar -xvzf hisat2/grch38_genome.tar.gz

# echo "Genome index downloaded and ready to use."

echo "Aligning reads."

# Alineament de les seqüències
hisat2 -q --rna-strandness R \ # El protocol es 'strandness'
-x hisat2/grch38/genome \ # Path a l'índex
-U data/demo_trimmed.fastq \ # Short reads
| samtools sort -o hisat2/demo_trimmed.bam # Ordenem i convertim en bam

echo "Alignement finished."

# Pas 3: quantificació de l'expressió gènica (taula de recompte)

# echo "Downloading annotation file."

# # Descàrrega del fitxer d'anotacions
# wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz \
# -P data/

# # Descompressió del fitxer
# gunzip data/Homo_sapiens.GRCh38.107.gtf.gz

# echo "Annotation file downloaded. Quantifying features."

# Taula de recompte amb featureCounts
featureCounts -s 2 \
-a data/Homo_sapiens.GRCh38.107.gtf \
-o quants/demo_featurecounts.txt \
data/demo_trimmed.bam

echo "featureCounts finished."

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
