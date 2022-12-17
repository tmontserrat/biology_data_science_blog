#!/bin/bash -li
shopt -s expand_aliases

SECONDS=0

#### WORKING ON wt_H3K4me3_read1.fastq and wt_H3K4me3_read2.fastq ####

### QUALITY CONTROL ###

# Move to working directory
cd ..

# Quality control
echo ""
fastqc data/wt_H3K4me3_read1.fastq data/wt_H3K4me3_read2.fastq -o quality/

# Remove low quality bases
trimmomatic PE -threads 1 -phred33 \
data/wt_H3K4me3_read1.fastq data/wt_H3K4me3_read2.fastq \
data/wt_H3K4me3_read1_trimmed_paired.fastq data/wt_H3K4me3_read1_trimmed_unpaired.fastq \
data/wt_H3K4me3_read2_trimmed_paired.fastq data/wt_H3K4me3_read2_trimmed_unpaired.fastq \
TRAILING:20 MINLEN:20

# Quality control post-trimming
fastqc data/wt_H3K4me3_read1_trimmed_paired.fastq data/wt_H3K4me3_read2_trimmed_paired.fastq -o quality/


### MAPPING READS ###

# Download reference genome
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz  -O reference/mm10.fa.gz 
gunzip reference/mm10.fa.gz

# Building the index
bowtie2-build reference/mm10.fa reference/mm10

# Mapping the reads
# bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]

echo "Mapping reads to reference genome."

bowtie2 -p 1 -q \
-x reference/mm10 \
-1 data/wt_H3K4me3_read1_trimmed.fastq -2 data/wt_H3K4me3_read2_trimmed.fastq \
-S alignment/wt_H3K4me3.sam

echo "Converting to bam, sorting and building an index."

# Convert to bam and sort
samtools view -b alignment/wt_H3K4me3.sam | samtools sort -o alignment/wt_H3K4me3.bam

# Index of the bam file
samtools index -b alignment/wt_H3K4me3.bam

# Index all the bam files we have
ls alignment/bam_files/*.bam | xargs -n 1 -P 1 samtools index -b

### ChIP-SEQ QUALITY CONTROL ###

# Compute read coverages over a large number of regions
echo "Calculating coverage per region."
multiBamSummary bins --bamfiles alignment/bam_files/*.bam \
--binSize 1000 \
--distanceBetweenBins 500 \
--region chrX \
--outFileName quality/readCounts.npz \
--outRawCounts quality/readCounts.tab

# Correlation between samples (using deepTools)
echo "Plotting correlation among samples."
plotCorrelation --corData quality/readCounts.npz \
--corMethod pearson \
--whatToPlot heatmap \
--plotFile quality/samples_correlation.pdf

# IP strength estimation
INPUT1="alignment/bam_files/wt_input_rep1.bam"
INPUT2="alignment/bam_files/wt_input_rep2.bam"
SAMPLES1=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep1")
SAMPLES2=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep2")

echo "Creating fingerprint plot for replicates 1."
plotFingerprint \
-b ${SAMPLES1} ${INPUT1} \
--labels CTCF H3K27me3 H3K4me3 Input \
--region chrX \
--numberOfSamples 10000 \
--plotFile quality/fingerprint_rep1.pdf

echo "Creating fingerprint plot for replicates 2."
plotFingerprint \
-b ${SAMPLES2} ${INPUT2} \
--labels CTCF H3K27me3 H3K4me3 Input \
--region chrX \
--numberOfSamples 10000 \
--plotFile quality/fingerprint_rep2.pdf


### NORMALIZATION ###

# Estimation of the sequencing depth
SAMPLES=$(ls alignment/bam_files/*.bam)
for i in $SAMPLES
do 
    samtools idxstats $i > alignment/sequencing_depth/`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`.txt
done

# Normalize by sequencing depth
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
GENOMESIZE=2308125349
for i in $SAMPLES
do
    bamCoverage --bam $i \
    --outFileName alignment/normalized_coverage/`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`.bedgraph \
    --outFileFormat bedgraph --normalizeUsing RPGC \
    --effectiveGenomeSize $GENOMESIZE \
    --binSize 25 \
    --region chrX
done

# Check regions with the highest coverage
sort -k4nr alignment/normalized_coverage/wt_H3K4me3_rep1.bedgraph | head

# Generate bigWig files
for i in $SAMPLES
do
    bamCoverage --bam $i \
    --outFileName alignment/normalized_coverage/`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`.bw \
    --outFileFormat bigwig --normalizeUsing RPGC \
    --effectiveGenomeSize $GENOMESIZE \
    --binSize 25 \
    --region chrX
done

# Input normalization
INPUT1="alignment/bam_files/wt_input_rep1.bam"
INPUT2="alignment/bam_files/wt_input_rep2.bam"
SAMPLES1=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep1")
SAMPLES2=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep2")

# bedgraph
echo "bedgraph"
echo "Starting Rep 1"
for i in $SAMPLES1
do
    echo " "
    echo $i
    bamCompare -b1 $i -b2 $INPUT1 \
    --binSize 50 \
    --operation log2 \
    --outFileName alignment/input_normalized_coverage/`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`.bedgraph \
    --outFileFormat bedgraph \
    --region chrX
done

echo "Rep 1 done. Starting Rep 2."

for i in $SAMPLES2
do
    echo " "
    echo $i
    bamCompare -b1 $i -b2 $INPUT2 \
    --binSize 50 \
    --operation log2 \
    --outFileName alignment/input_normalized_coverage/`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`.bedgraph \
    --outFileFormat bedgraph \
    --region chrX
done

# bigwig
echo "bigwig"
echo "Starting Rep 1"
for i in $SAMPLES1
do
    echo " "
    echo $i
    bamCompare -b1 $i -b2 $INPUT1 \
    --binSize 50 \
    --operation log2 \
    --outFileName alignment/input_normalized_coverage/`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`.bigwig \
    --outFileFormat bigwig \
    --region chrX
done

echo "Rep 1 done. Starting Rep 2."

for i in $SAMPLES2
do
    echo " "
    echo $i
    bamCompare -b1 $i -b2 $INPUT2 \
    --binSize 50 \
    --operation log2 \
    --outFileName alignment/input_normalized_coverage/`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`.bigwig \
    --outFileFormat bigwig \
    --region chrX
done

Check regions with the highest coverage
sort -k4gr alignment/input_normalized_coverage/wt_H3K4me3_rep1.bedgraph | head


### PEAK CALLING ###
SAMPLES1narrow=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep1" | grep -v "H3K27me")
SAMPLES1broad=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep1" | grep "H3K27me")
SAMPLES2narrow=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep2" | grep -v "H3K27me")
SAMPLES2broad=$(ls alignment/bam_files/*.bam | grep -v "input" | grep "rep2" | grep "H3K27me")

for i in $SAMPLES1narrow
do
    echo " "
    echo $i
    macs2 callpeak --treatment $i \
    --control $INPUT1 \
    --format BAMPE \
    --gsize mm \
    --name macs/macs_`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`
done

for i in $SAMPLES1broad
do
    echo " "
    echo $i
    macs2 callpeak --treatment $i \
    --control $INPUT1 \
    --format BAMPE \
    --gsize mm \
    --broad \
    --name macs/macs_`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`
done

echo " "
echo "Rep 1 done. Starting rep2."
echo " "

for i in $SAMPLES2narrow
do
    echo " "
    echo $i
    macs2 callpeak --treatment $i \
    --control $INPUT2 \
    --format BAMPE \
    --gsize mm \
    --name macs/macs_`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`
done

for i in $SAMPLES2broad
do
    echo " "
    echo $i
    macs2 callpeak --treatment $i \
    --control $INPUT2 \
    --format BAMPE \
    --gsize mm \
    --broad \
    --name macs/macs_`echo $i | cut -d "/" -f 3 | cut -d "." -f 1`
done

Check fold change of the peak located at chrX151519887-151522945
cat macs/macs_wt_H3K4me3_rep1_peaks.narrowPeak | grep 151519886 | cut -f 1,2,3,7

# Check number of peaks in rep 1 and rep 2
echo "Number of peaks rep 1:"
cat macs/macs_wt_H3K4me3_rep1_peaks.narrowPeak | wc -l
echo "Number of peaks rep 2:"
cat macs/macs_wt_H3K4me3_rep2_peaks.narrowPeak | wc -l

# Filter out those peaks with a fold change less than 50
awk '$7 > 50 {print $0}' macs/macs_wt_H3K4me3_rep1_peaks.narrowPeak | sort -k 4 -gr | head


### PLOTTING THE SIGNAL BETWEEN SAMPLES ###

# Concatenate the peaks from CTCG and H3K4me3 rep1
NARROWPEAKS=$(ls macs/*.narrowPeak)

echo ${NARROWPEAKS}

cat ${NARROWPEAKS} \
| cut -f 1-5 \
> results/concatenated_peaks/concatenated_narrowpeaks.bed

# Sort by chromosome and then by start position in ascending order
sortBed -i results/concatenated_peaks/concatenated_narrowpeaks.bed > results/concatenated_peaks/concatenated_narrowpeaks_sorted.bed

# Merge the regions of both files
mergeBed -i results/concatenated_peaks/concatenated_narrowpeaks_sorted.bed > results/merged_peaks/merged_narrowpeaks.bed

# Compute the scores for each region
SAMPLES=$(ls alignment/input_normalized_coverage/*.bigwig | grep -v "H3K27me3")

echo ${SAMPLES}

computeMatrix reference-point \
--regionsFileName results/merged_peaks/merged_narrowpeaks.bed \
--scoreFileName ${SAMPLES} \
--outFileName results/region_scores_matrix/matrix_narrowpeaks.gz \
--referencePoint center \
--upstream 3000 \
--downstream 3000

echo "Creating heatmap."

# Plot the heatmap
plotHeatmap -m results/region_scores_matrix/matrix_narrowpeaks.gz \
--outFileName results/heatmap/heatmap_narrowpeaks.pdf \
--refPointLabel center \
--legendLocation upper-left \
--kmeans 2 \
--yMax 7


### CLOSEST GENES ###

## Find shared peaks ##

# For H3K4me3
bedtools intersect \
-a macs/macs_wt_H3K4me3_rep1_peaks.narrowPeak \
-b macs/macs_wt_H3K4me3_rep2_peaks.narrowPeak \
> results/treatments_peaks/H3K4me3_overlaps.bed

# For H3K27me3
bedtools intersect \
-a macs/macs_wt_H3K27me3_rep1_peaks.broadPeak \
-b macs/macs_wt_H3K27me3_rep2_peaks.broadPeak \
> results/treatments_peaks/H3K27me3_overlaps.bed

# For CTCF
bedtools intersect \
-a macs/macs_wt_CTCF_rep1_peaks.narrowPeak \
-b macs/macs_wt_CTCF_rep2_peaks.narrowPeak \
> results/treatments_peaks/CTCF_overlaps.bed

## Download and unzip the gtf file ##
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz  -O annotations/mm10.refGene.gtf.gz 
gunzip annotations/mm10.refGene.gtf

# Convert gtf to bed ##

# Filter to keep chromosome X
awk '{ if ($1 == "chrX") { print } }' annotations/mm10.refGene.gtf \
> annotations/mm10.refGene.chrX.gtf

gtf2bed < annotations/mm10.refGene.chrX.gtf > annotations/mm10.refGene.chrX.bed

# Get closest genes ##

# Sort peaks
sort -k 1,1 -k2,2n results/treatments_peaks/CTCF_overlaps.bed \
> results/treatments_peaks/CTCF_overlaps_sorted.bed
sort -k 1,1 -k2,2n results/treatments_peaks/H3K27me3_overlaps.bed \
> results/treatments_peaks/H3K27me3_overlaps_sorted.bed
sort -k 1,1 -k2,2n results/treatments_peaks/H3K4me3_overlaps.bed \
> results/treatments_peaks/H3K4me3_overlaps_sorted.bed

# Change chromosome name from chrX to X
awk 'BEGIN { OFS = "\t" } $1="X"' annotations/mm10.refGene.chrX.bed \
| cut -f 1-9 \
> annotations/mm10.refGene.X.bed

# Find closest/overlapping genes from each peak
bedtools closest -io \
-a results/treatments_peaks/H3K4me3_overlaps_sorted.bed \
-b annotations/mm10.refGene.X.bed \
| cut -f 14 | sort | uniq \
> results/closest_genes/closest_genes_H3K4me3.bed

bedtools closest \
-a results/treatments_peaks/H3K27me3_overlaps_sorted.bed \
-b annotations/mm10.refGene.X.bed \
| cut -f 13 | sort | uniq \
> results/closest_genes/closest_genes_H3K27me3.bed

bedtools closest -io \
-a results/treatments_peaks/CTCF_overlaps_sorted.bed \
-b annotations/mm10.refGene.X.bed \
| cut -f 14 | sort | uniq \
> results/closest_genes/closest_genes_CTCF.bed

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."