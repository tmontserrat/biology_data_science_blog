#!/bin/bash

# Ref.: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/index.html
# Ref.: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3#Sec6

### PIPELINE TO PERFORM ATAC SEQ ANALYSIS ###

SECONDS=0
SAMPLE=CRR134566

# Change to the working directory
cd ..

### QUALITY CONTROL ###

# fastqc on the sample
echo "fastqc on ${SAMPLE}."
fastqc data/${SAMPLE}/${SAMPLE}_f1.fq.gz data/${SAMPLE}/${SAMPLE}_r2.fq.gz \
-o quality/fastqc/

echo "Running multiqc." 
multiqc quality/fastqc/ -o quality/multiqc/

# TrimGalore for adapter trimming
echo "TrimGalore on ${SAMPLE}."
trim_galore \
--phred33 \
--quality 20 \
--paired \
--cores 6 \
--output_dir trimmed_data/${SAMPLE}/ \
data/${SAMPLE}/${SAMPLE}_f1.fq.gz data/${SAMPLE}/${SAMPLE}_r2.fq.gz

# fastqc on trimmed fastq
echo "fastqc on ${SAMPLE} trimmed."
fastqc trimmed_data/${SAMPLE}/${SAMPLE}_f1_val_1.fq.gz trimmed_data/${SAMPLE}/${SAMPLE}_r2_val_2.fq.gz \
-o quality/fastqc_trimmed/

##########################




### REFERENCE GENOME ###

# Download
cd alignment/reference/
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/	Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz

# Unzip the file
gunzip Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz

# Create genome index
bowtie2-build alignment/reference/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa alignment/reference/tair10

###################




# ### ALIGNMENT ###

# Align fastq file, sort, outputs a bam file and create an index
echo "Aligning ${SAMPLE} using bowtie2."
bowtie2 -p 6 -q \
--very-sensitive \
-X 500 \
-x alignment/reference/tair10 \
-1 trimmed_data/${SAMPLE}/${SAMPLE}_f1_val_1.fq.gz \
-2 trimmed_data/${SAMPLE}/${SAMPLE}_r2_val_2.fq.gz \
2> quality/alignment/bowtie2/${SAMPLE}_bowtie2_stats.txt \
| samtools view -u - \
| samtools sort - > alignment/${SAMPLE}/${SAMPLE}_sort.bam

# Create index of the bam file
samtools index -b alignment/${SAMPLE}/${SAMPLE}_sort.bam

###############################




### FILTERING ALIGNED READS ###

# Filtering out organelles reads
echo "Filtering out reads aligned to mitochondrial and plastid genome for ${SAMPLE}."
samtools view -h alignment/${SAMPLE}/${SAMPLE}_sort.bam \
| grep -v "Mt" \
| grep -v "Pt" \
| samtools view -b - > alignment/${SAMPLE}/${SAMPLE}_sort_no_org.bam
# samtools index -b alignment/${SAMPLE}/${SAMPLE}_sort_no_org.bam

# Filtering out reads on blacklisted regions
echo "Filtering out reads aligned to blacklisted regions for ${SAMPLE}."
bedtools intersect \
-v -abam alignment/${SAMPLE}/${SAMPLE}_sort_no_org.bam \
-b blacklist_regions/blacklist.bed \
| samtools view -b - > alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black.bam
samtools index -b alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black.bam

# Remove duplicates
echo "Removing duplicated reads for ${SAMPLE}."
java -jar ${HOME}/picard/picard.jar MarkDuplicates \
--INPUT alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black.bam \
--OUTPUT alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup.bam \
--METRICS_FILE quality/alignment/duplicates/${SAMPLE}_duplicates.txt \
--REMOVE_DUPLICATES true

# Remove non-unique alignments
echo "Removing non-unique alignments for ${SAMPLE}."
samtools view -h -b -q 30 alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup.bam \
| samtools view -b - > alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup_q.bam
samtools index -b alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup_q.bam

# Computing quality metrics
echo "Running samtools stats and idxstats for ${SAMPLE}."
samtools stats alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup_q.bam > quality/alignment/stats/${SAMPLE}_sort_no_org_black_dup_q_bam_stats.txt
samtools idxstats alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup_q.bam > quality/alignment/seq_depth/${SAMPLE}_sort_no_org_black_dup_q_bam_seq_depth.txt

# Running mulitqc
echo "Running multiqc for ${SAMPLE}."
multiqc \
quality/fastqc/ \
quality/fastqc_trimmed \
quality/alignment/bowtie2 \
quality/alignment/duplicates \
quality/alignment/seq_depth \
quality/alignment/stats \
-o quality/multiqc/

####################




### SPECIFIC QC FOR ATAC-SEQ DATA ###

# Fragment length distribution
echo "Ploting fragment length distribution for ${SAMPLE}".
java -Xmx16G -jar ${HOME}/picard/picard.jar CollectInsertSizeMetrics \
-I alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup_q.bam \
-O quality/fragment_size_distribution/${SAMPLE}_frag_length.stats \
-H quality/fragment_size_distribution/${SAMPLE}_frag_length.pdf

# Perform quality control using ATACseqQC R package
echo "Running R script to use ATACseqQC package on ${SAMPLE}".
Rscript --vanilla scripts/atac_qc.R ${SAMPLE}

######################




## PEAK CALLING ###

echo "Calling peaks for ${SAMPLE} using macs2 callpeak."
macs2 callpeak --broad -t alignment/${SAMPLE}/${SAMPLE}_sort_no_org_black_dup_q.bam \
-n ${SAMPLE}_macs2 -f BAMPE \
-g 119481543 -q 0.05 --nomodel --keep-dup all \
--outdir results/peaks/




echo "Done."

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
