#!/bin/bash

SECONDS=0

# Change to working directory
cd ..

# Step 1: quality control
for i in SRR8435265 SRR8435266 SRR8435267 SRR8435268 SRR8435269 SRR8435270 SRR8435271 SRR8435272 SRR8435273 SRR8435274 SRR8435275 SRR8435276 SRR8435277 SRR8435278 SRR8435279 SRR8435280 SRR8435281 SRR8435282 SRR8435283 SRR8435284 SRR8435285 SRR8435286 SRR8435287 SRR8435288 SRR8435289 SRR8435290 SRR8435291 SRR8435292 SRR8435293 SRR8435294 SRR8435295 SRR8435296 SRR8435297 SRR8435298 SRR8435299 SRR8435300 SRR8435301 SRR8435302 SRR8435303 SRR8435304 SRR8435305 SRR8435306 SRR8435307 SRR8435308 SRR8435309 SRR8435310 SRR8435311 SRR8435312;
do
    echo "fastqc on ${i}"
    fastqc data/${i}/${i}_1.fastq.gz data/${i}/${i}_2.fastq.gz -o quality/fastqc/
done

echo "Quality control on raw data completed."

# multiqc quality/fastqc/ -o quality/multiqc/

# Step 2: alignment
for i in SRR8435265 SRR8435266 SRR8435267 SRR8435268 SRR8435269 SRR8435270 SRR8435271 SRR8435272 SRR8435273 SRR8435274 SRR8435275 SRR8435276 SRR8435277 SRR8435278 SRR8435279 SRR8435280 SRR8435281 SRR8435282 SRR8435283 SRR8435284 SRR8435285 SRR8435286 SRR8435287 SRR8435288 SRR8435289 SRR8435290 SRR8435291 SRR8435292 SRR8435293 SRR8435294 SRR8435295 SRR8435296 SRR8435297 SRR8435298 SRR8435299 SRR8435300 SRR8435301 SRR8435302 SRR8435303 SRR8435304 SRR8435305 SRR8435306 SRR8435307 SRR8435308 SRR8435309 SRR8435310 SRR8435311 SRR8435312;
do
    echo "Aligning ${i}."
    hisat2 -x ../hisat2/grch38/genome \
    -1 data/${i}/${i}_1.fastq.gz -2 data/${i}/${i}_2.fastq.gz \
    -t \
    -p 6 \
    --rna-strandness RF \
    | samtools sort -o alignment/${i}.bam
done

# Step 3: alignment statistics
for i in SRR8435265 SRR8435266 SRR8435267 SRR8435268 SRR8435269 SRR8435270 SRR8435271 SRR8435272 SRR8435273 SRR8435274 SRR8435275 SRR8435276 SRR8435277 SRR8435278 SRR8435279 SRR8435280 SRR8435281 SRR8435282 SRR8435283 SRR8435284 SRR8435285 SRR8435286 SRR8435287 SRR8435288 SRR8435289 SRR8435290 SRR8435291 SRR8435292 SRR8435293 SRR8435294 SRR8435295 SRR8435296 SRR8435297 SRR8435298 SRR8435299 SRR8435300 SRR8435301 SRR8435302 SRR8435303 SRR8435304 SRR8435305 SRR8435306 SRR8435307 SRR8435308 SRR8435309 SRR8435310 SRR8435311 SRR8435312;
do

    echo "samtools stats on ${i}."
    samtools stats alignment/${i}.bam > quality/alignment/${i}.txt

done

# multiqc quality/fastqc/ quality/alignment/ -o multiqc_quality/

# Step 4: quantification

# Sanity check chromosome names
# head ../gtf/Homo_sapiens.GRCh38.107.gtf | grep -v "#" | cut -f1,3,4,5 
# samtools index alignment/SRR8435265.bam
# samtools idxstats alignment/SRR8435265.bam | cut -f1-3 | head -n5

# Get the count matrix
echo "Quantifying gene expression with featureCounts."
featureCounts \
-T 6 \
-p \
-a ../gtf/Homo_sapiens.GRCh38.107.gtf \
-o quantification/count_matrix.txt \
-t exon \
-s 2 \
alignment/*.bam

# Remove not interested columns
cat quantification/count_matrix.txt | cut -f1,7-55 | \
grep -v "#" > quantification/genes_count_matrix.txt

# Add gene expression matrix to multiqc report
multiqc quality/fastqc/ quality/alignment/ quantification/ -o quality/multiqc/

echo "Done."

duration=$SECONDS

echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
