### SCRIPT TO INSPECT QUALITY OF THE ATAC-SEQ DATA USING ATACseqQC ###
# Ref.: https://nbis-workshop-epigenomics.readthedocs.io/en/latest/index.html


### LOADING LIBRARIES AND DATA ###

# Load libraries
suppressMessages(library("ATACseqQC"))
suppressMessages(library("BSgenome.Athaliana.TAIR.TAIR9"))
suppressMessages(library("TxDb.Athaliana.BioMart.plantsmart28"))
suppressMessages(library("ChIPpeakAnno"))
suppressMessages(library("Rsamtools"))

wd <- "/media/tomas/DATA/biodatascience/atac_seq/"

# Path to bam file
bamFileLabels <- "CRR134566"
bamFileLabels <- commandArgs(trailingOnly=TRUE)
print(bamFileLabels)
bamFile <- paste0(wd, "alignment/", bamFileLabels, "/", bamFileLabels, "_sort_no_org_black_dup_q.bam")




### STARTING ANALYSIS ###

# # Collect library statistics
# bam_qc <- bamQC(bamfile = bamFile, outPath = NULL)

# # Output to save files
outPath <- paste0(wd, "quality/ATACseqQC/", bamFileLabels)
# dir.create(outPath)

### Collect information on which BAM tags are present in out bam file
# Generate all possible tags
possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))

# Collect those present in our bam file
print("Collecting tags.")
bamTop100 <- scanBam(
  file = BamFile(bamFile, yieldSize = 100),
  param = ScanBamParam(tag = possibleTag)
)[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100) > 0]

#######


# Let's work just with autosome 
seqlev <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")

# Create granges object for each chromosome
which <- as(seqinfo(Athaliana)[seqlev], "GRanges")

# Change style of chromosome names to match the ones in the bam file
newStyle <- mapSeqlevels(seqlevels(which), "NCBI")
which <- renameSeqlevels(which, newStyle)

# Create an object with genomic alignments
print("Creating GAlignmentsList object.")
gal <- readBamFile(
  bamFile = bamFile, tag = tags, which = which, asMates = TRUE, bigFile = TRUE
)

### SHIFT ALIGNMENTS AND SPLIT THEM BY LENGTHS ###

print("Shifting alignment.")

# Output path
shiftedBamFile <- paste0(wd, "alignment/",  bamFileLabels, "/", bamFileLabels, "_shifted.bam")

# Shift alignments +4 and -5
gal_shifted_gr <- shiftGAlignmentsList(gal, outbam = shiftedBamFile)
saveRDS(
  gal_shifted_gr,
  file = paste0(wd, "r_objects/gal_shifted_gr_", bamFileLabels, ".RDS"),
  ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL
)
# gal_shifted_gr <- readRDS(file = paste0(wd, "r_objects/gal_shifted_gr_", bamFileLabels,".RDS"))


### Extract genomic locations of TSS
# Transcripts
txs <- transcripts(TxDb.Athaliana.BioMart.plantsmart28)
txs <- txs[seqnames(txs) %in% newStyle]
# genome <- Athaliana
# newStyleGenome <- mapSeqlevels(seqlevels(genome), "NCBI")
# genome <- renameSeqlevels(genome, newStyleGenome)

# Retrieve TSS of each transcript
TSS <- promoters(txs, upstream = 0, downstream = 1)
TSS <- unique(TSS)
export(TSS, paste0(outPath, "/tss.bed"))


### SIGNAL IN NFR AND MONONUCLEOSOME FRACTIONS

print("Spliting NFR and mono-, di- tri- nucleosome fraction.")
# outPath <- "."
# Split alignements by length
objs <- splitGAlignmentsByCut(
  gal_shifted_gr, 
  # txs = txs, genome = genome, 
  outPath = outPath
)
saveRDS(
  objs, file = paste0(wd, "r_objects/objs_", bamFileLabels, ".RDS")
)
objs <- readRDS(file = paste0(wd, "r_objects/objs_", bamFileLabels,".RDS"))


# Path to split bam files by length
bamFiles <- paste0(
  outPath, 
  c("/NucleosomeFree.bam", 
    "/mononucleosome.bam",
    "/dinucleosome.bam",
    "/trinucleosome.bam")
)

# Claculate signal around TSS and log transform
# librarySize <- estLibSize(bamFiles)
librarySize <- c(
  length(objs$NucleosomeFree)/2,
  length(objs$mononucleosome)/2,
  length(objs$dinucleosome)/2,
  length(objs$trinucleosome)/2
)
names(librarySize) <- c("NucleosomeFree", "mononucleosome", "dinucleosome", "trinucleosome")
print("Library size of each fraction:")
print(librarySize)

NTILE <- 101L
ups <- 1010L
dws <- 1010L

print("Computing TSS enrichment and creating plots.")
sigs <- enrichedFragments(
  gal = objs[c("NucleosomeFree", "mononucleosome", "dinucleosome", "trinucleosome")],
  # bamfiles = bamFiles,
  TSS = TSS,
  librarySize = librarySize,
  seqlev = "1",
  upstream = ups,
  downstream = dws,
  n.tile = NTILE,
  TSS.filter=0.5
)

sigs_logs2 <- lapply(sigs, function(.ele) log2(.ele + 1))


pdf(paste0(wd,  "quality/ATACseqQC/", bamFileLabels, "/", bamFileLabels, "_heatmap_splitbam.pdf"))
featureAlignedHeatmap(
  cvglists = sigs_logs2, 
  feature.gr = reCenterPeaks(TSS, width=ups+dws),
  zeroAt=.5,
  n.tile=NTILE, 
  color = colorRampPalette(c("white", "blue"))(50)
)
dev.off()

# Signal at TSS
out <- featureAlignedDistribution(
  cvglists = sigs,
  feature.gr = reCenterPeaks(TSS, width=ups+dws),
  zeroAt=0.5,
  n.tile=NTILE, type="l",
  ylab="Averaged coverage"
)

# Rescale the nucleosome-free and nucleosome signals to 0-1 for plotting
range01 <- function(x) {
  (x-min(x))/(max(x)-min(x))
}
out <- apply(out, 2, range01)

# Plot it
pdf(paste0(wd,  "quality/ATACseqQC/", bamFileLabels, "/", bamFileLabels, "_TSS_profile_splitbam.pdf"))
matplot(out, type="l", xaxt="n",
        xlab="Position (bp)",
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1,
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()




print("Done.")