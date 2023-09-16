### DGE ANALYSIS USING LIMMA-VOOM ###


##################################
### LOADING LIBRARIES AND DATA ###
##################################

# Load libraries
library(DESeq2)
library(limma)
library(edgeR)
library(RNAseqQC)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Glimma)

# Function to convert to factor
to_factor <- function(dataframe, variables_to_convert) {
  for (variable in variables_to_convert) {
    dataframe[, variable] <-  factor(dataframe[[variable]])
  }
  return(dataframe)
}

mergeReplicatesCounts <- function(counts, metadata, replicates.variable, variables.to.remove, group) {
  
  # Create a DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ 1
  )
  
  # Variance stabilizing transformation to moderate variance accross the mean
  vst <- vst(dds, blind = TRUE)
  
  # Plot to inspect the stabilization of the variance
  vst.plot <- mean_sd_plot(vst)
  vst_mat <- assay(vst)
  
  # Compute PCA from the VST matrix
  pca <- prcomp(t(vst_mat))
  
  # Extract loads
  loads <- round(pca$sdev^2/sum(pca$sdev^2) * 100, 1)
  
  # Prepare labels
  xlab <- c(paste("PC1", loads[1], "%"))
  ylab <- c(paste("PC2", loads[2], "%"))
  
  
  # Create data frame with metadata and PC1 and PC2 values for input to ggplot
  df <- cbind(metadata, pca$x)
  
  # Create PCA plot
  pca.plot <- ggplot(df, aes_string(
    x="PC1", y="PC2", 
    color = replicates.variable,
    shape = group),
    size = 6) + geom_point() +
    geom_text_repel(aes_string(label = replicates.variable), max.overlaps = 200) +
    xlab(xlab) + 
    ylab(ylab)
  
  # Hierarchical clustering
  # Compute pairwise correlation values
  vst_cor <- cor(vst_mat)
  clust.plot <- pheatmap(vst_cor)
  
  # Collapse technical replicates
  dds.collapsed <- collapseReplicates(
    object = dds,
    groupby = dds[[replicates.variable]],
    run = dds$run
  )
  
  # Extract new counts matrix
  new.counts <- counts(dds.collapsed, normalized = FALSE)
  
  # Extract new metadata
  new.metadata <- as.data.frame(colData(dds.collapsed))
  
  # Remove not needed columns
  new.metadata[, variables.to.remove] <- NULL
  
  return(list(
    "counts" = new.counts, "metadata" = new.metadata, 
    "vst.plot" = vst.plot, "pca.plot" = pca.plot, 
    "clust.plot" = clust.plot
    )
  )
  
}

# Load data
counts <- read.table(
  file = "../quantification/genes_count_matrix.txt",
  header = TRUE
)

# Gene id's to rownames
rownames(counts) <- counts$Geneid

# Remove first column
counts$Geneid <- NULL

# Empty vector to store the name of the samples
sample_names <- c()

# Recover sample names
for (i in 1:length(colnames(counts))) {
  sample_names[i] <- strsplit(colnames(counts)[i], split = "[.]")[[1]][2]
}

# Change sample names
colnames(counts) <- sample_names
saveRDS(counts, file = "../quantification/count_matrix.RDS")

# Load metadata
metadata <- read.table(
  file = "../quantification/metadata.txt",
  header = TRUE
)
metadata <- to_factor(metadata, colnames(metadata))


# Sample column to rownames
rownames(metadata) <- metadata$run

# Check that sample names match in the count matrix and in the metadata
all(colnames(counts) %in% rownames(metadata))
all(colnames(counts) == rownames(metadata))

# Merge technical replicates counts
collapsed.data <- mergeReplicatesCounts(
  counts = counts,
  metadata = metadata,
  replicates.variable = "sample",
  variables.to.remove = c("run", "experiment_accession"),
  group = "group"
)

# Extract counts matrix and metadata
counts.collapsed <- collapsed.data$counts
metadata.collapsed <- collapsed.data$metadata

# Sanity check
all(colnames(counts.collapsed) %in% rownames(metadata.collapsed))
all(colnames(counts.collapsed) == rownames(metadata.collapsed))




##########################
### DATA PREPROCESSING ###
##########################

# Create a DGEList object
dge.list <- DGEList(
  counts = counts.collapsed, 
  samples = metadata.collapsed, 
  genes = rownames(counts.collapsed)
)
colnames(dge.list$genes) <- "ENSEMBL"

bla <-  mapIds(
  x = org.Hs.eg.db,
  keys = dge.list$genes$ENSEMBL, 
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "asNA"
)


bla2 <-  mapIds(
  x = org.Hs.eg.db,
  keys = dge.list$genes$ENSEMBL, 
  column = "SYMBOL",
  keytype = "ENSEMBL"
)

# Add gene symbol to the gene information
dge.list$genes$SYMBOL <- mapIds(
  x = org.Hs.eg.db,
  keys = dge.list$genes$ENSEMBL, 
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# CPM transformation
cpm <- cpm(dge.list)

# Log2 values
lcpm <- cpm(dge.list, log = TRUE)

# Exploration of expression ranges
L <- mean(dge.list$samples$lib.size) / 1000000
M <- median(dge.list$samples$lib.size) / 1000000

c(L, M)

summary(lcpm)

# Number of genes with 0 counts in the 12 samples
table(rowSums(counts.collapsed == 0) == 12)

# Remove genes that are lowly expressed
keep.exprs <- filterByExpr(dge.list, group = "group")
dge.list <- dge.list[keep.exprs, keep.lib.size = FALSE]

# Distribution of counts before and after filtering
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(dge.list)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(dge.list), text.col=col, bty="n")
lcpm <- cpm(dge.list, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(dge.list), text.col=col, bty="n")

# Normalize gene expression distributions
dge.list <- calcNormFactors(dge.list, method = "TMM")

# Compare unnormalized and normalized expression values
dge.list.unnormalized <- dge.list
dge.list.unnormalized$samples$norm.factors <- 1

par(mfrow=c(1,2))
lcpm <- cpm(dge.list.unnormalized, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

lcpm <- cpm(dge.list, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

# Unsupervised clustering of samples
lcpm <- cpm(dge.list, log=TRUE)
par(mfrow=c(1,2))
col.group <- dge.list$samples$group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
groups <- dge.list$samples$group
col.group <- as.character(col.group)
col.plate <- dge.list$samples$plate
plate <- dge.list$samples$plate
levels(col.plate) <-  brewer.pal(nlevels(col.plate), "Set2")
col.plate <- as.character(col.plate)
plotMDS(lcpm, labels=groups, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=plate, col=col.plate, dim=c(3,4))
title(main="B. Pairs")

# Get normalized gene expression
lcpm <- cpm(dge.list, log=TRUE)

# Compute PCA from the log-cpm matrix
pca <- prcomp(t(lcpm))

# Extract loads
loads <- round(pca$sdev^2/sum(pca$sdev^2) * 100, 1)

# Create data frame with metadata and PC1 and PC2 values for input to ggplot
df <- cbind(dge.list$samples, pca$x)

# Create PCA plot
pca.plot <- ggplot(
  data = df, 
  mapping = aes(x=PC1, y=PC2, color = group),
  size = 6
) + geom_point() +
  geom_text_repel(aes_string(label = plate), max.overlaps = 200) +
  xlab(c(paste("PC1", loads[1], "%"))) + 
  ylab(c(paste("PC2", loads[2], "%")))
pca.plot

pca.plot <- ggplot(
  data = df, 
  mapping = aes(x=PC1, y=PC6, color = group),
  size = 6
) + geom_point() +
  geom_text_repel(aes_string(label = plate), max.overlaps = 200) +
  xlab(c(paste("PC1", loads[1], "%"))) + 
  ylab(c(paste("PC6", loads[6], "%")))
pca.plot


pca.plot2 <- ggplot(df, aes(
  x=PC3, y=PC4, 
  color = treatment),
  size = 6) + geom_point() +
  geom_text_repel(aes(label = treatment), max.overlaps = 200) +
  xlab(c(paste("PC3", loads[1], "%"))) + 
  ylab(c(paste("PC4", loads[2], "%")))
pca.plot2

# Hierarchical clustering
# Compute pairwise correlation values
lcpm_cor <- cor(lcpm)
clust.plot <- pheatmap(lcpm_cor)




########################################
### DIFFERENTIAL EXPRESSION ANALYSIS ###
########################################

# Extract the experimental group of each sample
group <- factor(dge.list$samples$group)

# Create design matrix and contrasts
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Set up contrasts
## Questions of interest:
## 1. Which genes respond to the treatment in EMT cell lines?
## 2. Which genes respond differently in EMT cell lines compared to non-EMT cell lines?
contr.matrix <- makeContrasts(
  EMTsal.vs.EMTmock = EMTsal - EMTmock,
  diff = (EMTsal-EMTmock) - (nonEMTsal - nonEMTmock),
  levels = design
)

# Remove heteroscedasticity from count data
# Ref.: https://support.bioconductor.org/p/59700/
# Ref.: https://support.bioconductor.org/p/114663/
par(mfrow=c(1,1))
voom.pre <- voom(counts = dge.list, design = design, plot = TRUE)
corfit <- duplicateCorrelation(
  object = voom.pre,
  design = design,
  block = voom.pre$targets$plate
)
voom <- voom(
  counts = dge.list,
  design = design,
  block = dge.list$samples$plate,
  correlation = corfit$consensus.correlation,
  plot = TRUE
)
corfit <- duplicateCorrelation(voom, design, block = voom$targets$plate)

# Fit model and make the contrasts for each gene
vfit <- lmFit(
  object = voom, 
  design = design,
  block = voom$targets$plate,
  correlation = corfit$consensus.correlation
)

# Compute estimated coefficients and standard errors for the contrasts
vfit <- contrasts.fit(fit = vfit, contrasts = contr.matrix)

# Apply empirical Bayes smoothing to the standard errors (compute moderated t statistic)
efit <- eBayes(vfit)

# Final model mean-variance trend
plotSA(efit, main="Final model: Mean-variance trend")

# Examine the number of DE genes
summary(decideTests(object = efit, adjust.method = "BH", p.value = 0.05))
top.efit.EMTsal.vs.EMTmock <- topTable(
  efit,
  number = nrow(efit),
  coef = "EMTsal.vs.EMTmock",
  adjust = "fdr"
)

top.efit.diff <- topTable(
  efit,
  number = nrow(efit),
  coef = "diff",
  adjust = "fdr"
)


# Consider minimum logFC for significance
tfit <- treat(vfit, lfc = 1)
summary(decideTests(object = tfit, adjust.method = "BH", p.value = 0.05))

# Extract results
top.tfit.EMTsal.vs.EMTmock <- topTreat(fit = tfit, coef = "EMTsal.vs.EMTmock", n = Inf)
top.tfit.EMTsal.vs.EMTmock <- top.tfit.EMTsal.vs.EMTmock %>% 
  filter(adj.P.Val < 0.05) %>% 
  arrange(desc(logFC))
rownames(top.tfit.EMTsal.vs.EMTmock) <- c(1:nrow(top.tfit.EMTsal.vs.EMTmock))

top.tfit.diff <- topTreat(tfit, coef = "diff", n = Inf)
top.tfit.diff <- top.diff %>% 
  filter(adj.P.Val < 0.05) %>% 
  arrange(desc(abs(logFC)))


# Visualize results
dt <- decideTests(tfit)
glimmaMA(
  x = tfit,
  counts = voom$E,
  transform.counts = "none",
  groups = voom$targets$group,
  coef = "EMTsal.vs.EMTmock",
  main = "EMTsal vs EMTmock",
  status = dt[, "EMTsal.vs.EMTmock"]
)

glimmaMA(
  x = tfit,
  counts = voom$E,
  transform.counts = "none",
  groups = voom$targets$group,
  coef = "diff",
  main = "Genes responding differently depending on the cell line",
  status = dt[, "diff"]
)




#############################
### FUNCTIONAL ENRICHMENT ###
#############################

library(clusterProfiler)

genes.emt <- top.tfit.EMTsal.vs.EMTmock$ENSEMBL
genes.non.emt <- top.tfit.nonEMTsal.vs.nonEMTmock$ENSEMBL
genes.diff <- top.tfit.diff$ENSEMBL

ego.emt <- enrichGO(
  gene = genes.emt,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego.emt.results <- ego.emt@result %>% 
  filter(p.adjust < 0.05 & qvalue < 0.05)
View(ego.emt.results)

ego.non.emt <- enrichGO(
  gene = genes.non.emt,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego.non.emt.results <- ego.non.emt@result %>% 
  filter(p.adjust < 0.05 & qvalue < 0.05)
View(ego.non.emt.results)



ego.diff <- enrichGO(
  gene = genes.diff,
  keyType = "ENSEMBL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego.diff.results <- ego.diff@result %>% 
  filter(p.adjust < 0.05 & qvalue < 0.05)
View(ego.diff.results)

save.image(file = "limma_voom_workflow.RData")










