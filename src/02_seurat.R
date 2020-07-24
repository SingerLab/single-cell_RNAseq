## following Seurat tutorial
## libraries
library(Seurat)
library(dplyr)

pbmc.data <- Read10X(data.dir = "DD1592_01/outs/filtered_feature_bc_matrix/")

( dense.size <- object.size(x = as.matrix(x = pbmc.data)) )
( sparse.size <- object.size(x = pbmc.data) )

dense.size / sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
##pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, 
##    project = "DD1592_10x_02")
## save.image("dd1592_10x_02.seurat.rda")
load("dd1592_10x_02.seurat.rda")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

## plots

pdf("Rplots.pdf")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
dev.off()


# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc,
                    subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf),
                    high.thresholds = c(8000, 0.05))


# By default, we employ a global-scaling normalization method “LogNormalize”
# that normalizes the gene expression measurements for each cell by the
# total expression, multiplies this by a scale factor (10,000 by default),
# and log-transforms the result.

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)


pdf("Rplots2.pdf", width = 12, height = 8)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()

save.image("dd1592_10x_02.seurat.rda")

## requires internet connection 
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"), do.par = FALSE, num.cores = 1)

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)

pdf("Rplots3.pdf")
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
dev.off()

# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)

PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, num.genes = 40, do.balanced = TRUE, label.columns = TRUE)

PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
         label.columns = FALSE, use.full = FALSE)

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = TRUE)

pdf("Rplots4.pdf", width = 24, height = 24)
JackStrawPlot(object = pbmc, PCs = 1:20)
PCElbowPlot(object = pbmc)
dev.off()


# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:16,
                     resolution = 0.8, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)

# While we do provide function-specific printing functions, the more general
# function to print calculation parameters is PrintCalcParams().

# Run Non-linear dimensional reduction (tSNE)
# Seurat continues to use tSNE as a powerful tool to visualize and explore
# these datasets. While we no longer advise clustering directly on tSNE
# components, cells within the graph-based clusters determined above should
# co-localize on the tSNE plot. This is because the tSNE aims to place 
# cells with similar local neighborhoods in high-dimensional space together
# in low-dimensional space. As input to the tSNE, we suggest using the same
# PCs as input to the clustering analysis, although computing the tSNE
# based on scaled gene expression is also supported using the genes.use
# argument.

pbmc <- RunTSNE(object = pbmc, dims.use = 1:16, do.fast = TRUE)

# note that you can set do.label=T to help label individual clusters
pdf("tSNE_DD1592_dim16_res0.8.pdf", width = 298/25.4, height = 210/25.4)
TSNEPlot(object = pbmc)
dev.off()

## save
saveRDS(pbmc, file = "DD1592_02_10x.seurat_tutorial.rds")


## Finding differentially expressed genes (cluster biomarkers)
# Seurat can help you find markers that define clusters via differential
# expression. By default, it identifes positive and negative markers of
# a single cluster (specified in ident.1), compared to all other cells.
# FindAllMarkers automates this process for all clusters, but you can
# also test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a gene to be detected at a minimum
# percentage in either of the two groups of cells, and the thresh.test
# argument requires a gene to be differentially expressed (on average)
# by some amount between the two groups. You can set both of these to
# 0, but with a dramatic increase in time - since this will test a
# large number of genes that are unlikely to be highly discriminatory.
# As another option to speed up these computations, max.cells.per.ident
# can be set. This will downsample each identity class to have no more
# cells than whatever this is set to. While there is generally going
# to be a loss in power, the speed increases can be significiant and
# the most highly differentially expressed genes will likely still
# rise to the top.

# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = 0:4, 
    min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
    thresh.use = 0.25)
