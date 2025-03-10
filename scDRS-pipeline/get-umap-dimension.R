### scDRS-score-tsne-visualization.R ##############################################################
# PURPOSE: use the tsne score to plot visualize tSNE embeddings 

### PREAMBLE ######################################################################################
# load libraries:
library(Seurat);
library(ggplot2);
library(cowplot);

# define specified paths:
data.dir <- '/u/project/geschwind/lixinzhe/data/';
session.save.path <- '/u/project/pasaniuc/lixinzhe/session-info/';
system.date <- Sys.Date();
save.path <- '/u/project/pasaniuc/lixinzhe/R_saves/';
code.path <- '/u/home/l/lixinzhe/project-github/methylation-RNA-xinzhe-rotation/code/';

# load data:
combined.obj <- readRDS(paste0(data.dir,'2023-02-17-combined-obj-preprocessed.rds'));

# set seed:
set.seed(103);

### DATA PROCESSING ###############################################################################
# grab out the count matrix:
count <- combined.obj@assays$RNA@counts;

# make a new seurat object:
seurat.count <- CreateSeuratObject(
    counts = count,
    project = "CreateSeuratObject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0
    );

# normalize:
seurat.count <- NormalizeData(
    seurat.count,
    normalization.method = "LogNormalize",
    scale.factor = 10000
    );

# Find variable features:
seurat.count <- FindVariableFeatures(seurat.count, selection.method = "vst", nfeatures = 8000);

# scale data:
seurat.count <- ScaleData(seurat.count, features = rownames(seurat.count));

# dimension reduction:
seurat.count <- RunPCA(seurat.count, features = VariableFeatures(object = seurat.count));
seurat.count <- RunUMAP(seurat.count, dims = 1:20);
seurat.count <- RunTSNE(seurat.count, dims = 1:20);
seurat.count$umap_0 <- seurat.count@reductions$umap@cell.embeddings[, 'UMAP_1'];
seurat.count$umap_1 <- seurat.count@reductions$umap@cell.embeddings[, 'UMAP_2'];

### VISUALIZATION UMAP LATENT EMBEDDING ###########################################################
stopifnot(colnames(seurat.count) == rownames(seurat.count@reductions$umap@cell.embeddings));
stopifnot(colnames(seurat.count) == rownames(combined.obj@meta.data));
stopifnot(identical(combined.obj@assays$RNA@counts, seurat.count@assays$RNA@counts));

# plot.df <- data.frame(
#     umap_1 = combined.obj$umap@cell.embeddings[, 'UMAP_1'],
#     umap_2 = combined.obj$umap@cell.embeddings[, 'UMAP_2'],
#     label = combined.obj$combined.clust.subclass,
#     score = risk.score[rownames(combined.obj@meta.data), 'norm_score'],
#     diagnosis = combined.obj$Diagnosis
#     );

# prepare meta data:
plot.df <- data.frame(
    umap_1 = seurat.count$umap_0,
    umap_2 = seurat.count$umap_1,
    label = combined.obj$combined.clust.subclass,
    diagnosis = combined.obj$Diagnosis
    );

# write out the meta data:
meta.out.path <- paste0()
write.table(
    plot.df,
    file = paste0(data.dir, 'bimodal-asd-visualization-meta.txt'),
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE,
    sep = '\t'
    );
