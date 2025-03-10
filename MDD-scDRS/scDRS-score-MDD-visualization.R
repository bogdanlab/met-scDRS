### scDRS-score-MDD-visualization.R ###############################################################
# purpose: visualize the latent embedding of MDD scRNA seq and color based on umap and scDRS

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
plotting.path <- '/u/project/pasaniuc/lixinzhe/plot/';

# load data:
combined.obj <- readRDS(paste0(save.path,'GSE144136-seurat-meta.rds'));

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

# plot.df <- data.frame(
#     umap_1 = combined.obj$umap@cell.embeddings[, 'UMAP_1'],
#     umap_2 = combined.obj$umap@cell.embeddings[, 'UMAP_2'],
#     label = combined.obj$combined.clust.subclass,
#     score = risk.score[rownames(combined.obj@meta.data), 'norm_score'],
#     diagnosis = combined.obj$Diagnosis
#     );

# crate the plotting data frame:
plot.df <- data.frame(
    umap_1 = seurat.count$umap_0,
    umap_2 = seurat.count$umap_1,
    label = combined.obj$cell.type
    );

# create the umap plot:
cell.type.umap <- ggplot(
    data = plot.df,
    aes(x = umap_1, y = umap_2, color = label)) +
    geom_point(size = 1) +
    theme_classic() +
    xlab('umap_1') +
    ylab('umap_2') +
    theme(text = element_text(size = 20)) +
    labs(fill = 'cell type')

plot.name <- paste0(
    plotting.path,
    system.date,
    '-UMAP-cell-type-MDD-scRNA.png'
    );

png(
    filename = plot.name,
    width = 15,
    height = 10,
    units = 'in',
    res = 400
    );
print(cell.type.umap);
dev.off();

### SAVE ##########################################################################################
combined.obj$umap_0 <- seurat.count@reductions$umap@cell.embeddings[, 'UMAP_1'];
combined.obj$umap_1 <- seurat.count@reductions$umap@cell.embeddings[, 'UMAP_2'];
mdd.meta <- combined.obj@meta.data;

# write the mdd meta to disk:
write.table(
    x = mdd.meta,
    sep = '\t',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE,
    file = paste0(data.dir, system.date, '-MDD-GSE144136-metadata.txt')
    );