### cell-type-umap.R ##############################################################################
# purpose: plot out the umap that is colored by cell type:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(circlize);
library(Seurat);

# define paths:
date <- Sys.Date()
meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv"
plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-GSE215353-cell-type-umap.png')

# load in the data:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

### VISUALIZE the UMAP plot #######################################################################
# create plot df:
plot.df <- data.frame(
    UMAP_1 = meta$UMAP_1,
    UMAP_2 = meta$UMAP_2,
    cell_type = as.factor(meta$X_MajorType)
    );
rownames(plot.df) = rownames(meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE215353 UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 20)) +
    theme(legend.position="none")
gplot.label <- LabelClusters(plot = gplot, id = 'cell_type', col = 'black', size = 5)

png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot.label)
dev.off();

# plot out the spatial region for the L2/3 neurons:
l23.exc <- rownames(meta)[meta$X_MajorType == 'L2/3-IT'];
l23.meta <- meta[l23.exc, ]

plot.df <- data.frame(
    UMAP_1 = l23.meta$UMAP_1,
    UMAP_2 = l23.meta$UMAP_2,
    tissue = as.factor(l23.meta$tissue),
    cell_type = as.factor(l23.meta$X_MajorType)
    );
rownames(plot.df) = rownames(l23.meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = tissue)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE215353 UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 20))

plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-GSE215353-tissue-with-legend-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = tissue)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE215353 UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(text = element_text(size = 20)) +
    theme(legend.position="none")

plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-GSE215353-tissue-no-legend-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();
