### cell-type-umap.R ##############################################################################
# purpose: plot out the umap that is colored by cell type:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(circlize);
library(Seurat);

# define paths:
date <- Sys.Date()
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv'
plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", date, '-GSE215353-cell-type-umap.pdf')

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

pdf(
    file = plot.path,
    width = 14,
    height = 14
    );
print(gplot.label)
dev.off();
