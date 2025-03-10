### subtype-analysis.R ############################################################################
# purpose: investigate how our data varies in subtype categories:

### PREAMBLE ######################################################################################
# load in libraries:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);
library(Seurat);

# load in data:
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv'
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

subtype_meta <- data.table::fread(
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/table-S5-cell-meta.csv',
    sep = ',',
    )
subtype_meta <- data.frame(subtype_meta, row.names = 1)

# read in the results:
trait.score <- read.table(
    file = '/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/PASS_MDD_Howard2019.score.gz',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
    )
trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr')
system.date <- Sys.Date()
p.cutoff <- 0.1

# load in all the traits z scores:
# read data:
score.files <- list.files('/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/', pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# read into the empty list:
for (result in score.files) {
    risk.score[[result]] <- read.table(
        file = result,
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );
    }

# simplify list names:
list.names <- gsub('/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/', '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

# load in the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

### PROCESS #######################################################################################
# merge the subtype information onto the actual meta:
meta$subtype <- subtype_meta[rownames(meta), 'SubType']

# make a umap on the L2/3 neurons colored by subtype:
# plot out the spatial region for the L2/3 neurons:
l23.exc <- rownames(meta)[meta$X_MajorType == 'L2/3-IT'];
l23.meta <- meta[l23.exc, ]
l23.meta$subtype <- paste0('L2/3 subtype ', as.numeric(as.factor(l23.meta$subtype)))

plot.df <- data.frame(
    UMAP_1 = l23.meta$UMAP_1,
    UMAP_2 = l23.meta$UMAP_2,
    tissue = as.factor(l23.meta$tissue),
    cell_type = as.factor(l23.meta$X_MajorType),
    met_scdrs = trait.score[l23.exc, 'zscore'],
    fdr = trait.score[l23.exc, 'fdr'],
    sub_type = as.factor(l23.meta$subtype)
    );
rownames(plot.df) = rownames(l23.meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = sub_type)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE215353 L2/3-IT subtypes UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(legend.position="none") +
    theme(text = element_text(size = 20))

gplot.label <- LabelClusters(plot = gplot, id = 'sub_type', col = 'black', size = 5)
plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", system.date, '-GSE215353-l23-subtype-with-legend-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot.label)
dev.off();

# also draw out the met-scDRS:
significant.cell <- rownames(plot.df)[plot.df$fdr < p.cutoff]
insignificant.cell <- setdiff(rownames(plot.df), significant.cell)

gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(data = plot.df[insignificant.cell, ], colour = 'grey') +
    geom_point(data = plot.df[significant.cell, ], aes(colour = met_scdrs)) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
    theme_classic() +
    ggtitle('L2/3-IT MDD Howard 2019') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    labs(color = "met-scDRS") +
    theme(text = element_text(size = 20))

# draw out the plot:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
output.path <- paste0(output.dir, system.date, '-MDD-scDRS-score-umap-l23-only.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

### VISUALIZATION - proportion heatmap ############################################################
# for each of the result, find out the number of significant cells in each cell type:
group.index <- 'subtype'
cell.type <- unique(l23.meta[, group.index]);
significance.matrix <- matrix(NA, nrow = length(risk.score), ncol = length(cell.type));
rownames(significance.matrix) <- names(risk.score);
colnames(significance.matrix) <- cell.type;
trait.class <- trait.info$Category[match(rownames(significance.matrix), trait.info$Trait_Identifier)];

# create the traits by cell type matrix:
plot.type <- 'proportion'
for (result in names(risk.score)) {
    # find the index for significant cells:
    significant.cell <- rownames(risk.score[[result]])[
        p.adjust(risk.score[[result]]$pval, method = 'fdr') < p.cutoff
        ];

    # for each cell type, check number of significant cells are part of that cell type:
    for (type in cell.type) {
        # locate cell id that belong in the cell type:
        cell.type.cell <- rownames(l23.meta)[l23.meta[, group.index] %in% type];

        # find the number of cells that are in each of the cell type category:
        if (plot.type == 'count') {
            significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell);
            } else {
                significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell) /
                    length(cell.type.cell);
            }
        }
    }

# draw out the proportion heatmap:
# select traits to plot:
publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

# reformat traits:
publication.traits <- gsub('PASS_', '', publication.traits)
rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
cell.type.order <- paste0('L2/3 subtype ', seq(1, ncol(significance.matrix)))

# make color function
col.fun <- colorRamp2(
    c(
        0,
        1
        ),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(
    at = c(
        0,
        1
        )
    );

# create heatmap:
plot <- Heatmap(
    as.matrix(significance.matrix)[publication.traits, cell.type.order],
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = cell.type.order,
    width = unit(10 * length(cell.type.order),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-subtype-brain-traits-proportion.png')
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

### Repeat with L5-IT neurons #####################################################################
# plot out the spatial region for the L2/3 neurons:
l5.exc <- rownames(meta)[meta$X_MajorType == 'L5-IT'];
l5.meta <- meta[l5.exc, ]
l5.meta$subtype <- paste0('L5 subtype ', as.numeric(as.factor(l5.meta$subtype)))

plot.df <- data.frame(
    UMAP_1 = l5.meta$UMAP_1,
    UMAP_2 = l5.meta$UMAP_2,
    tissue = as.factor(l5.meta$tissue),
    cell_type = as.factor(l5.meta$X_MajorType),
    met_scdrs = trait.score[l5.exc, 'zscore'],
    fdr = trait.score[l5.exc, 'fdr'],
    sub_type = as.factor(l5.meta$subtype)
    );
rownames(plot.df) = rownames(l5.meta);

# Create plot:
gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2, color = sub_type)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE215353 L5-IT subtypes UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP1') +
    ylab('UMAP2') +
    theme(legend.position="none") +
    theme(text = element_text(size = 20))

gplot.label <- LabelClusters(plot = gplot, id = 'sub_type', col = 'black', size = 5)
plot.path <- paste0("/u/home/l/lixinzhe/project-geschwind/plot/", system.date, '-GSE215353-l5-subtype-with-legend-umap.png')
png(
    filename = plot.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot.label)
dev.off();

# also draw out the met-scDRS:
significant.cell <- rownames(plot.df)[plot.df$fdr < p.cutoff]
insignificant.cell <- setdiff(rownames(plot.df), significant.cell)

gplot <- ggplot(plot.df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(data = plot.df[insignificant.cell, ], colour = 'grey') +
    geom_point(data = plot.df[significant.cell, ], aes(colour = met_scdrs)) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
    theme_classic() +
    ggtitle('L5-IT MDD Howard 2019') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    labs(color = "met-scDRS") +
    theme(text = element_text(size = 20))

# draw out the plot:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
output.path <- paste0(output.dir, system.date, '-MDD-scDRS-score-umap-l5-only.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

### proportion plot for l5-IT #####################################################################
# for each of the result, find out the number of significant cells in each cell type:
group.index <- 'subtype'
cell.type <- unique(l5.meta[, group.index]);
significance.matrix <- matrix(NA, nrow = length(risk.score), ncol = length(cell.type));
rownames(significance.matrix) <- names(risk.score);
colnames(significance.matrix) <- cell.type;
trait.class <- trait.info$Category[match(rownames(significance.matrix), trait.info$Trait_Identifier)];

# create the traits by cell type matrix:
plot.type <- 'proportion'
for (result in names(risk.score)) {
    # find the index for significant cells:
    significant.cell <- rownames(risk.score[[result]])[
        p.adjust(risk.score[[result]]$pval, method = 'fdr') < p.cutoff
        ];

    # for each cell type, check number of significant cells are part of that cell type:
    for (type in cell.type) {
        # locate cell id that belong in the cell type:
        cell.type.cell <- rownames(l5.meta)[l5.meta[, group.index] %in% type];

        # find the number of cells that are in each of the cell type category:
        if (plot.type == 'count') {
            significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell);
            } else {
                significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell) /
                    length(cell.type.cell);
            }
        }
    }

# draw out the proportion heatmap:
# select traits to plot:
publication.traits <- rownames(significance.matrix)[trait.class == 'brain'];

# reformat traits:
publication.traits <- gsub('PASS_', '', publication.traits)
rownames(significance.matrix) <- gsub('PASS_', '', rownames(significance.matrix))
cell.type.order <- paste0('L5 subtype ', seq(1, ncol(significance.matrix)))

# make color function
col.fun <- colorRamp2(
    c(
        0,
        1
        ),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(
    at = c(
        0,
        1
        )
    );

# create heatmap:
plot <- Heatmap(
    as.matrix(significance.matrix)[publication.traits, cell.type.order],
    name = 'Sig. cells',
    col = col.fun,
    rect_gp = gpar(col = "black", lwd = 2),
    row_order = publication.traits,
    column_order = cell.type.order,
    width = unit(10 * length(cell.type.order),"mm"),
    height = unit(10 * length(publication.traits),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    heatmap_legend_param = heatmap.legend.param
    );
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

# use the measured width and height for drawing:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-l5-IT-subtype-brain-traits-proportion.png')
png(
    filename = output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
dev.off();

### chisquare test between subtype and region #####################################################
# for l23:
contingency.table <- table(l23.meta$subtype, l23.meta$tissue)
print(chisq.test(contingency.table, simulate.p.value = TRUE, B = 200000))

# for l5:
contingency.table <- table(l5.meta$subtype, l5.meta$tissue)
print(chisq.test(contingency.table, simulate.p.value = TRUE, B = 200000))