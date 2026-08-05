### GSE215353-CpG-MDD-umap.R ######################################################################
# purpose: plot the umap for gse215353 on CpG methylation modality:

### PREAMBLE ######################################################################################
# load in packages:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);

# load in data:
# meta = read.table('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/meta_data_50k.tsv', sep = '\t', row.names = 1)
meta = read.table("/u/home/l/lixinzhe/project-geschwind/port/scratch/met_scdrs_dev/joint_umap_coords_281146.csv", sep = ',', row.names = 1, header = TRUE)

# read in the results:
trait.score <- read.table(
    file = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/ges215353_50k/qc_regression_gene_body_CGN/full/PASS_MDD_Howard2019.score.gz',
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
    )
rownames(trait.score) = gsub('.allc.tsv.gz', '', rownames(trait.score))
trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr')
system.date <- Sys.Date()


### VISUALIZATION #################################################################################
p.cutoff <- 0.1
significant.cell <- rownames(trait.score)[trait.score$fdr < p.cutoff]
insignificant.cell <- setdiff(rownames(trait.score), significant.cell)

shared_cells = intersect(rownames(meta), rownames(trait.score))
plot.df <- trait.score[shared_cells, ]
plot.df$umap1 <- meta[shared_cells, 'X_umap']
plot.df$umap2 <- meta[shared_cells, 'Y_umap']

gplot <- ggplot(plot.df, aes(x = umap1, y = umap2)) +
    geom_point(data = plot.df[insignificant.cell, ], colour = 'grey') +
    geom_point(data = plot.df[significant.cell, ], aes(colour = zscore)) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
    theme_classic() +
    ggtitle('MDD Howard 2019') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    labs(color = "met-scDRS") +
    theme(text = element_text(size = 20))

# draw out the plot:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
output.path <- paste0(output.dir, system.date, '-MDD-CpG-scDRS-score-umap.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();
