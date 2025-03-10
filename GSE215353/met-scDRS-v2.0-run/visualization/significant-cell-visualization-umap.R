### significant-cells-umap-visualization-script.R #################################################
# purpose: for a directory containing scDRS scores, 
# visualize the number of cells that are significant across all traits
# also visualize the proportion of cell types that are significant given a meta file

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    significant-cells-umap-visualization-script.R [--dir <scdrs> --meta_data <meta> --xaxis <xaxis> --yaxis <yaxis> --out <output> --cutoff <p>]

Options:
    --dir directory path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score
    --xaxis column name in meta that encode coordinate for dimensionality reduction 1
    --yaxis column name in meta that encode coordinate for dimensionality reduction 2
    --cutoff p value cutoff that user specifies
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.directory <- opts$dir;
xaxis <- opts$xaxis;
yaxis <- opts$yaxis;
output.dir <- opts$out;
p.cutoff <- as.numeric(opts$cutoff);
system.date <- Sys.Date();

# function testing:
# meta.data.path <- "/u/project/geschwind/lixinzhe/data/bimodal-asd-visualization-meta.txt"
# scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/fraction-ASD-methyl/mch/met-scDRS-v1-run/LXZ-74-traits/"
# xaxis <- 'umap_1'
# yaxis <- 'umap_2'
# group.index <- 'label'
# p.cutoff <- 0.1
# output.path <- '/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/test-plot.png'
# trait <- 'UKB_460K.cov_EDU_YEARS'

# load libraries:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);

# read data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

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
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

# print some results:
cat('number of risk scores loaded:', length(risk.score), '\n')

# subset meta data to the same set as cell:
stopifnot(rownames(risk.score[[result]]) %in% meta)

### visualize on the effect size on umap ##########################################################
for (trait in names(risk.score)){
    # first subset to the trait score of interest:
    trait.score <- risk.score[[trait]];
    trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr');
    significant.cell <- rownames(trait.score)[trait.score$fdr < p.cutoff];
    insignificant.cell <- setdiff(rownames(trait.score), significant.cell);

    # next, we will plot out the umap:
    plot.df <- trait.score;
    plot.df$umap1 <- meta[rownames(plot.df), xaxis];
    plot.df$umap2 <- meta[rownames(plot.df), yaxis];

    if (p.cutoff < 1) {
        gplot <- ggplot(plot.df, aes(x = umap1, y = umap2)) +
            geom_point(data = plot.df[insignificant.cell, ], colour = 'grey') +
            geom_point(data = plot.df[significant.cell, ], aes(colour = zscore)) +
            scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
            theme_classic() +
            ggtitle(gsub('_', ' ',trait)) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab(gsub('_', ' ',xaxis)) +
            ylab(gsub('_', ' ',yaxis)) +
            theme(text = element_text(size = 20))
        }
    if (p.cutoff == 1) {
        gplot <- ggplot(plot.df, aes(x = umap1, y = umap2)) +
            geom_point(data = plot.df[insignificant.cell, ], colour = 'grey') +
            geom_point(data = plot.df[significant.cell, ], aes(colour = zscore)) +
            scale_color_gradient(low = "#2c7bb6", high = "#d7191c") +
            theme_classic() +
            ggtitle(gsub('_', ' ',trait)) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab(gsub('_', ' ',xaxis)) +
            ylab(gsub('_', ' ',yaxis)) +
            theme(text = element_text(size = 20))
        }

    # draw out the plot:
    output.path <- paste0(output.dir, system.date, '-', trait, '-scDRS-score-umap.png')
    png(
        filename = output.path,
        width = 14,
        height = 14,
        units = 'in',
        res = 400
        );
    print(gplot)
    dev.off();
}
