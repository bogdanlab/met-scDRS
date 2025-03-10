system.date <- Sys.Date();
meta.data.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx";
mch.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE132489-full/";
mcg.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mcg/GSE132489-full/";
p.cutoff <- 0.1;
directories <- c(mch.directory, mcg.directory);
names(directories) <- c('mch', 'mcg');

# load libraries:
library(ggplot2);
library(ComplexHeatmap);
library(circlize);
library(readxl);
library(data.table);
library(Seurat);

# load meta data:
meta <- read_excel(meta.data.path, skip = 14);
meta <- as.data.frame(meta);
rownames(meta) <- meta[, '...1'];
colnames(meta)[1] <- 'cell';
meta <- meta[meta[, 'Pass QC'], ];

### Load in the scDRS run #########################################################################
### LEVEL 1:
# loop over modality for enrichment plot:
output.dir <- '/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/';
for (modality in c('mch')){

    scDRS.directory <- directories[modality];

    score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
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
    list.names <- gsub(scDRS.directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;

    for (trait in names(risk.score)){
        # first subset to the trait score of interest:
        trait.score <- risk.score[[trait]];
        trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr');
        significant.cell <- rownames(trait.score)[trait.score$fdr < p.cutoff];
        insignificant.cell <- setdiff(rownames(trait.score), significant.cell);

        # next, we will plot out the umap:
        plot.df <- trait.score;
        level <- 'L1'
        xaxis <- paste0(level, 'UMAP_0');
        yaxis <- paste0(level, 'UMAP_1');
        plot.df$umap1 <- meta[rownames(plot.df), xaxis];
        plot.df$umap2 <- meta[rownames(plot.df), yaxis];

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

        # draw out the plot:
        output.path <- paste0(output.dir, system.date, '-', modality, '-', trait, '-level-',level, '-scDRS-score-umap.png')
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
    }

# plot a umap 0 and 1 plot that looks at the cell class level:
plot.df$major_type <- meta[rownames(plot.df), 'MajorType'];
plot.df$class <- meta[rownames(plot.df), 'CellClass'];
gplot <- ggplot(plot.df, aes(x = umap1, y = umap2, color = class)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE132489 UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('L1_UMAP_1') +
    ylab('L1_UMAP_1') +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE132489-umap-class.png');
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

plot.df$cell_type <- as.factor(plot.df$major_type)
gplot <- ggplot(plot.df, aes(x = umap1, y = umap2, color = cell_type)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE132489 UMAP') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('L1_UMAP_1') +
    ylab('L1_UMAP_1') +
    theme(text = element_text(size = 20)) +
    theme(legend.position="none")
gplot.label <- LabelClusters(plot = gplot, id = 'cell_type', col = 'black', size = 5)

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE132489-umap-cell-type.png');
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot.label)
dev.off();

### Seperate out the inhibitory neurons and excitatory neurons and plot umap just with those cells:
# plot only excitatory neurons:
subset.meta <- meta[rownames(plot.df), ];
excitatory.cell <- rownames(subset.meta)[subset.meta$CellClass == 'Exc'];

gplot <- ggplot(subset.meta[excitatory.cell, ], aes(x = L1UMAP_0, y = L1UMAP_1, color = MajorType)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE132489 UMAP Excitatory Neurons') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('L1_UMAP_1') +
    ylab('L1_UMAP_1') +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE132489-level-1-excitatory-major-type.png');
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

# next for the inhibitory neurons, plot the the level one umap:
subset.meta <- meta[rownames(plot.df), ];
inhibitory.cell <- rownames(subset.meta)[subset.meta$CellClass == 'Inh'];

gplot <- ggplot(subset.meta[inhibitory.cell, ], aes(x = L1UMAP_0, y = L1UMAP_1, color = MajorType)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE132489 UMAP Inhibitory Neurons') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('L1_UMAP_1') +
    ylab('L1_UMAP_1') +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE132489-level-1-inh-major-type.png');
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

### LEVEL 1:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
# for diseases:
for (modality in c('mch')){

    scDRS.directory <- directories[modality];

    score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
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
    list.names <- gsub(scDRS.directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;

    for (trait in c("PASS_MDD_Howard2019", "PASS_BIP_Mullins2021")){
        # first subset to the trait score of interest:
        trait.score <- risk.score[[trait]];
        trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr');
        significant.cell <- rownames(trait.score)[trait.score$fdr < p.cutoff];
        insignificant.cell <- setdiff(rownames(trait.score), significant.cell);

        # next, we will plot out the umap:
        plot.df <- trait.score;
        level <- 'L1'
        xaxis <- paste0(level, 'UMAP_0');
        yaxis <- paste0(level, 'UMAP_1');
        plot.df$umap1 <- meta[rownames(plot.df), xaxis];
        plot.df$umap2 <- meta[rownames(plot.df), yaxis];

        insig.neuron.plot <- intersect(insignificant.cell, inhibitory.cell);
        sig.neuron.plot <- intersect(significant.cell, inhibitory.cell)

        gplot <- ggplot(plot.df[inhibitory.cell, ], aes(x = umap1, y = umap2)) +
            geom_point(data = plot.df[insig.neuron.plot, ], colour = 'grey') +
            geom_point(data = plot.df[sig.neuron.plot, ], aes(colour = zscore)) +
            scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
            theme_classic() +
            ggtitle(gsub('_', ' ',trait)) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab(gsub('_', ' ',xaxis)) +
            ylab(gsub('_', ' ',yaxis)) +
            theme(text = element_text(size = 20))

        # draw out the plot:
        output.path <- paste0(output.dir, system.date, '-', modality, '-', trait, '-level-',level, '-scDRS-score-umap.png')
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
    }


# for diseases:
for (modality in c('mch')){

    scDRS.directory <- directories[modality];

    score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
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
    list.names <- gsub(scDRS.directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;

    for (trait in c("PASS_MDD_Howard2019", "PASS_BIP_Mullins2021")){
        # first subset to the trait score of interest:
        trait.score <- risk.score[[trait]];
        trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr');
        significant.cell <- rownames(trait.score)[trait.score$fdr < p.cutoff];
        insignificant.cell <- setdiff(rownames(trait.score), significant.cell);

        # next, we will plot out the umap:
        plot.df <- trait.score;
        level <- 'L1'
        xaxis <- paste0(level, 'UMAP_0');
        yaxis <- paste0(level, 'UMAP_1');
        plot.df$umap1 <- meta[rownames(plot.df), xaxis];
        plot.df$umap2 <- meta[rownames(plot.df), yaxis];

        insig.neuron.plot <- intersect(insignificant.cell, excitatory.cell);
        sig.neuron.plot <- intersect(significant.cell, excitatory.cell)

        gplot <- ggplot(plot.df[excitatory.cell, ], aes(x = umap1, y = umap2)) +
            geom_point(data = plot.df[insig.neuron.plot, ], colour = 'grey') +
            geom_point(data = plot.df[sig.neuron.plot, ], aes(colour = zscore)) +
            scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
            theme_classic() +
            ggtitle(gsub('_', ' ',trait)) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab(gsub('_', ' ',xaxis)) +
            ylab(gsub('_', ' ',yaxis)) +
            theme(text = element_text(size = 20))

        # draw out the plot:
        output.path <- paste0(output.dir, system.date, '-', modality, '-', trait, '-level-',level, '-exc-scDRS-score-umap.png')
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
    }


### LEVEL 2:
# next for the inhibitory neurons, plot the the level two umap:
gplot <- ggplot(subset.meta[excitatory.cell, ], aes(x = L2UMAP_0, y = L2UMAP_1, color = MajorType)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE132489 UMAP Excitatory Neurons') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('L2_UMAP_1') +
    ylab('L2_UMAP_1') +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE132489-level-2-excitatory-major-type.png');
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

gplot <- ggplot(subset.meta[inhibitory.cell, ], aes(x = L2UMAP_0, y = L2UMAP_1, color = MajorType)) +
    geom_point() +
    theme_classic() +
    ggtitle('GSE132489 UMAP Inhibitory Neurons') +
    theme(plot.title = element_text(hjust=0.5)) +
    xlab('L2_UMAP_1') +
    ylab('L2_UMAP_1') +
    theme(text = element_text(size = 20))

output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE132489-level-2-inh-major-type.png');
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

# for diseases:
for (modality in c('mch')){

    scDRS.directory <- directories[modality];

    score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
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
    list.names <- gsub(scDRS.directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;

    for (trait in c("PASS_MDD_Howard2019", "PASS_BIP_Mullins2021")){
        # first subset to the trait score of interest:
        trait.score <- risk.score[[trait]];
        trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr');
        significant.cell <- rownames(trait.score)[trait.score$fdr < p.cutoff];
        insignificant.cell <- setdiff(rownames(trait.score), significant.cell);

        # next, we will plot out the umap:
        plot.df <- trait.score;
        level <- 'L2'
        xaxis <- paste0(level, 'UMAP_0');
        yaxis <- paste0(level, 'UMAP_1');
        plot.df$umap1 <- meta[rownames(plot.df), xaxis];
        plot.df$umap2 <- meta[rownames(plot.df), yaxis];

        insig.neuron.plot <- intersect(insignificant.cell, inhibitory.cell);
        sig.neuron.plot <- intersect(significant.cell, inhibitory.cell)

        gplot <- ggplot(plot.df[inhibitory.cell, ], aes(x = umap1, y = umap2)) +
            geom_point(data = plot.df[insig.neuron.plot, ], colour = 'grey') +
            geom_point(data = plot.df[sig.neuron.plot, ], aes(colour = zscore)) +
            scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
            theme_classic() +
            ggtitle(gsub('_', ' ',trait)) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab(gsub('_', ' ',xaxis)) +
            ylab(gsub('_', ' ',yaxis)) +
            theme(text = element_text(size = 20))

        # draw out the plot:
        output.path <- paste0(output.dir, system.date, '-', modality, '-', trait, '-level-',level, '-inh-scDRS-score-umap.png')
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
    }


# for diseases:
for (modality in c('mch')){

    scDRS.directory <- directories[modality];

    score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
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
    list.names <- gsub(scDRS.directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;

    for (trait in c("PASS_MDD_Howard2019", "PASS_BIP_Mullins2021")){
        # first subset to the trait score of interest:
        trait.score <- risk.score[[trait]];
        trait.score$fdr <- p.adjust(trait.score$pval, method = 'fdr');
        significant.cell <- rownames(trait.score)[trait.score$fdr < p.cutoff];
        insignificant.cell <- setdiff(rownames(trait.score), significant.cell);

        # next, we will plot out the umap:
        plot.df <- trait.score;
        level <- 'L2'
        xaxis <- paste0(level, 'UMAP_0');
        yaxis <- paste0(level, 'UMAP_1');
        plot.df$umap1 <- meta[rownames(plot.df), xaxis];
        plot.df$umap2 <- meta[rownames(plot.df), yaxis];

        insig.neuron.plot <- intersect(insignificant.cell, excitatory.cell);
        sig.neuron.plot <- intersect(significant.cell, excitatory.cell)

        gplot <- ggplot(plot.df[excitatory.cell, ], aes(x = umap1, y = umap2)) +
            geom_point(data = plot.df[insig.neuron.plot, ], colour = 'grey') +
            geom_point(data = plot.df[sig.neuron.plot, ], aes(colour = zscore)) +
            scale_color_gradient(low = "#fee0d2", high = "#de2d26") +
            theme_classic() +
            ggtitle(gsub('_', ' ',trait)) +
            theme(plot.title = element_text(hjust=0.5)) +
            xlab(gsub('_', ' ',xaxis)) +
            ylab(gsub('_', ' ',yaxis)) +
            theme(text = element_text(size = 20))

        # draw out the plot:
        output.path <- paste0(output.dir, system.date, '-', modality, '-', trait, '-level-',level, '-exc-scDRS-score-umap.png')
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
    }

