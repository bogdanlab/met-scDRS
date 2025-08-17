### met-scDRS-score-visualization-region.R ########################################################
# purpose: visualize the met-scDRS score gradient with respect to the region gradient

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(dplyr)
library(data.table)

# define parameters:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date();

# load in data:
# load in meta data:
meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv"
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# load in the scDRS score:
scDRS.file <- "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/PASS_MDD_Howard2019.full_score.gz"
# load in the data
full.score <- fread(
    file = scDRS.file,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE,
    data.table = FALSE
    );
rownames(full.score) <- full.score$cell

### ANALYSIS ######################################################################################
# get the set of cells that we wish to plot:
all(rownames(meta) == rownames(full.score))
meta$score = full.score$zscore
cell_type_tissue_average = meta %>% 
    group_by(X_MajorType, tissue) %>% 
    summarize(average = mean(score), sd = sd(score), cell_num = length(score)) %>% data.frame()
filtered = cell_type_tissue_average[cell_type_tissue_average$cell_num > 50, ]

# lets look at the 
l23.cell <- rownames(meta)[meta$X_MajorType == "L2/3-IT"]
prefrontal.cell <- rownames(meta)[meta$tissue == 'prefrontal cortex']
a46.cell <- rownames(meta)[meta$tissue == 'Brodmann (1909) area 46']
itg.cell <- rownames(meta)[meta$tissue == 'inferior temporal gyrus']
mtg.cell <- rownames(meta)[meta$tissue == 'middle temporal gyrus']
visual.cortex.cell <- rownames(meta)[meta$tissue == 'primary visual cortex']

# first we wish to plot the middle and inferior temporal gyrus met-scDRS score for significant cell
full.score$fdr <- p.adjust(full.score$pval, method = 'fdr')
significant.cell <- rownames(full.score)[full.score$fdr < 1]
insignificant.cell <- setdiff(rownames(full.score), significant.cell)

# generate plot df:
full.df <- data.frame(
    umap1 = meta$UMAP_1,
    umap2 = meta$UMAP_2,
    score = full.score[rownames(meta), 'zscore'],
    tissue = meta$tissue,
    cell_type = meta$X_MajorType
    )
rownames(full.df) <- rownames(meta)

# subset to the set of cells we wish to plot:
plot.cell <- c(
    intersect(l23.cell, a46.cell),
    intersect(l23.cell, prefrontal.cell),
    intersect(l23.cell, itg.cell),
    intersect(l23.cell, mtg.cell),
    intersect(l23.cell, visual.cortex.cell)
    )
plot.cell <- intersect(plot.cell, significant.cell)

# lets plot this:
plot.df <- full.df[plot.cell, ]
plot.df$tissue <- factor(
    plot.df$tissue,
    levels = c(
        'primary visual cortex',
        'prefrontal cortex',
        'Brodmann (1909) area 46',
        'inferior temporal gyrus',
        'middle temporal gyrus')
        )
gplot <- ggplot(plot.df, aes(x = tissue, y = score, fill = tissue)) +
    geom_boxplot() +
    scale_fill_manual(
        values = c(
            'primary visual cortex' = '#fdc086',
            'Brodmann (1909) area 46' = '#fc8d62',
            'prefrontal cortex' = '#8da0cb',
            'inferior temporal gyrus' = '#66c2a5',
            'middle temporal gyrus' = '#e78ac3'
            )
        ) +
    xlab('tissue') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20))

output.path <- paste0(output.dir, system.date, '-MDD-l23-lg-A46-itg-mtg-scDRS-score-boxplot.png')
png(
    filename = output.path,
    width = 8,
    height = 7,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

###########################################################################################
######                                    for MGC                                    ######
###########################################################################################
# lets look at the 
mgc.cell <- rownames(meta)[meta$X_MajorType == "MGC"]
prefrontal.cell <- rownames(meta)[meta$tissue == 'prefrontal cortex']
a46.cell <- rownames(meta)[meta$tissue == 'Brodmann (1909) area 46']
itg.cell <- rownames(meta)[meta$tissue == 'inferior temporal gyrus']
mtg.cell <- rownames(meta)[meta$tissue == 'middle temporal gyrus']
visual.cortex.cell <- rownames(meta)[meta$tissue == 'primary visual cortex']

# first we wish to plot the middle and inferior temporal gyrus met-scDRS score for significant cell
full.score$fdr <- p.adjust(full.score$pval, method = 'fdr')
significant.cell <- rownames(full.score)[full.score$fdr < 1]
insignificant.cell <- setdiff(rownames(full.score), significant.cell)

# generate plot df:
full.df <- data.frame(
    umap1 = meta$UMAP_1,
    umap2 = meta$UMAP_2,
    score = full.score[rownames(meta), 'zscore'],
    tissue = meta$tissue,
    cell_type = meta$X_MajorType
    )
rownames(full.df) <- rownames(meta)

# subset to the set of cells we wish to plot:
plot.cell <- c(
    intersect(mgc.cell, a46.cell),
    intersect(mgc.cell, prefrontal.cell),
    intersect(mgc.cell, itg.cell),
    intersect(mgc.cell, mtg.cell),
    intersect(mgc.cell, visual.cortex.cell)
    )
plot.cell <- intersect(plot.cell, significant.cell)

# lets plot this:
plot.df <- full.df[plot.cell, ]
plot.df$tissue <- factor(
    plot.df$tissue,
    levels = c(
        'primary visual cortex',
        'prefrontal cortex',
        'Brodmann (1909) area 46',
        'inferior temporal gyrus',
        'middle temporal gyrus')
        )
gplot <- ggplot(plot.df, aes(x = tissue, y = score, fill = tissue)) +
    geom_boxplot() +
    scale_fill_manual(
        values = c(
            'primary visual cortex' = '#fdc086',
            'Brodmann (1909) area 46' = '#fc8d62',
            'prefrontal cortex' = '#8da0cb',
            'inferior temporal gyrus' = '#66c2a5',
            'middle temporal gyrus' = '#e78ac3'
            )
        ) +
    xlab('tissue') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20))

output.path <- paste0(output.dir, system.date, '-MDD-MGC-lg-A46-itg-mtg-scDRS-score-boxplot.png')
png(
    filename = output.path,
    width = 8,
    height = 7,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

###########################################################################################
######                                        For OPC                                ######
###########################################################################################
# for OPC:
# lets look at the 
opc.cell <- rownames(meta)[meta$X_MajorType == "OPC"]
prefrontal.cell <- rownames(meta)[meta$tissue == 'prefrontal cortex']
a46.cell <- rownames(meta)[meta$tissue == 'Brodmann (1909) area 46']
itg.cell <- rownames(meta)[meta$tissue == 'inferior temporal gyrus']
mtg.cell <- rownames(meta)[meta$tissue == 'middle temporal gyrus']
visual.cortex.cell <- rownames(meta)[meta$tissue == 'primary visual cortex']

# first we wish to plot the middle and inferior temporal gyrus met-scDRS score for significant cell
full.score$fdr <- p.adjust(full.score$pval, method = 'fdr')
significant.cell <- rownames(full.score)[full.score$fdr < 1]
insignificant.cell <- setdiff(rownames(full.score), significant.cell)

# generate plot df:
full.df <- data.frame(
    umap1 = meta$UMAP_1,
    umap2 = meta$UMAP_2,
    score = full.score[rownames(meta), 'zscore'],
    tissue = meta$tissue,
    cell_type = meta$X_MajorType
    )
rownames(full.df) <- rownames(meta)

# subset to the set of cells we wish to plot:
plot.cell <- c(
    intersect(opc.cell, a46.cell),
    intersect(opc.cell, prefrontal.cell),
    intersect(opc.cell, itg.cell),
    intersect(opc.cell, mtg.cell),
    intersect(opc.cell, visual.cortex.cell)
    )
plot.cell <- intersect(plot.cell, significant.cell)

# lets plot this:
plot.df <- full.df[plot.cell, ]
plot.df$tissue <- factor(
    plot.df$tissue,
    levels = c(
        'primary visual cortex',
        'prefrontal cortex',
        'Brodmann (1909) area 46',
        'inferior temporal gyrus',
        'middle temporal gyrus')
        )
gplot <- ggplot(plot.df, aes(x = tissue, y = score, fill = tissue)) +
    geom_boxplot() +
    scale_fill_manual(
        values = c(
            'primary visual cortex' = '#fdc086',
            'Brodmann (1909) area 46' = '#fc8d62',
            'prefrontal cortex' = '#8da0cb',
            'inferior temporal gyrus' = '#66c2a5',
            'middle temporal gyrus' = '#e78ac3'
            )
        ) +
    xlab('tissue') +
    ylab('met-scDRS') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(text = element_text(size = 20))

output.path <- paste0(output.dir, system.date, '-MDD-opc-lg-A46-itg-mtg-scDRS-score-boxplot.png')
png(
    filename = output.path,
    width = 8,
    height = 7,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();
