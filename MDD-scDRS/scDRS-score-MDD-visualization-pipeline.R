### scDRS-score-MDD-visualization-pipeline.R ######################################################
# purpose: a script that plot the umap visualization for command line scDRS outputs

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

# gather arguments:
args <- commandArgs(trailingOnly = TRUE);
sc.drs.path <- args[1];
trait.name <- args[2];
plotting.path <- args[3];

# load in the visualization data:
mdd.meta <- read.table(
    sep = '\t',
    header = TRUE,
    row.names = 1,
    file = paste0(data.dir, '2023-05-07-MDD-GSE144136-metadata.txt')
    );

# load in the scDRS score data points:
risk.score <- read.table(
    file = sc.drs.path,
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
    );

### VISUALIZATION UMAP LATENT EMBEDDING ###########################################################
stopifnot(rownames(risk.score) == rownames(mdd.meta));

plot.df <- data.frame(
    umap_1 = mdd.meta$umap_0,
    umap_2 = mdd.meta$umap_1,
    label = mdd.meta$cell.type,
    score = risk.score[, 'norm_score'],
    neg_log10_fdr = -log(p.adjust(risk.score$pval, method = 'fdr')),
    diagnosis = mdd.meta$diagnosis,
    project = 'MDD-GSE144136'
    );
rownames(plot.df) <- rownames(mdd.meta);

# plot out scatter plot:
cell.type.umap <- ggplot(
    data = plot.df,
    aes(x = umap_1, y = umap_2, color = label)) +
    geom_point(size = 1) +
    theme_classic() +
    xlab('umap_1') +
    ylab('umap_2') +
    theme(text = element_text(size = 20)) +
    guides(color = guide_legend(ncol = 1)) +
    labs(fill = 'cell type')

score.umap <- ggplot(
    data = plot.df,
    aes(x = umap_1, y = umap_2, color = score)) +
    scale_color_gradientn(colours = c("#91bfdb" , "#ffffbf", "#fc8d59")) +
    geom_point(size = 1) +
    theme_classic() +
    xlab('umap_1') +
    ylab('umap_2') +
    theme(text = element_text(size = 20)) +
    labs(fill = 'normalized score')

# visualize the score p value:
score.pval.umap <- ggplot(
    data = plot.df,
    aes(x = umap_1, y = umap_2, color = neg_log10_fdr)) +
    scale_color_gradientn(colours = c("#91bfdb" , "#ffffbf", "#fc8d59")) +
    geom_point(size = 1) +
    theme_classic() +
    xlab('umap_1') +
    ylab('umap_2') +
    theme(text = element_text(size = 20)) +
    labs(fill = '-log10(fdr)')

plot.name <- paste0(
    plotting.path,
    system.date,
    '-',
    trait.name,
    '-UMAP-cell-type-scDRS-score.png'
    );

png(
    filename = plot.name,
    width = 30,
    height = 10,
    units = 'in',
    res = 400
    );
plot.list <- list(cell.type.umap, score.umap, score.pval.umap);
print(plot_grid(plotlist = plot.list, nrow = 1));
dev.off();

### VISUALIZE DIAGNOSIS SCORE DISTRIBUTION ########################################################
# make a violin plot for diseased and control cell population:
violin.plot <- ggplot(
    data = plot.df,
    aes(x = diagnosis, y = score, color = diagnosis)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width=0.1) +
    # facet_wrap(~ donor) +
    theme_classic() +
    xlab('diagnosis') +
    ylab('normalized score') +
    theme(text = element_text(size = 20)) +
    labs(fill = 'Diagnosis') +
    ggtitle('score distribution by control and disease status')

# now output the plot:
plot.name <- paste0(
    plotting.path,
    system.date,
    '-',
    trait.name,
    '-disease-score-distribution.png'
    );

png(
    filename = plot.name,
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );

print(violin.plot);
dev.off();

# make a violoin plot for diseased and control cell population by cell type:
cell.type.violin.plot <- ggplot(
    data = plot.df,
    aes(x = diagnosis, y = score, color = diagnosis)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width=0.1) +
    facet_wrap(~ label) +
    theme_classic() +
    xlab('diagnosis') +
    ylab('normalized score') +
    theme(text = element_text(size = 20)) +
    labs(fill = 'Diagnosis') +
    ggtitle('Score distribution by control and disease status per cell type')

# now output the plot:
plot.name <- paste0(
    plotting.path,
    system.date,
    '-',
    trait.name,
    '-per-cell-type-disease-score-distribution.png'
    );

png(
    filename = plot.name,
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );

print(cell.type.violin.plot);
dev.off();

# also plot out the violin plots that doesn't split by diagnosis:
cell.type.violin.plot <- ggplot(
    data = plot.df,
    aes(x = project, y = score)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width=0.1) +
    facet_wrap(~ label) +
    theme_classic() +
    xlab('') +
    ylab('normalized score') +
    theme(text = element_text(size = 20)) +
    ggtitle('Score distribution by cell type')

# now output the plot:
plot.name <- paste0(
    plotting.path,
    system.date,
    '-',
    trait.name,
    '-score-distribution-without-diagnosis-per-cell-type.png'
    );

png(
    filename = plot.name,
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );

print(cell.type.violin.plot);
dev.off();

### ANALYSIS - STATISTICAL T TEST #################################################################
score.hypothesis.test <- function(stat.df) {
    # initiate empty lists:
    cell.types <- unique(stat.df$label);
    cell.type.ttest <- vector('list', length = length(cell.types));
    names(cell.type.ttest) <- cell.types;
    diagnosis.test <- cell.type.ttest;

    # initiate summary table:
    columns <- c('from_zero', 'diagnosis')
    summary.table <- data.frame(matrix(NA, nrow = length(cell.types), ncol = length(columns)));
    rownames(summary.table) <- cell.types;
    colnames(summary.table) <- columns;

    # grab out the control index:
    control.cells <- rownames(stat.df)[stat.df$diagnosis == 'Control'];
    disease.cells <- rownames(stat.df)[stat.df$diagnosis == 'Suicide'];

    for (cell.type in cell.types) {
        # grab out the index:
        cell.type.index <- stat.df$label == cell.type;
        cell.type.index <- rownames(stat.df)[cell.type.index];

        control.index <- intersect(control.cells, cell.type.index);
        disease.index <- intersect(disease.cells, cell.type.index);

        # add check:
        stopifnot(stat.df[disease.index, 'diagnosis'] == 'Suicide');
        stopifnot(stat.df[control.index, 'diagnosis'] == 'Control');
        stopifnot(stat.df[cell.type.index, 'label'] == cell.type);

        # perform statistical test:
        cell.type.ttest[[cell.type]] <- t.test(x = stat.df[cell.type.index, 'score']);
        diagnosis.test[[cell.type]] <- t.test(
            x = stat.df[control.index, 'score'],
            y = stat.df[disease.index, 'score']
            );

        # summarize into a matrix form:
        summary.table[cell.type, 'from_zero'] <- cell.type.ttest[[cell.type]]$p.value
        summary.table[cell.type, 'diagnosis'] <- diagnosis.test[[cell.type]]$p.value
        }

    # also compute global significance:
    global.significance <- t.test(
        x = stat.df[control.cells, 'score'],
        y = stat.df[disease.cells, 'score']
        );
    global.result <- c(NA, global.significance$p.value);
    full.rownames <- c(rownames(summary.table), 'all_cells');
    summary.table <- rbind(summary.table, global.result);
    rownames(summary.table) <- full.rownames;

    # return the summary table:
    return(summary.table);
    }

score.ttest.summary <- score.hypothesis.test(plot.df);
table.name <- paste0(
    plotting.path,
    system.date,
    '-',
    trait.name,
    '-score-distribution-t-test.txt'
    );

write.table(
    x = score.ttest.summary,
    file = table.name,
    sep ='\t',
    quote = FALSE,
    col.names = TRUE,
    row.names = TRUE
    );

### SESSION INFO ##################################################################################
writeLines(
    text = capture.output(sessionInfo()),
    con = paste0(session.save.path, system.date, '-sessionInfo.txt')
    );