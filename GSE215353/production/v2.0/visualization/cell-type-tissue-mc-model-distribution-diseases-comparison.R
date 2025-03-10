### cell-type-tissue-mc-model-distribution-diseases-comparison.R ##################################
# purpose: quantify cell type specific dissection region signal across disorders

### PREAMBLE ######################################################################################
# load in packages:
require(docopt)
require(circlize)
require(ggplot2)
require(tidyverse)
require(data.table)
require(dplyr)
require(ComplexHeatmap)

# define parameters:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date()
zscore.count <- 2 # will count how many tissue, disection pattern with > zscore.count

# specify data path:
result.dir <- "/u/scratch/l/lixinzhe/tmp-file/region-cell-type-analysis-sig-only/"
disease.dir <- list.dirs(result.dir, full.names = FALSE, recursive = FALSE)

# load in the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

# make a progress bar:
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(disease.dir),
    clear = FALSE
    );

### ANALYSIS ######################################################################################
# find cell-type : tissue pairs that are with z score > 3 in each diseases

# create some placeholders:
disease.z.count <- rep(NA, times = length(disease.dir))
names(disease.z.count) <- disease.dir
z.score.storage <- vector('list', length = length(disease.dir))
names(z.score.storage) <- disease.dir

# for each disease: load in data and compute z score
for (disease in disease.dir){
    ## first load in the data
    # get the mc model result directories:
    model.result.dir <- paste0(result.dir, disease, '/')
    mc.results <- list.files(model.result.dir, pattern = '*-model-summary.csv');

    # load in all the cell type : tissue pair from the disease:
    result.collection <- NULL;
    for (result.file in mc.results) {
        pair.index <- gsub('-MC-linear-model-summary.csv', '', result.file)
        result.collection[[pair.index]] <- data.frame(
                fread(
                file = paste0(model.result.dir, result.file),
                sep = ',',
                data.table = FALSE
                ),
            row.names = 1
            )
        }
    
    ## next compute z score:
    # create a placeholder:
    permuted.value <- rep(NA, length = length(result.collection))
    names(permuted.value) <- names(result.collection)
    beta <- zscore <- permuted.value

    # compute the p value and z score:
    for (result in names(result.collection)) {
        foreground.effect <- result.collection[[result]]['foreground', 'beta']
        if (is.na(foreground.effect)) {
            permuted.value[result] <- NA;
            beta[result] <- NA;
            } else {
                numerator <- 1 + sum(abs(result.collection[[result]][, 'beta']) > abs(foreground.effect))
                denominator <- nrow(result.collection[[result]])
                pval <- numerator / denominator
                permuted.value[result] <- pval
                beta[result] <- result.collection[[result]]["foreground", 'beta']

                # also compute the z score using (foreground - background mean) / background sd
                # grep out the control mc draws:
                control.index <- grep('ctrl_norm_score_*', rownames(result.collection[[result]]))
                z.mean <- mean(result.collection[[result]][control.index, 'beta'])
                z.sd <- sd(result.collection[[result]][control.index, 'beta'])
                zscore[result] <- (foreground.effect - z.mean) / z.sd
            }
        }

    # count how many z score result that is over zscore.count:
    disease.z.count[disease] <- sum(zscore > zscore.count, na.rm = TRUE)

    # store the z score:
    z.score.storage[[disease]] <- zscore

    # update progress bar:
    pb$tick()
    }

### VISUALIZE #####################################################################################
# bar plot on the top 3 diseases for each trait categories:
plot.df <- data.frame(
    disease = names(disease.z.count),
    count = disease.z.count,
    category = trait.info$Category[match(names(disease.z.count), trait.info$Trait_Identifier)]
    )
plot.df <- plot.df %>% group_by(category) %>% summarize(average = mean(count), sd = sd(count));
plot.df <- as.data.frame(plot.df)

# make visualization that give me cell type by trait:
gplot <- ggplot(
    plot.df,
    aes(x = category, y = average)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width=.2,
                    position=position_dodge(.9)) +
    theme_classic() +
    xlab('Traits categories') +
    ylab('average of tissue:cell type pairs with high z score') +
    theme(legend.position = "none") +
    theme(text = element_text(size = 30))

output.path <- paste0(output.dir, system.date, '-average of tissue cell type pairs with high z score in MC models.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

### ANALYSIS ######################################################################################
# we look at the distribution for brain traits:
brain.traits <- trait.info$Trait_Identifier[trait.info$Category == 'brain']
non.brain.traits <- setdiff(trait.info$Trait_Identifier, brain.traits)
brain.distribution <- NULL
for (disease in brain.traits) {
    brain.distribution <- c(
        brain.distribution, 
        names(na.omit(z.score.storage[[disease]][z.score.storage[[disease]] > zscore.count]))
        )
    cat(disease, 'have high z score in ', names(na.omit(z.score.storage[[disease]][z.score.storage[[disease]] > zscore.count])), '\n')
    }

# look at the distribution for non brain traits:
nonbrain.distribution <- NULL
for (disease in non.brain.traits) {
    nonbrain.distribution <- c(
        nonbrain.distribution, 
        names(na.omit(z.score.storage[[disease]][z.score.storage[[disease]] > zscore.count]))
        )
    }

# create a table:
sort(table(brain.distribution))
sort(table(nonbrain.distribution))

# make a heatmap to investigate if there are any differences:
## first we will make one for the non brain traits:

# grab out the issue and region names:
cell.type <- unique(sapply(strsplit(names(zscore), ":"), '[[', 1));
tissue <- unique(sapply(strsplit(names(zscore), ":"), '[[', 2))

# create place holder for heatmap:
record.matrix <- matrix(NA, nrow = length(tissue), ncol = length(cell.type))
rownames(record.matrix) <- tissue
colnames(record.matrix) <- cell.type
brain.record.matrix <- record.matrix

### for non brain traits:
cell.type.tissue.pair.index <- data.frame(
    cell_type = sapply(strsplit(names(table(nonbrain.distribution)), ":"), '[[', 1),
    tissue = sapply(strsplit(names(table(nonbrain.distribution)), ":"), '[[', 2),
    count = table(nonbrain.distribution)
    )
stopifnot(!duplicated(cell.type.tissue.pair.index$count.nonbrain.distribution))

# find the set of tissue and cell type pair that has less than the z score:
zeros <- NULL
for (disease in non.brain.traits) {
    zeros <- c(
        zeros, 
        names(na.omit(z.score.storage[[disease]][z.score.storage[[disease]] <= zscore.count]))
        )
    }
zeros <- unique(zeros)

# seperate out the z score and index:
zero.index <- data.frame(
    cell_type = sapply(strsplit(names(table(zeros)), ":"), '[[', 1),
    tissue = sapply(strsplit(names(table(zeros)), ":"), '[[', 2),
    count = 0
    )

# fill in the zeros to place holder:
for (line in seq(1, nrow(zero.index))) {
    cell.type <- zero.index[line, 'cell_type']
    tissue <- zero.index[line, 'tissue']
    count = zero.index[line, 'count']
    record.matrix[tissue, cell.type] <- count
    }

# fill in the real counts
for (line in seq(1, nrow(cell.type.tissue.pair.index))) {
    cell.type <- cell.type.tissue.pair.index[line, 'cell_type']
    tissue <- cell.type.tissue.pair.index[line, 'tissue']
    count = cell.type.tissue.pair.index[line, 'count.Freq']
    record.matrix[tissue, cell.type] <- count
    }

## for brain traits:
# find the set of tissue and cell type pair that has less than the z score:
zeros <- NULL
for (disease in brain.traits) {
    zeros <- c(
        zeros, 
        names(na.omit(z.score.storage[[disease]][z.score.storage[[disease]] <= zscore.count]))
        )
    }
zeros <- unique(zeros)

# find the set of tissue and cell type pair that has greater than the z score:
cell.type.tissue.pair.index <- data.frame(
    cell_type = sapply(strsplit(names(table(brain.distribution)), ":"), '[[', 1),
    tissue = sapply(strsplit(names(table(brain.distribution)), ":"), '[[', 2),
    count = table(brain.distribution)
    )
stopifnot(!duplicated(cell.type.tissue.pair.index$count.nonbrain.distribution))

# fill in the cell type tissue pairs that has 0 z score > zscore.count
zero.index <- data.frame(
    cell_type = sapply(strsplit(names(table(zeros)), ":"), '[[', 1),
    tissue = sapply(strsplit(names(table(zeros)), ":"), '[[', 2),
    count = 0
    )

# fill in the place holder:
for (line in seq(1, nrow(zero.index))) {
    cell.type <- zero.index[line, 'cell_type']
    tissue <- zero.index[line, 'tissue']
    count = zero.index[line, 'count']
    brain.record.matrix[tissue, cell.type] <- count
    }

for (line in seq(1, nrow(cell.type.tissue.pair.index))) {
    cell.type <- cell.type.tissue.pair.index[line, 'cell_type']
    tissue <- cell.type.tissue.pair.index[line, 'tissue']
    count = cell.type.tissue.pair.index[line, 'count.Freq']
    brain.record.matrix[tissue, cell.type] <- count
    }

### VISUALIZE #####################################################################################
# lets create a heatmap:
cell.type <- c(
    "L2-3-IT",
    "L4-IT",
    "L5-ET",
    "L5-IT",
    "L5-6-NP",
    "L6-CT",
    "L6-IT",
    "L6-IT-Car3",
    "L6b",
    "Amy-Exc",
    "CA1",
    "CA3",
    "DG",
    "HIP-Misc1",
    "HIP-Misc2",
    "CB",
    "Chd7",
    "Foxp2",
    "MSN-D1",
    "MSN-D2",
    "PKJ",
    "PN",
    "Lamp5",
    "Lamp5-Lhx6",
    "Pvalb",
    "Pvalb-ChC",
    "Sncg",
    "Sst",
    "SubCtx-Cplx",
    "THM-Exc",
    "THM-Inh",
    "THM-MB",
    "Vip",
    "ASC",
    "EC",
    "MGC",
    "ODC",
    "OPC",
    "PC",
    "VLMC"
    );

column.split = c(
    rep('Excitatory', 15),
    rep('Inhibitory', 18),
    rep('Others', 7)
    )

# residualize the two heatmap:
residual.table <- brain.record.matrix / length(brain.traits) - record.matrix / length(non.brain.traits)
plot.list <- list(record.matrix, brain.record.matrix, residual.table)
names(plot.list) <- c('nonbrain', 'brain', 'residual')

# make a heatmap:
for (mode in c('nonbrain', 'brain', 'residual')){
    plot.df <- plot.list[[mode]]
    plot.df <- plot.df[, cell.type]

    # create color function
    col.fun <- colorRamp2(
        c(
            0,
            max(plot.df, na.rm = TRUE)
            ),
        c('white', '#ca0020')
        );
    heatmap.legend.param <- list(
        at = c(
            0,
            round(
                max(plot.df, na.rm = TRUE),
                digit = 2 # round to 2 digit
                )
            )
        );
    
    if (mode == 'residual') {
        color.bar.name <- 'residuals'

        # create a new color function if we are plotting the residualized heatmap:
        col.fun <- colorRamp2(
            c(min(plot.df, na.rm = TRUE), 0, max(plot.df, na.rm = TRUE)),
            c('#0571b0', 'white', '#ca0020')
            )
        
        # create a new heatmap legend:
        heatmap.legend.param <- list(
            at = c(
                round(min(plot.df, na.rm = TRUE), digit = 2),
                0,
                round(max(plot.df, na.rm = TRUE), digit = 2)
                )
            )
        } else {
            color.bar.name <- paste0('number of \n tissue : cell type \n pairs with \n zscore > ', zscore.count)
        }

    # Now create the heatmap:
    plot <- Heatmap(
        as.matrix(plot.df),
        name = color.bar.name,
        col = col.fun,
        rect_gp = gpar(col = "black", lwd = 2),
        # row_order = publication.traits,
        column_order = colnames(plot.df),
        width = unit(10 * ncol(plot.df),"mm"),
        height = unit(10 * nrow(plot.df),"mm"),
        column_names_gp = grid::gpar(fontsize = 15),
        row_names_gp = grid::gpar(fontsize = 15),
        # row_split = row.split,
        cluster_rows = FALSE,
        column_split = column.split,
        heatmap_legend_param = heatmap.legend.param
        );
    plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

    # measure the size of the heatmap:
    heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
    heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

    # use the measured width and height for drawing:
    plot.path <- paste0(output.dir, system.date, '-region-cell-type-zscore-count-', mode, '.png')
    png(
        filename = plot.path,
        width = heatmap.width,
        height = heatmap.height,
        units = 'in',
        res = 400
        );
    draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));
    dev.off();
    }
