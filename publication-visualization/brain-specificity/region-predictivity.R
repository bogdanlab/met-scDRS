### region-predictivity.R #########################################################################
# purpose: plot out the MC z score of predictivity across traits 

### PREAMBLE ######################################################################################
# load libraries:
require(docopt)
require(circlize);
require(ComplexHeatmap);
require(ggplot2);
require(tidyverse);
require(ggrepel);
require(data.table);

# define parameters:
output.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date();

# Extract the arguments
args <- commandArgs(trailingOnly = TRUE)
model.result.dir <- args[1]

# load in the data:
start <- Sys.time()
model.result.dir <- "/u/scratch/l/lixinzhe/tmp-file/mc_test/"
traits <- list.dirs(model.result.dir, full.names = FALSE, recursive = FALSE)

# create a progress bar:
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(traits),
    clear = FALSE
    );

trait.collection <- NULL;
for (trait in traits) {
    trait.result.dir <- paste0(model.result.dir, trait, '/')
    mc.results <- list.files(trait.result.dir, pattern = '*-model-summary.csv')
    for (result.file in mc.results) {
        pair.index <- gsub('-MC-linear-model-summary.csv', '', result.file)
        trait.collection[[trait]][[pair.index]] <- data.frame(
            fread(
                file = paste0(trait.result.dir, result.file),
                sep = ',',
                data.table = FALSE
                ),
            row.names = 1
            )
        }
    pb$tick()
}

# load in the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

### COMPUTE PERMUTED VALUE ########################################################################
# compute the permuted value:
cat('compute permuted value: \n')

# create a placeholder:
zscore.result <- vector('list', length = length(traits))
names(zscore.result) <- traits
permuted.result <- zscore.result
beta.result <- zscore.result

for (trait in traits){
    # create place holder:
    result.collection <- trait.collection[[trait]]
    permuted.value <- rep(NA, length = length(result.collection))
    names(permuted.value) <- names(result.collection)
    beta <- zscore <- permuted.value

    # compute the p values:
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
    # record the z score and the p values:
    zscore.result[[trait]] <- zscore
    permuted.result[[trait]] <- permuted.value
    beta.result[[trait]] <- beta
}

# shape up our results for output in matrix format:
zscore <- Reduce('rbind', zscore.result)
pval <- Reduce('rbind', permuted.result)
beta <- Reduce('rbind', beta.result)
rownames(zscore) <- rownames(pval) <- rownames(beta) <- traits

### VISUALIZATION #################################################################################
# remove the columns that contain only NAs:
na.mask <- apply(zscore, 2, FUN = function(x) sum(is.na(x)))
na.mask <- na.mask == length(traits)
zscore <- zscore[, !na.mask]

# now plot out the distribution:
brain.index <- trait.info$Trait_Identifier[trait.info$Category == 'brain']
brain.trait.zscore <- Reduce('c', zscore[brain.index, ])
nonbrain.zscore <- Reduce('c', zscore[setdiff(rownames(zscore), brain.index), ])

# create plot df:
brain.df <- data.frame(score = brain.trait.zscore, class = 'brain')
nonbrain.df <- data.frame(score = nonbrain.zscore, class = 'non-brain')
plot.df <- rbind(brain.df, nonbrain.df)
brain.index <- plot.df$class == 'brain'
nonbrain.index <- plot.df$class == 'non-brain'

# create the visualization:
gplot <- ggplot(plot.df, aes(score)) + 
    geom_histogram(data = plot.df[nonbrain.index, ], fill = "#bdbdbd", alpha = 0.5, bins = 50) +
    geom_histogram(data = plot.df[brain.index, ], fill = "#de2d26", alpha = 0.5, bins = 50) + 
    theme_classic() +
    xlab('region trait prediction zscore') +
    ylab('frequency') +
    theme(text = element_text(size = 20))

plot.name <- paste0(output.path, system.date, '-brain-vs-nonbrain-spatial-zscore.png')
png(
    filename = plot.name,
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot);
dev.off();
