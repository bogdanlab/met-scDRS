### brain-region-linear-model.R ###################################################################
# purpose: Build linear model disease z score with brain region as explainatory variable in each
# of the cell type

### PREAMBLE ######################################################################################
# gather the different regions
require(docopt)
require(circlize);
require(tidyverse);

'Usage:
    brain-region-linear-model [--score_dir <scdrs> --meta_data <meta> --field1 <group1> --field2 <group2> --p_cutoff <cutoff> --num_cutoff <num_cutoff> --min_cell_num <min_num> --out <output>]

Options:
    --score_dir path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field1 categorical variable in meta data
    --field2 second categorical variable in meta, will test for significance of each level of field 2 within each level of field 1
    --p_cutoff FDR cutoff for cell to be included in this analysis for each disease [default: 1] 
    --num_cutoff minimum number of significant cells required for each disease to have for the disease to be included in this analysis [default: 0]
    --min_cell_num minimum number of cells to build linear model, minimum should be 2 [default: 2]
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.directory <- opts$score_dir;
group1.index <- opts$field1;
group2.index <- opts$field2;
p.cutoff <- as.numeric(opts$p_cutoff);
cell.cutoff <- as.numeric(opts$num_cutoff);
output.path <- opts$out;
system.date <- Sys.Date();
min.cell.num <- as.numeric(opts$min_cell_num);

# for code testing:
# meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv"
# scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/"
# group1.index <- "X_MajorType"
# group2.index <- "tissue"
# p.cutoff <- 1;
# cell.cutoff <- 0;
# output.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'
# min.cell.num <- 100
# top.visualize <- 10;

# read data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
all.risk.score <- vector('list', length = length(score.files));
names(all.risk.score) <- score.files;

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# make sure that the group index is in the region:
stopifnot(group1.index %in% colnames(meta));

# read into the empty list:
relavent.disease <- rep(FALSE, length = length(score.files));
names(relavent.disease) <- score.files

# create progress bar:
cat('loading in trait scores \n');
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(score.files),
    clear = FALSE
    );

# load in disease score:
for (result in score.files) {
    all.risk.score[[result]] <- read.table(
        file = result,
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );

    # make sure that the risk score have the same set of cells as meta in the same order:
    stopifnot(rownames(all.risk.score[[result]]) %in% rownames(meta));

    # filter to significant only cells:
    fdr <- p.adjust(all.risk.score[[result]]$pval, method = 'fdr');
    significant.cell <- fdr <= p.cutoff;
    
    # filter disease if it has < X significant cells:
    if (sum(significant.cell) > cell.cutoff) {
        relavent.disease[result] <- TRUE;
        all.risk.score[[result]] <- all.risk.score[[result]][significant.cell, ];
        }
    
    pb$tick()
    }

# filter the diseases:
risk.score <- all.risk.score[relavent.disease];

# simplify list names:
list.names <- gsub(scDRS.directory, '', names(risk.score));
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

### MODEL BUILDING ################################################################################
# Next start building the model:
# constructing linear model between each level of field 2 meta against risk score, within each level
# of field 1 meta.
group1.levels <- unique(meta[, group1.index]);
group2.levels <- unique(meta[, group2.index]);

# create progress bar:
total.operations <- length(risk.score) * length(group1.levels) * length(group2.levels);
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = total.operations,
    clear = FALSE
    );

# create place holder:
effect.collection <- pval.collection <- vector('list', length = length(risk.score));
names(effect.collection) <- names(pval.collection) <- names(risk.score);
for (disease in names(risk.score)){
    # create place holder:
    pval.collection[[disease]] <- matrix(NA, ncol = length(group2.levels), nrow = length(group1.levels));
    colnames(pval.collection[[disease]]) <- group2.levels;
    rownames(pval.collection[[disease]]) <- group1.levels;
    effect.collection[[disease]] <- pval.collection[[disease]];

    for (group1.level in group1.levels){
        # obtain the group1 level cells:
        group1.level.cell <- rownames(meta)[meta[, group1.index] == group1.level]

        for(group2.level in group2.levels){
            # obtain the group2 level cells:
            group2.level.cell <- rownames(meta)[meta[, group2.index] == group2.level]
            group12.level.cell <- intersect(group1.level.cell, group2.level.cell);
            group12.level.cell <- intersect(group12.level.cell, rownames(risk.score[[disease]]))

            # create a design matrix:
            design <- rep(0, length = length(group1.level.cell));
            names(design) <- group1.level.cell;
            design[group12.level.cell] <- 1;

            # grab the response vector:
            design.df <- data.frame(
                zscore = risk.score[[disease]][group1.level.cell, 'zscore'],
                covariate = design
                );
            
            # compute model:
            if (length(group12.level.cell) > min.cell.num) {
                result <- lm(zscore ~ covariate, data = design.df)
                pval <- summary(result)$coefficients['covariate', 'Pr(>|t|)'];
                effect <- result$coefficients['covariate'];
            } else {
                pval <- NA;
                effect <- NA;
            }

            # save the p value:
            pval.collection[[disease]][group1.level, group2.level] <- pval;
            effect.collection[[disease]][group1.level, group2.level] <- effect;

            # report progress:
            pb$tick()
        }
    }
}

### VISUALIZATION #################################################################################
# we wish to first flatten our matrixes of p values:
convert_long <- function(collection, val_to = 'pval') {
    # create a placeholder:
    collection_long <- vector('list', length = length(collection))
    names(collection_long) <- names(collection)

    # for each of the trait, convert it into a long format:
    for (trait in names(collection)){
        # convert each of the disease p value from a heatmap to a long format:
        collection[[trait]] <- data.frame(collection[[trait]])
        cell_type <- rownames(collection[[trait]])
        collection[[trait]] <- cbind(cell_type, collection[[trait]])
        collection_long[[trait]] <- data.frame(
            pivot_longer(
                collection[[trait]],
                cols = 2:ncol(collection[[trait]]),
                names_to = 'tissue',
                values_to = val_to
                )
            );
        collection_long[[trait]]$disease = trait;
    }

    # concatenate across traits:
    concat_long <- Reduce('rbind', collection_long)
    return(concat_long)
}

# convert the data into a long format:
long.p.val <- convert_long(pval.collection, val_to = 'pval')
long.effect <- convert_long(effect.collection, val_to = 'effect')

# adjust the p value:
correction.factor <- Reduce('sum', lapply(pval.collection, FUN = function(x) sum(!is.na(x))));
long.p.val$bf_p <- p.adjust(long.p.val$pval, method = 'bonferroni', n = correction.factor);
long.p.val$effect <- long.effect$effect

# write out results:
write.table(
    x = long.p.val,
    file = paste0(output.path, system.date, '-effect-size-',group1.index, '-', group2.index, '-disease-table.csv'),
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
    );
