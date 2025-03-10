### monte-carlo-linear-model.R ##########################################################
# purpose: use the monte carlo draws to compute test statistics of group1 specific
# group 2 effect toward disease score prediction

### PREAMBLE ############################################################################
# gather the different regions
require(docopt)
require(circlize);
require(tidyverse);
require(data.table);

'Usage:
    monte-carlo-linear-model.R [--score_file <scdrs> --meta_data <meta> --field1 <group1> --field2 <group2> --p_cutoff <cutoff> --min_sig_cell <min_sig_cell> --min_cell_num <min_num> --out <output>]

Options:
    --score_file path to full scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field1 categorical variable in meta data
    --field2 second categorical variable in meta, will test for significance of each level of field 2 within each level of field 1
    --p_cutoff FDR cutoff for cell to be included in this analysis for each disease [default: 1] 
    --min_sig_cell minimum number of significant cells required for each disease to have for the disease to be included in this analysis [default: 0]
    --min_cell_num minimum number of cells to build linear model, minimum is at least 2 [default: 2]
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.file <- opts$score_file;
group1.index <- opts$field1;
group2.index <- opts$field2;
p.cutoff <- as.numeric(opts$p_cutoff);
min.significant.cell <- as.numeric(opts$min_sig_cell);
output.path <- opts$out;
system.date <- Sys.Date();
min.cell.num <- as.numeric(opts$min_cell_num);

# for code testing:
# meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv"
# scDRS.file <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/PASS_MDD_Howard2019.full_score.gz"
# group1.index <- "X_MajorType"
# group2.index <- "tissue"
# p.cutoff <- 0.1;
# min.significant.cell <- 100;
# output.path <- '/u/scratch/l/lixinzhe/tmp-file/test/'
# min.cell.num <- 100

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# make sure that the group index is in the region:
stopifnot(group1.index %in% colnames(meta));

# load in the data
full.score <- fread(
    file = scDRS.file,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE,
    data.table = FALSE
    );
rownames(full.score) <- full.score$cell

# use the p.cutoff flag to filter for the set of cells that we care about:
if (p.cutoff < 1) {
    # filter for the cells that are significant:
    full.score$fdr <- p.adjust(full.score$pval, method = 'fdr')
    full.score.origin <- full.score
    significant.cells <- full.score$fdr <= p.cutoff
    full.score <- full.score.origin[significant.cells, ]
    num.sig.cell.data <- nrow(full.score)
    } else {
        num.sig.cell.data <- nrow(full.score)
    }

### MODEL BUILDING ################################################################################
# Next start building the model:
# constructing linear model between each level of field 2 meta against risk score, within each level
# of field 1 meta.
if (num.sig.cell.data >= min.significant.cell) {
    group1.levels <- unique(meta[, group1.index]);
    group2.levels <- unique(meta[, group2.index]);
    control.score.index <- colnames(full.score)[grep('ctrl_norm_score', colnames(full.score))];
    cat('loaded', length(control.score.index), 'control scores \n');
    group12.pairs <- expand.grid(group1.levels, group2.levels, stringsAsFactors = FALSE);
    group12.pairs <- paste0(group12.pairs$Var1, ':', group12.pairs$Var2);

    # create progress bar:
    total.operations <- (1 + length(control.score.index)) * length(group1.levels) * length(group2.levels);
    pb <- progress::progress_bar$new(
        format = "[:bar] (:current/:total)",
        total = total.operations,
        clear = FALSE
        );

    # create place holder:
    result.collection <- vector('list', length = length(group12.pairs));
    names(result.collection) <- group12.pairs;
    for(group12.pair in group12.pairs){
        result.collection[[group12.pair]] <- data.frame(matrix(NA, ncol = 4, nrow = 1 + length(control.score.index)))
        rownames(result.collection[[group12.pair]]) <- c('foreground', control.score.index)
        colnames(result.collection[[group12.pair]]) <- c('group1', 'group2', 'beta', 'p')

    }

    for (group1.level in group1.levels){
        # get the group1 level index:
        group1.level.cell <- rownames(meta)[meta[, group1.index] == group1.level]
        group1.level.cell <- intersect(group1.level.cell, rownames(full.score))
        for (group2.level in group2.levels){
            # get group2 level index:
            pair.index <- paste0(group1.level, ':', group2.level);
            stopifnot(pair.index %in% names(result.collection));
            group2.level.cell <- rownames(meta)[meta[, group2.index] == group2.level]
            group2.level.cell <- intersect(group2.level.cell, rownames(full.score))
            group12.level.cell <- intersect(group1.level.cell, group2.level.cell);

            # assign group labels into the result matrix:
            result.collection[[pair.index]]$group1 <- group1.level;
            result.collection[[pair.index]]$group2 <- group2.level;

            # if there are more than X number of group12 level cells, we will construct model:
            if (length(group12.level.cell) > min.cell.num) {
                # create a design matrix:
                design <- rep(0, length = length(group1.level.cell));
                names(design) <- group1.level.cell;
                design[group12.level.cell] <- 1;
                
                # for each MC instance, grab out the design matrix：
                # grab the foreground response vector:
                design.df <- data.frame(
                    zscore = full.score[group1.level.cell, 'norm_score'],
                    covariate = design
                    );
                # compute model:
                result <- lm(zscore ~ covariate, data = design.df)
                pval <- summary(result)$coefficients['covariate', 'Pr(>|t|)'];
                effect <- result$coefficients['covariate'];

                # record the computed p value and effect:
                result.collection[[pair.index]]['foreground', 'beta'] <- effect;
                result.collection[[pair.index]]['foreground', 'p'] <- pval;
                pb$tick()
                
                # also compute the MC instances:
                for (control_score in control.score.index){
                    # write out the design data frame index:
                    design.df <- data.frame(
                        zscore = full.score[group1.level.cell, control_score],
                        covariate = design
                        );
                    
                    # compute model:
                    result <- lm(zscore ~ covariate, data = design.df)
                    pval <- summary(result)$coefficients['covariate', 'Pr(>|t|)'];
                    effect <- result$coefficients['covariate'];

                    # record the computed p value and effect:
                    result.collection[[pair.index]][control_score, 'beta'] <- effect;
                    result.collection[[pair.index]][control_score, 'p'] <- pval;
                    pb$tick()
                    }
                } else {
                    pval <- NA;
                    effect <- NA;
                    # record the computed p value and effect:
                    result.collection[[pair.index]]['foreground', 'beta'] <- effect;
                    result.collection[[pair.index]]['foreground', 'p'] <- pval;
                    pb$tick()

                    # also put the NAs into the control scores:
                    for (control_score in control.score.index){
                        # record the computed p value and effect:
                        result.collection[[pair.index]][control_score, 'beta'] <- effect;
                        result.collection[[pair.index]][control_score, 'p'] <- pval;
                        pb$tick()
                        }
                }
            }
        }

    for (pair.index in names(result.collection)) {
        index <- rownames(result.collection[[pair.index]]);
        result.collection[[pair.index]] <- cbind(index, result.collection[[pair.index]]);
        pair <- gsub(' ','-', gsub('/', '-', pair.index));
        fwrite(
            x = result.collection[[pair.index]],
            sep = ',',
            file = paste0(output.path, pair, '-MC-linear-model-summary.csv'),
            row.names = FALSE,
            col.names = TRUE
            )
        }

}
