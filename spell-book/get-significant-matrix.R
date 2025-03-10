get.significant.matrix <- function(risk.score.collection, meta.data, group.index, plot.type = 'proportion') {
    # create place holder for the group index:
    feature <- unique(meta.data[, group.index]);
    significance.matrix <- matrix(NA, nrow = length(risk.score.collection), ncol = length(feature));
    rownames(significance.matrix) <- names(risk.score.collection);
    colnames(significance.matrix) <- feature;

    # create the traits by cell type matrix:
    for (result in names(risk.score.collection)) {
        # find the index for significant cells:
        significant.cell <- rownames(risk.score.collection[[result]])[
            p.adjust(risk.score.collection[[result]]$pval, method = 'fdr') < p.cutoff
            ];

        # for each cell type, check number of significant cells are part of that cell type:
        for (feature.level in feature) {
            # locate cell id that belong in the cell type:
            feature.level.cell <- rownames(meta.data)[meta.data[, group.index] %in% feature.level];

            # find the number of cells that are in each of the cell type category:
            if (plot.type == 'count') {
                significance.matrix[result, feature.level] <- sum(significant.cell %in% feature.level.cell);
                } else {
                    significance.matrix[result, feature.level] <- sum(significant.cell %in% feature.level.cell) /
                        length(feature.level.cell);
                }
            }
        }
    return(significance.matrix);
    }

### UNIT TESTING ##################################################################################
# test.score <- risk.score[1:2]
# test.score[[1]]$pval = c(rep(1,15000),rep(0, 15000))
# test.score[[2]]$pval = c(rep(1,15000),rep(0, 15000))
# test.meta.data <- meta[rownames(test.score[[1]]),]
# test.meta.data$CellClass <- c(rep('Exc', 10000), rep('Inh', 10000), rep('NonN', 10000))

# sig.matrix.proportion <- lapply(
#     X = c('CellClass', 'CellClass'),
#     get.significant.matrix,
#     risk.score.collection = test.score,
#     meta.data = test.meta.data,
#     plot.type = 'proportion'
#     )

# since the test case have 30000 cells, the expected output :
#                           Exc Inh NonN
# PASS_ADHD_Demontis2018       0 0.5    1
# PASS_Alzheimers_Jansen2019   0 0.5    1
# is returned

