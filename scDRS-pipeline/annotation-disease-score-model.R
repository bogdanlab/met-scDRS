### annotation-disease-score-model.R ##############################################################
# purpose: predict if annotation 1 or the interaction effect of annotation 2 with 1 is significant
# output: 
# case 1 (if you only have one annotation)
# 1) the p value for annotation 1 in predicting met-scDRS risk score

# case 2 (if you have two annotation)
# 1) the p value for annotation 2 on predicting significant cells' met-scDRS risk score 
# within each level of annotation 1
# 2) the expected risk score for each level of annotation 2 for each level in annotation 1

### PREAMBLE ######################################################################################
# gather the different regions
require(docopt)
require(circlize);
require(ComplexHeatmap)

'Usage:
    annotation-disease-score-model [--score_dir <scdrs> --meta_data <meta> --field1 <group1> --field2 <group2> --p_cutoff <cutoff> --num_cutoff <num_cutoff> --out <output>]

Options:
    --score_dir path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field1 categorical variable in meta data
    --field2 second categorical variable in meta, will build an interaction with field 1 if exists. LRT will test for interaction significance [default: empty]
    --p_cutoff FDR cutoff for cell to be included in this analysis for each disease [default: 0.1] 
    --num_cutoff minimum number of significant cells required for each disease to have for the disease to be included in this analysis [default: 0]
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

# for code testing:
# meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv"
# scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/"
# group1.index <- "X_MajorType"
# group2.index <- "tissue"
# p.cutoff <- 0.1;
# cell.cutoff <- 0;
# output.path <- '/u/home/l/lixinzhe/project-geschwind/plot/'

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
    }

# filter the diseases:
risk.score <- all.risk.score[relavent.disease];

# simplify list names:
list.names <- gsub(scDRS.directory, '', names(risk.score));
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

### BUILD MODEL ###################################################################################
# if second term does not exist, we will simply test for covariates being significant or not:
# using likelihood ratio test between full model against model with only intercept term
if (group2.index == 'empty') {
    # create list for housing the full model:
    model.collection <- vector('list', length = length(risk.score));
    names(model.collection) <- names(risk.score);
    pval.collection <- rep(0, length = length(risk.score));
    names(pval.collection) <- names(risk.score);

    # create progress bar:
    pb <- progress::progress_bar$new(
        format = "[:bar] (:current/:total)",
        total = length(risk.score),
        clear = FALSE
        );

    # for each of the disease, we run the model:
    cat('Starting likelihood ratio test: \n')
    for (result in names(risk.score)){
        # Create formula
        formula <- as.formula(paste("~", group1.index))
        covariate <- model.matrix(object = formula, data = meta[rownames(risk.score[[result]]), ])
        covariate <- covariate[, -1];
        
        # call linear model:
        full.model <- lm(risk.score[[result]][, 'zscore'] ~ covariate);
        small.model <- lm(risk.score[[result]][, 'zscore'] ~ 1)

        # perform likelihood ratio test:
        pval <- as.numeric(anova(small.model, full.model)[2, 'Pr(>F)']);
        pval.collection[result] <- pval;
        model.collection[[result]] <- full.model;

        # update progress bar:
        pb$tick()
        }
    }

# if field 2 is not empty, we will be subsetting to each level in first term
# and then compute association with second term:
if (group2.index != 'empty') {
    # create progress bar:
    pb <- progress::progress_bar$new(
        format = "[:bar] (:current/:total)",
        total = length(risk.score) * length(unique(meta[, group1.index])),
        clear = FALSE
        );
  
    # create list for housing the full model:
    model.collection <- vector('list', length = length(risk.score));
    names(model.collection) <- names(risk.score);
    pval.collection <- model.collection;

    # initiate the list for each teirs of disease:
    for (result in names(risk.score)){
        model.collection[[result]] <- vector('list', length = length(unique(meta[, group1.index])));
        names(model.collection[[result]]) <- unique(meta[, group1.index]);
        pval.collection[[result]] <- rep(NA, length = length(unique(meta[, group1.index])));
        names(pval.collection[[result]]) <- unique(meta[, group1.index]);
        }

    # subset to each groups and start computation:
    for (result in names(risk.score)){
        for(group1.level in unique(meta[, group1.index])){
            # group 1 level samples:
            cell.in.group.1.level <- rownames(meta)[meta[, group1.index] == group1.level]
            significant.cell <- rownames(risk.score[[result]]);
            cell.in.scope <- intersect(cell.in.group.1.level, significant.cell);

            # perform check:
            if (length(cell.in.scope) > 1 + length(unique(meta[, group2.index]))) {
                if (length(unique(meta[cell.in.scope, group2.index])) > 1) {
                    # create formula
                    formula <- as.formula(paste("~", group2.index))
                    covariate <- model.matrix(object = formula, data = meta[cell.in.scope, ])
                    covariate <- covariate[, -1];
                    
                    # create models:
                    full.model <- lm(risk.score[[result]][cell.in.scope, 'zscore'] ~ covariate);
                    small.model <- lm(risk.score[[result]][cell.in.scope, 'zscore'] ~ 1)

                    # perform likelihood ratio test:
                    pval <- as.numeric(anova(small.model, full.model)[2, 'Pr(>F)']);
                    pval.collection[[result]][group1.level] <- pval;
                    model.collection[[result]][[group1.level]] <- full.model;
                    }
                }
            pb$tick()
            }
        }
    }

### OUTPUT ########################################################################################
# first aggregate the output into a matrix of p values:
annotation2.significance <- data.frame(Reduce('rbind', pval.collection));
rownames(annotation2.significance) <- names(pval.collection);

# output table:
write.table(
    annotation2.significance,
    file = paste0(output.path, '-', group2.index, '-pvalue-in-predicting-zscore-for-each-level-in-', group1.index, '.csv'),
    sep = ',',
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
    )

### VISUALIZE #####################################################################################
# p value correction:
significance.matrix <- -log10(annotation2.significance + 10e-324);

# make visualization:
# design olor function:
cutoff <- -log10(0.05 / sum(!is.na(significance.matrix)));

col.fun <- colorRamp2(
    c(
        0,
        max(significance.matrix, na.rm = TRUE)
        ),
    c('white', '#de2d26')
    );
heatmap.legend.param <- list(
    title_position = 'topcenter',
    direction = 'horizontal',
    at = c(
        0,
        signif(max(significance.matrix, na.rm = TRUE), 3)
        )
    );

# plot out the heatmap:
plot <- Heatmap(
    as.matrix(significance.matrix),
    name = '-log10(P + 10e-324)',
    col = col.fun,
    na_col = 'gray',
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "black", lwd = 2),
    cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(as.matrix(significance.matrix)[i, j]) && as.matrix(significance.matrix)[i, j] > cutoff) {
                grid.text("*", x, y, gp = gpar(col = "black", fontsize = 16))
            }
        },
    width = unit(10 * ncol(significance.matrix),"mm"),
    height = unit(10 * nrow(significance.matrix),"mm"),
    column_names_gp = grid::gpar(fontsize = 15),
    row_names_gp = grid::gpar(fontsize = 15),
    heatmap_legend_param = heatmap.legend.param
    );

# automatically detect plot size and make plot:
plot.size <- draw(plot, heatmap_legend_side = 'left', padding = unit(c(10, 10, 10, 70), "mm"));

# measure the size of the heatmap:
heatmap.width <- convertX(ComplexHeatmap:::width(plot.size), "inch", valueOnly = TRUE);
heatmap.height <- convertY(ComplexHeatmap:::height(plot.size), "inch", valueOnly = TRUE)

png(
    filename = paste0(output.path, '-', group2.index, '-pvalue-in-predicting-zscore-for-each-level-in-', group1.index, '.png'),
    width = heatmap.width ,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, padding = unit(c(10, 10, 10, 70), "mm"), heatmap_legend_side = "bottom");
dev.off();
