### group-specific-MC-LRT-test.R ##################################################################
# purpose: an analysis on how to do group specific MC test with respect to an annotation 

### PREAMBLE ######################################################################################
# load in specified requirements:
require(docopt);
require(dplyr);
require(data.table);
require(progress);
require(ComplexHeatmap);
require(circlize);

'Usage:
    annotation-disease-score-model [--score_dir <scdrs> --meta_data <meta> --field1 <group1> --field2 <group2> --out <output>]

Options:
    --score_dir path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field1 categorical variable in in meta data that you would like to subset into each level
    --field2 categorical variable in meta which you would like to aggregate the scDRS Z score in within the level of field 1
    --disease which of the trait you would want to look at
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.directory <- opts$score_dir;
group1.index <- opts$field1;
group2.index <- opts$field2;
output.path <- opts$out;
system.date <- Sys.Date();

# for function testing:
meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv"
scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/"
group1.index <- "X_MajorType"
group2.index <- "X_Region"
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-GSE215353-within-cell-type-MC-test')

# read data:
score.files <- list.files(scDRS.directory, pattern = '\\.full_score.gz', full.names = TRUE);
all.risk.score <- vector('list', length = length(score.files));
names(all.risk.score) <- score.files;

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# make sure the group index is in the meta data:
stopifnot(group1.index %in% colnames(meta));
stopifnot(group2.index %in% colnames(meta));

# create progress bar:
cat('start reading in the case and control scores:\n')
# create progress bar:
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(score.files),
    clear = FALSE
    );

for (result in score.files) {
    data <- fread(
        file = result,
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE,
        data.table = FALSE
        );
    all.risk.score[[result]] <- data.frame(data, row.names = 1);
    pb$tick()
    }

# simplify list names:
list.names <- gsub(scDRS.directory, '', names(all.risk.score));
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);
names(all.risk.score) <- list.names

### MODEL #########################################################################################
# create a progress bar:
cat('Computing test statistics for each disease \n')
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(all.risk.score),
    clear = FALSE
    );

# for each of the level in the group 1 variable, 
group1.level <- unique(meta[, group1.index]);
group2.level <- unique(meta[, group2.index]);

# create place holder for our process:
f.collection <- vector('list', length = length(list.names));
names(f.collection) <- names(all.risk.score)

# start computation for each of the disease
for (disease in names(all.risk.score)){
    # subset within each of the disease:
    risk.score <- all.risk.score[[disease]];

    # loop over each level of the first variable:
    for (level in group1.level) {
        # subset to the cells within the categorical variable 1:
        cells.in.level <- rownames(meta)[meta[, group1.index] == level];
        cells.in.level <- intersect(cells.in.level, rownames(all.risk.score[[disease]]))

        # perform check:
        # 1) more cells than variables
        # 2) more more than one level in group2
        if (length(cells.in.level) > length(group2.level) & length(unique(meta[cells.in.level, group2.index] > 1))) {
            # create formula:
            formula <- as.formula(paste("~", group2.index))
            covariate <- model.matrix(object = formula, data = meta[cells.in.level, ])
            covariate <- covariate[, -1];

            # create models:
            full.model <- lm(risk.score[cells.in.level, 'norm_score'] ~ covariate);
            small.model <- lm(risk.score[cells.in.level, 'norm_score'] ~ 1)

            # perform likelihood ratio test:
            f.stat <- as.numeric(anova(small.model, full.model)[2, 'F']);
            f.collection[[disease]][level] <- f.stat;
        }
    }
    pb$tick()
}

### NULL DISTRIBUTION #############################################################################
# find out how much models we need to compute:
control.scores <- colnames(risk.score)[grep('ctrl_norm_score_', colnames(risk.score))]
job.num <- sum(sapply(f.collection, length)) * length(control.scores)

# Compute null distribution:
cat('Computing the background distriution: \n')
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = job.num,
    clear = FALSE
    );

# create place holder for our process:
null.collection <- vector('list', length = length(list.names));
names(null.collection) <- names(all.risk.score)

# start computation for each of the disease
for (disease in names(all.risk.score)){
    # subset within each of the disease:
    risk.score <- all.risk.score[[disease]];
    null.collection[[disease]] <- vector('list', length = length(group1.level));
    names(null.collection[[disease]]) <- group1.level;

    # loop over each level of the first variable:
    for (level in group1.level) {
        null.collection[[result]][[level]] <- rep(0, length = length(control.scores));
        names(null.collection[[result]][[level]]) <- control.scores;

        # subset to the cells within the categorical variable 1:
        cells.in.level <- rownames(meta)[meta[, group1.index] == level];
        cells.in.level <- intersect(cells.in.level, rownames(all.risk.score[[disease]]));

        # perform check:
        # 1) more cells than variables
        # 2) more more than one level in group2
        if (length(cells.in.level) > length(group2.level) & length(unique(meta[cells.in.level, group2.index] > 1))) {
            # create formula:
            formula <- as.formula(paste("~", group2.index))
            covariate <- model.matrix(object = formula, data = meta[cells.in.level, ])
            covariate <- covariate[, -1];

            # create models:
            for (control.norm.score in control.scores){
                full.model <- lm(risk.score[cells.in.level, control.norm.score] ~ covariate);
                small.model <- lm(risk.score[cells.in.level, control.norm.score] ~ 1)

                # perform likelihood ratio test:
                f.stat <- as.numeric(anova(small.model, full.model)[2, 'F']);
                null.collection[[disease]][[level]][control.norm.score] <- f.stat;
                pb$tick()
                }
            }
        }
    }

### OUTPUT ########################################################################################
# get the p values for each of the disease for each of the cell type:
mc.matrix <- matrix(NA, ncol = length(f.collection), nrow = length(f.collection[[1]]))
colnames(mc.matrix) <- names(f.collection)
rownames(mc.matrix) <- names(f.collection[[1]])

for (disease in names(f.collection)) {
    for (level in names(f.collection[[disease]])) {
        # get the actual test statistics:
        mc.p <- (1 + sum(f.collection[[disease]][level] >= null.collection[[disease]][[level]])) / 
            (length(null.collection[[disease]][[level]]) + 1)

        # write data into place holder:
        mc.matrix[level, disease] <- mc.p
    }
}
# save this mc.matrix:
csv.output.path <- paste0(output.path, '.csv');
write.table(
    mc.matrix,
    file = csv.output.path,
    sep = ',',
    row.names = TRUE,
    col.names = TRUE
    )

# visualize this matrix:
significance.matrix <- t(-log10(mc.matrix));
# make visualization:
# design olor function:
cutoff <- -log10(0.05);

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
    name = '-log10(P)',
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

heatmap.output.path <- paste0(output.path, '.png')
png(
    filename = heatmap.output.path,
    width = heatmap.width,
    height = heatmap.height,
    units = 'in',
    res = 400
    );
draw(plot, padding = unit(c(10, 10, 10, 70), "mm"), heatmap_legend_side = "bottom");
dev.off();
