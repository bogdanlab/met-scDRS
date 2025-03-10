### cov-generate.R ################################################################################
# purpose: output a covaraite file for running met-scDRS with:

### PREAMBLE ######################################################################################
# load in file:
### covariate visualization.R #####################################################################
# purpose: validate our signal by making sure it is not driven by some covaraites:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(cowplot)

# specify paths:
scDRS.directory <- '/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/';
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data_with_rowSum.csv';
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date()

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

subset_meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv'
    )

# read trait info:
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

# load in data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# create progress bar:
cat('loading scores \n')
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(score.files),
    clear = FALSE
    );

# read into the empty list:
for (result in score.files) {
    risk.score[[result]] <- fread(
        file = result,
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE
        );
    risk.score[[result]] <- data.frame(risk.score[[result]], row.names = 1)
    pb$tick()
    }

# simplify list names:
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

# make a cov file:
covariates <- c('donor_id', 'batch', 'rowSum', 'X_pool')
cov <- meta[rownames(risk.score[[1]]), covariates]
cov$batch <- paste0('batch', cov$batch)
cov$X_pool <- paste0('pool', cov$X_pool)
cov$rowSum <- scale(cov$rowSum)

# prepare design matrix for regression:
design.matrix <- model.matrix(formula(~ 0 + X_pool + batch + donor_id + rowSum), data = cov);

# test for colinearity:
linear.combination.check <- alias(meta$X_mCGFrac ~ design.matrix);
alias.features <- gsub('design.matrix','',rownames(linear.combination.check$Complete));
full.rank.design <- design.matrix[, setdiff(colnames(design.matrix), alias.features)];

# return user feed back:
warning(paste0('alias: ', alias.features, ' detected and removed from design matrix \n'));
cat('alias checked, adding in intercept term!', '\n');

# add back the constant term:
const <- 1;
full.rank.design <- cbind(const, full.rank.design); 
cell <- rownames(full.rank.design)

# check if all columns are numeric:
stopifnot(is.numeric(full.rank.design));
stopifnot(!is.na(full.rank.design));

# return the covariate matrix:
output.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-covarite.cov"
full.rank.design <- cbind(cell, full.rank.design)
write.table(
    full.rank.design,
    file = output.path,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
    );

subset_path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-covarite-subset.cov"
write.table(
    full.rank.design[rownames(subset_meta), ],
    file = subset_path,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
    );