### get-covariates.R ##############################################################################
# PURPOSE: Prepare meta data into covariates matrix for scDRS:
# NOTE: RData and rds file are assumed to be Seurat objects!

### PREAMBLE ######################################################################################
require(Seurat);
require(sceasy);

# gather user inputs:
args <- commandArgs(trailingOnly = TRUE);

# try loading obj:
seurat.path <- args[1];
file.type <- gsub('.*\\.', '', seurat.path);
seurat.obj <- switch(
    file.type,
    'rds' = readRDS(seurat.path),
    'RData' = print('please supply rds file!')
    );

# get meta data:
meta.data <- seurat.obj@meta.data;

# get the list of covariates user wanted to correct for:
covariate <- args[2];

# get the output path:
output.path <- args[3];

### CONSTRUCT DESIGN MATRIX #######################################################################
# prepare design matrix for regression:
design.matrix <- model.matrix(formula(covariate), data = meta.data);

# test for colinearity:
linear.combination.check <- alias(seurat.obj@assays$RNA@data[1,] ~ design.matrix);
alias.features <- gsub('design.matrix','',rownames(linear.combination.check$Complete));
full.rank.design <- design.matrix[, setdiff(colnames(design.matrix), alias.features)];

# return user feed back:
warning(paste0('alias: ', alias.features, ' detected and removed from design matrix \n'));
cat('alias checked, adding in intercept term!', '\n');

# add back the constant term:
const <- 1;
full.rank.design <- cbind(const, full.rank.design); 

# check if all columns are numeric:
stopifnot(is.numeric(full.rank.design));

# return the covariate matrix:
write.table(
    full.rank.design,
    file = output.path,
    sep = '\t',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );
