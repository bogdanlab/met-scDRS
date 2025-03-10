### gene-checking.R ###############################################################################
# purpose: check the genes that are present between the production data and the subsetted data

### PREAMBLE ######################################################################################
# load in library:
library('data.table')
library('ggplot2')
library('ggvenn')

# load in the data:
subset.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mch.csv'
production.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/processed-unique-mch.csv'

# load in subset data
subset <- fread(
    input = subset.path,
    sep = ',',
    nrows = 2,
    header = FALSE,
    data.table = FALSE
    )

# load in production data:
production <- fread(
    input = subset.path,
    sep = ',',
    nrows = 2,
    header = FALSE,
    data.table = FALSE
    )

### CHECKING ######################################################################################
# check the number of genes:
stopifnot(ncol(subset) == ncol(production))
cat('number of genes in subset data produced by R: ', ncol(subset), '\n')
cat('number of genes in production data produced by python: ', ncol(production), '\n')
