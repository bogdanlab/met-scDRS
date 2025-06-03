### csv-h5ad-conversion.R #########################################################################
# purpose: convert an csv file into a h5ad file 

### PREAMBLE ######################################################################################
# load in the libraries:
library(tidyverse);
library(Seurat);
library(sceasy);
library(data.table);
library(readxl);

# define the input and its help page:
require(docopt)
'Usage:
    csv-h5ad-conversion.R [--data_matrix <input> --output_file <group>]

Options:
    --data_matrix csv file path that houses processed data matrix
    --output_file output csv
]' -> doc

# collect user input: 
opts <- docopt(doc)
csv.path <- opts$data_matrix;
output.path <- opts$output_file;

# for function testing:
# csv.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v1.2/processed-mch-gene-name.csv';
# output.path <- '/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/processed-met-scDRS-mch.csv';

### PROCESSING ####################################################################################
# grab out the processed data:
fraction <- fread(
    file = csv.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    data.table = TRUE,
    nThread = 1
    );
cat('loaded fraction! \n');
print(dim(fraction));

# extract out the rownames:
cell.names <- fraction$cell;
fraction[, cell := NULL];

# find genes that has low variance:
fraction.variance <- fraction[, lapply(.SD, var)];
fraction.variance <- unlist(fraction.variance);
low.variance.cutoff <- quantile(fraction.variance, prob = 0.05);
low.var.genes <- names(fraction.variance)[fraction.variance < low.variance.cutoff];

# remove genes that have low variance:
low.var.genes.index <- match(low.var.genes, colnames(fraction));
set(fraction, j = low.var.genes.index, value = NULL);

### CONVERSION ####################################################################################
# check the result is still bound between 1 and 0:
stopifnot(fraction >= 0 & fraction <= 1);
cols.to.update <- colnames(fraction);
fraction[, (cols.to.update) := lapply(.SD, function(x) 1 - x), .SDcols = cols.to.update]
stopifnot(fraction >= 0 & fraction <= 1);

# add back the cell id to the data:
column.order <- c('cell', colnames(fraction));
fraction[, cell := cell.names];
setcolorder(fraction, column.order);

# write the table:
write.table(
    fraction,
    file = output.path,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = ','
    );
