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
    --output_file h5ad output path
]' -> doc

# collect user input: 
opts <- docopt(doc)
csv.path <- opts$data_matrix;
h5ad.path <- opts$output_file;

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
fraction <- fraction %>%
    column_to_rownames('cell');

cat('dimension of data loaded for conversion:\n');
print(dim(fraction));

# remove all genes that has zero variance:
zero.var.genes <- apply(fraction, 2, sd);
zero.var.genes <- names(zero.var.genes)[zero.var.genes == 0]
fraction <- fraction[, setdiff(colnames(fraction), zero.var.genes)]

### CONVERSION ####################################################################################
# check the result is still bound between 1 and 0:
stopifnot(fraction >= 0 & fraction <= 1);
fraction <- t(fraction); # transposed so the cells in column, genes in row

# put the result into Seurat object:
fraction.seurat <- CreateSeuratObject(
    counts = fraction,
    min.cells = 0,
    min.features = 0
    );

# printing the first 5 cells and first 5 genes of the :
cat('printing first 5 cells and genes: \n');
print(fraction[1:5,1:5]);

# convert the seurat object into the h5ad format:
convertFormat(
    fraction.seurat,
    from = "seurat",
    to = "anndata",
    outFile = h5ad.path
    );
