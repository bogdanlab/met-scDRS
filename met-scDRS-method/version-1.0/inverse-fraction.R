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
# csv.path <- '/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.csv';
# output.path <- '/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv';

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
fraction.variance <- apply(fraction, 2, sd);
low.variance.cutoff <- quantile(fraction.variance, prob = 0.05);
low.var.genes <- names(fraction.variance)[fraction.variance < low.variance.cutoff];
fraction <- fraction[, setdiff(colnames(fraction), low.var.genes)];

### CONVERSION ####################################################################################
# check the result is still bound between 1 and 0:
stopifnot(fraction >= 0 & fraction <= 1);
inverse.fraction <- 1 - fraction;
stopifnot(inverse.fraction >= 0 & inverse.fraction <= 1);

# add back the cell id to the data:
cell <- rownames(inverse.fraction);
inverse.fraction <- cbind(cell, inverse.fraction);

# write the table:
write.table(
    inverse.fraction,
    file = output.path,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = ','
    );
