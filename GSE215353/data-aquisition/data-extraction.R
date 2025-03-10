### data-extraction.R #############################################################################
# purpose: extract the data out of the rds object into csv file

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    data-extraction.R [--rds <rds_path> --percentage <percent> --out <output_path> ]

Options:
    --rds_path path to which the rds object will be read in
    --percentage the percent of cells that you would like to subset to
    --output_path file path at which the fraction should be output to
]' -> doc

### DATA LOADING ##################################################################################
# collect user input: 
opts <- docopt(doc)
rds.path <- opts$rds_path;
output.path <- opts$output_path;
subset.percent <- as.numeric(opts$percentage);
system.date <- Sys.Date();
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/';

# load in the gene name converter function:
source(paste0(function.path, 'gene-name-converter.R'));

# script testing:
# rds.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mch/940105e9-6d36-4a70-9712-de3a808c4a09.rds"
# output.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/"
# subset.percent <- 0.2;

# load libraries:
library(Seurat);
library(data.table);

# set a seed:
set.seed(103);

# load in the RDS data:
methylation <- readRDS(rds.path);

# extract meta data:
meta.data <- methylation@meta.data;

# prepare the fraction for output:
fraction <- methylation@assays$RNA@data;

# remove the seurat object from memory to save space
rm(list=c('methylation'));
gc();

# get the cell id to subset to:
cell.number <- round(subset.percent * ncol(fraction));
subset.cell <- sample(colnames(fraction), cell.number, replace = FALSE);
fraction <- fraction[, subset.cell];
fraction <- data.frame(t(fraction));
cell <- rownames(fraction);
fraction <- cbind(cell, fraction);

### OUTPUT ########################################################################################
# sample name:
sample.name <- gsub('.rds', '', gsub('.*/', '', rds.path))

# write data:
cat('writing fraction for ', sample.name, '\n')
fwrite(
    fraction,
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    file = paste0(output.path, sample.name, '-fraction.csv'),
    quote = 'auto',
    nThread = 1
    );

fwrite(
    meta.data[subset.cell, ],
    sep = ',',
    row.names = TRUE,
    col.names = TRUE,
    file = paste0(output.path, sample.name, '-meta.csv'),
    quote = 'auto',
    nThread = 1
    );
