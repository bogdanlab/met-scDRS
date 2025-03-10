### add-annotation-h5ad.R #########################################################################
# purpose: add cell level meta to an h5ad object
# this function is meant to be called from command line

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    significant-cells-visualization-script.R [--fraction <fraction> --meta_data <meta> --field <group> --out <output>]

Options:
    --fraction the csv file that we would need to add the cell wise annotation
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field name which column need to be added to the h5ad adata.obs.columns ; can be "none"
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
fraction.path <- opts$fraction;
group.index <- opts$field;
output.path <- opts$out;
system.date <- Sys.Date();

# for testing:
# fraction.path <- "/u/project/geschwind/lixinzhe/data/2024_01_29_Phase3data/processed/version-1.2/processed-met-scDRS-mch-ASD-bimodal-methylation.csv"
# meta.data.path <- "/u/project/geschwind/lixinzhe/data/2024_01_29_Phase3data/processed/version-1.2/bimodal-meta-data.csv" 
# group.index <- "joint_subclass"
# output.path <- '/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/test-scripting.h5ad'

# load in the libraries:
library(tidyverse);
library(Seurat);
library(sceasy);
library(data.table);
library(readxl);

# grab out the processed data:
fraction <- fread(
    file = fraction.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    data.table = TRUE,
    nThread = 1
    );
fraction <- fraction %>%
    column_to_rownames('cell');

cat('fraction loaded \n')

# read in the meta data:
meta <- read.table(
    sep = ',',
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );
cat('meta loaded \n')

### CONVERT DATA ##################################################################################
# make sure that the cell id are identical:
stopifnot(rownames(meta) == rownames(fraction))

# first make a collection of Seurat object:
cat('creating seurat object from fraction \n')
if (group.index == 'none'){
    fraction.seurat <- CreateSeuratObject(
        counts = t(fraction),
        meta.data = meta,
        min.cells = 0,
        min.features = 0
        );
    } else{
        fraction.seurat <- CreateSeuratObject(
            counts = t(fraction),
            min.cells = 0,
            min.features = 0
            );
        annotation = meta[, group.index];
        fraction.seurat@meta.data <- cbind(fraction.seurat@meta.data, annotation)
    }
cat('fraction prepared in seurat format, converting now \n')

# remove the fraction to save memory:
rm(list = c('fraction'));
gc()

# convert the seurat object into the h5ad format:
convertFormat(
    fraction.seurat,
    from = "seurat",
    to = "anndata",
    outFile = output.path
    );
cat('conversion completed! \n')
