### mcg-processing.R ##############################################################################
# purpose: processing script for GSE215353 mcg data

### PREAMBLE ######################################################################################
# load in the libraries:
library(Seurat);
library(sceasy);
library(data.table);
library(tidyverse);

# define specified paths:
data.dir <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/';
session.save.path <- '/u/project/pasaniuc/lixinzhe/session-info/';
system.date <- Sys.Date();
save.path <- '/u/project/pasaniuc/lixinzhe/R_saves/';
code.path <- '/u/home/l/lixinzhe/project-github/methylation-RNA-xinzhe-rotation/code/';
plotting.path <- '/u/project/pasaniuc/lixinzhe/plot/';
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/'

# define path for fractioned methylation data:
mcg.fraction.path <- paste0(data.dir, 'mcg/GSE215353-concat.csv');
meta.path <- paste0(data.dir, 'B06_metadata_clusters_obs.csv');

# load in methylation fractions:
mcg.fraction <- fread(
    file = mcg.fraction.path,
    sep = ',',
    header = TRUE,
    data.table = TRUE,
    nThread = 1
    );

# load in the meta data:
mcg.meta <- fread(
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv',
    sep = ',',
    header = TRUE,
    data.table = FALSE,
    nThread = 1
    );

# load in the gene name converter function:
source(paste0(function.path, 'gene-name-converter.R'));

### PROCESSING ####################################################################################
# first assign rownames:
cell.names <- mcg.fraction$cell;
mcg.fraction[, cell := NULL];

# convert the gene names:
hgnc.symbol <- gene.name.converter(vec = colnames(mcg.fraction), from = 'unversioned');
confusing.genes <- which(hgnc.symbol == '' | is.na(hgnc.symbol) | duplicated(hgnc.symbol));

# remove the genes that does not have a gene name from biomart:
set(mcg.fraction, j = confusing.genes, value = NULL);

# rename the mcg genes to the converted names:
hgnc.symbol <- gene.name.converter(vec = colnames(mcg.fraction), from = 'unversioned');
colnames(mcg.fraction) <- hgnc.symbol;

### OUTPUT ########################################################################################
# add back the cell index into the data table:
column.order <- c('cell', colnames(mcg.fraction));
mcg.fraction[, cell := cell.names];
setcolorder(mcg.fraction, column.order);

# Finally, output the data:
fwrite(
    mcg.fraction,
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    file = paste0(data.dir, '/v1.2/processed-mcg-gene-name.csv'),
    quote = 'auto',
    nThread = 1
    );
