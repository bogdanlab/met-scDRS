### mch-processing.R ##############################################################################
# purpose: processing script for GSE215353 mch data

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
mch.fraction.path <- paste0(data.dir, 'mch/GSE215353-concat.csv');
mcg.fraction.path <- paste0(data.dir, 'mcg/GSE215353-concat.csv');
meta.path <- paste0(data.dir, 'B06_metadata_clusters_obs.csv');

# load in methylation fractions:
mch.fraction <- fread(
    file = mch.fraction.path,
    sep = ',',
    header = TRUE,
    data.table = TRUE,
    nThread = 1
    );

# load in the meta data:
mch.meta <- fread(
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv',
    sep = ',',
    header = TRUE,
    data.table = FALSE,
    nThread = 1
    );

# load in the gene name converter function:
source(paste0(function.path, 'gene-name-converter.R'));

### PROCESSING ####################################################################################
# first assign rownames:
cell.names <- mch.fraction$cell;
mch.fraction[, cell := NULL];

# convert the gene names:
hgnc.symbol <- gene.name.converter(vec = colnames(mch.fraction), from = 'unversioned');
confusing.genes <- which(hgnc.symbol == '' | is.na(hgnc.symbol) | duplicated(hgnc.symbol));

# remove the genes that does not have a gene name from biomart:
set(mch.fraction, j = confusing.genes, value = NULL);

# rename the mch genes to the converted names:
hgnc.symbol <- gene.name.converter(vec = colnames(mch.fraction), from = 'unversioned');
colnames(mch.fraction) <- hgnc.symbol;

### OUTPUT ########################################################################################
# add back the cell index into the data table:
column.order <- c('cell', colnames(mch.fraction));
mch.fraction[, cell := cell.names];
setcolorder(mch.fraction, column.order);

# Finally, output the data:
fwrite(
    mch.fraction,
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    file = paste0(data.dir, '/v1.2/processed-mch-gene-name.csv'),
    quote = 'auto',
    nThread = 1
    );
