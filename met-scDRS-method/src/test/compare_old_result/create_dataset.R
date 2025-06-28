### gse132489_small_sample.R ######################################################################
# create a smaller sample data from the subsetted csv data for faster development:
library(tidyverse);
library(Seurat);
library(sceasy);
library(data.table);
library(readxl);

csv.path <-  "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/subset_randomized_mch_gene_fraction_copy.csv"
h5ad.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/subset_randomized_mch_gene_fraction_copy.h5ad"

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
