### seurat-h5ad.R #################################################################################
# PURPOSE: convert custom seurat object into a h5ad object that can be used by the scDRS:

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
    'RData' = load(seurat.path)
    );

# check if the seurat object is successfully loaded:
if(is.null(seurat.obj)) {
    stop('Seurat object is not loaded, terminating!')
    }

# get the output path:
if (!is.character(args[2])) {
    stop('Second argument should be path for outputing the magam stats.');
    }   else {
    output.path <- args[2];
    }
### CONVERSION ####################################################################################
# grab out the count matrix:
count <- seurat.obj@assays$RNA@counts;

# make a new seurat object:
seurat.count <- CreateSeuratObject(
    counts = count,
    project = "CreateSeuratObject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0
    );

# make sure the seurat count matrixes are identical:
stopifnot(identical(count, seurat.count@assays$RNA@counts));

# now convert it to a count:
convertFormat(
    seurat.count,
    from = "seurat",
    to = "anndata",
    outFile = output.path 
    );
