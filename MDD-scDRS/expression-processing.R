### expression-processing.R #################################################################################
# PURPOSE: convert expression matrix into a h5ad object that can be used by the scDRS:

### PREAMBLE ######################################################################################
# load in the libraries:
library(Seurat);
library(sceasy);

# define specified paths:
data.dir <- '/u/project/geschwind/lixinzhe/data/GSE144136-MDD/';
session.save.path <- '/u/project/pasaniuc/lixinzhe/session-info/';
system.date <- Sys.Date();
save.path <- '/u/project/pasaniuc/lixinzhe/R_saves/';
code.path <- '/u/home/l/lixinzhe/project-github/methylation-RNA-xinzhe-rotation/code/';
plotting.path <- '/u/project/pasaniuc/lixinzhe/plot/';
h5ad.output.path <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/MDD-RNA.h5ad';
covariate.output.path <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/MDD-covariates.txt';

# load in data:
mdd.rna <- ReadMtx(
    mtx = paste0(data.dir, 'GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz'),
    cells = paste0(data.dir, 'GSE144136_CellNames.csv.gz'),
    features = paste0(data.dir, 'GSE144136_GeneNames.csv.gz'),
    cell.column = 2,
    feature.column = 2,
    cell.sep = ',',
    feature.sep = ',',
    skip.cell = 1,
    skip.feature = 1,
    mtx.transpose = FALSE,
    unique.features = FALSE,
    strip.suffix = FALSE
    );

# load in the remaining list of cells after filtering:
filtered.cell.path <- list.files(
    path = '/u/project/geschwind/lixinzhe/data/GSE144136-MDD',
    pattern = '.*_filtered_cells.csv.gz',
    full.names = TRUE
    );
filter.label <- gsub('.*/GSE144136_' , '', filtered.cell.path);
filter.label <- gsub('_filtered_cells.csv.gz', '', filter.label);
names(filtered.cell.path) <- filter.label;

# load in the cell names:
filter.cell <- vector('list', length = length(filter.label));
names(filter.cell) <- filter.label;
for (cell.type in filter.label) {
    filter.cell[[cell.type]] <- read.table(
        file = filtered.cell.path[cell.type],
        sep = ',',
        row.names = 1,
        header = TRUE
        )[, 1];
    }
cell.keep <- Reduce('union', filter.cell);

### PREPARE EXPRESSION ############################################################################
# grab out the cell type information of each cell:
cell.type <- sapply(strsplit(colnames(mdd.rna), '\\.'), '[[', 1);

# grab out the donor information:
remaining.info <- sapply(strsplit(colnames(mdd.rna), '\\.'), '[[', 2);
splitted.string <- strsplit(remaining.info, '_');

# grab out the rest of the covariates:
covariate.name <- c('donor', 'diagnosis', 'batch', 'barcode');
meta.data <- matrix(
    unlist(splitted.string),
    ncol = length(covariate.name),
    byrow = T
    );
colnames(meta.data) <- covariate.name;

# prepare the meta data for the seurat object:
meta.data <- data.frame(meta.data);
meta.data$cell.type <- cell.type;
rownames(meta.data) <- colnames(mdd.rna);

# make seurat object:
seurat.count <- CreateSeuratObject(
    counts = mdd.rna[, cell.keep],
    project = "CreateSeuratObject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0
    );

# output the h5ad file:
convertFormat(
    seurat.count,
    from = "seurat",
    to = "anndata",
    outFile = h5ad.output.path
    );

# output the seurat object with the meta:
seurat.mdd.rna <- CreateSeuratObject(
    counts = mdd.rna[, cell.keep],
    project = "CreateSeuratObject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0,
    meta.data = meta.data[cell.keep, ]
    );
saveRDS(seurat.mdd.rna, paste0(save.path, 'GSE144136-seurat-meta.rds'));

### PREPARE COVARIATES ############################################################################
# prepare design matrix for regression:
design.matrix <- model.matrix(
    ~ 0 +
    batch,
    data = meta.data
    );

# test for colinearity:
linear.combination.check <- alias(mdd.rna[1,] ~ design.matrix);
alias.features <- gsub('design.matrix','',rownames(linear.combination.check$Complete));
full.rank.design <- design.matrix[, setdiff(colnames(design.matrix), alias.features)];

# add back the constant term:
const <- 1;
full.rank.design <- cbind(const, full.rank.design); 

# output the full rank matrix:
stopifnot(is.numeric(full.rank.design));

# return the covariate matrix:
write.table(
    full.rank.design,
    file = covariate.output.path,
    sep = '\t',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );
