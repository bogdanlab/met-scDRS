### preprocess.R ##################################################################################
# purpose: preprocessing smartseq dataset:

### PREAMBLE ######################################################################################
# load in libraries:
library(Seurat)
library(data.table)
library(sceasy)

# define paths:
expression_path <- '/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/matrix.csv'
meta_path <- '/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/metadata.csv'
tsne_path <- '/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/tsne.csv'

# load in data:
aibs_expr <- fread(expression_path, sep = ',', data.table = FALSE)
meta <- fread(meta_path, sep = ',', data.table = FALSE)
tsne <- fread(tsne_path, sep = ',', data.table = FALSE)

### PROCESS #######################################################################################
# get common cells that have tsne, meta, and expr:
rownames(aibs_expr) <- aibs_expr[, 'sample_name']
rownames(meta) <- meta[, 'sample_name']
rownames(tsne) <- tsne[, 'sample_name']
aibs_expr <- aibs_expr[, -match('sample_name', colnames(aibs_expr))]
cell_order <- intersect(rownames(aibs_expr), intersect(meta$sample_name, tsne$sample_name))

# re-align dataset:
aibs_expr <- aibs_expr[cell_order, ]
meta <- meta[cell_order, ]
tsne <- tsne[cell_order, ]

# bind the tsne into meta:
meta <- cbind(meta, tsne)

# obtain design matrix:
design_matrix <- model.matrix(
    ~ 0 +
    donor_sex_label +
    external_donor_name_label,
    data = meta
    );

# test for colinearity:
linear_combination_check <- alias(aibs_expr[,1] ~ design_matrix);
alias_features <- gsub('design_matrix','',rownames(linear_combination_check$Complete));
full_rank_design <- design_matrix[, setdiff(colnames(design_matrix), alias_features)];

# add back the constant term:
const <- 1;
full_rank_design <- cbind(const, full_rank_design); 
cell <- rownames(full_rank_design)
full_rank_design <- cbind(cell, full_rank_design)

# return the covariate matrix:
cov_out <- '/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/expr_covariate.cov'
write.table(
    full_rank_design,
    file = cov_out,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
    );

# next output h5ad file:
seurat_count <- CreateSeuratObject(
    counts = t(aibs_expr), # transpose to cells in column
    project = "CreateSeuratObject",
    assay = "RNA",
    min.cells = 0,
    min.features = 0
    );

# add in checks:
stopifnot(colnames(seurat_count) == rownames(full_rank_design));

# output the h5ad file:
h5ad_out <- '/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/expr.h5ad'
convertFormat(
    seurat_count,
    from = "seurat",
    to = "anndata",
    outFile = h5ad_out
    );

# output meta file:
meta_out <- '/u/home/l/lixinzhe/project-geschwind/data/AIBS_human_smartseq/meta.csv'
write.table(
    meta,
    file = meta_out,
    sep = ',',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
    );