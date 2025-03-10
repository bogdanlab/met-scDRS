### data-format-convert.R ######################################################################
# purpose: output methylation data as a h5ad object for scDRS

### PREAMBLE ######################################################################################
# load in the libraries:
library(tidyverse);
library(Seurat);
library(sceasy);
library(data.table);
library(readxl);

# define specified paths:
data.dir <- '/u/project/geschwind/lixinzhe/data/';
session.save.path <- '/u/project/pasaniuc/lixinzhe/session-info/';
system.date <- Sys.Date();
save.path <- '/u/project/pasaniuc/lixinzhe/R_saves/';
plotting.path <- '/u/home/l/lixinzhe/project-geschwind/plot/';
h5ad.mch.path <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/subset-GSE132489-mch.h5ad';
h5ad.mcg.path <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/subset-GSE132489-mcg.h5ad';
covariate.prefix <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/subset-GSE132489-covariates-'

# define data input paths:
mch.path <- "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/subset_randomized_mch_gene_fraction.csv"
mcg.path <- "/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/subset_randomized_mcg_gene_fraction.csv"
meta.data.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx"

# define input looping structure:
data.path <- c(mcg.path, mch.path);
output.path <- c(h5ad.mcg.path, h5ad.mch.path);
names(data.path) <- names(output.path) <- c('mcg', 'mch');

# load meta data:
snmc.meta <- read_excel(meta.data.path, skip = 14);
snmc.meta <- data.frame(snmc.meta, row.names = 1);
high.quality.cell <- rownames(snmc.meta)[snmc.meta$Pass.QC]

### PREPARE COVARIATES ############################################################################
# prepare empty list to house the nfeatrue calcualted:
nfeature <- vector('list', length = length(data.path));
names(nfeature) <- c('mcg', 'mch');

# for each modality, convert data:
for(modality in c('mcg', 'mch')) {
    # first we will load in the mcg results:
    cat('processing: ', modality, '\n');
    cat('data loaded from: ', data.path[modality], '\n');
    fraction <- fread(
        file = data.path[modality],
        header = TRUE,
        sep = ',',
        stringsAsFactors = FALSE,
        data.table = TRUE,
        nThread = 1
        );
    fraction <- fraction %>%
        column_to_rownames('cell');

    # quick data check:
    stopifnot(rownames(fraction) %in% high.quality.cell);

    ### CONVERSION ####################################################################################
    # get the set of data to regress out:
    reversed.fraction <- t(1 - fraction);

    # check the result is still bound between 1 and 0:
    stopifnot(reversed.fraction >= 0 & reversed.fraction <= 1);

    # put the result into Seurat object:
    fraction.seurat <- CreateSeuratObject(
        counts = reversed.fraction,
        min.cells = 0,
        min.features = 0
        );
    
    # prepare a design matrix:
    design.matrix <- model.matrix(
        ~ 0 +
        Sample,
        data = snmc.meta[colnames(reversed.fraction),]
        );
    convertFormat(
        fraction.seurat,
        from = "seurat",
        to = "anndata",
        outFile = output.path[modality]
        );
    
    # prepare covariates:
    # add the n_feature to covariates:
    nfeature[[modality]] <- fraction.seurat$nFeature_RNA;
    stopifnot(rownames(fraction.seurat@meta.data) == rownames(design.matrix))
    design.matrix <- cbind(fraction.seurat$nFeature_RNA, design.matrix);
    colnames(design.matrix)[1] <- paste0('nfeature_', modality)

    # test for colinearity:
    linear.combination.check <- alias(reversed.fraction[1,] ~ design.matrix);
    alias.features <- gsub('design.matrix','',rownames(linear.combination.check$Complete));
    full.rank.design <- design.matrix[, setdiff(colnames(design.matrix), alias.features)];

    # add in the constant and rownames into the matrix:
    const <- 1;
    cell <- rownames(full.rank.design);
    full.rank.design <- cbind(const, full.rank.design); 
    full.rank.design <- cbind(cell, full.rank.design); 

    covariate.output.path <- paste0(covariate.prefix, modality, '.txt');
    write.table(
        full.rank.design,
        file = covariate.output.path,
        sep = '\t',
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
        );

    # print progress:
    cat(modality, 'has been processed, data outoutted to: ', output.path[modality], '\n')
    }



