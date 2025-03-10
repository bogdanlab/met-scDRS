### meta-processing.R #############################################################################
# purpose: concatenate the meta data, make a plotting data frame if needed

### PREAMBLE ######################################################################################
# load in libraries:
library(data.table);

# spcifify meta data file directories:
mch.meta.dir <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mch/'
mcg.meta.dir <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/mcg/'

# load in the files:
mch.metas <- list.files(mch.meta.dir, pattern = '-meta');
mcg.metas <- list.files(mcg.meta.dir, pattern = '-meta');

# initiate the meta file placeholder
mch.meta.file <- vector('list', length = length(mch.metas));
names(mch.meta.file) <- mch.metas;

# initiate the meta file place holder:
mcg.meta.file <- vector('list', length = length(mcg.metas));
names(mcg.meta.file) <- mcg.metas;

# for each of the meta file, read it in and piece them together:
for (file in mch.metas) {
    meta.file <- fread(
        file = paste0(mch.meta.dir, file),
        sep = ',',
        header = TRUE,
        data.table = FALSE,
        nThread = 1
        )
    mch.meta.file[[file]] <- meta.file;
    }

for (file in mcg.metas) {
    meta.file <- fread(
        file = paste0(mcg.meta.dir, file),
        sep = ',',
        header = TRUE,
        data.table = FALSE,
        nThread = 1
        )
    mcg.meta.file[[file]] <- meta.file;
    }

# load in the embedding:
# first grab out the embedding for cells that are within the mch:
files.to.read <- list.files('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mch/');
embedding <- vector('list', length = length(files.to.read));
names(embedding) <- files.to.read;

# read in the rds file:
for (rds.file in files.to.read) {
    methylation<- readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mch/', rds.file));
    embedding[[rds.file]] <- methylation@reductions$umap@cell.embeddings;
    rm(list = c('methylation'))
    gc()
    }

# concatenate into one file:
mch.umap <- data.frame(Reduce('rbind', embedding));

# next for mcg:
files.to.read <- list.files('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mcg/');
embedding <- vector('list', length = length(files.to.read));
names(embedding) <- files.to.read;

# read in the rds file:
for (rds.file in files.to.read) {
    methylation<- readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/mcg/', rds.file))
    embedding[[rds.file]] <- methylation@reductions$umap@cell.embeddings;
    rm(list = c('methylation'))
    gc()
    }

# concatenate into one file:
mcg.umap <- data.frame(Reduce('rbind', embedding));

### PROCESSING ####################################################################################
# combine the metas:
mcg.meta <- data.frame(Reduce('rbind', mcg.meta.file), row.names = 1);
mch.meta <- data.frame(Reduce('rbind', mch.meta.file), row.names = 1);

# next add the UMAP into the meta data:
mcg.meta <- cbind(mcg.meta, mcg.umap[rownames(mcg.meta), ])
mch.meta <- cbind(mch.meta, mch.umap[rownames(mch.meta), ])

### OUTPUT ########################################################################################
# grab out the cell name and put it into the matrix for output:
cell <- row.names(mch.umap);
mch.umap <- cbind(cell, mch.umap);
cell <- row.names(mcg.umap);
mcg.umap <- cbind(cell, mcg.umap);

# write out the embedding space:
fwrite(
    mch.umap,
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/full-mch-umap.csv',
    row.names = FALSE,
    col.names = TRUE,
    quote = 'auto',
    nThread = 1
    )
fwrite(
    mcg.umap,
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/full-mcg-umap.csv',
    row.names = FALSE,
    col.names = TRUE,
    quote = 'auto',
    nThread = 1
    )

# write out the table:
cell <- row.names(mch.meta);
mch.meta <- cbind(cell, mch.meta);
fwrite(
    mch.meta,
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv',
    quote = 'auto',
    nThread = 1
    );

# for mcg
cell <- row.names(mcg.meta);
mcg.meta <- cbind(cell, mcg.meta)
fwrite(
    mcg.meta,
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    file = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mcg-meta.csv',
    quote = 'auto',
    nThread = 1
    );
