### data-QC.R #####################################################################################
# purpose: for the concatenated methylation coverage and counts
# conduct QC based on meta data and curate the data itself and a more usable set of meta data

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    data-QC.R [--count <count_matrix> --coverage <coverage_matrix> --meta <meta_file> --output <output_dir>]

Options:
    --count count matrix in csv style (first column = rownames)
    --coverage coverage matrix in csv style
    --meta path to meta data on cells associated in the dataset
    --output directory to output file
]' -> doc

# collect user input: 
opts <- docopt(doc)
count.path <- opts$count;
coverage.path <- opts$coverage;
meta.data.path <- opts$meta;
output.path <- opts$output;
system.date <- Sys.Date();
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/';

# define the chunking parameters:
lines.num <- 114706;
chunk.size <- 5000;
n.chunks <- floor(lines.num / chunk.size);

# load in the gene name converter function:
source(paste0(function.path, 'gene-name-converter.R'));

# script testing:
# count.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mch_count.csv"
# coverage.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mch_coverage.csv"
# meta.data.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx"
# output.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/chunk/mch_gene_qced_fraction"

# load libraries:
library(readxl);
library(data.table);

# load meta data:
snmc.meta <- read_excel(meta.data.path, skip = 14);
snmc.meta <- as.data.frame(snmc.meta);
rownames(snmc.meta) <- snmc.meta[, '...1'];

### DATA HANDLING ###########################################################################################
# load the first small chunks of the data for colnames and such:
count <- fread(
    file = count.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    nrow = 30,
    data.table = FALSE,
    nThread = 1
    );
coverage <- fread(
    file = coverage.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    nrow = 30,
    data.table = FALSE,
    nThread = 1
    );
stopifnot(colnames(count) == colnames(coverage));

# genes:
gene <- colnames(count);

# grab out the set of genes that can be mapped to hgnc symbols:
hgnc.symbol <- gene.name.converter(
    vec = gene,
    from = 'versioned',
    mart.tsv = '/u/project/geschwind/lixinzhe/data/2023-08-09-mart-grcm39-gene-names-dictionary.txt'
    );
confusing.genes <- which(hgnc.symbol == '' | is.na(hgnc.symbol) | duplicated(hgnc.symbol));
not.confusing.genes <- hgnc.symbol[-confusing.genes];
hgnc.symbol[confusing.genes] <- gene[confusing.genes];

# now grab out the set of ids that passed the QC:
quality.cell <- rownames(snmc.meta)[snmc.meta[, 'Pass QC']];

# initiate placeholder for count the processed cells:
processed.num <- rep(-1, times = n.chunks);

# process the data one chunk at a time:
for (chunk in 0:n.chunks) {
    cat('processing chunk', chunk + 1, '\n');

    # load count data:
    count <- fread(
        file = count.path,
        header = 'auto',
        sep = ',',
        stringsAsFactors = FALSE,
        nrow = chunk.size,
        skip = chunk * chunk.size,
        data.table = FALSE,
        nThread = 1
        );

    # load coverage data:
    coverage <- fread(
        file = coverage.path,
        header = 'auto',
        sep = ',',
        stringsAsFactors = FALSE,
        nrow = chunk.size,
        skip = chunk * chunk.size,
        data.table = FALSE,
        nThread = 1
        );

    # rename the column names and rownames:
    colnames(count) <- colnames(coverage) <- hgnc.symbol;
    rownames(count) <- count[, 'cell'];
    rownames(coverage) <- coverage[, 'cell'];
    stopifnot(rownames(coverage) == rownames(count));

    # remove the cell field from the data:
    count <- count[, -match('cell', colnames(count))];
    coverage <- coverage[, -match('cell', colnames(coverage))];

    # get the fraction:
    quality.cell.in.slice <- intersect(quality.cell, rownames(count));
    fraction <- count[quality.cell.in.slice, not.confusing.genes] / 
        coverage[quality.cell.in.slice, not.confusing.genes];
    
    # perform data check:
    stopifnot(is.nan(fraction[coverage[quality.cell.in.slice, not.confusing.genes] == 0]));

    # refill cells with 0 coverage in fraction to 0:
    fraction[coverage[quality.cell.in.slice, not.confusing.genes] == 0] <- 0;

    # check again if there are any NaNs:
    stopifnot(sum(apply(coverage[quality.cell.in.slice, not.confusing.genes], 2, is.nan)) == 0);

    # save the number of cells processed in this chunk:
    processed.num[chunk + 1] <- nrow(fraction);

    # write the data out for every chunks:
    fwrite(
        fraction,
        file = paste0(output.path, '-chunk', chunk, '.csv'),
        sep = ',',
        col.names = TRUE,
        row.names = TRUE,
        quote = 'auto'
        );
    }

# print out the total number of cells processed:
cat('processed', sum(processed.num), 'cells \n');
cat('total high QC cells', length(quality.cell), 'cells \n');
