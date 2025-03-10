### data-validation.R #############################################################################
# purpose: spot check computed fraction

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    data-validation.R [--count <count_matrix> --coverage <coverage_matrix> --fraction <fraction_matrix> --meta <meta_file> --lines <lines_check>]

Options:
    --count count matrix in csv style (first column = rownames)
    --coverage coverage matrix in csv style
    --fraction fraction matrix computed from count / coverage
    --meta path to meta data on cells associated in the dataset
    --lines number of lines to read for checking the data
]' -> doc

### DATA LOADING ##################################################################################
# collect user input: 
opts <- docopt(doc)
count.path <- opts$count;
coverage.path <- opts$coverage;
meta.data.path <- opts$meta;
fraction.path <- opts$fraction;
lines <- as.numeric(opts$lines);
system.date <- Sys.Date();
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/';

# load in the gene name converter function:
source(paste0(function.path, 'gene-name-converter.R'));

# script testing:
# count.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mcg_count.csv"
# coverage.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mcg_coverage.csv"
# meta.data.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx"
# fraction.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/unique_mcg_gene_fraction.csv"
# lines <- 1000;

# load libraries:
library(readxl);
library(data.table);

# load meta data:
snmc.meta <- read_excel(meta.data.path, skip = 14);
snmc.meta <- as.data.frame(snmc.meta);
rownames(snmc.meta) <- snmc.meta[, '...1'];

# load in the data:
count <- fread(
    file = count.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    nrow = lines,
    data.table = FALSE,
    nThread = 1
    );
coverage <- fread(
    file = coverage.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    nrow = lines,
    data.table = FALSE,
    nThread = 1
    );
fraction <- fread(
    file = fraction.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    nrow = lines,
    data.table = FALSE,
    nThread = 1
    );

### DATA FORMATING ################################################################################
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

# rename the column names and rownames:
colnames(count) <- colnames(coverage) <- hgnc.symbol;
rownames(count) <- count[, 'cell'];
rownames(coverage) <- coverage[, 'cell'];
stopifnot(rownames(coverage) == rownames(count));
rownames(fraction) <- fraction[, 'cell'];

# remove the cell field from the data:
count <- count[, -match('cell', colnames(count))];
coverage <- coverage[, -match('cell', colnames(coverage))];
fraction <- fraction[, -match('cell', colnames(fraction))]

# subset the count data, coverage data and fraction data to the same set of rows and columns:
# now grab out the set of ids in fraction:
quality.cell <- intersect(rownames(fraction), rownames(count))
count <- count[quality.cell, not.confusing.genes];
coverage <- coverage[quality.cell, not.confusing.genes];
fraction <- fraction[quality.cell, ];

### DATA CHECKING #################################################################################
# check lossless gene names:
stopifnot(colnames(coverage) == colnames(fraction));
stopifnot(colnames(count) == colnames(fraction));

# generate random integers:
check.num <- 100;
rows.index <- sample(1:length(quality.cell), size = check.num, replace = TRUE);
columns.index <- sample(1:ncol(count), size = check.num, replace = TRUE);

for (index in seq(1, check.num)) {
    # get the indexes:
    count.cell <- count[rows.index[index], columns.index[index]];
    coverage.cell <- coverage[rows.index[index], columns.index[index]];
    fraction.cell <- fraction[rows.index[index], columns.index[index]];
    
    # spot check:
    if (coverage.cell == 0){
        check.cell <- 0;
    } else {
        check.cell <- count.cell / coverage.cell;
    }
    stopifnot(signif(check.cell, digits = 10) == signif(fraction.cell, 10));
}

