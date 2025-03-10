### high-variance-gene-permutation.R ############################################################
# purpose: the script extracts out genes with high variance genes and combine those genes
# with permuted random gwas weights

### PREAMBLE ######################################################################################
# load in the libraries:
library(tidyverse);
library(Seurat);
library(sceasy);
library(data.table);
library(readxl);

# load in user specified arguments:
# define the input and its help page:
require(docopt)
'Usage:
    random-genes-permutation.R [--data_matrix <input> --gs_file <gene_sets> --gene_number <number> --quantile <quantile> --replication <num_rep> --output_dir <output>]

Options:
    --data_matrix csv file path that houses processed data matrix
    --gs_file full gs file where we can choose a random gene set
    --gene_number number of genes to permute
    --replication number of replication on simulation would you like
    --quantile quantile for defining high variance genes [default: 0.75]
    --output_dir directory to output the processed gene set file
]' -> doc

# collect user input: 
opts <- docopt(doc)
csv.path <- opts$data_matrix;
gs.path <- opts$gs_file;
gene.number <- as.numeric(opts$gene_number);
replication <- as.numeric(opts$replication);
quantile.cutoff <- as.numeric(opts$quantile);
output.dir <- opts$output_dir;
seed <- seq(1, replication);

# for script testing:
# csv.path <- '/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv'
# gs.path <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs'
# gene.number <- 100;
# output.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/'
# replication <- 3;
# seed <- seq(1, replication);

### DATA LADING ###################################################################################
# load in a function used to decompose gs file:
gs.decomposer <- function(gs) {
    # split by , and get the genes:
    gs.split <- unlist(strsplit(gs, ','));
    gs.genes <- gsub(':.*', '', gs.split);
    gs.weight <- gsub('.*:', '', gs.split);
    gs.weight <- as.numeric(gs.weight);
    names(gs.weight) <- gs.genes;

    # return the gs.weight:
    return(gs.weight)
    }

# load in the data:
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

# load in gs file
trait.gs <- read.table(
    file = gs.path,
    sep = '\t',
    header = TRUE
    );
# split gene sets:
trait.gene.set <- lapply(trait.gs$GENESET, gs.decomposer);
names(trait.gene.set) <- trait.gs$TRAIT;

### PERMUTATION ###################################################################################
# compute the top variance:
feature.variance <- apply(fraction, 2, var);
top.variance.gene <- feature.variance[feature.variance > quantile(feature.variance, probs = quantile.cutoff)];

for (simulation in seed) {
    set.seed(simulation + gene.number);

    # grab random genes from the data:
    random.genes <- sample(names(top.variance.gene), size = gene.number, replace = FALSE);

    # grab random trait:
    random.trait <- sample(trait.gs$TRAIT, size = 1);
    random.weight <- trait.gene.set[[random.trait]];
    random.weight <- sort(random.weight, decreasing = TRUE);

    # grab out the top X weight:
    random.weight.subset <- head(random.weight, n = gene.number)

    # permute the subsetted weight:
    randomized.weight <- sample(random.weight.subset);
    names(randomized.weight) <- random.genes;

    # make randomized weight into a gs format:
    random.gs <- paste0(random.genes, ':', randomized.weight);
    random.gs <- paste0(random.gs, collapse = ',');

    ### DATA CHECK ################################################################################
    # decompose the gs:
    decomposed.random.gs <- gs.decomposer(random.gs);
    
    # check that the genes weight that were in the gs actually have the highest gs weight:
    random.trait.top.weight.gene <- tail(sort(trait.gene.set[[random.trait]]), gene.number)
    stopifnot(random.trait.top.weight.gene %in% decomposed.random.gs);
    stopifnot(length(random.trait.top.weight.gene) == gene.number);

    # check that the genes that were selected actually have the highest variance:
    sorted.variance <- sort(feature.variance, decreasing = TRUE);
    gene.only <- names(decomposed.random.gs);
    stopifnot(gene.only %in% names(sorted.variance[1:ceiling(ncol(fraction) * (1 - quantile.cutoff))]));
    stopifnot(length(gene.only) %in% gene.number);

    # ensure that the genes are not sorted in the output:
    stopifnot(is.unsorted(decomposed.random.gs));

    ### OUTPUT ####################################################################################
    output.path <- paste0(output.dir, 'seed-', simulation, '-high-variance-genes-permuted-weights-', gene.number, '-genes.gs');
    gs.df <- data.frame(
        TRAIT = paste0('NULL_SIM_TOP_VARIANCE_', gene.number, '_genes_seed_', simulation),
        GENESET = random.gs
        );

    write.table(
        x = gs.df,
        file = output.path,
        quote = FALSE,
        sep = '\t',
        row.names = FALSE,
        col.names = TRUE
        );
    
    # log information about the file:
    # Open an existing log file for appending
    log.path <- gsub('\\.gs', '\\.log', output.path);
    log.file <- file(log.path, open = "a");

    # Use cat to append to the log file
    cat('seed: ', simulation, '\n', file = log.file)
    cat('random trait selected: ', random.trait, '\n', file = log.file);
    cat('random gene set outputted to: \n', output.path, '\n', file = log.file);

    # Close the log file when you're done
    close(log.file);
}








