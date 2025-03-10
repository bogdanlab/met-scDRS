### get-magma.stats.R #############################################################################
# PURPOSE: this is a part of the pipeline script of scDRS for generating z score statistics from 
# MAGMA outputs using custom GWAS

### PREAMBLE ######################################################################################
# gather arguments that were passed from command line:
# The first argument should be the output of the magma
# The second argument should be the ncbi file that used in magma to generate the summary
# the third argument should be the output path

# gather arguments from command line:
args <- commandArgs(trailingOnly = TRUE);

# check argument and assign input path variable:
if (!is.character(args[1])) {
    stop('First argument should be path for reading in magma file.');
    }   else {
    # read in the table of the magma output:
    magma.out.path <- args[1];
    magma.summary <- read.table(
        file = magma.out.path,
        row.names = 1,
        header = TRUE,
        stringsAsFactors = FALSE
        );
    }

# check the argument and assign ncbi loc path variable:
if (!is.character(args[2])) {
    stop('Second argument should be path for ncbi.gene.loc file.');
    }  else {
    # read in the ncbi location file:
    ncbi.loc.path <- args[2];
    ncbi.loc <- read.table(
        file = ncbi.loc.path,
        sep = '\t',
        header = FALSE,
        stringsAsFactors = FALSE
        );
    colnames(ncbi.loc) <- c('entrez', 'chr', 'start', 'end', 'strand', 'hgnc');
    }

# check the argument and assign output path variable
if (!is.character(args[3])) {
    stop('Third argument should be path for outputing the magam stats.');
    }   else {
    output.path <- args[3];
    }

# match the entrez id from the input:
gs.munge.input <- data.frame(
    gene = ncbi.loc$hgnc[match(rownames(magma.summary), ncbi.loc$entrez)],
    trait = magma.summary$ZSTAT
    );

# output the result:
cat('gs-munge zscore file prepared:\n', output.path, '\n');
write.table(
    x = gs.munge.input,
    file = output.path,
    quote = FALSE,
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE
    );