### get-cell-type-proportion.R ####################################################################
# PURPOSE: for scDRS significant cells, find out what proportion of them are in what cell type

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    get-cell-type-proportion.R [--score <scdrs> --meta_data <meta> --field <group> --out <output> --cutoff <p> --min <min.cell>]

Options:
    --score path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field name in meta which you would like to compute proportion of significant cells in
    --cutoff p value cutoff that user specifies
    --min minimum amount of cells that need to be significant for using it for proportion analysis
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.score.path <- opts$score;
group.index <- opts$field;
output.path <- opts$out;
p.cutoff <- as.numeric(opts$cutoff);
min.significant.cell <- as.numeric(opts$min);
system.date <- Sys.Date();

# load in the user input:
meta <- read.table(
    sep = '\t',
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# load in the scDRS score data points:
risk.score <- read.table(
    file = scDRS.score.path,
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
    );

# perform checks to make sure ordering of cell names are the same:
stopifnot(rownames(meta) == rownames(risk.score));

### ANALYSIS - PROPORTION #########################################################################
# find the set of scDRS significant cells:
risk.score$fdr.p <- p.adjust(risk.score$pval, method = 'fdr');
significant.cell <- rownames(risk.score)[risk.score$fdr.p < p.cutoff];

# add a if loop to check if there are actually significant cells:
if (length(significant.cell) > min.significant.cell) {
    # find the cell type identity for those significant cells:
    cell.type.proportion <- table(meta[significant.cell, group.index]) / length(significant.cell) * 100;
    cell.type.proportion <- data.frame(sort(cell.type.proportion, decreasing = TRUE));
    colnames(cell.type.proportion) <- c('group', 'percentage');
    cell.type.proportion$cell.number <- sort(table(meta[significant.cell, group.index]), decreasing = TRUE);

    # output the proportion:
    cat('% significant cells in each cell type outputted to: ', output.path, '\n');
    write.table(
        cell.type.proportion,
        file = output.path,
        sep = '\t',
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
        );
    } else {
    cell.type.proportion <- data.frame(
        group = unique(meta[,group.index]),
        percentage = 0,
        cell.number = 0
        );
    cat('% significant cells in each cell type outputted to: ', output.path, '\n');
    write.table(
        cell.type.proportion,
        file = output.path,
        sep = '\t',
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
        );
    }
