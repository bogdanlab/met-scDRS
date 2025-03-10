### gs-decomposer.R ###############################################################################
# purpose: decompose the gene set file and split it into individual genes:

# INPUT:
# gs: the GENESET column of the gs table

# OUTPUT:
# gs.weight: a named numbers where number are the gs weights and the names are the genes

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
