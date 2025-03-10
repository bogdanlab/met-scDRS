### gene-name-converter.R #########################################################################
# purpose: convert gene ids from between forms:

# INPUT:
# vec: character vectors that contains the gene name:
# from: what is the source gene name, e.g.: ensembl gene id? or ensembl transcript id?

gene.name.converter <- function(
    vec,
    from = 'versioned',
    mart.tsv = '/u/project/geschwind/lixinzhe/data/2023-08-09-mart-grch38-gene-names-dictionary.txt',
    human.mouse.tsv = '/u/home/l/lixinzhe/project-geschwind/data/2023-10-15-mouse-human-biomart-grch38-p14.txt') {
    # read in the dictionary:
    mart <- read.table(
        file = mart.tsv,
        sep = '\t',
        header = TRUE
        );
    human.mouse.mart <- read.table(
        file = human.mouse.tsv,
        sep = '\t',
        header = TRUE
        );

    if (from == 'versioned') {
        # if we need to convert from transcript id, strip the version number:
        stripped <- gsub('-.*', '' ,gsub('\\.', '-', vec));

        # next grab out from our dictionary:
        gene.name <- mart$Gene.name[match(stripped, mart$Gene.stable.ID)];
        
        # return the result:
        return(gene.name);
        
        } else if (from == 'unversioned') {
        # directly convert gene names:
        gene.name <- mart$Gene.name[match(vec, mart$Gene.stable.ID)]

        # return result:
        return(gene.name);

        } else if (from == 'mouse_to_human') {
            gene.name <- human.mouse.mart$Gene.name[match(vec, human.mouse.mart$Mouse.gene.name)];
            return(gene.name);
        } else {
        # if the from statement is not recognized:
        cat('supported argument: versioned, unversioned, mosue_to_human');
        break
        }

    }
