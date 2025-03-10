expr.ds.cor <- function(score = NULL, expression = NULL, gene.set = NULL, cor.method = 'spearman') {
    # INPUT:
    # score: the computed disease score 
    # expression: the expression matrix either RNA or methylation expression
    
    # OUTPUT:
    # ranked genes with highest correlation between the two:

    # perform checks that the traits are the same in score and gs_file:
    stopifnot(all(names(score) %in% names(gene.set)));
    stopifnot(length(score) == length(gene.set));

    # initiate:
    result.collection <- vector('list', length = length(score));
    names(result.collection) <- names(score);

    # to limit the correlation test space, we would only test the correlation between
    # gene expression and score for genes that are present in the gs file for each disease
    for (disease in names(score)) {
        # find the set of genes that we will do correlation for:
        gene.interest <- intersect(names(gene.set[[disease]]), colnames(expression));
        
        # perform checks that the cell id are the same order between expression and score:
        stopifnot(rownames(score[[disease]]) == rownames(expression));

        # for every genes in the gene interest list, compute correlation:
        score.correlation <- sapply(
            gene.interest,
            FUN = function(gene) cor(expression[, gene], score[[disease]]$zscore, method = cor.method)
            )
        
        result.collection[[disease]] <- score.correlation;
    }
    return(result.collection);
}

# fake.score <- NULL
# fake.score[['diseaseA']] = data.frame(zscore = rnorm(1000), id = paste0('fake', seq(1, 1000)), row.names = 2)
# fake.score[['diseaseB']] = data.frame(zscore = rnorm(1000), id = paste0('fake', seq(1, 1000)), row.names = 2)
# fake.score[['diseaseC']] = data.frame(zscore = rnorm(1000), id = paste0('fake', seq(1, 1000)), row.names = 2)

# fake.expr <- data.frame(
#     geneA = fake.score[['diseaseA']],
#     geneB = fake.score[['diseaseB']],
#     geneC = rnorm(1000)
#     )
# colnames(fake.expr) <- c('geneA', 'geneB', 'geneC')
# rownames(fake.expr) <- paste0('fake', seq(1, 1000))

# fake.gene.set <- NULL;
# fake.gene.set[['diseaseA']] = c(1)
# names(fake.gene.set[['diseaseA']]) = 'geneA'
# fake.gene.set[['diseaseB']] = c(1)
# names(fake.gene.set[['diseaseB']]) = 'geneB'
# fake.gene.set[['diseaseC']] = c(1)
# names(fake.gene.set[['diseaseC']]) = 'geneC'

# test.case <- expr.ds.cor(fake.score, fake.expr, fake.gene.set)

# # return: 
# $diseaseA
# geneA 
#     1 

# $diseaseB
# geneB 
#     1 

# $diseaseC
#      geneC 
# -0.0392494