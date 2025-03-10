gene.ontology.caller <- function(x, lambda = 'lambda.min', background = NULL, terms, visualize = FALSE, ...) {
    # INPUT:
    # x: can be either glmnet model or character vector of hgnc symbol
    # lambda: the model to use if x is glmnet model
    # if vector of charactets, the lambda argument is ignored
    # terms: loaded gmt argument used in enricher TERM2ID argument
    # background: hgnc symbol for background set of enrichment analysis
    # visualize: turn on to return a visualizable result

    # OUTPUT:
    # enrichment.result: results from the enricher call
    
    require(glmnet);
    require(org.Hs.eg.db);

    # first process background set:
    # set up conditionals:
    if (FALSE == is.null(background)) {
        hgnc.background <- background
        }
    else if (is.null(background) & 'cv.glmnet' == class(x)) {
        cat('background unspecified, using genes in glmnet model \n');
        hgnc.background <- rownames(coef(x));
        }
    else if(is.null(background) & 'cv.glmnet' != class(x)) {
        stop('background specificaton required, terminating! \n')
        }

    # first check if the x is just a list of characters:
    if (is.character(x)) { 
        hgnc.foreground <- x;
        }
    
    else {
        # extract the non zero coefficients from the inputed glmnet object:
        hgnc.foreground <- get.nonzero.coef(model = x, lambda = lambda, ranked = TRUE);
        hgnc.background <- rownames(coef(x));
        }

    # convert the gene names from foreground into entrez id:
    entrez.foreground <- mapIds(org.Hs.eg.db, hgnc.foreground, 'ENTREZID', 'SYMBOL');
    entrez.foreground <- na.omit(entrez.foreground);

    # convert the gene names from background into entrez id:
    entrez.background <- mapIds(org.Hs.eg.db, hgnc.background, 'ENTREZID', 'SYMBOL');
    entrez.background <- na.omit(entrez.background);

    # perform check for redundancies:
    stopifnot(FALSE == duplicated(entrez.foreground));
    stopifnot(FALSE == duplicated(entrez.background));

    # compute gene ontology enrichment:
    ontology.enrichment <- enricher(
        gene = entrez.foreground,
        pAdjustMethod = "fdr",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        universe = entrez.background,
        TERM2GENE = terms,
        ...
        );
    if (visualize) {
        enrichment.result <- ontology.enrichment;
    } else {
        enrichment.result <- ontology.enrichment@result;
    }

    # return the result:
    return(enrichment.result);
    }
