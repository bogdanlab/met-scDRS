### go_pathway_counter.R ##########################################################################
pathway_count = function(go.result, trait.info){
    # go.result should be a list that contains all readable clusterprofiler output
    # trait info should have trait identifier and trait category as columns
    sig.pathways <- NULL
    for (disease in names(go.result)){
        result <- go.result[[disease]]@result
        sig.pathways[[disease]] <- rownames(result)[result$p.adjust < 0.05] 
    }
    
    brain.trait <- trait.info$Trait_Identifier[trait.info$Category == 'brain'];
    nonbrain.trait <- trait.info$Trait_Identifier[trait.info$Category != 'brain'];

    # find the set of union pathways that are significant in all brain traits:
    union.pathways <- Reduce('union', sig.pathways)

    # initiate result matrix:
    pathways.counter <- matrix(0, nrow = length(brain.trait), ncol = length(union.pathways))
    rownames(pathways.counter) <- brain.trait
    colnames(pathways.counter) <- union.pathways

    # also initiate result matrix for non brain traits;
    nonbrain.counter <- matrix(0, nrow = length(nonbrain.trait), ncol = length(union.pathways))
    rownames(nonbrain.counter) <- nonbrain.trait
    colnames(nonbrain.counter) <- union.pathways

    pathways.pval <- matrix(NA, nrow = length(names(go.result)), ncol = length(union.pathways))
    rownames(pathways.pval) <-  names(go.result)
    colnames(pathways.pval) <- union.pathways

    for (disease in brain.trait) {
        pathways.counter[disease, sig.pathways[[disease]]] <- 1
    }

    for (disease in nonbrain.trait) {
        nonbrain.counter[disease, sig.pathways[[disease]]] <- 1
    }

    cat('average number of significant pathways in non brain traits:', mean(rowSums(nonbrain.counter)), '\n')
    return(list(pathways.counter, nonbrain.counter))
}