### jaccard-computer.R ############################################################################
# purpose: compute jaccard index between features given a list of features:

# input:
# feature list: list of characters

# output:
# jaccard.index: the matrix of jaccard index with dimension feature * feature

jaccard.computer <- function(feature.list) {
    # create a jaccard matrix:
    jaccard.matrix <- matrix(NA, ncol = length(feature.list), nrow = length(feature.list));

    # for each features i,j: compute the jaccard index between the two sets of features:
    for (feature_i in seq(1, length(feature.list))) {
        for (feature_j in seq(1, length(feature.list))) {
            # find the length of the union set:
            union.set <- union(feature.list[[feature_i]], feature.list[[feature_j]])
            union.length <- length(union.set);

            # find the lenght of the intersection set:
            intersection.set <- intersect(feature.list[[feature_i]], feature.list[[feature_j]])
            intersection.length <- length(intersection.set);
            
            # fill in the jaccard matrix:
            jaccard.matrix[feature_i, feature_j] <- intersection.length / union.length;
            }
        }
    
    # return the jaccard matrix:
    return(jaccard.matrix);
    }
