### pca-reconstruction ############################################################################
# purpose: a function to reconstruct original data from pca:
# input:
# score: the matrix of the $x from prcomp function
# loading: the matrix of $rotation from prcomp function
# rescale: should the reconstruction be rescaled based on the mean and sd
# mean: mean used to scale the input pca
# sd: sd used to scale the unput pca

# output:
# reconstruction: the reconstructed matrix with n x p dimension

pca.reconstruction <- function(score, loading, rescale = FALSE, mean, sd) {
  reconstruction <- score %*% t(loading)
  if (rescale){
    for(feature in seq(1, ncol(reconstruction))) {
      reconstruction[, feature] <- reconstruction[, feature] * sd[feature] + mean[feature]
      }
  }
  return(reconstruction)
}