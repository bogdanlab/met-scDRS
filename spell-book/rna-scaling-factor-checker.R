### rna-scaling-factor-checker.R ##################################################################
# purpose: find the normalization factor of the RNA sequencing data from Seurat:

rna.scaling.factor.checker <- function(rna) {
    # compute the scaling factor:
    computed.scaling.factor <- colSums(exp(rna) - 1);
    
    # summarize if all cells have the same scaling factor:
    scaling.factor.summary <- table(computed.scaling.factor);
    return(scaling.factor.summary);
}