### simulation-checker ############################################################################
# purpose: script to make sure that all of the perturbations have different sets of cell and genes

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2);
library(dplyr);

# also specify the perturbation record paths:
overlap.perturbation <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/fixed-overlap/'
effect.perturbation <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/fixed-effect/'

# bundle the result directories:
perturbation.paths <- c(overlap.perturbation, effect.perturbation);
names(perturbation.paths) <- c('overlap', 'effect');

# input the system date:
system.date <- Sys.Date();
significance.cutoff <- 0.1;

### DATA LOADING ##################################################################################
for (simulation.mode  in names(perturbation.paths)){
    # initiate empty perturbation list for each of the simulation mode:
    perturbation <- NULL;

    # load in the perturbation within that simulation mode:
    causal.directory <- perturbation.paths[simulation.mode];
    perturbation.files <- list.files(causal.directory, pattern = '\\.csv', full.names = TRUE);
    for(perturbation.meta in perturbation.files) {        
        perturbation[[perturbation.meta]] <- read.table(
            file = perturbation.meta,
            sep = ',',
            header = TRUE,
            stringsAsFactors = FALSE
            );

        }

    # ensure non of the perturbation is identical to each other:    
    for (first.instance in seq(1, length(perturbation))){
        for (second.instance in seq(1, length(perturbation))){
            if (first.instance != second.instance) {
                stopifnot(
                    FALSE == identical(
                        perturbation[[first.instance]],
                        perturbation[[second.instance]])
                    )
                }
            }
        }
    }
