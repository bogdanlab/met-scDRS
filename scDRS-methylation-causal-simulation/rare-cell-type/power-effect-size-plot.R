### power effect size plot.R ######################################################################
# purpose: plot the effect size of the perturbation vs power from the causal simulation:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2);
library(dplyr);

# specify different directories:
cell.type.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/result/rare-cell-type/"

# also specify the perturbation record paths:
cell.type.perturbation <- '/u/scratch/l/lixinzhe/tmp-file/causal-simulation/rare-cell-type/'

# bundle the result directories:
result.directories <- c(cell.type.directory);
names(result.directories) <- c('cell_type');

# input the system date:
system.date <- Sys.Date();
significance.cutoff <- 0.1;

### DATA LOADING ##################################################################################
for(simulation.mode in 'cell_type') {
    scDRS.directory <- result.directories[simulation.mode];
    score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
    risk.score <- vector('list', length = length(score.files));
    names(risk.score) <- score.files;

    # read into the empty list:
    for (result in score.files) {
        risk.score[[result]] <- read.table(
            file = result,
            sep = '\t',
            header = TRUE,
            row.names = 1,
            stringsAsFactors = FALSE
            );
        }

    # simplify list names:
    list.names <- gsub(scDRS.directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.h5ad.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;

    # similarly, load in the the causal genes, causal cells:
    perturb.csv <- paste0(list.names, '.csv');
    causal.perturbation.path <- paste0(cell.type.perturbation, perturb.csv);
    perturbation <- vector('list', length = length(score.files));
    names(perturbation) <- causal.perturbation.path;

    # read into the empty perturbation list:
    for(perturbation.meta in causal.perturbation.path) {
        perturbation[[perturbation.meta]] <- read.table(
            file = perturbation.meta,
            sep = ',',
            header = TRUE,
            stringsAsFactors = FALSE
            );
        }

    # simplify list names:
    list.names <- gsub('.*/', '', causal.perturbation.path);
    list.names <- gsub('\\.csv', '', list.names);
    names(perturbation) <- list.names;

    # run a check to make sure the names are same between perturbation and score files:
    stopifnot(names(perturbation) == names(risk.score));

    # initiate place holder:
    power <- rep(NA, length(perturbation));
    names(power) <- names(perturbation);

    # for each of the perturbed, mark its overlap, seed, effect size:
    seed <- gsub('-cell-type.*', '', names(perturbation));
    cell.type <- gsub('seed-.*-', '', gsub('-down-sample.*', '', names(perturbation)));
    down.sample <- as.numeric(gsub('-proportion-.*', '', gsub('.*-down-sample-', '', names(perturbation))));
    proportion <- as.numeric(gsub('-effect-.*', '', gsub('.*proportion-', '', names(perturbation))));
    effect.size <- as.numeric(gsub('-overlap-.*', '', gsub('.*effect-', '', names(perturbation))));
    overlap <- as.numeric(gsub('-causal-.*', '', gsub('.*overlap-', '', names(perturbation))));
    
    # next, for each of the perturbation result, find the power:
    for(perturbation.meta in names(perturbation)) {
        # extract cells from the result:
        risk.score[[perturbation.meta]]$fdr <- p.adjust(risk.score[[perturbation.meta]]$pval, method = 'fdr');
        significant.cell <- rownames(risk.score[[perturbation.meta]])[risk.score[[perturbation.meta]]$fdr < significance.cutoff];

        # get perturbed cells:
        perturbed.cell <- na.omit(perturbation[[perturbation.meta]]$causal.cells);

        # calculate power:
        power[perturbation.meta] <- length(intersect(significant.cell, perturbed.cell)) / length(perturbed.cell);
        
        }
    
    # prepare the plotting dataframe:
    summary.data.frame <- data.frame(
        power = power,
        effect = effect.size,
        cell_type = cell.type,
        proportion = proportion,
        down_sample = down.sample,
        mode = paste0('fixed_', simulation.mode)
        )

    # aggregate the mean and sd from simmulation:
    plot.df <- summary.data.frame %>% group_by(down_sample, cell_type) %>% summarize(average = mean(power), sd = sd(power));
    plot.df <- as.data.frame(plot.df);

    # for data checking:
    # mean(summary.data.frame[summary.data.frame$effect == 1.001,'power']) # same as dplyr result
    # sd(summary.data.frame[summary.data.frame$effect == 1.001,'power']) # same as dplyr result

    ### VISUALIZE POWER ###########################################################################
    # visualize the power:
    output.path <- paste0(
        '/u/home/l/lixinzhe/project-geschwind/plot/',
        system.date,
        '-',
        simulation.mode,
        '-power-effect-plot.png'
        );
    gplot <- ggplot(
        plot.df,
        aes(x = down_sample, y = average, col = cell_type)) +
        geom_point(position = position_dodge(width = 0.05), size = 3) +
        geom_pointrange(aes(ymin = average - 2 * sd, ymax = average + 2 * sd), position = position_dodge(width = 0.05)) + 
        theme_classic() +
        ylab('power') +
        xlab('perturbation effect') +
        theme(text = element_text(size = 20))

    png(
        filename = output.path,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
        );
    print(gplot);
    dev.off();

}
