### quantile-quantile-plot-visualization.R ########################################################
# purpose: make visualization script that look at the theoratical p value's distribution to the 
# observed p values distribution. Theoratical p value is presumed to be following normal distribution

### PREAMBLE ######################################################################################
# specify different directories:
variance.directory <- '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-75percentile/high_variance/result/'
random.directory <- '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-95percentile/random/result/'
abundance.directory <- '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-75percentile/hypomethylation/result/'

# bundle them together:
scDRS.directories <- c(variance.directory, random.directory, abundance.directory);
names(scDRS.directories) <- c('variance', 'random', 'abundance');

# load libraries:
require(ggplot2);
system.date <- Sys.Date();

# read data:
for (simulation.mode in names(scDRS.directories)) {
    scDRS.directory <- scDRS.directories[simulation.mode];
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
        # risk.score[[result]]$pval = runif(nrow(risk.score[[1]]), min = 1/(nrow(risk.score[[1]])));
        # risk.score[[result]]$pval = seq(1/10000, 1, by= 0.0001)
        }

    # simplify list names:
    list.names <- gsub(scDRS.directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;

    ### GET PLOTTING MATRIX ###########################################################################
    # split the data by groups of 100 genes, 500 genes, and 1000 genes:
    one.hundred.set <- grep('NULL_SIM_.*_100_genes_seed_[0-9]', names(risk.score));
    one.hundred.set <- risk.score[one.hundred.set];
    five.hundred.set <- grep('NULL_SIM_.*_500_genes_seed_[0-9]', names(risk.score));
    five.hundred.set <- risk.score[five.hundred.set];
    one.thousand.set <- grep('NULL_SIM_.*_1000_genes_seed_[0-9]', names(risk.score))
    one.thousand.set <- risk.score[one.thousand.set];
    all.set <- list(one.hundred.set, five.hundred.set, one.thousand.set);
    names(all.set) <- c('one_hundred', 'five_hundred', 'one_thousand');

    # initiate the ranked.p.distribution:
    ranked.p.distribution <- vector('list', length = length(all.set));
    names(ranked.p.distribution) <- names(all.set);
    for (gene.number in names(all.set)){
        ranked.p.distribution[[gene.number]] <- vector('list', length = nrow(risk.score[[1]]));
        names(ranked.p.distribution[[gene.number]]) <- paste0('rank', seq(1, nrow(risk.score[[1]])));

        for (rank.cell in names(ranked.p.distribution[[gene.number]])){
            ranked.p.distribution[[gene.number]][[rank.cell]] <- vector('list', length = length(all.set[[gene.number]]));
            names(ranked.p.distribution[[gene.number]][[rank.cell]]) <- names(all.set[[gene.number]]);
        }
    }

    # now we will fill in the nested list:
    for (gene.number in names(all.set)){
        plotting.set <- all.set[[gene.number]];
        ordered.simulation <- vector('list', length = length(plotting.set));
        names(ordered.simulation) <- names(plotting.set);

        for (simulation in names(plotting.set)) {

            # first we will need to order the cells by p values in each seed for all sets:
            ordering.index <- order(plotting.set[[simulation]]$pval);
            ordered.simulation[[simulation]] <- plotting.set[[simulation]][ordering.index, ];
            
            # next we will grab out the cell ranked cell from the data:
            for (rank.cell in seq(1, nrow(ordered.simulation[[simulation]]))){
                ranked.p.distribution[[gene.number]][[rank.cell]][[simulation]] <- ordered.simulation[[simulation]][rank.cell, 'pval']
                }
            }
        }

    # lets collapse the innest layer by summarization:
    for (gene.number in names(ranked.p.distribution)) {
        for (rank.cell in names(ranked.p.distribution[[gene.number]])) {
            simulation.results <- Reduce('c', ranked.p.distribution[[gene.number]][[rank.cell]])
            plot.df <- data.frame(
                feature_number = gene.number,
                rank = rank.cell,
                mean = -log10(mean((simulation.results))),
                sd = sd(simulation.results),
                theoratical = -log10(as.numeric(gsub('rank', '', rank.cell)) / (length(ranked.p.distribution[[gene.number]]) + 1))
                );

            ranked.p.distribution[[gene.number]][[rank.cell]] <- plot.df;
            }
        ranked.p.distribution[[gene.number]] <- Reduce('rbind', ranked.p.distribution[[gene.number]]);
        }
    plot.df <- Reduce('rbind', ranked.p.distribution);

    ### VISUALIZATION #################################################################################
    # visaulize the theoratical vs real p value without the legend first
    output.path <- paste0(
        '/u/home/l/lixinzhe/project-geschwind/plot/',
        system.date,
        '-',
        simulation.mode,
        '-75quantile-pvalue-pvalue-plot.png'
        );
    gplot <- ggplot(
        plot.df,
        aes(x = theoratical, y = mean, group = feature_number, color = feature_number)) +
        geom_pointrange(aes(ymin = mean - 2 * sd, ymax= mean + 2 * sd)) + 
        theme_classic() +
        ylab(expression(-log[10]("observed p value"))) +
        xlab(expression(-log[10]("theoratical p value"))) +
        geom_abline(slope = 1, intercept = 0) +
        scale_color_manual(
            values = c("one_hundred" = "#fbb4ae", "five_hundred" = "#b3cde3", "one_thousand" = "#ccebc5"),
            labels = c("100 genes", "500 genes", "1000 genes")
            ) +
        labs(color = '# permuted genes') +
        theme(text = element_text(size = 20)) +
        theme(legend.position = "none")

    png(
        filename = output.path,
        width = 5,
        height = 5,
        units = 'in',
        res = 400
        );
    print(gplot);
    dev.off();

    # make a duplicated plot but with the legend
    output.path <- paste0(
        '/u/home/l/lixinzhe/project-geschwind/plot/',
        system.date,
        '-',
        simulation.mode,
        '-75quantile-pvalue-pvalue-plot-with-legend.png'
        );
    gplot <- ggplot(
        plot.df,
        aes(x = theoratical, y = mean, group = feature_number, color = feature_number)) +
        geom_pointrange(aes(ymin = mean - 2 * sd, ymax= mean + 2 * sd)) + 
        theme_classic() +
        ylab(expression(-log[10]("observed p value"))) +
        xlab(expression(-log[10]("theoratical p value"))) +
        geom_abline(slope = 1, intercept = 0) +
        scale_color_manual(
            values = c("one_hundred" = "#fbb4ae", "five_hundred" = "#b3cde3", "one_thousand" = "#ccebc5"),
            labels = c("100 genes", "500 genes", "1000 genes")
            ) +
        labs(color = '# permuted genes') +
        theme(text = element_text(size = 20))

    png(
        filename = output.path,
        width = 5,
        height = 5,
        units = 'in',
        res = 400
        );
    print(gplot);
    dev.off();
    cat('the maximum standard deviation in calibration simulation is:', max(plot.df$sd), 'for simulation mode:', simulation.mode, '\n')
    }

#the maximum standard deviation in calibration simulation is: 0.003677928 for simulation mode: variance 
#the maximum standard deviation in calibration simulation is: 0.003895855 for simulation mode: random 
#the maximum standard deviation in calibration simulation is: 0.008768931 for simulation mode: abundance 