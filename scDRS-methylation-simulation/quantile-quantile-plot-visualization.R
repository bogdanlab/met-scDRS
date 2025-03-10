### quantile-quantile-plot-visualization.R ########################################################
# purpose: make visualization script that look at the theoratical p value's distribution to the 
# observed p values distribution. Theoratical p value is presumed to be following normal distribution

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    significant-cells-visualization-script.R [--dir <scdrs> --meta_data <meta> --field <group> --out <output> --cutoff <p> --plot_type <count>]

Options:
    --dir directory path to scDRS score file (first column = rownames)
    --out path to output file
]' -> doc

# collect user input: 
opts <- docopt(doc)
scDRS.directory <- opts$dir;
output.path <- opts$out;

# function testing:
# scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/random-simulation-output/TOP_VARIANCE/"
# output.path <- '/u/home/l/lixinzhe/project-geschwind/plot/10-05-2023-random-simulation-qq-plot.png'

# load libraries:
require(ggplot2);

# read data:
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
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

### GET PLOTTING MATRIX ###########################################################################
# split the data by groups of 100 genes, 500 genes, and 1000 genes:
one.hundred.set <- grep('NULL_SIM_TOP_VARIANCE_100_genes_seed_[0-9]', names(risk.score));
one.hundred.set <- risk.score[one.hundred.set];
five.hundred.set <- grep('NULL_SIM_TOP_VARIANCE_500_genes_seed_[0-9]', names(risk.score));
five.hundred.set <- risk.score[five.hundred.set];
one.thousand.set <- grep('NULL_SIM_TOP_VARIANCE_1000_genes_seed_[0-9]', names(risk.score))
one.thousand.set <- risk.score[one.thousand.set];
all.set <- list(one.hundred.set, five.hundred.set, one.thousand.set);
names(all.set) <- c('one_hundred', 'five_hundred', 'one_thousand');

# initiate place holder:
theoratical.collection <- observed.collection <- vector('list', length = length(all.set));
names(theoratical.collection) <- names(observed.collection) <- names(all.set);
plot.df.collection <- observed.collection;

# now we will fill in the nested list:
for (gene.number in names(all.set)){
    plotting.set <- all.set[[gene.number]];
    
    # initiate place holder:
    ordered.simulation <- vector('list', length = length(plotting.set));
    names(ordered.simulation) <- names(plotting.set);
    theoratical.collection[[gene.number]] <- observed.collection[[gene.number]] <- ordered.simulation;
    
    # loop over the list ot give us the theoratical and observerd quantile
    for (simulation in names(plotting.set)) {

        # first we will need to order the cells by p values in each seed for all sets:
        ordering.index <- order(plotting.set[[simulation]]$pval);
        ordered.simulation[[simulation]] <- plotting.set[[simulation]][ordering.index, ];
        ordered.simulation[[simulation]]$theoratical <- seq(1/nrow(ordered.simulation[[simulation]]), 1, length = nrow(ordered.simulation[[simulation]]))

        # next we can compute the quantile for each of the data:
        theoratical.quantile <- quantile(
            x = -log10(ordered.simulation[[simulation]]$theoratical),
            probs = seq(0, 1, length.out = 25)
            )
        observed.quantile <- quantile(
            x = -log10(ordered.simulation[[simulation]]$pval),
            probs = seq(0, 1, length.out = 25)
            )

        # save the quantile information:
        theoratical.collection[[gene.number]][[simulation]] <- theoratical.quantile;
        observed.collection[[gene.number]][[simulation]] <- observed.quantile;
        }
    
    # concatenate the simulations:
    theoratical.collection[[gene.number]] <- Reduce('rbind', theoratical.collection[[gene.number]]);
    observed.collection[[gene.number]] <- Reduce('rbind', observed.collection[[gene.number]]);

    # calculate the mean and sd for plotting:
    quantile.mean <- colMeans(observed.collection[[gene.number]]);
    quantile.sd <- apply(observed.collection[[gene.number]], 2, sd);
    plot.df <- data.frame(cbind(quantile.mean, quantile.sd));
    plot.df$theoratical <- theoratical.quantile;
    plot.df$gene.number <- gene.number;
    plot.df.collection[[gene.number]] <- plot.df;
    }
plot.df <- Reduce('rbind', plot.df.collection);

### VISUALIZATION #################################################################################
gplot <- ggplot(
    plot.df,
    aes(x = theoratical, y = quantile.mean, group = gene.number, color = gene.number)) +
    geom_pointrange(aes(ymin = quantile.mean - quantile.sd, ymax= quantile.mean + quantile.sd)) + 
    theme_classic() +
    xlab('theoratical -log10 P quantiles') +
    ylab('observed -log10 P quantiles') +
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
