### significance-bar-plot.R #######################################################################
# purpose: plot the number of significance cells in each of the simulation method:

### PREAMBLE ######################################################################################
# load in libraries
library(ggplot2);

# load in data:
top.variance.call <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation-95percentile/random-simulation-output/TOP_VARIANCE/';
top.fraction.call <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation-95percentile/random-simulation-output/TOP_FRACTION/';
random.call <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/null-simulation/random-simulation-output/RANDOM/';
system.date <- Sys.Date();

# read data:
result.collection <- NULL;
for (directory in c(top.variance.call, top.fraction.call, random.call)) {
    score.files <- list.files(directory, pattern = '\\.score.gz', full.names = TRUE);
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
    list.names <- gsub(directory, '', score.files);
    list.names <- gsub('/', '', list.names);
    list.names <- gsub('\\.score.gz', '', list.names);

    # rename the list names:
    names(risk.score) <- list.names;
    result.collection[[directory]] <- risk.score;
    }
names(result.collection) <- c('top.variance', 'top.fraction', 'random');

### SUMMARIZE DATA ################################################################################
# initiate a plotting dataframe:
summarized.data <- NULL;
for (permute.method in names(result.collection)){ 
    num.genes.seed <- names(result.collection[[permute.method]]);
    one.hundred.genes <- num.genes.seed[grep('_100_gene', num.genes.seed)];
    five.hundred.genes <- num.genes.seed[grep('_500_gene', num.genes.seed)];
    one.thousand.genes <- num.genes.seed[grep('_1000_gene', num.genes.seed)];
    genes.number.permuted <- list(one.hundred.genes, five.hundred.genes, one.thousand.genes);
    names(genes.number.permuted) <- c('one.hundred', 'five.hundred', 'one.thousand')

    # calculate the number of cells that exceed significance cutoff:
    for(num.genes in names(genes.number.permuted)) {
        permuted.sets <- genes.number.permuted[[num.genes]];
        significance.across.seed <- Reduce(
            'c',
            lapply(result.collection[[permute.method]][permuted.sets], FUN = function(x) 
                sum(p.adjust(x$pval, method = 'fdr') < 0.1)
                )
            );
        significance.average <- mean(significance.across.seed);
        significance.se <- sd(significance.across.seed);

        # put into plotting dataframe:
        summarized.data[[permute.method]][[num.genes]] <- data.frame(
            mode = permute.method,
            significance_average = significance.average,
            significance_se = significance.se,
            num_genes = num.genes
            );        
        }
    summarized.data[[permute.method]] <- Reduce('rbind', summarized.data[[permute.method]]);
    }
plot.df <- Reduce('rbind', summarized.data);

### VISUALIZE DATA ################################################################################
# plot the data as a bar chart:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, 'random-simulation-bar-plot-95quantile.png');

# ready the ggplot:
gplot <- ggplot(plot.df, aes(x = mode, y = significance_average, fill = num_genes)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(
        aes(
            ymin = significance_average - significance_se,
            ymax = significance_average + significance_se
            ),
        width=0.2,
        position = position_dodge(0.9)
        )  +
    theme_classic() +
    xlab('average number of significant cells') +
    ylab('simulation method') +
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
