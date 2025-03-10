### get-magma.stats.R #############################################################################
# PURPOSE: this is a part of the pipeline script of scDRS for generating z score statistics from 
# MAGMA outputs using custom GWAS

### PREAMBLE ######################################################################################
# gather arguments that were passed from command line:
# The first argument should be the output of the scDRS
# The second argument should be the meta data that contains the coloring label
# The third argument will be the column name for which the coloring label is
# The fourth argument should be the plot output:

# gather arguments from command line:
args <- commandArgs(trailingOnly = TRUE);

# read in the result table:
score.path <- args[1];
score <- read.table(
    file = score.path,
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
    );

# load in library:
require(ggplot2)

# get the label:
meta.data.path <- args[2];
if (meta.data.path != ''){
    meta.data <- read.table(
        file = meta.data.path,
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );
    }

# get the plot output path:
plot.name <- args[4];

### PLOT DISTRIBUTION #############################################################################
# prepare plotting dataframe:
color.index <- args[3];

if (color.index != ''){
    plot.df <- data.frame(
        p_value = score$pval,
        color_label = factor(meta.data[rownames(score), color.index])
        );
    
    gplot <- ggplot(plot.df, aes(x = p_value, color = color_label)) +
        geom_histogram(fill="white", alpha=0.5, position="identity") +
        theme_classic() +
        xlab('p value') +
        ylab('frequency') +
        theme(text = element_text(size = 20))

    } else {
        # grab out the plot df and make gplot:
        plot.df <- data.frame(p_value = score$pval);
        gplot <- ggplot(plot.df, aes(x = p_value)) + 
            geom_histogram(color="black", fill="white") +
            theme_classic() +
            xlab('p value') +
            ylab('frequency') +
            theme(text = element_text(size = 20))

        }

png(
    filename = plot.name,
    width = 10,
    height = 10,
    units = 'in',
    res = 400
    );
print(gplot);
dev.off();