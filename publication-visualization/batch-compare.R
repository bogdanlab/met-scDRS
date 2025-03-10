### batch-compare.R ###############################################################################
# purpose: compare the without batch correction and with batch correction on the identity of the cells:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(cowplot)

# load in the new propoortion:
scDRS.directory <- '/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-dev-with-cov/'
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;
system.date <- Sys.Date()
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'

# create progress bar:
cat('loading scores \n')
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(score.files),
    clear = FALSE
    );

for (result in score.files) {
    risk.score[[result]] <- fread(
        file = result,
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE
        );
    risk.score[[result]] <- data.frame(risk.score[[result]], row.names = 1)
    pb$tick()
    }
# simplify list names:
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(risk.score) <- list.names;

# load in the old propoortion:
scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/"
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
old.risk.score <- vector('list', length = length(score.files));
names(old.risk.score) <- score.files;

# create progress bar:
cat('loading scores \n')
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(old.risk.score),
    clear = FALSE
    );

for (result in score.files) {
    old.risk.score[[result]] <- fread(
        file = result,
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE
        );
    old.risk.score[[result]] <- data.frame(old.risk.score[[result]], row.names = 1)
    pb$tick()
    }

# simplify list names:
list.names <- gsub(scDRS.directory, '', score.files);
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

# rename the list names:
names(old.risk.score) <- list.names;

### PROCESS #######################################################################################
comparison_correlation <- rep(NA, length = length(risk.score))
names(comparison_correlation) <- names(risk.score)
for(disease in names(risk.score)) {
    new <- risk.score[[disease]]
    old <- old.risk.score[[disease]]
    comparison_correlation[disease] <- cor(new$zscore, old$zscore, method = 'spearman')
    }

# make a supplementary plot on the comparison for MDD:
plot_df <- data.frame(
    MDD_batch = old.risk.score[['PASS_MDD_Howard2019']]$zscore,
    MDD_correct = risk.score[['PASS_MDD_Howard2019']]$zscore
    )

gplot <- ggplot(plot_df, aes(x = MDD_batch, y = MDD_correct)) +
    geom_point() + 
    labs(
        title = 'met-scDRS pre and post batch comparison',
        x = "without batch correction",
        y = "with batch correction"
        ) + 
    theme_classic() +
    geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
    theme(
    plot.title = element_text(hjust = 0.5)  # Center align the title
    )

# Create a combined plot using cowplot
output.path <- paste0(output.dir, system.date, '-pre-post-batch-correction-plot.png')
png(
    filename = output.path,
    width = 5,
    height = 5,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();

# # create a box plot on the pre and post batch correction:
# gplot <- ggplot(data.frame(comparison_correlation), aes(x = "", y = comparison_correlation)) +
#     geom_violin(fill = 'cornflowerblue', alpha = 0.5) +
#     scale_y_continuous(limits = c(0, 1)) +
#     labs(
#         title = "pre and post batch correction correlation across disorders",
#         y = "spearman correlation coefficients"
#         ) +
#     theme_classic() +
#     theme(
#         plot.title = element_text(hjust = 0.5)  # Center align the title
#     )
# # Create a combined plot using cowplot
# output.path <- paste0(output.dir, system.date, '-pre-post-batch-correction-box-plot.png')
# png(
#     filename = output.path,
#     width = 5,
#     height = 5,
#     units = 'in',
#     res = 400
#     );
# print(gplot)
# dev.off();
