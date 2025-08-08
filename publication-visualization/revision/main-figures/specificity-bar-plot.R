### specificity-bar-plot.R ########################################################################
# purpose: plot out bar chart that show case our results are brain specific on GSE215353

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)

# specify paths:
scDRS.directory <- '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/';
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv';
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date()

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# read trait info:
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

# load in data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# create progress bar:
cat('loading scores \n')
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(score.files),
    clear = FALSE
    );

# read into the empty list:
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

### DISTILL INFO ##################################################################################
# for each of the risk score, calculate the proportion of significant cells for each of the cell type
# for each of the result, find out the number of significant cells in each cell type:
group.index <- 'X_CellClass'
p.cutoff <- 0.1;
plot.type <- 'proportion'
cell.type <- unique(meta[, group.index]);
significance.matrix <- matrix(NA, nrow = length(risk.score), ncol = length(cell.type));
rownames(significance.matrix) <- names(risk.score);
colnames(significance.matrix) <- cell.type;

# create the traits by cell type matrix:
for (result in names(risk.score)) {
    # find the index for significant cells:
    significant.cell <- rownames(risk.score[[result]])[
        p.adjust(risk.score[[result]]$pval, method = 'fdr') < p.cutoff
        ];

    # for each cell type, check number of significant cells are part of that cell type:
    for (type in cell.type) {
        # locate cell id that belong in the cell type:
        cell.type.cell <- rownames(meta)[meta[, group.index] %in% type];

        # find the number of cells that are in each of the cell type category:
        if (plot.type == 'count') {
            significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell);
            } else {
                significance.matrix[result, type] <- sum(significant.cell %in% cell.type.cell) /
                    length(cell.type.cell);
            }
        }
    }

# make the significant proportion into a long format for plotting:
plot.df <- melt(significance.matrix);
colnames(plot.df) <- c('disease', 'cell_type', 'proportion');
plot.df$trait_class <- trait.info$Category[match(plot.df$disease, trait.info$Trait_Identifier)]

# summarize the result base on the trait class:
plot.df$proportion <- plot.df$proportion * 100
plot.df <- plot.df %>% group_by(cell_type, trait_class) %>% summarize(average = mean(proportion), sd = sd(proportion));
plot.df <- as.data.frame(plot.df)

### VISUALIZATION #################################################################################
# make visualization that give me cell type by trait:
gplot <- ggplot(
    plot.df,
    aes(x = trait_class, y = average, fill = cell_type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.2,
                    position=position_dodge(.9)) +
    theme_classic() +
    scale_fill_manual(values = c(
        "Non-neuronal Cells" = "#b3e2cd",
        "Inhibitory and Subcortical Neurons" = "#cbd5e8",
        'Excitatory Neurons' = '#fdcdac')
        ) +
    xlab('Traits categories') +
    ylab('Average % significant cells in cell class') +
    labs(fill = "cell class") +
    theme(legend.position = "none") +
    theme(text = element_text(size = 30))

output.path <- paste0(output.dir, system.date, '-scDRS-cell-class-by-disease-category-proportion.png')
png(
    filename = output.path,
    width = 14,
    height = 14,
    units = 'in',
    res = 400
    );
print(gplot)
dev.off();
