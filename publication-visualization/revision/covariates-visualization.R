### covariate visualization.R #####################################################################
# purpose: validate our signal by making sure it is not driven by some covaraites:

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(cowplot)

# specify paths:
scDRS.directory <- "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/"
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data_with_rowSum_mch.csv';
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

# also load in the significant matrix:
proportion.matrix <- read.table(
    file = '/u/home/l/lixinzhe/project-geschwind/plot/2025-08-01-revision-GSE215353-production-fraction-mch-cell-type-significance-proportion.csv',
    sep = ',',
    row.names = 1,
    header = TRUE,
    check.names = FALSE
    )

### PROCESS #######################################################################################
# compute the correlation beteen mch expression and proportion of significance:
# on all cells:
summarized.mch <- meta %>% group_by(X_MajorType) %>% summarize(avg_mch = mean(rowSum_centered)) %>% data.frame(row.names = 1)
# on excitatory neurons:
cor.result <- matrix(NA, nrow = nrow(proportion.matrix), ncol = 2)
colnames(cor.result) <- c('rho', 'pval')
rownames(cor.result) <- rownames(proportion.matrix)
for (disease in rownames(proportion.matrix)){
    summarized.mch$proportion <- t(proportion.matrix[disease, rownames(summarized.mch)])
    cor.result[disease, 'pval']<- cor.test(summarized.mch$proportion, summarized.mch$avg_mch, method = 'spearman')$p.value
    cor.result[disease, 'rho'] <- cor.test(summarized.mch$proportion, summarized.mch$avg_mch, method = 'spearman')$estimate
    }

print(cor.result)
number_of_significant_global_correlation = sum(na.omit(p.adjust(cor.result[,'pval'], method = 'fdr')) < 0.05)
cat('number of significant correlation between mch and proportion of significant cell type: ', number_of_significant_global_correlation, '\n')

# visualize the number of significant cells in each of the covariates:
for (result in names(risk.score)) {
    risk.score[[result]]$FDR <- p.adjust(risk.score[[result]]$pval, method = 'fdr')
    risk.score[[result]]$disease <- result
    risk.score[[result]]$number_significant <- sum(risk.score[[result]]$FDR < 0.1)
    risk.score[[result]]$donor <- meta$donor_id[match(rownames(risk.score[[result]]), rownames(meta))]
    risk.score[[result]]$batch <- meta$batch[match(rownames(risk.score[[result]]), rownames(meta))]
    }

plot.df.list <- gplot.list <- NULL
for (result in names(risk.score)){
    plot.df.list[[result]] <- risk.score[[result]] %>%
    group_by(donor) %>%
    summarize(count_FDR = sum(FDR < 0.1, na.rm = TRUE)) %>%
    ungroup()
    }

for (result in names(plot.df.list)) {
    gplot.list[[result]] <- ggplot(plot.df.list[[result]], aes(x = donor, y = count_FDR, fill = donor)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        labs(
            title = gsub('PASS_', '', gsub('UKB_460K.','', result)),
            x = "donor",
            y = "cell count"
            ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create a combined plot using cowplot
output.path <- paste0(output.dir, system.date, '-met-scDRS-sig-count-by-donor.png')
png(
    filename = output.path,
    width = 20,
    height = 60,
    units = 'in',
    res = 400
    );
combined_plot <- plot_grid(plotlist = gplot.list, ncol = 4)
print(combined_plot)
dev.off();

### also do this for batch:
plot.df.list <- gplot.list <- NULL
for (result in names(risk.score)){
    plot.df.list[[result]] <- risk.score[[result]] %>%
    group_by(batch) %>%
    summarize(count_FDR = sum(FDR < 0.1, na.rm = TRUE)) %>%
    ungroup()
    }

for (result in names(plot.df.list)) {
    gplot.list[[result]] <- ggplot(plot.df.list[[result]], aes(x = batch, y = count_FDR, fill = as.factor(batch))) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        labs(
            title = gsub('PASS_', '', gsub('UKB_460K.','', result)),
            x = "batch",
            y = "cell count"
            ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create a combined plot using cowplot
output.path <- paste0(output.dir, system.date, '-met-scDRS-sig-count-by-batch.png')
png(
    filename = output.path,
    width = 20,
    height = 60,
    units = 'in',
    res = 400
    );
combined_plot <- plot_grid(plotlist = gplot.list, ncol = 4)
print(combined_plot)
dev.off();
