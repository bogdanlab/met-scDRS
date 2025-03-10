### sig-cells-cell-class-cell-type-breakdown.R ####################################################
# purpose: allow us to find the population structure of significant cells

### PREAMBLE ######################################################################################
# load in the libraries:
library(ggplot2)
library(webr)
library(dplyr)

# load in the data:
scDRS.directory <- '/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-full/'
disease <- 'PASS_MDD_Howard2019.score.gz'
disease.score <- data.table::fread(
    file = paste0(scDRS.directory, disease),
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE
    )
disease.score <- data.frame(disease.score, row.names = 1)

# load in meta:
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv';
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );
system.date <- Sys.Date()

### VISUALIZATION #################################################################################
# extract out the set of cells that are signifnicant:
disease.score$fdr <- p.adjust(disease.score$pval, method = 'fdr')
significant.cell <- rownames(disease.score)[disease.score$fdr < 0.1]

# prepare plot df:
significant.info <- data.frame(
    cell = significant.cell,
    class = meta[match(significant.cell, rownames(meta)), 'X_CellClass'],
    type = meta[match(significant.cell, rownames(meta)), 'X_MajorType'],
    frequency = 1
    )

plot.df = as.data.frame(significant.info %>% group_by(class, type) %>% summarise(n = sum(frequency)))
stopifnot(sum(plot.df$n) == length(significant.cell))
rownames(plot.df) <- paste0(plot.df$class, '-', plot.df$type)

# write the table:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-PASS_MDD_Howard2019-sig-cell-type-structure.txt')
write.table(
    plot.df,
    file = output.path,
    sep = ',',
    quote = FALSE,
)

# from the plot df, retain only the top 3 cell type of each class in number
top3.type <- as.data.frame(plot.df %>% group_by(class) %>% top_n(3, n))
rownames(top3.type) <- paste0(top3.type$class, '-', top3.type$type)

# set the cells that are not in the top 3 of each class as 'other'
small.types <- setdiff(rownames(plot.df), rownames(top3.type))
plot.df[small.types, 'type'] <- 'else'
plot.df = as.data.frame(plot.df %>% group_by(class, type) %>% summarise(n = sum(n)))
stopifnot(sum(plot.df$n) == length(significant.cell))

# write the table:
output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-PASS_MDD_Howard2019-sig-cell-type-structure-aggregate.txt')
write.table(
    plot.df,
    file = output.path,
    sep = ',',
    quote = FALSE,
)

# plotting transferred to Excel and powerpoint