### met-scDRS-score-visualization-region.R ########################################################
# purpose: visualize the met-scDRS score gradient with respect to the region gradient

### PREAMBLE ######################################################################################
# load in libraries:
library(ggplot2)
library(dplyr)
library(data.table)

# define parameters:
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date();

# load in data:
# load in meta data:
meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv"
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# load in the scDRS score:
scDRS.file <- "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/PASS_MDD_Howard2019.full_score.gz"
# load in the data
full.score <- fread(
    file = scDRS.file,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE,
    data.table = FALSE
    );
rownames(full.score) <- full.score$cell

###########################################################################################
######                        MDD significant by cell class                          ######
###########################################################################################
# first we wish to plot the middle and inferior temporal gyrus met-scDRS score for significant cell
full.score$fdr <- p.adjust(full.score$pval, method = 'fdr')
significant.cell <- rownames(full.score)[full.score$fdr < 0.1]
insignificant.cell <- setdiff(rownames(full.score), significant.cell)

# get the set of cells that we wish to plot:
all(rownames(meta) == rownames(full.score))
meta$score = full.score$zscore

# this result is for making the donut plot:
significant_cell_class = table(meta[significant.cell, 'X_CellClass'])

# drop last to remove grouping of X_MajorType byt keeping X_CellClass
meta[significant.cell, ] %>%
  group_by(X_CellClass, X_MajorType) %>%
  summarize(number_cells = n(), .groups = "drop_last") %>%
  arrange(X_CellClass, desc(number_cells)) %>%
  as.data.frame()
