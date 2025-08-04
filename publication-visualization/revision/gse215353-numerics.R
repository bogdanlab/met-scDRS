### gse215353-numerics.R ##########################################################################
# purpose: get the numbers for publication

### PREAMBLE ######################################################################################
# load libraries:
require(ggplot2);
require(ComplexHeatmap);
require(circlize);
library(dplyr)

# read data:
scDRS.directory <- "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/"
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
risk.score <- vector('list', length = length(score.files));
names(risk.score) <- score.files;

# read meta:
meta.data.path = '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv'
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

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

###########################################################################################
######                                         MDD                                   ######
###########################################################################################
# mdd result:
mdd_result = risk.score[['PASS_MDD_Howard2019']]
mdd_result$FDR = p.adjust(mdd_result$pval, method = 'fdr')
mdd_result$significant = mdd_result$FDR < 0.1

# mdd result:
stopifnot(all(rownames(mdd_result) == rownames(meta)))
meta$significant = mdd_result$significant

# summarize:
meta = cbind(mdd_result, meta)

# split groups:
meta %>% group_by(X_MajorType) %>% group_split() -> grouped_list
names(grouped_list) = sapply(grouped_list, FUN=function(x) unique(x$X_MajorType))

# for each group, identify gradient with respect to tissue region:
for (cell_type in names(grouped_list)){
    print(grouped_list[[cell_type]] %>% group_by(X_Region) %>% summarize(mean(norm_score)))

}



meta %>% group_by(X_MajorType) %>% summarize(prop= (mean(significant)))

for 