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
###  Met-scDRS identifies disease relevant cells in human methylation atlas for MDD     ###
###########################################################################################
# mdd result:
mdd_result = risk.score[['PASS_MDD_Howard2019']]
mdd_result$FDR = p.adjust(mdd_result$pval, method = 'fdr')
mdd_result$significant = mdd_result$FDR < 0.1

cat('We identified a total of', sum(mdd_result$significant), 'cells with significant met-scDRS \n')

# mdd result:
stopifnot(all(rownames(mdd_result) == rownames(meta)))

# summarize:
meta = cbind(mdd_result, meta)

# split groups:
meta %>% group_by(X_CellClass) %>% group_split() -> grouped_list
names(grouped_list) = sapply(grouped_list, FUN=function(x) unique(x$X_CellClass))

# for each group, identify gradient with respect to tissue region:
for (cell_class in names(grouped_list)){
    print(grouped_list[[cell_class]] %>% summarize(sum(significant)))
    percent_significant = as.numeric(grouped_list[[cell_class]] %>% summarize(sum(significant)) / sum(mdd_result$significant))
    cat(percent_significant * 100, '% of significant cells in cell class', cell_class, '\n')
}

# cell type:
meta %>% group_by(X_MajorType) %>% group_split() -> grouped_list
names(grouped_list) = sapply(grouped_list, FUN=function(x) unique(x$X_MajorType))

# initiate result_df
result_df = data.frame(matrix(NA, ncol=2, nrow = length(grouped_list)))
rownames(result_df) = names(grouped_list)
colnames(result_df) = c('num_sig', 'prop_sig')

for (cell_type in names(grouped_list)){
    result_df[cell_type, 'num_sig']= (grouped_list[[cell_type]] %>% summarize(sum(significant)))
    result_df[cell_type, 'prop_sig'] = as.numeric(grouped_list[[cell_type]] %>% summarize(sum(significant)) / nrow(grouped_list[[cell_type]]) ) * 100
}
result_df = result_df[order(result_df$num_sig),]
tail(result_df, 3)

write.table(
    result_df,
    file = '/u/home/l/lixinzhe/project-geschwind/port/met_scdrs_supp_table/supplementary_table_1.csv',
    sep = ','
    )

# add the cell class to the table as well:
result_df$cell_class = meta$X_CellClass[match(rownames(result_df), meta$X_MajorType)]

###########################################################################################
######                          load in the average mch                              ######
###########################################################################################
