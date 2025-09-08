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
######          divergent and convergent on cell types with brain traits             ######
###########################################################################################
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

# for each of the trait, get fdr corrected p values:
trait.info$number_significant = NA
for (result in names(risk.score)) {
    # get the adj_pval:
    adj_pval = p.adjust(risk.score[[result]][, 'pval'], method = 'fdr')
    risk.score[[result]]$adj_pval = adj_pval
    
    # identify the adjusted pval;
    significant_cell = rownames(risk.score[[result]])[risk.score[[result]]$adj_pval < 0.1]
    trait.info[match(result, trait.info$Trait_Identifier), 'number_significant'] = length(significant_cell)
    }

# summarize by category:
trait.info %>% group_by(Category) %>% summarize(median = median(number_significant), sd = sd(number_significant))

###########################################################################################
######                                     Region stuff                              ######
###########################################################################################
script="/u/home/l/lixinzhe/project-github/met-scDRS/publication-visualization/revision/main-figures/common-axis-cell-type-region-predictiveness-heatmap.R"

###########################################################################################
######                                    GO analysis                                ######
###########################################################################################
library(clusterProfiler);
source('/u/home/l/lixinzhe/project-github/met-scDRS/spell-book/go_pathway_counter.R')

# read in the correlation results:
score.mch.cor = readRDS('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/score_mch_correlation.rds')
print('for MDD, here is the most correlated genes methylation to met-scDRS')
tail(sort(score.mch.cor[[1]]))

# load in the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

# read in the go terms:
go_results = list.files('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', pattern = '*_readable_go_results.rds')
go_terms = NULL
trait.info$sig_pathway = NA
rownames(trait.info) = trait.info$Trait_Identifier

for (result in go_results){
    file_path = paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', result)
    # load in the GO result:
    trait = gsub('_readable_go_results.rds', '', result)
    go_terms[[trait]] = readRDS(file_path)

    # count the number of significant terms in each trait:
    trait.info[trait, 'sig_pathway'] = sum(go_terms[[trait]]@result$p.adjust < 0.05)

}
print(trait.info %>% group_by(Category) %>% summarise(average = mean(sig_pathway)))

pathway_count_mat = pathway_count(go_terms, trait.info)

# common pathways among brain traits:
common_paths = tail(names(sort(colSums(pathway_count_mat[[1]]))))
print(tail(sort(colSums(pathway_count_mat[[1]]))))

# grab out the more specific pathways:
specific = colnames(pathway_count_mat[[1]])[4 == colSums(pathway_count_mat[[1]])]
print('more specific pathways and diseases')
print(pathway_count_mat[[1]][, specific])

print('for postsynaptic specialization:')
print(pathway_count_mat[[1]][,'GOCC_POSTSYNAPTIC_SPECIALIZATION'])

print('central nervous system developmental pathway:')
print(pathway_count_mat[[1]][,'GOBP_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT'])


bipolar_significant_pathways = colnames(pathway_count_mat[[1]])[pathway_count_mat[[1]]['PASS_BIP_Mullins2021',] > 0]
cat(paste0('significant pathways in bipolar: ', bipolar_significant_pathways), '\n')

# grab out the ion channel related pathways:
print('channel related:')
print(pathway_count_mat[[1]][,grep('CHANNEL', colnames(pathway_count_mat[[1]]))])

# output the pathway count matrix for both brain and non brain as our results:
write.table(pathway_count_mat[[1]], sep = '\t', file = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-brain-trait-signficant-pathways.txt'))
write.table(pathway_count_mat[[2]], sep = '\t', file = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-non-brain-trait-signficant-pathways.txt'))


###########################################################################################
######                              robustness analysis                              ######
###########################################################################################
# load in the null tests:
p_matrix = read.table('/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse215353_full/arcsine-null/PASS_MDD_Howard2019_calibration_p_metric.txt', sep = '\t', header = T)
bad_nulls_num = sum(p_matrix$ctrl_norm_adjusted < 0.1)
cat(bad_nulls_num, 'nulls are significantly not normal in MDD \n')

# load in the mouse MDD result:
mouse_mdd_scdrs = read.table(
    file = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges132489_full/mean_var_length_arcsine/PASS_MDD_Howard2019.score.gz',
    sep = '\t',
    header = TRUE,
    row.names = 1
    )
mouse_mdd_scdrs$fdr = p.adjust(mouse_mdd_scdrs$pval, method = 'fdr')
mouse_meta = read.table(
    file = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/processed-full/all_meta_with_rowSum.csv",
    sep = ',',
    row.names = 1,
    header = TRUE
    )

# look at the number of significant cells by cell type:
significant_cells = rownames(mouse_mdd_scdrs)[mouse_mdd_scdrs$fdr < 0.1]
table(mouse_meta[significant_cells, 'MajorType']) / length(significant_cells) * 100

### null-significance.R ###########################################################################
arcsine_dir = "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/arcsine/null_distribution/"
no_norm_dir = "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/untransformed/null_distribution/"
logit_dir = "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/logit/null_distribution/"
library_size = "/u/home/l/lixinzhe/project-geschwind/plot/met_scdrs_revision/gse132489_30K/normalization/library/null_distribution/"

# directories:
dirs = c(arcsine_dir, no_norm_dir, logit_dir, library_size)
names(dirs) = c('arcsine', 'unnormalized', 'logit', 'library')

# for each of the directory, grab out the p values:
score_collection = NULL
for (normalization in names(dirs)){
    trait_files = list.files(dirs[normalization], pattern = '.txt')
    
    # for each of the traits, load in the summary:
    for(trait in trait_files){
        trait_name = gsub('_calibration_p_metric.txt', '', trait)
        score_collection[[normalization]][[trait_name]] = read.table(file = paste0(dirs[normalization], trait), sep = '\t', header = TRUE, row.names = 1)
        }
    }

summary_fx = function(trait_list){
    collector = rep(NA, length = length(trait_list))
    names(collector) = names(trait_list)
    for (trait in names(trait_list)){
        summary_matrix = trait_list[[trait]]
        significant_num = sum(summary_matrix[, 'ctrl_norm_adjusted'] < 0.05)
        collector[trait] = significant_num
    }
    return(collector)
}

# for each of the normalization schemes, get the summary:
num_sig_across_normalization = lapply(score_collection, FUN=function(x) summary_fx(x))

# for each of the normalization schemes, find the average number:
print('average number of significant cells for normalizations across 75 traits')
print(lapply(num_sig_across_normalization, mean))

print('standard deviation for number of significant cells for normalization across 75 traits')
print(lapply(num_sig_across_normalization, sd))

