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
