### gse215353_normality check #####################################################################
# purpose: check the first/second moment, Third and Fourth moment distributions

### PREAMBLE ######################################################################################
# load in the data:
library(data.table)
library(ggplot2)
library(ggrepel)

# get the trait info:
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);


# define a function to compute moments for us:
moments_compute = function(observation){
    moment1 = mean(observation)
    moment2 = mean((observation - moment1)^2)
    moment3 = mean((observation - moment1)^3) / (moment2^(3/2))
    moment4 = mean((observation - moment1)^4) / (moment2^2)
    
    collector = c(moment1, moment2, moment3, moment4)
    names(collector) = c('mean', 'variance', 'skewness', 'kurtosis')
    return(collector)
    }


nulls_dir = "/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/"
traits = list.files(nulls_dir, pattern = '.full_score.gz')

# load in the meta:
meta = fread('/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv', sep = ',', header = TRUE, data.table = FALSE)
rownames(meta) = meta[, 'cell']

# grab out the brain traits:
brain_traits = trait.info[trait.info$Category=='brain', 1]

trait_nulls = NULL
for (brain_trait in brain_traits){
    trait = traits[grep(brain_trait, traits)]
    trait_name = gsub('.full_score.gz', '', trait)
    met_scdrs = fread(file = paste0(nulls_dir, trait), sep = '\t', header = TRUE)
    
    # get the null's columns:
    null_ix = grep('ctrl_raw_score', colnames(met_scdrs))
    
    # compute the first four moments:
    moments = t(apply(met_scdrs[, ..null_ix], 1, moments_compute))
    # rename the rows:
    rownames(moments) <- met_scdrs[, cell]
    moments = data.frame(moments)
    
    # grab out the cell type label:
    moments$cell_type = meta[match(rownames(moments), rownames(meta)), '_MajorType']
    
    # Compute average location for labels
    labels_dt <- aggregate(
        cbind(mean, variance, skewness, kurtosis) ~ cell_type,
        data = moments, 
        FUN = mean
        )
    
    # plot the skewness vs kurtosis:
    p2 <- ggplot(moments, aes(x = skewness, y = kurtosis, color = cell_type)) +
    geom_point(size = 2, alpha = 0.5, show.legend = FALSE) +
    geom_text_repel(
        max.overlaps = Inf,
        data = labels_dt,
        aes(x = skewness, y = kurtosis, label = cell_type),  # <- explicitly map x,y
        color = "black", size = 4, inherit.aes = FALSE
    ) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5)) +
    theme(text = element_text(size = 20)) +
    labs(title = trait_name, x = "skewness", y = "kurtosis")

    plot_path = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', Sys.Date(), '-', trait_name, '-moments-of-nulls.png')
    png(
        filename = plot_path,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
        );
    print(p2);
    dev.off();

}


