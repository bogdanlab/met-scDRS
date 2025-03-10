### GO-distribution.R #############################################################################
# purpose: plot out the distribution of gene ontology terms p values

### PREAMBLE ######################################################################################
# load in the different libraries:
# load libraries:
require(ggplot2)
require(tidyverse)
require(data.table)
require(ggvenn)
library(clusterProfiler);
library(enrichplot);

# load in the gene ontology result:
go.result <- readRDS("/u/project/pasaniuc/lixinzhe/R_saves/2024-07-08-75-traits-GO-result.rds")
trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
output.dir <- '/u/home/l/lixinzhe/project-geschwind/plot/'
system.date <- Sys.Date()
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv';

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# read trait info:
trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

### ANALYSIS: find the set of pathways that are commonly flagged by brain traits ##################
# extract the significant pathways:
sig.pathways <- NULL
for (disease in names(go.result)){
    result <- go.result[[disease]]@result
    sig.pathways[[disease]] <- rownames(result)[result$p.adjust < 0.05] 
}

# find out for brain trait, what is the common pathways:
# first find out waht traits are brain traits:
brain.trait <- trait.info$Trait_Identifier[trait.info$Category == 'brain'];
nonbrain.trait <- trait.info$Trait_Identifier[trait.info$Category != 'brain'];

# find the set of union pathways that are significant in all brain traits:
union.pathways <- Reduce('union', sig.pathways)

# initiate result matrix:
pathways.counter <- matrix(0, nrow = length(brain.trait), ncol = length(union.pathways))
rownames(pathways.counter) <- brain.trait
colnames(pathways.counter) <- union.pathways

# also initiate result matrix for non brain traits;
nonbrain.counter <- matrix(0, nrow = length(nonbrain.trait), ncol = length(union.pathways))
rownames(nonbrain.counter) <- nonbrain.trait
colnames(nonbrain.counter) <- union.pathways

pathways.pval <- matrix(NA, nrow = length(names(go.result)), ncol = length(union.pathways))
rownames(pathways.pval) <-  names(go.result)
colnames(pathways.pval) <- union.pathways

for (disease in brain.trait) {
    pathways.counter[disease, sig.pathways[[disease]]] <- 1
}

for (disease in nonbrain.trait) {
    nonbrain.counter[disease, sig.pathways[[disease]]] <- 1
}

# print the average number significant pathways in non brain traits that are significant:
cat('average number of significant pathways in non brain traits:', mean(rowSums(nonbrain.counter)), '\n')

for (disease in names(go.result)){
    pathways.pval[disease, sig.pathways[[disease]]] <- go.result[[disease]]@result[sig.pathways[[disease]], 'p.adjust']
    
    # visualize the gene ontology for MDD:
    dot.plot <- dotplot(
        go.result[[disease]],
        showCategory=10
        ) + 
        ggtitle(paste0("gene ontology for ", disease)) +
        theme_classic() +
        theme(plot.title = element_text(hjust=0.5)) +
        theme(text = element_text(size = 20))

    
    # plot out the umap:
    dot.path <- '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/dot/'
    plot.name <- paste0(dot.path, system.date, '-', disease, '-gene-ontology-enrichment-dotplot.png')
    png(
        filename = plot.name,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
        );
    print(dot.plot);
    dev.off();
}

# write the table:
write.table(
    pathways.counter, 
    file = paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, 'significant-pathways-counter.csv'),
    sep = ',',
    quote = FALSE
    )
