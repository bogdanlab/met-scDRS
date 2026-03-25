### GO-enrichment.R ###############################################################################
# purpose: compute the gene ontology enrichment for GSE215353 data

### PREAMBLE ######################################################################################
# load in the libraries:
library(Seurat);
library(sceasy);
library(ggplot2);
library(data.table);
library(ggpubr);
library(clusterProfiler);
library(enrichplot);
library(cowplot);
library(grid);

# define number of top selected genes for computing the go enrichments:
ontology.gene.num <- 100;

# define specified paths:
data.dir <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/';
session.save.path <- '/u/project/pasaniuc/lixinzhe/session-info/';
system.date <- Sys.Date();
save.path <- '/u/project/pasaniuc/lixinzhe/R_saves/';
code.path <- '/u/home/l/lixinzhe/project-github/methylation-RNA-xinzhe-rotation/code/';
plotting.path <- '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/';
function.path <- '/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/'

# define the gene set path:
gs.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/';

# read in the gmt file downloaded from the GSEA website for up to data enrichment estimation:
gsea.c5 <- read.gmt('/u/project/geschwind/lixinzhe/data/c5.all.v2023.1.Hs.entrez.gmt');

# load in the methylation fraction:
meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv'

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# load in the function:
source(paste0(function.path, 'gene-name-converter.R'));
source(paste0(function.path, 'gs-decomposer.R'));
source(paste0(function.path, 'score-loader.R'));
source(paste0(function.path, 'disease-score-expr-correlation.R'))
source(paste0(function.path, 'gene-ontology-caller.R'))

# load in the gene set from the 74 traits gs file:
trait.gs <- read.table(
    file = paste0(gs.dir, 'magma_10kb_top1000_zscore.75_traits.rv1.gs'),
    sep = '\t',
    header = TRUE
    );

# split gene sets:
trait.gene.set <- lapply(trait.gs$GENESET, gs.decomposer);
names(trait.gene.set) <- trait.gs$TRAIT;

# load in the scores:
score <- score.loader('/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/')

# load in the rds:
#initiate place holder:
traits.interest <- names(score)

for (disease in traits.interest){
    # read in the correlate and read in 
    score.mch.cor = readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_score_mch_correlation.rds'))
    readable.result = readRDS(paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_readable_go_results.rds'))
    
    dot.plot <- dotplot(
        readable.result,
        showCategory = 5
    ) + 
        ggtitle(gsub("PASS ", "", gsub("_", ' ', disease))) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5),
            text = element_text(size = 25)
        )

    # extract legend
    legend <- cowplot::get_legend(dot.plot)

    # save legend only
    legend.path <- '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/dot/'
    legend.name <- paste0(
        legend.path,
        system.date, '-', disease, '-gene-ontology-enrichment-dotplot-legend.png'
    )

    png(
        filename = legend.name,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
    )
    grid::grid.newpage()
    grid::grid.draw(legend)
    dev.off()

    # remove legend from main plot
    dot.plot.nolegend <- dot.plot + theme(legend.position = "none")

    # save plot only
    plot.name <- paste0(
        legend.path,
        system.date, '-', disease, '-gene-ontology-enrichment-dotplot.png'
    )

    png(
        filename = plot.name,
        width = 10,
        height = 5,
        units = 'in',
        res = 400
    )
    print(dot.plot.nolegend)
    dev.off()
    
    ### NETWORK PLOT ###
    # also visualize the network plot:
    cor.gene <- head(sort(score.mch.cor[[disease]], decreasing = TRUE), ontology.gene.num);
    network.plot <- cnetplot(
        readable.result,
        categorySize = "pvalue",
        color.params = list(foldChange = cor.gene, edge = TRUE),
        circular = TRUE) + 
        ggtitle(disease) +
        theme_classic() +
        theme(text = element_text(size = 25)) +
        theme(plot.title = element_text(hjust=0.5)) +
        theme(legend.text = element_text(size = 12))
        # theme(legend.position = "none")

    # grab out the legend:
    network.path <- '/u/scratch/l/lixinzhe/tmp-file/tmp-plot/network/'
    plot.name <- paste0(network.path, system.date, '-', disease, '-gene-ontology-enrichment-network-legend.png')
    png(
        filename = plot.name,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
        );

    legend <- cowplot::get_legend(network.plot);
    grid.newpage()
    grid.draw(legend)
    dev.off()

    # plot out the network only:
    plot.name <- paste0(network.path, system.date, '-', disease, '-gene-ontology-enrichment-network.png')
    network.plot <- network.plot + theme(legend.position = "none")
    png(
        filename = plot.name,
        width = 10,
        height = 10,
        units = 'in',
        res = 400
        );
    print(network.plot);
    dev.off();
}

# ###########################################################################################
# ######          visualize the number of significant pathways per trait               ######
# ###########################################################################################
# # load in the trait info:
# trait.info.path <- '/u/home/l/lixinzhe/project-geschwind/data/tait-classification.txt';
# trait.info <- read.table(file = trait.info.path, sep = '\t', header = TRUE);

# # read in the go terms:
# go_results = list.files('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', pattern = '*_readable_go_results.rds')
# go_terms = NULL
# trait.info$sig_pathway = NA
# rownames(trait.info) = trait.info$Trait_Identifier

# for (result in go_results){
#     file_path = paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', result)
#     # load in the GO result:
#     trait = gsub('_readable_go_results.rds', '', result)
#     go_terms[[trait]] = readRDS(file_path)

#     # count the number of significant terms in each trait:
#     trait.info[trait, 'sig_pathway'] = sum(go_terms[[trait]]@result$p.adjust < 0.05)

# }

# # get the bar chart:
# plot_df = trait.info %>% group_by(Category) %>% summarise(median = median(sig_pathway))

# ggplot(plot_df, aes(x = Category, y = median)) +
#   geom_col() +
#   theme_classic() +
#   labs(x = "Category", y = "Median") +
#   theme(
#     text = element_text(size = 16),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
  
# png(
#     filename = plot.name,
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 400
#     );
# print(network.plot);
# dev.off();