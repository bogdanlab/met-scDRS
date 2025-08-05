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

# specify the data path:
mch.fraction.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/met_scdrs_processed-mch-v1_1_1_rc1.csv"

# define the gene set path:
gs.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/';

# read in the gmt file downloaded from the GSEA website for up to data enrichment estimation:
gsea.c5 <- read.gmt('/u/project/geschwind/lixinzhe/data/c5.all.v2023.1.Hs.entrez.gmt');

# load in the methylation fraction:
# load in methylation fractions:
mch.fraction <- fread(
    file = mch.fraction.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    data.table = FALSE,
    nThread = 1
    );
rownames(mch.fraction) <- mch.fraction$cell

meta.data.path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv'

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# cehck if all the cels in the data is in the meta:
stopifnot(rownames(mch.fraction) %in% rownames(meta))

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

# compute the expression and score correlation:
score.mch.cor <- expr.ds.cor(
    score = score['PASS_MDD_Howard2019'],
    expression = mch.fraction,
    gene.set = trait.gene.set['PASS_MDD_Howard2019']
    );

### ANALYSIS: GENE ONTOLOGY #######################################################################
traits.interest <- names(score.mch.cor)
top.genes <- vector('list', length = length(traits.interest));
names(top.genes) <- traits.interest;
gene.ontology.result <- top.genes;

for (disease in traits.interest) {
    # extract top genes that are used for gene ontology:
    top.genes[[disease]] <- names(head(sort(score.mch.cor[[disease]], decreasing = TRUE), ontology.gene.num));
    genes = top.genes[[disease]];

    # get the background genes:
    bg.gene <- intersect(names(trait.gene.set[[disease]]), colnames(mch.fraction));

    # call GO enrichment analysis:
    gene.ontology.result[[disease]] <- gene.ontology.caller(
        x = genes,
        background = bg.gene,
        terms = gsea.c5,
        visualize = TRUE
        );

    # visualize the gene ontology for MDD:
    dot.plot <- dotplot(
        gene.ontology.result[[disease]],
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

network.traits <- c(
    'PASS_MDD_Howard2019'
    );
for (disease in network.traits){
    # also visualize the network plot:
    cor.gene <- head(sort(score.mch.cor[[disease]], decreasing = TRUE), ontology.gene.num);
    readable.result <- setReadable(gene.ontology.result[[disease]], 'org.Hs.eg.db', 'ENTREZID')
    network.plot <- cnetplot(
        readable.result,
        categorySize = "pvalue",
        color.params = list(foldChange = cor.gene, edge = TRUE),
        circular = TRUE) + 
        ggtitle(disease) +
        theme_classic() +
        theme(text = element_text(size = 20)) +
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

saveRDS(score.mch.cor, '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/score_mch_correlation.rds')
saveRDS(readable.result, '/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/ontology_readable_result.rds')
