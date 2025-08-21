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

### ANALYSIS: GENE ONTOLOGY #######################################################################
#initiate place holder:
traits.interest <- names(score)
top.genes <- vector('list', length = length(traits.interest));
names(top.genes) <- traits.interest;

# get progress bar:
pb <- progress::progress_bar$new(
    format = "[:bar] (:current/:total)",
    total = length(traits.interest),
    clear = FALSE
    );

# for each disease call ontology caller:        
for (disease in traits.interest) {
    score.mch.cor <- expr.ds.cor(
        score = score[disease],
        expression = mch.fraction,
        gene.set = trait.gene.set[disease]
        );
    
    # extract top genes that are used for gene ontology:
    top.genes[[disease]] <- names(head(sort(score.mch.cor[[disease]], decreasing = TRUE), ontology.gene.num));
    genes = top.genes[[disease]];

    # get the background genes:
    bg.gene <- intersect(names(trait.gene.set[[disease]]), colnames(mch.fraction));

    # call GO enrichment analysis:
    gene.ontology.result <- gene.ontology.caller(
        x = genes,
        background = bg.gene,
        terms = gsea.c5,
        visualize = TRUE
        );
    readable.result <- setReadable(gene.ontology.result, 'org.Hs.eg.db', 'ENTREZID')
    
    # save the results:
    saveRDS(score.mch.cor, paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_score_mch_correlation.rds'))
    saveRDS(readable.result, paste0('/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/GO/', disease, '_readable_go_results.rds'))
    pb$tick()
    }

