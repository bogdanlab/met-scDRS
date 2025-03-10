### gene-set-differences-magam-kangcheng.R ########################################################
# purpose: quantify the differences between the magma gene set that I have vs the gs file 
# kang cheng published

### PREAMBLE ######################################################################################
# load in Kangcheng's gs file:
kangcheng.gs.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/';
magma.gs.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/magma-out/'
plotting.path <- '/u/project/pasaniuc/lixinzhe/plot/';
system.date <- Sys.Date();
session.save.path <- '/u/project/pasaniuc/lixinzhe/session-info/';

# load library:
library('ggvenn');

# read in the 74 traits published gs:
kangcheng.gs <- read.table(
    file = paste0(kangcheng.gs.dir, 'magma_10kb_top1000_zscore.74_traits.rv1.gs'),
    sep = '\t',
    header = TRUE
    );

# read in the MDD magma gs:
mdd.gs <- read.table(
    file = paste0(magma.gs.dir, 'MDD-GWAS-grch37-munge-gs-output.gs'),
    sep = '\t',
    header = TRUE
    );

# make a function that decomposes the geneset line output:
gs.decomposer <- function(gs) {
    # split by , and get the genes:
    gs.split <- unlist(strsplit(gs, ','));
    gs.genes <- gsub(':.*', '', gs.split);
    gs.weight <- gsub('.*:', '', gs.split);
    gs.weight <- as.numeric(gs.weight);
    names(gs.weight) <- gs.genes;

    # return the gs.weight:
    return(gs.weight)
    }

### ANALYSIS: JACCARD INDEX #######################################################################
# grab out KC's MDD genes:
kc.mdd.gs <- kangcheng.gs$GENESET[kangcheng.gs$TRAIT == 'PASS_MDD_Howard2019'];
kc.mdd.genes <- gs.decomposer(kc.mdd.gs);

# decompose magma MDD gs genes:
magma.mdd.genes <- gs.decomposer(mdd.gs$GENESET);

# compute the jaccard index between my mdd genes and kangcheng's mdd genes:
common.mdd.genes <- intersect(
    names(kc.mdd.genes),
    names(magma.mdd.genes)
    );
union.mdd.genes <- union(
    names(kc.mdd.genes),
    names(magma.mdd.genes)
    );

mdd.jaccard.index <- length(common.mdd.genes) / length(union.mdd.genes) * 100; # 16.75%

# find the correlation between common weights:
spearman.rho <- cor(kc.mdd.genes[common.mdd.genes], magma.mdd.genes[common.mdd.genes], method = 'spearman');
print(spearman.rho);

### VISUALIZATION #################################################################################
plot.df <- list(
    publication = names(kc.mdd.genes),
    my_result = names(magma.mdd.genes)
    );
venn.diagram <- ggvenn(
    plot.df, 
    fill_color = c("#0073C2FF", "#CD534CFF"),
    stroke_size = 0.5,
    set_name_size = 4
    );

plot.name <- paste0(
    plotting.path,
    system.date,
    '-MDD-Howard-gs-similarity-venn-diagram.png'
    );

png(
    filename = plot.name,
    width = 20,
    height = 10,
    units = 'in',
    res = 400
    );
print(venn.diagram);
dev.off();

### SESSION INFO ##################################################################################
writeLines(
    text = capture.output(sessionInfo()),
    con = paste0(session.save.path, system.date, '-sessionInfo.txt')
    );