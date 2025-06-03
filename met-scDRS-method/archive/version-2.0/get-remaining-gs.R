### get-remaining-gs.R ############################################################################
# purpose: compare with the output folder, find the set of gene set that has not finished computing
# and then resubmit those jobs

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    get-remaininig-gs.R [--scDRS_dir <scDRS_dir> --gs_file <gs_file> --output_gs <output_gs>]

Options:
    --scDRS_dir current directories that houses the unfinished computed scores
    --gs_file gene set file that were originally used to compute the scores in the folder
    --output_gs path to file that contains the uncomputed gene sets
]' -> doc

# collect user input: 
opts <- docopt(doc)
scDRS.dir <- opts$scDRS_dir;
gs.path <- opts$gs_file;
output.gs <- opts$output_gs;

# load in deocmposer function:
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

# for function testing:
# scDRS.dir <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE132489-full/"
# gs.path <- "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs"
# output.gs <- "/u/scratch/l/lixinzhe/tmp-file/test/remaining-gs-file.gs"

### GET REMAINING GENESET #########################################################################
# get the files in the scDRS folder:
compiled <- list.files(scDRS.dir, pattern = '\\.score.gz');
finished <- gsub('\\.score.gz', '', compiled)

# load in the gs file:
trait.gs <- read.table(
    file = gs.path,
    sep = '\t',
    header = TRUE
    );

# get the gs:
rownames(trait.gs) <- trait.gs$TRAIT;
remaining.gs <- trait.gs[setdiff(trait.gs$TRAIT, finished), ]

# output the remaining gs:
write.table(
    remaining.gs,
    file = output.gs,
    sep = '\t',
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
    );
