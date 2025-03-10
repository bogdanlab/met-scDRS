### parallel-splitter.R ###########################################################################
# purpose: split a input gene set into multiple different gene set for parallel compute

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    get-remaininig-gs.R [--scDRS_dir <scDRS_dir> --gs_file <gs_file> --output_gs <output_gs>]

Options:
    --gs_file gene set file that were originally used to compute the scores in the folder
    --output_gs path to file that contains the uncomputed gene sets
]' -> doc

# collect user input: 
opts <- docopt(doc)
scDRS.dir <- opts$scDRS_dir;
gs.path <- opts$gs_file;
output.gs <- opts$output_gs;

# for function test:
# gs.path <- "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.75_traits.rv1.gs"
# output.gs <- "/u/scratch/l/lixinzhe/tmp-file/job-array/GSE215353-mch-with-cov-dev.gs"

# load in deocmposer function:
source('/u/home/l/lixinzhe/project-github/spells/scDDS/gs-out.R')
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

# make another function:
split_list <- function(input_list, split_num = 40) {
    # Determine the length of each split
    split_size <- length(input_list) %/% split_num
    remainder <- length(input_list) %% split_num

    # Initialize the result list
    result <- vector("list", split_num)
    start <- 1

    for (i in 1:split_num) {
    # Calculate the end index for each part
    end <- start + split_size + ifelse(i <= remainder, 1, 0) - 1
    result[[i]] <- input_list[start:end]
    start <- end + 1
    }
  
    return(result)
    }

# read in the gs list;
trait.gs <- read.table(
    file = gs.path,
    sep = '\t',
    header = TRUE
    );

# split gene sets:
trait.gene.set <- lapply(trait.gs$GENESET, gs.decomposer);
names(trait.gene.set) <- trait.gs$TRAIT;
splitted_gs <- split_list(trait.gene.set)

# output the splitted gene set into the gs files:
# output them into collection:
for (job in seq(1, length(splitted_gs))) {
    gs.out(out_list = splitted_gs[[job]], path = paste0(output.gs, job))
}