### purturbation.R ################################################################################
# purpose: perturb data for causal simulation 

### PREAMBLE ######################################################################################
# load in the libraries:
library(tidyverse);
library(Seurat);
library(sceasy);
library(data.table);
library(readxl);
library(biomaRt);

# define the input and its help page:
require(docopt)
'Usage:
    perturbation.R [--data_matrix <csv> --gs_file <gs_file> --trait <trait> --perturbation_effect_start <start> --perturbation_effect_end <end> --effect_step <effect_step>  --gene_number <gene_number> --cell_number <cell_number> --overlap_start <overlap_start> --overlap_end <overlap_end> --overlap_step <overlap_step> --replication <replications> --output_dir <dir>]

Options:
    --data_matrix csv file path that houses processed data matrix
    --perturbation_effect_start start point of perturbation effect [default: 1.25]
    --perturbation_effect_end end point of perturbation effect [default: 1.25]
    --effect_step step at which the effect size perturbation should increase [default: 0.01]
    --gs_file gene set file
    --trait trait from gene set file to use for simulation [default: "UKB_460K.body_HEIGHTz"]
    --gene_number gene number [default: 1000]
    --cell_number cell number [default: 500]
    --overlap_start starting point of overlap [default: 0.25]
    --overlap_end ending point of overlap [default: 0.25]
    --overlap_step step at which the overlap perturbation should increase [default: 0.1]
    --replication number of simulation to replicate [default: 100]
    --output_dir output dir
]' -> doc

# collect user input: 
opts <- docopt(doc)
csv.path <- opts$data_matrix;
gs.path <- opts$gs_file;
trait <- opts$trait;
perturbation.effect.start <- as.numeric(opts$perturbation_effect_start);
perturbation.effect.stop <- as.numeric(opts$perturbation_effect_end);
effect.step <- as.numeric(opts$effect_step);
perturbation.effect <- seq(perturbation.effect.start, perturbation.effect.stop, by = effect.step);
causal.gene.number <- as.numeric(opts$gene_number);
causal.cell.number <- as.numeric(opts$cell_number);
causal.overlap.start <- as.numeric(opts$overlap_start);
causal.overlap.end <- as.numeric(opts$overlap_end);
overlap.step <- as.numeric(opts$overlap_step)
causal.overlap <- seq(causal.overlap.start, causal.overlap.end, by = overlap.step);
simulation.replication <- as.numeric(opts$replication);
output.directory <- opts$output_dir;

# for function testing:
# csv.path <- '/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv';
# gs.path <- "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs";
# trait <- 'UKB_460K.body_HEIGHTz';
# perturbation.effect <- 1.25;
# causal.gene.number <- 1000;
# causal.cell.number <- 500;
# causal.overlap <- 0.25;
# simulation.replication <- 1;
# output.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/causal-simulation/fraction/";

# load in function for converting gene names between mouse to human:
source('/u/home/l/lixinzhe/project-github/scDRS-applications/spell-book/gene-name-converter.R');

### PROCESSING ####################################################################################
# grab out the processed data:
fraction <- fread(
    file = csv.path,
    header = TRUE,
    sep = ',',
    stringsAsFactors = FALSE,
    data.table = TRUE,
    nThread = 1
    );
fraction <- fraction %>%
    column_to_rownames('cell');
cat('fraction loaded\n');

# load in the gene set:
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

# load in gs file
trait.gs <- read.table(
    file = gs.path,
    sep = '\t',
    header = TRUE
    );
# split gene sets:
trait.gene.set <- lapply(trait.gs$GENESET, gs.decomposer);
names(trait.gene.set) <- trait.gs$TRAIT;

# get the gene name conversions:
human.gene.name <- gene.name.converter(colnames(fraction), from = 'mouse_to_human');

# print out the input:
cat(
    'input parameters: \n',
    'causal gene number: ', causal.gene.number, '\n',
    'causal cell number: ', causal.cell.number, '\n',
    'trait: ', trait, '\n',
    'causal effect size: ', perturbation.effect, '\n',
    'causal overlap: ', causal.overlap, '\n'
    );

### PERTURB DATA ##################################################################################
for (seed in seq(1, simulation.replication)){
    set.seed(seed);
    for (effect.size in perturbation.effect){
        for (overlap in causal.overlap) {

            # file check module:
            # define the file path that we should check if existed:
            h5ad.path <- paste0(
                output.directory,
                'seed-',
                seed,
                '-effect-',
                effect.size,
                '-overlap-',
                overlap,
                '-causal-simulation.h5ad'
                );
            # perform file check:
            if (file.exists(h5ad.path)) {
                # if the file already existed, print a message and do nothing:
                cat('the perturbation file already exist! Skipping perturbation \n')
                } else {

            # if the file does not exist yet:
            # print out the details about the current iteration:
            cat('replication = ', seed, '\n');
            cat('effect size perturbing', effect.size, '\n');

            # get causal cells:
            causal.cell <- sample(rownames(fraction), causal.cell.number);

            # get causal genes:
            # first we will have to get the mouse genes that exist in gene set file:
            gwas.overlap.gene <- sample(
                intersect(names(trait.gene.set[[trait]]), na.omit(human.gene.name)),
                size = length(names(trait.gene.set[[trait]])) * overlap
                );
            causal.gene.from.gs <- colnames(fraction)[match(gwas.overlap.gene, human.gene.name)];

            # next we will have to randomly get the rest of the genes:
            random.causal.gene <- sample(
                setdiff(colnames(fraction), causal.gene.from.gs),
                size = causal.gene.number - length(gwas.overlap.gene)
                );
            causal.gene <- c(causal.gene.from.gs, random.causal.gene);

            # perturb the fraction matrix by causal genes and causal cells:
            perturbed.fraction <- fraction;
            perturbed.fraction[causal.cell, causal.gene] <- fraction[causal.cell, causal.gene] * effect.size;

            # cap the perturbed fraction to be bound between 0 and 1:
            perturbed.fraction[perturbed.fraction >= 1] = 1

            ### CONVERSION ####################################################################################
            # check the result is still bound between 1 and 0:
            stopifnot(perturbed.fraction >= 0 & perturbed.fraction <= 1);
            perturbed.fraction <- t(perturbed.fraction); # transposed so the cells in column, genes in row

            # put the result into Seurat object:
            fraction.seurat <- CreateSeuratObject(
                counts = perturbed.fraction,
                min.cells = 0,
                min.features = 0
                );

            # printing the first 5 cells and first 5 genes of the :
            cat('printing first 5 cells and genes: \n');
            print(perturbed.fraction[1:5,1:5]);

            # convert the seurat object into the h5ad format:
            h5ad.path <- paste0(
                output.directory,
                'seed-',
                seed,
                '-effect-',
                effect.size,
                '-overlap-',
                overlap,
                '-causal-simulation.h5ad'
                );
            convertFormat(
                fraction.seurat,
                from = "seurat",
                to = "anndata",
                outFile = h5ad.path
                );
            
            # should also output a log file that states what are the permuted info:
            perturbation.log <- data.frame(
                causal.genes = causal.gene,
                causal.cells = causal.cell,
                effect.size = effect.size,
                overlap = overlap
                );
            perturbation.log$causal.cells[duplicated(perturbation.log$causal.cells)] <- NA;

            # output the perturbation file:
            file.path <- gsub('h5ad', 'csv', h5ad.path);
            write.table(perturbation.log, file = file.path, sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE);

            # tidy up:
            rm(list = c('fraction.seurat', 'perturbed.fraction'));
            gc();
            }
        }
    }
}