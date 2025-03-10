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
    perturbation.R [--data_matrix <csv> --cell_type_meta <cell_type_metadata> --gs_file <gs_file> --trait <trait> --perturbation_effect_start <start> --perturbation_effect_end <end> --effect_step <effect_step> --gene_number <gene_number> --cell_type <cell_type> --proportion_to_perturb <proportion> --down_sample_to <down_sample_to> --down_sample_rate <down_sample_rate> --overlap_start <overlap_start> --overlap_end <overlap_end> --overlap_step <overlap_step> --replication <replications> --output_dir <dir>]

Options:
    --data_matrix csv file path that houses processed data matrix
    --cell_type_meta csv file path that have cell type information for each cells
    --perturbation_effect_start start point of perturbation effect [default: 1.25]
    --perturbation_effect_end end point of perturbation effect [default: 1.25]
    --effect_step step at which the effect size perturbation should increase [default: 0.01]
    --gs_file gene set file
    --trait trait from gene set file to use for simulation [default: "UKB_460K.body_HEIGHTz"]
    --gene_number gene number [default: 1000]
    --cell_type cell type to perturb
    --proportion_to_perturb [default: 0.1]
    --down_sample_to minima for downsampling [default: 0.2]
    --down_sample_rate [default: 0.2]
    --overlap_start starting point of overlap [default: 0.25]
    --overlap_end ending point of overlap [default: 0.25]
    --overlap_step step at which the overlap perturbation should increase [default: 0.1]
    --replication number of simulation to replicate [default: 100]
    --output_dir output dir
]' -> doc

# collect user input: 
opts <- docopt(doc)
csv.path <- opts$data_matrix;
meta.path <- opts$cell_type_meta;
gs.path <- opts$gs_file;
trait <- opts$trait;
perturbation.effect.start <- as.numeric(opts$perturbation_effect_start);
perturbation.effect.stop <- as.numeric(opts$perturbation_effect_end);
effect.step <- as.numeric(opts$effect_step);
perturbation.effect <- seq(perturbation.effect.start, perturbation.effect.stop, by = effect.step);
causal.gene.number <- as.numeric(opts$gene_number);
causal.cell.type <- as.character(opts$cell_type);
causal.proportion <- as.numeric(opts$proportion_to_perturb);
down.sample.to <- as.numeric(opts$down_sample_to);
down.sample.rate <- as.numeric(opts$down_sample_rate);
down.sample.rate <- seq(down.sample.to, 1, by = down.sample.rate);
causal.overlap.start <- as.numeric(opts$overlap_start);
causal.overlap.end <- as.numeric(opts$overlap_end);
overlap.step <- as.numeric(opts$overlap_step)
causal.overlap <- seq(causal.overlap.start, causal.overlap.end, by = overlap.step);
simulation.replication <- as.numeric(opts$replication);
output.directory <- opts$output_dir;

# for function testing:
# csv.path <- '/u/project/pasaniuc/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/inverted-simulation-subset-GSE132489-mch.csv';
# meta.path <- '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/cell_type_meta_only.csv';
# gs.path <- "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs";
# trait <- 'UKB_460K.body_HEIGHTz';
# perturbation.effect.start <- 1.05;
# perturbation.effect.stop <- 1.05;
# effect.step <- 0.01;
# perturbation.effect <- seq(perturbation.effect.start, perturbation.effect.stop, by = effect.step);
# causal.gene.number <- 1000;
# causal.cell.type <- 'VLMC';
# causal.proportion <- 0.25;
# down.sample.to <- 0.2;
# down.sample.rate <- 0.2;
# down.sample.rate <- seq(down.sample.to, 1, by = down.sample.rate);
# causal.overlap.start <- 0.25;
# causal.overlap.end <- 0.25;
# overlap.step <- 0.1;
# causal.overlap <- seq(causal.overlap.start, causal.overlap.end, by = overlap.step);
# simulation.replication <- 1;
# output.directory <- "/u/scratch/l/lixinzhe/tmp-file/causal-simulation/rare-cell-type/";

# printing out the opts loaded:
cat('### PRINTING OPTIONS LOADED ###\n')
print(opts);

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

# remove all genes that has zero variance:
zero.var.genes <- apply(fraction, 2, sd);
zero.var.genes <- names(zero.var.genes)[zero.var.genes == 0]
fraction <- fraction[, setdiff(colnames(fraction), zero.var.genes)]

# load in the meta data:
cell.type.meta <- read.table(file = meta.path, sep = ',', header = TRUE);
rownames(cell.type.meta) <- cell.type.meta$cell;

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
    'causal cell type: ', causal.cell.type, '\n',
    'causal proportion within cell type to perturb:', causal.proportion, '\n',
    'down sampling rate:', down.sample.rate, '\n',
    'trait: ', trait, '\n',
    'causal effect size: ', perturbation.effect, '\n',
    'causal overlap: ', causal.overlap, '\n'
    );

### PERTURB DATA ##################################################################################
# first we will find cells that are in our dataset that belong to the causal cell type:
cell.in.causal.cell.type <- rownames(cell.type.meta)[cell.type.meta[, 2] == causal.cell.type];
cell.in.causal.cell.type <- intersect(cell.in.causal.cell.type, rownames(fraction));

# make a copy of original fraction so we can operate on it for every downsampling iteration:
original.fraction <- fraction;
total.step <- simulation.replication * length(perturbation.effect) * length(causal.overlap) * length(down.sample.rate)
current.step <- 0;

start.time <- Sys.time();
for (seed in seq(1, simulation.replication)){
    set.seed(seed);
    for (effect.size in perturbation.effect){
        for (overlap in causal.overlap) {
            for (missingness in down.sample.rate) {
                # count the step:
                current.step <- current.step + 1;

                # print out the details about the current iteration:
                cat('replication = ', seed, '\n');
                cat('effect size perturbing', effect.size, '\n');
                
                # down sampling target cell type:
                down.sample.number <- round(length(cell.in.causal.cell.type) * missingness);
                cells.retain <- sample(cell.in.causal.cell.type, down.sample.number);
                cells.to.remove <- setdiff(cell.in.causal.cell.type, cells.retain);
                
                # remove the cells to we want to downsample from the causal cell type from fraction:
                fraction <- original.fraction[setdiff(rownames(original.fraction), cells.to.remove), ]

                # get causal cells:
                causal.cell.number <- round(length(cells.retain) * causal.proportion);
                causal.cell <- sample(cells.retain, size = causal.cell.number);

                # get causal genes:
                # first we will have to get the mouse genes that exist in gene set file:
                gwas.overlap.gene <- sample(
                    intersect(names(trait.gene.set[[trait]]), na.omit(human.gene.name)),
                    size = length(names(trait.gene.set[[trait]])) * causal.overlap
                    );
                causal.gene.from.gs <- colnames(fraction)[match(gwas.overlap.gene, human.gene.name)];

                # next we will have to randomly get the rest of the genes:
                random.causal.gene <- sample(
                    setdiff(colnames(fraction), causal.gene.from.gs),
                    size = causal.gene.number - length(gwas.overlap.gene)
                    );
                causal.gene <- c(causal.gene.from.gs, random.causal.gene);

                # report the number of causal genes:
                cat('number of perturbed genes: ', length(causal.gene), '\n');
                cat('number of causal genes from geneset file: ', length(causal.gene.from.gs), '\n');
                cat('number of perturbed random genes: ', length(random.causal.gene), '\n');

                # lets print out some information about the number of removed cells, causal cells:
                cat('causal cell type: ', causal.cell.type, '\n');
                cat('number of cells in causal cell type: ', cell.in.causal.cell.type, '\n');
                cat('number of cells removed: ', length(cells.to.remove), '\n');
                cat('number of causal cells: ', causal.cell.number, '\n');

                # perturb the fraction matrix by causal genes and causal cells:
                perturbed.fraction <- fraction;
                perturbed.fraction[causal.cell, causal.gene] <- fraction[causal.cell, causal.gene] * effect.size;

                # cap the perturbed fraction to be bound between 0 and 1:
                perturbed.fraction[perturbed.fraction >= 1] = 1

                # add checks to make sure that all the cells that are perturbed are from the cell type:
                stopifnot(cell.type.meta[causal.cell, 2] == causal.cell.type);
                stopifnot(cell.type.meta[cells.to.remove, 2] == causal.cell.type);
                stopifnot(length(cells.to.remove) + length(cells.retain) == length(cell.in.causal.cell.type));

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
                    '-cell-type-',
                    causal.cell.type,
                    '-down-sample-',
                    missingness,
                    '-proportion-',
                    causal.proportion,
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
                length(causal.cell) <- length(causal.gene);
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

                # update progress bar:
                cat('\n', '\n', '\n');
                cat(current.step, 'out of ', total.step, ' simulation completed! \n');
                cat('\n', '\n', '\n');
                }
            }
        }
    }
end.time <- Sys.time();
print(end.time - start.time);