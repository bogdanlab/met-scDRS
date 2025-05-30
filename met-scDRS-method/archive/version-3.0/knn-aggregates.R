### knn-aggregates.R ##############################################################################
# purpose: extend the met-scDRS to CpG methylation by using knn as an aggregate:

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
require(progress)
require(pbmcapply)
require(FNN)
require(data.table)

'Usage:
    knn-aggregates.R [--meta_data <meta_path> --fraction <methylation_fraction> --dr1 <column_index_1> --dr2 <column_index_2> --k <neibour_number> --threads <threads> --output <output>]

Options:
    --meta_data file path for meta data that contains umap or tsne
    --fraction methylation fraction matrix to aggregate over
    --dr1 first column index for dimensionality reduction (e.g.: UMAP1)
    --dr2 second column index for dimensionality reduction (e.g.: UMAP2)
    --k number of neibours to aggregate over
    --threads number of threads for parrallel computation
    --output path aggregated output
]' -> doc

# collect user input: 
opts <- docopt(doc)
meta_data_path <- opts$meta_data
fraction <- opts$fraction
dr1 <- opts$dr1
dr2 <- opts$dr2
neibour_num <- as.numeric(opts$k)
threads <- as.numeric(opts$threads)
output <- opts$output

# print out what are the arguments:
print(str(opts))

# for test:
# meta_data_path <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/production/meta_data.csv'
# fraction <- '/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/v2.0/processed-met-scDRS-mcg.csv'
# dr1 <- 'UMAP_1'
# dr2 <- 'UMAP_2'
# neibour_num <- 10
# threads <- 4

# load in the data:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta_data_path
    );

# load in methylation fractions:
cat('reading data \n')
fraction <- fread(
    file = fraction,
    sep = ',',
    # nrows = 500,
    header = TRUE,
    data.table = TRUE
    );

### PROCESS #######################################################################################
# subset to dimensionally reduced features;
meta <- meta[fraction$cell, ]
features <- c(dr1, dr2)
knn_result <- get.knn(as.matrix(meta[, features]), k = neibour_num)
knn_index <- knn_result$nn.index
rownames(knn_index) <- rownames(fraction)

# make progress bar:
cat('aggregating with neibours \n')
# for each of the cell in the fraction, aggregate the five nearest neibours:
# note: byte wise identical to the non-multithreaded version
start <- Sys.time()
results <- pbmclapply(
    X = 1:nrow(fraction),
    FUN = function(cell) {
        neighbour <- knn_index[cell, ]
        index <- sort(c(cell, neighbour))
        colSums(as.matrix(fraction[index, -1, with = FALSE]))
        },
    mc.style = 'ETA',
    ignore.interactive = TRUE,
    mc.cores = threads
    )

# merge the data by rbind and add the cell column back in:
merged_data <- as.data.table(do.call(rbind, results))
cell <- fraction$cell
merged_data <- cbind(cell, merged_data)
end <- Sys.time()

# print out the ellapse time:
print(end - start)

# non multithreaded:
# creates a progress bar:
# pb <- progress::progress_bar$new(
#     format = "[:bar] :current/:total (:percent)",
#     total = nrow(fraction),
#     clear = FALSE
#     );
# merged_data <- data.table()
# for (cell in seq(1, nrow(fraction))){
#     # grab out the cells that are the cell's neibour
#     neighbour <- knn_index[cell, ]
#     index <- sort(c(cell, neighbour))
    
#     # Sum the cell and its neibours up:
#     col_sums <- as.list(rowSums(t(as.matrix(fraction[index, 2:ncol(fraction), with = FALSE]))))
#     merged_data <- rbind(merged_data, col_sums)

#     # update progress bar:
#     pb$tick()
# }
# merged_data <- cbind(fraction$cell, merged_data)

### OUTPUT ########################################################################################
# output the merged data:
fwrite(
    merged_data,
    sep = ',',
    row.names = FALSE,
    col.names = TRUE,
    file = output,
    quote = 'auto',
    nThread = 1
    );

### MAYBE HAVE A SECOND ALGO THAT DOES MERGE CELLS WITHIN CELL TYPE