### brain-region-score-aggregation.R ##############################################################
# purpose: aggregate the met-scDRS Z score for each dissection regions for cell type

### PREAMBLE ######################################################################################
# load in specified requirements:
require(docopt);
require(data.table);
require(dplyr);

'Usage:
    annotation-disease-score-model [--score_dir <scdrs> --meta_data <meta> --field1 <group1> --field2 <group2> --p_cutoff <cutoff> --out <output>]

Options:
    --score_dir path to scDRS score file (first column = rownames)
    --meta_data path to meta data on cells associated with the score (first column = rownames)
    --field1 categorical variable in in meta data that you would like to subset into each level
    --field2 categorical variable in meta which you would like to aggregate the scDRS Z score in within the level of field 1
    --p_cutoff FDR cutoff for cell to be included in the aggregation [default: 0.1]
    --out directory to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
meta.data.path <- opts$meta_data;
scDRS.directory <- opts$score_dir;
group1.index <- opts$field1;
group2.index <- opts$field2;
disease <- opts$disease;
p.cutoff <- as.numeric(opts$p_cutoff);
output.dir <- opts$out;
system.date <- Sys.Date();

# for function testing:
# meta.data.path <- "/u/home/l/lixinzhe/project-geschwind/data/GSE215353/processed/meta/subset-mch-meta.csv"
# scDRS.directory <- "/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v2.0/mch/GSE215353-subset/"
# group1.index <- "X_MajorType"
# group2.index <- "X_Region"
# disease <- "PASS_MDD_Howard2019"
# p.cutoff <- 0.1;
# output.path <- '/u/scratch/l/lixinzhe/tmp-file/MDD-aggregate.csv'

# read data:
score.files <- list.files(scDRS.directory, pattern = '\\.score.gz', full.names = TRUE);
all.risk.score <- vector('list', length = length(score.files));
names(all.risk.score) <- score.files;

# read meta:
meta <- read.csv(
    header = TRUE,
    row.names = 1,
    file = meta.data.path
    );

# make sure the group index is in the meta data:
stopifnot(group1.index %in% colnames(meta));
stopifnot(group2.index %in% colnames(meta));

for (result in score.files) {
    all.risk.score[[result]] <- fread(
        file = result,
        sep = '\t',
        header = TRUE,
        stringsAsFactors = FALSE
        );
    all.risk.score[[result]] <- data.frame(all.risk.score[[result]], row.names = 1)
    }

# simplify list names:
list.names <- gsub(scDRS.directory, '', names(all.risk.score));
list.names <- gsub('/', '', list.names);
list.names <- gsub('\\.score.gz', '', list.names);

### PROCESS #######################################################################################
# rename the list names:
names(all.risk.score) <- list.names;

# make sure the disease is within one of the score file:
for (disease in names(all.risk.score)){
    disease.risk.score <- all.risk.score[[disease]];

    # compute the FDR for the disease score:
    disease.risk.score$fdr <- p.adjust(disease.risk.score$pval, method = 'fdr');

    # get the set of cells that have the fdr:
    significant.cell <- rownames(disease.risk.score)[disease.risk.score$fdr <= p.cutoff];
    significant.score <- disease.risk.score[significant.cell, ];

    # put the group1 index and group2 index into the matirx:
    summarize.matrix <- cbind(significant.score, meta[significant.cell, group1.index]);
    summarize.matrix <- cbind(summarize.matrix, meta[significant.cell, group2.index]);
    colnames(summarize.matrix) <- c(colnames(significant.score), 'group1.index', 'group2.index');

    # summarize the data:
    plot.df <- summarize.matrix %>% group_by(group1.index, group2.index) %>% summarize(mean_score = mean(zscore), .groups = 'drop');
    plot.df <- as.data.frame(plot.df);

    # write out the result:
    output.path <- paste0(output.dir, system.date, '-', disease, '-', group1.index, '-', group2.index, 'summary.csv')
    write.table(plot.df, file = output.path, sep = ',')
    }
