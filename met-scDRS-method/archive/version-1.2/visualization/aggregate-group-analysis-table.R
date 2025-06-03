### aggregate-group-analysis-table.R ##############################################################
# purpose: aggregate group analysis into a table:

### PREAMBLE ######################################################################################
# specifiy input variables:
require(docopt)
'Usage:
    significant-cells-visualization-script.R [--result_dir <result> --field <field> --out <output>]

Options:
    --result_dir the directories where the result of the group level analysis are located
    --field which of the field would we like to aggregate, can be hetero_mcp or assoc_mcp [default: assoc_mcp]
    --out path to output file

]' -> doc

# collect user input: 
opts <- docopt(doc)
group.analysis.dir <- opts$result_dir;
group.index <- opts$field;
output.path <- opts$out;
system.date <- Sys.Date();

# for testing purposese:
# group.analysis.dir <- '/u/project/geschwind/lixinzhe/scDRS-output/met-scDRS-v1.2/mch/bimodality-QC-75traits-with-cov/heterogeneity-analysis/';
# group.index <- 'assoc_mcp';
# output.path <- paste0('/u/home/l/lixinzhe/project-geschwind/plot/', system.date, '-group-analysis-aggregated', group.index);

# load in the data:
group.result.file <- list.files(group.analysis.dir, pattern = 'scdrs_group.annotation');

# initiate place holder:
result.list <- vector('list', length = length(group.result.file));
names(result.list) <- group.result.file;
subset.list <- result.list;

### RESULT EXTRACTION #############################################################################
# extract out result:
for (file in group.result.file) {
    # for each of item in the result list, read in the result:
    result.file <- read.table(
        file = paste0(group.analysis.dir, file),
        sep = '\t',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );
    
    result.list[[file]] <- result.file;

    # now extract the data:
    subset.list[[file]] <- result.file[, group.index, drop = FALSE]
    }

# combind the results across traits:
aggregated.result <- Reduce('cbind', subset.list);

# apply fdr correction for each of column:
fdr.adjusted.result <- apply(aggregated.result, 2, p.adjust, method = 'fdr');
colnames(fdr.adjusted.result) <- gsub('.scdrs_group.annotation', '', group.result.file);

# output the data:
write.table(
    fdr.adjusted.result,
    file = output.path,
    sep = ',',
    quote = FALSE,
    row.names = TRUE,
    col.names = TRUE
    );
