### get_cell_type_meta_data.R #####################################################################
# purpose: take in the meta data and output a csv file where it only contain cell id in the first
# column and then the cell type information on the second column

### PREAMBLE ######################################################################################
# load in library:
library(readxl);

# specify input to the script:
require(docopt)
'Usage:
    get_cell_type_meta_data.R [--meta_csv <meta_csv> --meta_column <cell_type> --output <output_file>]

Options:
    --meta_csv path for meta data matrix
    --meta_column column name that encode the column that you wish to extract out
    --output output path to output to
]' -> doc

# gather arguments:
opts <- docopt(doc)
meta.data.path <- opts$meta_csv;
meta.column <- opts$meta_column;
output.path <- opts$output;

# for testing:
meta.data.path = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx"
meta.column <- "MajorType"
output.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/cell_type_meta_only.csv"

### EXTRACT INFO ##################################################################################
# get the information from meta:
meta <- read_excel(meta.data.path, skip = 14);
meta <- as.data.frame(meta);
rownames(meta) <- meta[, '...1'];
colnames(meta)[1] <- 'cell';
meta <- meta[meta[, 'Pass QC'], ];
meta.output <- meta[, c('cell', meta.column)];

# output the extracted meta:
write.table(
    meta.output,
    file = output.path,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = ','
    );
