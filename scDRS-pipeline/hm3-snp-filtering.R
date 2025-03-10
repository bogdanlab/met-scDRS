### hm3-snp-filtering.R ###########################################################################
# purpose: filter GWAS based on hapmap 3 snps

### PREAMBLE ######################################################################################
# define the input and its help page:
require(docopt)
'Usage:
    hm3-snp-filtering.R [--GWAS <gwas_file> --hm3_rsid_list <hm3> --column_use <rsid_column> --out <output>]

Options:
    --GWAS path to GWAS summary statistics
    --hm3_rsid_list a list of hm3 rsid to retain
    --column_use column that enclose rsid
    --out output of hm3 rsid filtered GWAS summary statistics
]' -> doc

# load in packages:
library(data.table);

# collect user input: 
opts <- docopt(doc)
gwas.path <- opts$GWAS;
hm3.path <- opts$hm3_rsid_list;
rsid.column <- opts$column_use;
output.path <- opts$out;

# test the function:
# gwas.path <- "/u/project/geschwind/lixinzhe/data/ASD-GWAS/iPSYCH-PGC_ASD_Nov2017";
# hm3.path <- "/u/project/pasaniuc/pasaniucdata/DATA/LDSC/1000G_EUR/1000G_hm3_noMHC.rsid";
# rsid.column <- "SNP"
# output.path <- "/u/scratch/l/lixinzhe/tmp-file/test-scripting-tmp/hm3-filtered-iPSYCH-PGC_ASD_Nov2017.txt"

# data loading:
gwas <- fread(file = gwas.path, sep = '\t', data.table = FALSE);
hm3.snps <- fread(file = hm3.path, sep = '\n', data.table = FALSE, header = FALSE);

# get the set of snps that are in the hm3:
cat('original number of snps: ', nrow(gwas), '\n');
common.snps <- intersect(gwas[, rsid.column], hm3.snps[, 1]);
cat('number of snps in GWAS present in hm3: ', length(common.snps), '\n');

# filter the snps based on hm3:
hm3.filtered.gwas <- gwas[gwas[, rsid.column] %in% common.snps, ];

# write the data:
fwrite(
    hm3.filtered.gwas,
    file =  output.path,
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
    );
