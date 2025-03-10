### meta-extraction.R #############################################################################
# purpose: grab out the meta data with the same set of cells as the mch:

### PREAMBLE ######################################################################################
# load in library:
library(readxl);

# load in meta data:
meta.data.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx";
meta <- read_excel(meta.data.path, skip = 14);
meta <- data.frame(meta, row.names = 1);

# load in the cells that are in mch:
mch.score.file <- "/u/project/geschwind/lixinzhe/scDRS-output/fraction-gse132489-methyl//mch/KC/all/UKB_460K.repro_MENOPAUSE_AGE.full_score.gz";
mch.score <- read.table(
    file = mch.score.file,
    sep = '\t',
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
    );

# subset to the same set of cells as in mch:
cell.kept <- rownames(mch.score);
subset.meta <- meta[cell.kept, ];

write.table(
    subset.meta,
    file = '/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch-30k-subset-meta.csv',
    sep = ',',
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
    )
