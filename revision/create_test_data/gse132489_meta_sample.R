### gse132489_meta_sample.R #######################################################################
# purpose: generate gse132489 meta data for covariate correction testing

### PROCESS ######################################################################################
library(data.table)
library(readxl)

# load in meta data:
meta.data.path <- "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx";
meta <- read_excel(meta.data.path, skip = 14);
meta <- as.data.frame(meta)
rownames(meta) <- meta[, '...1'];
colnames(meta)[1] <- 'cell';

# load in the simulation data frame
csv.path <- "/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.csv"
data = data.frame(fread(csv.path), row.names = 1)

# subset the meta data:
meta <- meta[rownames(data), ]

# get covariates design matrix:
covariates = c('CH_Frac', 'CG_Frac', 'Sample')
const = 1
cov_design = cbind(const, meta[, covariates])
cov_design_matrix = model.matrix(~1 + CH_Frac + CG_Frac + Sample, cov_design)
colnames(cov_design_matrix)[1] = 'const'

# add the cell name into the matrix:
cell = rownames(cov_design_matrix)
cov_design_matrix = cbind(cell, cov_design_matrix)

write.table(
    cov_design_matrix,
    sep = '\t',
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE,
    file = '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-covariate-file.tsv'
    )
