###########################################################################################
######                                       PURPOSE                                 ######
###########################################################################################
# Compare between 2 MDD GWAS
# one being PGC MDD 2025 , and another is the Howard et al MDD GWAS


###########################################################################################
######                                      PREAMBLE                                 ######
###########################################################################################
# load in the packages:
library(ggplot2)
library(data.table)
library(dplyr)

# load in the files
mdd_2025 = data.frame(
    fread(
        file = '/u/home/l/lixinzhe/project-geschwind/data/MDD-GWAS/met_scdrs/out/PGC_MDD_2025.score.gz',
        sep = '\t',
        data.table = FALSE,
        header = TRUE
        ),
    row.names = 1
    )

# load in the Howard et al GWAS:
mdd_howard = data.frame(
    fread(
        file = '/u/home/l/lixinzhe/project-cluo/result/met-scDRS/revision/v1.1/ges215353_full/mean_var_length_arcsine/PASS_MDD_Howard2019.score.gz',
        sep = '\t',
        data.table = FALSE,
        header = TRUE
        ),
    row.names = 1
    )

###########################################################################################
######                                PROCESS - COMPARE                              ######
###########################################################################################
# fdr compute:
mdd_2025$fdr = p.adjust(mdd_2025$pval, method = 'fdr')
mdd_howard$fdr = p.adjust(mdd_howard$pval, method = 'fdr')

# compute correlation between the two:
cor(mdd_2025$norm_score, mdd_howard$norm_score) # return 0.7
