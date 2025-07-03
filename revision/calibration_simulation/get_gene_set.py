### get_gene_set.py ###############################################################################
# purpose: get the gene set for calibration studies:

### PREAMBLE ######################################################################################
import met_scdrs
import polars as pl
import subprocess

###########################################################################################
######                                    preprocess                                 ######
###########################################################################################
# define parameters:
H5AD_FILE = '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/simulation-subset-GSE132489-mch.h5ad'
PREPROCESS_METHOD = 'inverse'
VARIANCE_CLIP = 5
TRANSFORMATION = 'arcsine'
WEIGHT_OPT = 'inv_std'
CTRL_MATCH_OPT = 'mean_var_length'
VERBOSE = True

# load in h5ad:
h5ad = met_scdrs.util.load_h5ad(H5AD_FILE)

# normalze h5ad file:
h5ad = met_scdrs.normalize(
    h5ad_obj = h5ad,
    method = PREPROCESS_METHOD,
    variance_clip = VARIANCE_CLIP,
    transformation = TRANSFORMATION,
    verbose = VERBOSE
)

met_scdrs.preprocess(
    h5ad,
    cov=None,
    n_mean_bin=10,
    n_var_bin=10,
    n_length_bin = 10,
    copy=False,
    weight_option=WEIGHT_OPT,
    ctrl_match_key=CTRL_MATCH_OPT,
    verbose = VERBOSE)

###########################################################################################
######                                    select genes                              ######
###########################################################################################
# h5ad to csv:
h5ad_pl = pl.DataFrame(
    h5ad.X,
    schema = h5ad.var_names.to_list()
).with_columns([
    pl.Series('cell', h5ad.obs_names.to_list())
])

# reformat columns:
h5ad_pl = h5ad_pl.select(['cell'] + h5ad_pl.columns[:-1])

# output:
h5ad_pl.write_csv('/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/normalized-preprocessed-simulation-subset-GSE132489-mch.csv')

# invoke subprocess:
for gene_num in [100, 500, 1000]:
    cmd = f"""\
    Rscript /u/home/l/lixinzhe/project-github/met-scDRS/scDRS-methylation-simulation/high-expression-gene-permutation.R \
        --data_matrix '/u/home/l/lixinzhe/project-geschwind/data/Liu_et_al_2021_methylation_gse132489/normalized-preprocessed-simulation-subset-GSE132489-mch.csv' \
        --gs_file "/u/project/geschwind/lixinzhe/scDRS-output/magma-out/Kangcheng-gs/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs" \
        --gene_number "{gene_num}" \
        --quantile "0.95" \
        --replication "100" \
        --output_dir "/u/home/l/lixinzhe/project-geschwind/port/scratch/revision/null-simulation-95percentile/" > /u/scratch/l/lixinzhe/tmp-file/tmp-simulation/gene_num_{gene_num}.logfile
    """
    print(f'getting gene set for null simulation with number of genes {gene_num}')
    subprocess.run(cmd, shell=True)

