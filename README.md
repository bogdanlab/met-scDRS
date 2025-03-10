# met-scDRS tutorial

met-scDRS is a statistical tool for finding disease associated cells and genes in single cell methylome data.
This tutorial will illustrate how we used met-scDRS to find:

- Major depression disorder (MDD) associated cells in GSE215353 methylome atlas
- prioritized gene set compared to GWAS MAGMA genes
- brain regions disease association heterogeneity

## Finding MDD associated cells in GSE215353 methylome atlas
### data download:
[Single-Cell DNA Methylation and 3D Genome Human Brain Atlas](https://cellxgene.cziscience.com/collections/fdebfda9-bb9a-4b4b-97e5-651097ea07b0)

### preprocessing:
filtered for low variance genes (5th percentile or lower)
```py
gene_variances = pd.Series(merged_adata.X.var(axis=0), index=merged_adata.var_names)
percentile_5th = gene_variances.quantile(0.05)
variance_mask = gene_variances >= percentile_5th
merged_adata = merged_adata[:, variance_mask]
```

as well as taking 1-raw fraction:
```py
merged_adata.X = 1 - merged_adata.X
```

outputting to h5ad file:
```py
merged_adata.write('processed-mch.h5ad')
```

full script documented in:
```sh
python inverse-fraction.py \
    --fraction_path ${methylation_csv} \
    --output_path ${processed_fraction}

```
options:
 - fraction_path: path to raw fraction csv that you would like to filter for low variance genes and inverse the fraction
 - output_path: path for writing the processed fraction csv

```sh
python csv-to-h5ad.py \
    --fraction_path ${processed_fraction}
    --meta_path ${meta_data_csv}
    --h5ad_path ${output_h5ad}
```
options:
 - fraction_path: path to processed fraction output of inverse-fraction.py
 - meta_data_csv: path to cells' meta data. e.g.: transcription batch, donor
 - h5ad_path: output path for converted h5ad

### Computation:
```sh
scdrs compute-score \
    --h5ad-file ${h5ad_path} \
    --h5ad-species ${species} \
    --gs-file ${gs_file} \
    --gs-species human \
    --out-folder ${output_dir} \
    --flag-filter-data False \
    --flag-raw-count False \
    --n-ctrl 1000 \
    --weight_opt "inv_std" \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True
```
options:
 - h5ad-file: the processed h5ad file for with inverse methylation fraction (1-raw fraction)
 - h5ad-species: which species is sequenced in h5ad, supports mouse or human
 - gs-file: gene set file that contain the disease putative genes [more info](https://martinjzhang.github.io/scDRS/file_format.html)
 - gs-species: which species is putative disease genes coded for, supports mouse or human
 - out-folder: folder for output, create subdirectories for each gene set enclosed in gs-file
 - n-ctrl: number of control gene sets to sample
 - flag-return-ctrl-raw-score: if you wish met-scDRS to return raw scores for sampled control genes
 - flag-return-ctrl-norm-score: if you wish met-scDRS to return normalized score for sampled control genes

## Gene set prioritization:
Compute the gene wise correlation between the processed methylation level and the met-scDRS risk score
[here](https://github.com/xinzhe-ucla/scDRS-applications/blob/main/code/GSE215353/production/v2.0/visualization/GO-enrichment.R) is how we computed the prioritized genes using met-scDRS risk scores and the processed methylation level. At its core, we use the spearman correlation between risk score and methylation level to prioritized genes


## License

MIT
