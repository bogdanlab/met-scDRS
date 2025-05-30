# met-scDRS tutorial

met-scDRS is a statistical tool for finding disease associated cells and genes in single cell methylome data.
This tutorial will illustrate how we used met-scDRS to find:

- Major depression disorder (MDD) associated cells in GSE215353 methylome atlas
- prioritized gene set compared to GWAS MAGMA genes
- brain regions disease association heterogeneity

## Finding MDD associated cells in GSE215353 methylome atlas
### data download:
[Single-Cell DNA Methylation and 3D Genome Human Brain Atlas](https://cellxgene.cziscience.com/collections/fdebfda9-bb9a-4b4b-97e5-651097ea07b0)

### Installation of met-scDRS
```bash
# Clone the repo
git clone git@github.com:bogdanlab/met-scDRS.git
cd met-scDRS/met-scDRS-method/V1.0.0-rc1

# (Optional) Create a new environment
conda create -n metscdrs-env python=3.10 -y
conda activate metscdrs-env

# Install in editable mode
pip install -e .
```

## License
MIT
