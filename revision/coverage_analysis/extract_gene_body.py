import re
import pandas as pd
from tqdm import tqdm

# load in the gtf file:
output_dir = '/u/home/l/lixinzhe/project-geschwind/port/scDRS/gtf_info/feature_out'
gene_body_out = f"{output_dir}/gene_body.bed"
promoter_out = f"{output_dir}/promoter_2kb_500bp.bed"
exon_out = f"{output_dir}/exon_merged_segments.bed"
intron_out = f"{output_dir}/intron_segments.bed"
gtf_file = '/u/home/l/lixinzhe/project-geschwind/port/scDRS/gtf_info/gencode.v50.primary_assembly.annotation.gtf'

###########################################################################################
######                                    Parse file                                 ######
###########################################################################################
# define a  helper:
def parse_attr(attr, key):
    for item in attr.strip().split(";"):
        item = item.strip()
        if not item:
            continue

        parts = item.split(None, 1)  # split on any whitespace, once
        if len(parts) != 2:
            continue

        k, v = parts

        if k == key:
            return v.strip().strip('"')

    return None

# load in the data:
genes = []
exons = []

with open(gtf_file) as f:
    total_lines = sum(1 for _ in f)

with open(gtf_file) as f:
    for line in tqdm(f, total = total_lines):
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        chrom, source, feature, start, end, score, strand, frame, attr = fields

        start = int(start)
        end = int(end)

        gene_id = parse_attr(attr, "gene_id")
        gene_name = parse_attr(attr, "gene_name")
        gene_type = parse_attr(attr, "gene_type")

        if gene_id is None:
            continue

        # Convert GTF to BED:
        # GTF: 1-based inclusive
        # BED: 0-based half-open
        bed_start = start - 1
        bed_end = end

        if feature == "gene":
            genes.append({
                "chrom": chrom,
                "start": bed_start,
                "end": bed_end,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
                "strand": strand,
            })

        elif feature == "exon":
            transcript_id = parse_attr(attr, "transcript_id")
            exon_number = parse_attr(attr, "exon_number")

            exons.append({
                "chrom": chrom,
                "start": bed_start,
                "end": bed_end,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
                "transcript_id": transcript_id,
                "exon_number": exon_number,
                "strand": strand,
            })

genes_df = pd.DataFrame(genes)
exons_df = pd.DataFrame(exons)

###########################################################################################
######                                    EXTRACT GENE                               ######
###########################################################################################
# figure out gene body region;
bed_df = genes_df.loc[
    genes_df.gene_type == 'protein_coding',
    ["chrom", "start", "end", "gene_name"]].copy()

bed_df = bed_df.groupby(["chrom", "gene_name"], as_index=False).agg(
    start = ('start', "min"),
    end = ('end', 'max')
    )[["chrom", "start", "end", "gene_name"]]

# filter to autosomal chromosome:
valid_chroms = [f"chr{i}" for i in range(1, 23)]

# remove genes that are not in autosomal chromosome:
bed_df_out = bed_df.loc[bed_df['chrom'].isin(valid_chroms), :]

# output:
bed_df_out.to_csv(
    gene_body_out,
    sep="\t",
    header=False,
    index=False
)