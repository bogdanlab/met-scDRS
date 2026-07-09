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

# only get the chromosome 1-22 and protein coding.
exons_df = exons_df[exons_df["chrom"].isin([f"chr{i}" for i in range(1, 23)])].copy()
exons_df = exons_df.loc[exons_df.gene_type == 'protein_coding', ]
genes_df = genes_df[genes_df["chrom"].isin([f"chr{i}" for i in range(1, 23)])].copy()
genes_df = genes_df.loc[genes_df.gene_type == 'protein_coding', ]

###########################################################################################
######                            Create promoter bed file                           ######
###########################################################################################
def get_promoter_bed(genes_df, tss_start_flank = 2000, tss_end_flank = 500):
    df = genes_df.copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    
    promoter_rows = []
    grouped = df.groupby(["gene_name", "chrom", "strand"], sort=False)
    
    for (gene_name, chrom, strand), sub in grouped:
        gene_id = sub["gene_id"].iloc[0]
        if strand == "+":
            # TSS is the most upstream/smallest start
            tss = sub["start"].min()
            promoter_start = tss - tss_start_flank
            promoter_end = tss + tss_end_flank

        elif strand == "-":
            # TSS is the largest end for negative strand
            tss = sub["end"].max()
            promoter_start = tss - tss_end_flank
            promoter_end = tss + tss_start_flank
        
        else:
            continue
        promoter_start = max(promoter_start, 0)
        promoter_rows.append({
            "chrom": chrom,
            "start": promoter_start,
            "end": promoter_end,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "strand": strand,
            "tss": tss,
        })
    promoter_bed = pd.DataFrame(promoter_rows)
    promoter_bed = promoter_bed.sort_values(
        ["chrom", "start", "end", "gene_name"]
    ).reset_index(drop=True)
    promoter_bed = promoter_bed[["chrom", "start", "end", "gene_name"]]
    return promoter_bed

promoter_bed = get_promoter_bed(genes_df)
promoter_bed.columns = ['#chr', 'start', 'end', 'gene']

###########################################################################################
######                                Create exon file                               ######
###########################################################################################
def get_exon_bed(df, exon_id_col="gene_exon_id", exon_order="genomic"):
    # safety:
    df = df.copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    
    # we first sort the df:
    df = df.drop_duplicates(["chrom", "start", "end", "gene_name", "strand"])
    df = df.sort_values(["gene_name", "chrom", "strand", "start", "end"])
    
    # next we can querry each row:
    merged_rows = []
    grouped = df.groupby(["gene_name", "chrom", "strand"], sort=False)
    
    for (gene_name, chrom, strand), sub in tqdm(
        grouped,
        total=grouped.ngroups,
        desc="Merging exon intervals"
    ):
        # getting first instance of the gene name:
        gene_id = sub["gene_id"].iloc[0]
        current_start = None
        current_end = None
        
        # for each exon line:
        for row in sub.itertuples(index=False):
            start = row.start
            end = row.end
            
            # initiate as None
            if current_start is None:
                current_start = start
                current_end = end
                continue
            
            # if the next exon start is more upstream than end, the window broadens
            if start <= current_end:
                current_end = max(current_end, end)
            
            # if next exon start is downstream than last exons end, it means disjoint exons
            else:
                merged_rows.append({
                    "chrom": chrom,
                    "start": current_start,
                    "end": current_end,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "strand": strand,
                })
                current_start = start
                current_end = end
        
        # for the last exon, append to merge:
        if current_start is not None:
            merged_rows.append({
                "chrom": chrom,
                "start": current_start,
                "end": current_end,
                "gene_id": gene_id,
                "gene_name": gene_name,
                "strand": strand,
            })
    
    # merge everything 
    merged = pd.DataFrame(merged_rows)
    
    # Assign exon numbers within each gene
    output_rows = []
    grouped_merged = merged.groupby(["gene_name", "chrom", "strand"], sort=False)
    for (gene_name, chrom, strand), sub in tqdm(
        grouped_merged,
        total=grouped_merged.ngroups,
        desc="Assigning exon IDs"
        ):
        strand = sub["strand"].iloc[0]
        
        if exon_order == "genomic":
            sub = sub.sort_values(["chrom", "start", "end"]).copy()
        elif exon_order == "transcript":
            if strand == "+":
                sub = sub.sort_values(["chrom", "start", "end"]).copy()
            elif strand == "-":
                sub = sub.sort_values(["chrom", "start", "end"], ascending=[True, False, False]).copy()
            else:
                sub = sub.sort_values(["chrom", "start", "end"]).copy()
        else:
            raise ValueError("exon_order must be either 'genomic' or 'transcript'")
        
        # after sorting, add in the id:
        sub["merged_exon_number"] = range(1, len(sub) + 1)
        sub[exon_id_col] = (
            sub["gene_name"].astype(str)
            + "_exon"
            + sub["merged_exon_number"].astype(str)
        )
        output_rows.append(sub)
    
    # concatenate all the exons:
    exon_union = pd.concat(output_rows, ignore_index=True)
    
    # Final BED sorting
    exon_union = exon_union.sort_values(
        ["chrom", "start", "end", "gene_name"]
    ).reset_index(drop=True)
    return exon_union

gene_name_check = (
    exons_df.groupby("gene_name")
    .agg(
        n_gene_ids=("gene_id", "nunique"),
        n_chroms=("chrom", "nunique"),
        n_strands=("strand", "nunique")
    )
    .query("n_gene_ids > 1 or n_chroms > 1 or n_strands > 1")
)
exon_union = get_exon_bed(exons_df)

# keep only the chromosome 1:22
exon_bed = exon_union[['chrom', 'start', 'end', 'gene_exon_id']]
exon_bed.columns = ['#chr', 'start', 'end', 'gene']

###########################################################################################
######                                    intron file                                ######
###########################################################################################
def get_intron_bed(
    genes_df,
    exon_union_df,
    intron_id_col = 'gene_intron_id',
    intron_order = 'genomic'
):
    # input parsing
    required_gene_cols = ["chrom", "start", "end", "gene_id", "gene_name", "strand"]
    required_exon_cols = ["chrom", "start", "end", "gene_id", "gene_name", "strand"]
    
    missing_gene = [c for c in required_gene_cols if c not in genes_df.columns]
    missing_exon = [c for c in required_exon_cols if c not in exon_union_df.columns]
    
    if missing_gene:
        raise ValueError(f"genes_df missing columns: {missing_gene}")
    if missing_exon:
        raise ValueError(f"exon_union_df missing columns: {missing_exon}")
    genes = genes_df[required_gene_cols].copy()
    exons = exon_union_df[required_exon_cols].copy()
    
    # data type management:
    genes["start"] = genes["start"].astype(int)
    genes["end"] = genes["end"].astype(int)
    exons["start"] = exons["start"].astype(int)
    exons["end"] = exons["end"].astype(int)
    
    genes = genes.drop_duplicates(["gene_name", "chrom", "strand"]).copy()
    exons = exons.drop_duplicates(["chrom", "start", "end", "gene_name", "strand"])
    
    # Sort
    genes = genes.sort_values(["gene_name", "chrom", "strand", "start", "end"])
    exons = exons.sort_values(["gene_name", "chrom", "strand", "start", "end"])
    
    intron_rows = []
    
    # Group exons for lookup
    exon_groups = {
        key: sub.sort_values(["start", "end"]).copy()
        for key, sub in exons.groupby(["gene_name", "chrom", "strand"], sort=False)
    }
    
    for gene in tqdm(genes.itertuples(index=False), total = len(genes), desc = 'Getting introns'):
        chrom = gene.chrom
        gene_start = gene.start
        gene_end = gene.end
        gene_id = gene.gene_id
        gene_name = gene.gene_name
        strand = gene.strand
        
        # get the key:
        key = (gene_name, chrom, strand)
        if key not in exon_groups:
            continue
        
        sub_exons = exon_groups[key]
        current_pos = None
        
        for exon in sub_exons.itertuples(index=False):
            exon_start = exon.start
            exon_end = exon.end
        
            # Gap before the start of the exon is the intron:
            if current_pos is not None:
                if current_pos < exon_start:
                    intron_rows.append({
                        "chrom": chrom,
                        "start": current_pos,
                        "end": exon_start,
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "strand": strand,
                    })
            else:
                current_pos = exon_end
            # move the current pos to the end of the exon:
            current_pos = max(current_pos, exon_end)
        
    # intron df:
    intron_df = pd.DataFrame(intron_rows)
    if intron_df.empty:
        return intron_df
    
    output_rows = []
    grouped_merged = intron_df.groupby(["gene_name", "chrom", "strand"], sort=False)
    for (gene_name, chrom, strand), sub in tqdm(
        grouped_merged,
        total=grouped_merged.ngroups,
        desc="Assigning intron IDs"
    ):
        strand = sub['strand'].iloc[0]
        if intron_order == "genomic":
            sub = sub.sort_values(["chrom", "start", "end"]).copy()
        elif intron_order == "transcript":
            if strand == "+":
                sub = sub.sort_values(["chrom", "start", "end"]).copy()
            elif strand == "-":
                sub = sub.sort_values(
                    ["chrom", "start", "end"],
                    ascending=[True, False, False]
                ).copy()
            else:
                sub = sub.sort_values(["chrom", "start", "end"]).copy()
        else:
            raise ValueError("intron_order must be either 'genomic' or 'transcript'")
        
        # add the intron id:
        sub["merged_intron_number"] = range(1, len(sub) + 1)
        sub[intron_id_col] = (
            sub["gene_name"].astype(str)
            + "_intron"
            + sub["merged_intron_number"].astype(str)
        )
        output_rows.append(sub)
    
    intron_df = pd.concat(output_rows, ignore_index=True)
    intron_df = intron_df.sort_values(
        ["chrom", "start", "end", "gene_name"]
    ).reset_index(drop=True)
    return intron_df

intron_df = get_intron_bed(genes_df, exon_union, intron_id_col = 'gene_intron_id', intron_order = 'genomic')
intron_bed = intron_df[['chrom', 'start', 'end', 'gene_intron_id']]
intron_bed.columns = ['#chr', 'start', 'end', 'gene']

###########################################################################################
######                                    OUTPUT                                     ######
###########################################################################################
# finalize the files:
print(promoter_bed["gene"].duplicated().sum())
print(exon_bed["gene"].duplicated().sum())
print(intron_bed["gene"].duplicated().sum())


# output these as bed file:
exon_bed.to_csv(
    exon_out,
    sep="\t",
    header=False,
    index=False
)

# output these as bed file:
intron_bed.to_csv(
    intron_out,
    sep="\t",
    header=False,
    index=False
)

# output promoter:
promoter_bed.to_csv(
    promoter_out,
    sep="\t",
    header=False,
    index=False
)