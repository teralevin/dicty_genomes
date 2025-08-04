#!/usr/bin/env python3
r"""
Convert BLAST hits of TEs into bedGraph repeat density format with sliding window bins.

Usage:
    python 05_generate_TE_bedgraph.py \
        -i path/to/blast_hits.csv \
        -o path/to/output.bedgraph \
        [-f TE_family] \
        [-b BIN_SIZE] \
        [--step_size STEP_SIZE] \
        [--mapping CONTIG_MAPPING_FILE]

Options:
    -i, --input       Path to input BLAST CSV file (required)
    -o, --output      Path to output bedGraph file (required)
    -f, --family      (Optional) Specific TE family name (from qseqid) for example "DIRS1"
    -b, --bin_size    (Optional) Bin size in bp (default: 10000)
    --step_size       (Optional) Step size for overlapping bins (default: 5000)
    --mapping         (Optional) Path to contig renaming file (tab-separated: old new) 
"""
import pandas as pd
import numpy as np
import math
import argparse

#Chromosome name normalization 
def normalize_chrom_names(df: pd.DataFrame, column):
    df[column] = (
        df[column]
        .str.replace(r"v[14]_f1_", "", regex=True)
        .str.replace(r"_p\d+$", "", regex=True)
        .str.replace(r"^([A-Za-z0-9]+)_\1_(contig|scaffold|chr)_", r"\1_\2_", regex=True)
    )
    return df

#Load BLAST CSV with your column structure 
def load_blast_te_hits(filepath: str) -> pd.DataFrame:
    cols = [
        "sseqid", "bitscore", "evalue", "qlen", "qstart", "qend",
        "slen", "sstart", "send", "sseq", "length", "qseqid"
    ]
    return pd.read_csv(filepath, names=cols)

#Filter by TE family (from qseqid) 
def filter_by_family(df: pd.DataFrame, family: str | None) -> pd.DataFrame:
    if family:
        df = df[df["qseqid"] == family]
        print(f" Filtered to {len(df)} hits from family: {family}")
    else:
        print(f" Using all TE hits: {len(df)} entries")
    return df

#Extract genomic coordinates
def get_te_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    df["start"] = df[["sstart", "send"]].min(axis=1)
    df["end"] = df[["sstart", "send"]].max(axis=1)
    return df

#Use slen (subject length) for bin generation 
def get_contig_lengths_from_slen(df: pd.DataFrame) -> dict:
    return df.groupby("sseqid")["slen"].max().to_dict()

#Build overlapping bins 
def build_bins(contig: str, length: int, bin_size: int, step_size: int) -> pd.DataFrame:
    starts = np.arange(0, length, step_size)
    return pd.DataFrame({
        "contig": contig,
        "start": starts,
        "end": starts + bin_size
    })

#Count TEs per overlapping bin
def count_tes_per_bin(df: pd.DataFrame, bin_size: int = 10000, step_size: int = 5000) -> pd.DataFrame:
    contig_lengths = get_contig_lengths_from_slen(df)
    all_bins = []

    for contig, sub_df in df.groupby("sseqid"):
        length = math.ceil(contig_lengths[contig] / step_size) * step_size
        bins = build_bins(contig, length, bin_size, step_size)
        bins["count"] = 0

        for _, row in sub_df.iterrows():
            overlap = (bins["start"] < row["end"]) & (bins["end"] > row["start"])
            bins.loc[overlap, "count"] += 1

        all_bins.append(bins)

    return pd.concat(all_bins)

# Apply optional contig renaming using a mapping file
def rename_chromosomes_in_bedgraph(df: pd.DataFrame, mapping_file: str) -> pd.DataFrame:
    mapping = pd.read_csv(mapping_file, sep=r"\s+", engine="python", header=None, names=["old_chrom", "new_chrom"])
    rename_dict = dict(zip(mapping["old_chrom"], mapping["new_chrom"]))
    df["contig"] = df["contig"].map(rename_dict).fillna(df["contig"])
    return df

#Save to bedGraph 
def save_bedgraph(df: pd.DataFrame, output_file: str):
    df[["contig", "start", "end", "count"]].to_csv(output_file, sep="\t", header=False, index=False)
    print(f" Saved bedGraph: {output_file}")


def main(input_csv: str, output_bedgraph: str, family: str | None, bin_size: int, mapping_file: str | None, step_size: int):
    df = load_blast_te_hits(input_csv)
    df = filter_by_family(df.drop_duplicates(), family)
    df = get_te_coordinates(df)
    df = normalize_chrom_names(df, "sseqid")
    binned_df = count_tes_per_bin(df, bin_size, step_size)

    if mapping_file:
        binned_df = rename_chromosomes_in_bedgraph(binned_df, mapping_file)

    save_bedgraph(binned_df, output_bedgraph)

#Command-line Interface 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert BLAST hits of TEs into bedGraph repeat density format with sliding window bins.")
    parser.add_argument("-i", "--input", required=True, help="Path to input BLAST CSV file")
    parser.add_argument("-o", "--output", required=True, help="Path to output bedGraph file")
    parser.add_argument("-f", "--family", default=None, help="(Optional) Specific TE family name (from qseqid)")
    parser.add_argument("-b", "--bin_size", type=int, default=10000, help="(Optional) Bin size in bp (default: 10000)")
    parser.add_argument("--step_size", type=int, default=5000, help="(Optional) Step size for overlapping bins (default: 5000)")
    parser.add_argument("--mapping", help="(Optional) Path to contig renaming file (tab-separated: old new)")

    args = parser.parse_args()
    main(args.input, args.output, args.family, args.bin_size, args.mapping, args.step_size)
