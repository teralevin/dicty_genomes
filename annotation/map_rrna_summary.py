#!/usr/bin/env python3
"""
map_rrna_summary.py — Map gene IDs to rRNA elements and summarize BLAST hits

USAGE
    # Map gene IDs to rRNA elements:
    map_rrna_summary.py map -s RRNA_HITS_SUMMARY.tsv -g gene_to_rrna_table.tsv -o mapped_rrna_hits.tsv

    # Summarize mapped rRNA hits per species:
    map_rrna_summary.py summarize -i mapped_rrna_hits.tsv -o rrna_per_species_summary.tsv

COMMANDS
    map        Map gene IDs to rRNA elements using a lookup table.
    summarize  Summarize rRNA hits per species from mapped data.

OPTIONS
    map
        -s, --summary      Input BLAST hits summary TSV
        -g, --gene_table   Gene-to-rRNA mapping TSV
        -o, --output       Output TSV with rRNA_element column
    summarize
        -i, --input        Mapped summary TSV
        -o, --output       Output species × rRNA_element summary TSV
"""

import pandas as pd
import argparse


# === Chromosome name normalization ===
def normalize_chrom_names(df, column):
    df[column] = (
        df[column]
        .str.replace(r"v[14]_f1_", "", regex=True)
        .str.replace(r"_p\d+$", "", regex=True)
        .str.replace(r"^([A-Za-z0-9]+)_\1_(contig|scaffold|chr)_", r"\1_\2_", regex=True)
    )
    return df

def map_rrna_elements(summary_file, gene_table_file, output_file):
    # Load the summary table
    summary_df = pd.read_csv(summary_file, sep="\t", header=None, names=[
        "species", "gene_id", "contig", "start", "end", "identity", "length"
    ])
    summary_df = normalize_chrom_names(summary_df, "contig")
    # Load the gene-to-rRNA-element mapping
    mapping_df = pd.read_csv(gene_table_file, sep="\t")
    gene_to_rrna = dict(zip(mapping_df["gene_id"], mapping_df["rRNA_element"]))

    # Map gene_id to rRNA_element
    summary_df["rRNA_element"] = summary_df["gene_id"].map(gene_to_rrna).fillna("NA")

    # Save the output
    summary_df.to_csv(output_file, sep="\t", index=False)
    print(f"Mapped rRNA elements written to: {output_file}")

def summarize_rrna_by_contig(mapped_summary_file, output_file):
    df = pd.read_csv(mapped_summary_file, sep="\t", header=0)

    # Deduplicate to one row per (species, contig, rRNA_element)
    dedup_df = df[["species", "contig", "rRNA_element"]].drop_duplicates()

    # Pivot table: species as rows, rRNA_element as columns, values = list of contigs
    pivot_df = (
        dedup_df.groupby(["species", "rRNA_element"])["contig"]
        .unique()
        .apply(lambda contigs: "[" + ", ".join(contigs) + "]")
        .unstack(fill_value="[]")
        .reindex(columns=[
            "17S_rRNA-1", "17S_rRNA-2", "5S_rRNA-1", "5S_rRNA-2",
            "5.8S_rRNA-1", "5.8S_rRNA-2", "26S_rRNA-1", "26S_rRNA-2"
        ])
        .reset_index()
    )

    pivot_df.to_csv(output_file, sep="\t", index=False)
    print(f"rRNA summary per species written to: {output_file}")


# === Command-line Interface ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map gene IDs to rRNA elements and summarize BLAST hits.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: map
    map_parser = subparsers.add_parser("map", help="Map gene_id to rRNA element")
    map_parser.add_argument("-s", "--summary", required=True, help="Input summary file")
    map_parser.add_argument("-g", "--gene_table", required=True, help="Gene table with gene_id to rRNA_element")
    map_parser.add_argument("-o", "--output", required=True, help="Output file for mapped summary")

    # Subcommand: summarize
    sum_parser = subparsers.add_parser("summarize", help="Summarize mapped rRNA hits by species")
    sum_parser.add_argument("-i", "--input", required=True, help="Mapped summary file with rRNA_element column")
    sum_parser.add_argument("-o", "--output", required=True, help="Output summary file (per species)")

    args = parser.parse_args()

    if args.command == "map":
        map_rrna_elements(args.summary, args.gene_table, args.output)
    elif args.command == "summarize":
        summarize_rrna_by_contig(args.input, args.output)
