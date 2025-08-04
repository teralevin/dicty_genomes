#!/usr/bin/env python3

"""
map_tblastn_results.py — Parse TBLASTN output, summarize proteins per contig, and map to gene names.

USAGE:
    map_tblastn_results.py -i INPUT_TBLASTN -o OUTPUT_MAPPING
        [-s SUMMARY_TSV -c CONTIG_LIST]
        [-m MAPPING_CSV -g OUTPUT_WITH_GENES]

FUNCTIONS:
    1. parse_tblastn:
       Maps protein_id → normalized contigs from raw TBLASTN output.

    2. summarize_proteins_per_contig:
       Summarizes protein counts per contig, including zero-hit contigs.

    3. map_proteins_to_genes:
       Annotates contigs with gene names using a protein_id → gene mapping.

OPTIONS:
    -i, --input          Raw TBLASTN output (format 6). **Required**.
    -o, --output         TSV mapping protein_id → contigs. **Required**.
    -s, --summary        Optional: Summary TSV of protein counts per contig.
    -c, --contig_list    TXT file of contigs to include in the summary.
    -m, --mapping_file   CSV mapping protein_id → gene names.
    -g, --output_with_genes TSV with gene names added to the summary.

EXAMPLES:
    # Parse TBLASTN to mapping:
    map_tblastn_results.py -i raw.tblastn.outfmt6 -o protein_to_contigs.tsv

    # Summarize per contig:
    map_tblastn_results.py -i raw.tblastn.outfmt6 -o protein_to_contigs.tsv -s contig_summary.tsv -c my_contigs.txt

    # Annotate summary with gene names:
    map_tblastn_results.py -i raw.tblastn.outfmt6 -o protein_to_contigs.tsv -s contig_summary.tsv -c my_contigs.txt -m protein_name_map.csv -g contig_summary_with_genes.tsv
"""

import argparse
import pandas as pd
from collections import defaultdict


# === Chromosome name normalization ===
def normalize_chrom_names(df, column):
    df[column] = (
        df[column]
        .str.replace(r"v[14]_f1_", "", regex=True)
        .str.replace(r"_p\d+$", "", regex=True)
        .str.replace(r"^([A-Za-z0-9]+)_\1_(contig|scaffold|chr)_", r"\1_\2_", regex=True)
    )
    return df

def parse_tblastn(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t", header=None, usecols=[0, 1],
                     names=["protein_id", "contig"])
    # Normalize contig names
    df = normalize_chrom_names(df, "contig")
    # Group contigs by protein ID
    mapping = defaultdict(set)
    for _, row in df.iterrows():
        mapping[row["protein_id"]].add(row["contig"])

    # Convert to DataFrame
    records = [{"protein_id": prot, "contigs": ", ".join(sorted(list(contigs)))} for prot, contigs in mapping.items()]
    result_df = pd.DataFrame(records).sort_values("protein_id")

    # Save output
    result_df.to_csv(output_file, sep="\t", index=False)
    print(f"Saved protein-to-contig mapping to: {output_file}")
import pandas as pd

def summarize_proteins_per_contig(blast_tsv, contig_list_txt, output_summary_tsv):
    """
    For each contig from the provided list, summarize how many mitochondrial proteins matched it
    using parsed tblastn output. Ensure all input contigs appear, even with 0 hits.
    """
    # Load TBLASTN results
    df = pd.read_csv(blast_tsv, sep="\t")
    if "protein_id" not in df.columns or "contigs" not in df.columns:
        raise ValueError("Parsed blast file must contain 'protein_id' and 'contig' columns.")

    # --- Expand comma-separated contigs ---
    df["contigs"] = df["contigs"].str.split(",").apply(lambda x: [c.strip() for c in x if c.strip()])
    df = df.explode("contigs").reset_index(drop=True)

    # Load and normalize contig list
    contig_df = pd.read_csv(contig_list_txt, header=None, names=["contig"])
    contig_df = normalize_chrom_names(contig_df, "contig")

    # --- Group tblastn hits: contig → list of proteins ---
    grouped = (
        df.groupby("contigs")["protein_id"]
        .agg(protein_count=lambda x: len(set(x)),
             proteins=lambda x: ",".join(sorted(set(x))))
        .reset_index()
    )
    grouped = grouped.rename(columns={"contigs": "contig"})

    # --- Merge with all expected contigs ---
    summary = pd.merge(contig_df, grouped, on="contig", how="left")

    # Fill contigs with no hits
    summary["protein_count"] = summary["protein_count"].fillna(0).astype(int)
    summary["proteins"] = summary["proteins"].fillna("")

    # Save
    summary.to_csv(output_summary_tsv, sep="\t", index=False)
    print(f" Complete summary written to {output_summary_tsv} with {len(summary)} contigs.")

def map_proteins_to_genes(protein_summary_tsv, mapping_file, output_with_genes_tsv):
    """
    Maps protein IDs in summary to gene names using a mapping file.
    """
    # Load input
    df = pd.read_csv(protein_summary_tsv, sep="\t")
    mapping = pd.read_csv(mapping_file)  # columns: protein_id, protein_name
    mapping_dict = dict(zip(mapping["protein_id"], mapping["protein_name"]))

    # Map proteins to gene names
    def map_genes(protein_str):
        ids = protein_str.split(",")
        return ",".join([mapping_dict.get(pid, "NA") for pid in ids])

    df["protein_names"] = df["proteins"].apply(map_genes)

    #drop the protein_id column
    df = df.drop(columns=["proteins"])
    # Save output
    df.to_csv(output_with_genes_tsv, sep="\t", index=False)
    print(f"Gene name-annotated summary written to {output_with_genes_tsv}")

# === CLI ===
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse TBLASTN results and summarize contigs per mitochondrial protein.")
    parser.add_argument("-i", "--input", required=True, help="TBLASTN input file (outfmt 6)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file mapping protein -> contigs")
    parser.add_argument("-s", "--summary", help="Output summary TSV with protein counts per contig")
    parser.add_argument("-c", "--contig_list", help="List of contigs to filter the summary")
    parser.add_argument("-m", "--mapping_file", help="Mapping file for protein IDs to gene names")
    parser.add_argument("-g", "--output_with_genes", help="Output summary with gene names")
    args = parser.parse_args()

    parse_tblastn(args.input, args.output)
    if args.summary and args.contig_list:
        summarize_proteins_per_contig(args.output, args.contig_list, args.summary)
    if args.mapping_file and args.output_with_genes:
        map_proteins_to_genes(args.summary, args.mapping_file, args.output_with_genes)
