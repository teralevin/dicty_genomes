#!/usr/bin/env python3
r"""Extract the Tgr locus data from GFF and FASTA files to run BLAST.

Usage:
    python 01_tgr_extraction.py \
        --tgr-file path/to/tgr.bedgraph.txt \
        --gff-file path/to/all_species_w_ref.gff \
        [--slide-window 100000] \
        [--sequences-dir ./sequences] \
        [--output-dir ./tgr_locus_sequences]
"""

import argparse
import os
from collections import defaultdict
import pandas as pd


def normalize_chrom_names(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """
    Normalize chromosome names by:
      - Removing 'v1_f1_' or 'v4_f1_'
      - Removing any '_p<digits>' suffix
      - Fixing duplicated prefixes like 'KGL29Av1_KGL29Av1_contig_50'
    """
    df[column] = (
        df[column]
        .str.replace(r"v[14]_f1_", "", regex=True)
        .str.replace(r"_p\d+$", "", regex=True)
        .str.replace(
            r"^([A-Za-z0-9]+)_\1_(contig|scaffold|chr)_",
            r"\1_\2_",
            regex=True,
        )
    )
    return df


def generate_tgr_locus_dict(tgr_file: str) -> dict[str, tuple[int, int]]:
    """
    Parses a BEDGRAPH‐style file to get Tgr locus coordinates per chromosome.
    """
    tgr_df = pd.read_csv(
        tgr_file,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "value"],
        skiprows=[0],
    )
    clean_tgr = normalize_chrom_names(tgr_df, "chrom")

    tgr_dict: dict[str, Tuple[int, int]] = {}
    for _, row in clean_tgr.iterrows():
        chrom = row["chrom"].strip()
        start, end = int(row["start"]), int(row["end"])
        tgr_dict[chrom] = (start, end)
        print(f"[TGR] {chrom}: {start}-{end}")
    return tgr_dict


def extract_tgr_locus(
    gff_file: str,
    tgr_locus_dictionary: dict[str, tuple[int, int]],
    slide_window: int,
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """
    From a GFF file, pull out all genes whose coordinates fall within each
    Tgr locus ± slide_window.
    """
    gff_df = pd.read_csv(
        gff_file,
        sep="\t",
        header=None,
        names=["chrom", "gene", "start", "end"],
    )
    clean_gff = normalize_chrom_names(gff_df, "chrom")

    result: dict[str, list[str]] = {}
    filtered_rows = []

    for chrom, (raw_start, raw_end) in tgr_locus_dictionary.items():
        chrom = chrom.strip()
        start, end = sorted((raw_start, raw_end))
        start = max(0, start - slide_window)
        end += slide_window

        df_sub = clean_gff[
            (clean_gff["chrom"].str.strip() == chrom)
            & (clean_gff["start"] >= start)
            & (clean_gff["end"] <= end)
        ]
        genes = df_sub["gene"].drop_duplicates().tolist()
        result[chrom] = genes
        filtered_rows.append(df_sub)
        print(f"[GFF] {chrom} → {len(genes)} genes in {start}-{end}")

    combined = pd.concat(filtered_rows, ignore_index=True)
    return result, combined


def extract_sequences(
    grouped: dict[str, dict[str, tuple[int, int]]],
    sequences_dir: str,
    out_dir: str,
):
    os.makedirs(out_dir, exist_ok=True)

    for prefix, chrom_map in grouped.items():
        fasta_name = prefix.lower()
        # special case for KGL29Av1 → kgl29a.faa
        if prefix == "KGL29Av1":
            fasta_name = fasta_name[:-2]
        fpath = os.path.join(sequences_dir, f"{fasta_name}.faa")
        if not os.path.exists(fpath):
            print(f"[WARN] FASTA not found: {fpath}")
            continue

        targets = {g for genes in chrom_map.values() for g in genes}
        out_faa = os.path.join(out_dir, f"{fasta_name}_tgr_locus.faa")

        with open(fpath) as inp, open(out_faa, "w") as outp:
            keep = False
            for line in inp:
                if line.startswith(">"):
                    fid = line[1:].split()[0].split("-")[0]
                    keep = fid in targets
                if keep:
                    outp.write(line)

        nseq = sum(1 for _ in open(out_faa) if _.startswith(">"))
        print(f"[OUT] {prefix}: {nseq} sequences → {out_faa}")


def main():
    p = argparse.ArgumentParser(
        description="Extract Tgr locus genes and sequences from GFF+FASTA"
    )
    p.add_argument(
        "--tgr-file",
        required=True,
        help="Bedgraph‐style file with Tgr loci (chrom, start, end, ...)",
    )
    p.add_argument(
        "--gff-file",
        required=True,
        help="GFF file containing all species annotations",
    )
    p.add_argument(
        "--slide-window",
        type=int,
        default=100000,
        help="Number of bp to extend each side of the locus",
    )
    p.add_argument(
        "--sequences-dir",
        default="./sequences",
        help="Directory containing per-species .faa FASTA files",
    )
    p.add_argument(
        "--output-dir",
        default="./tgr_locus_sequences",
        help="Where to write extracted GFFs and FASTAs",
    )

    args = p.parse_args()

    tgr_dict = generate_tgr_locus_dict(args.tgr_file)
    tgr_genes, combined_gff = extract_tgr_locus(
        args.gff_file, tgr_dict, args.slide_window
    )

    # write combined GFF
    os.makedirs(args.output_dir, exist_ok=True)
    gff_out = os.path.join(args.output_dir, "tgr_locus.gff")
    combined_gff.to_csv(gff_out, sep="\t", header=False, index=False)
    print(f"[SAVED] Combined GFF → {gff_out}")

    # group by prefix logic
    grouped = defaultdict(dict)
    for chrom, genes in tgr_genes.items():
        if "-" in chrom:
            pref = chrom.split("-")[0]
        elif "_" in chrom:
            pref = chrom.split("_")[0]
        elif "CM" in chrom:
            pref = "tnsc14"
        elif "chr" in chrom:
            pref = "ax4"
        else:
            pref = chrom
        grouped[pref][chrom] = genes

    extract_sequences(grouped, args.sequences_dir, args.output_dir)


if __name__ == "__main__":
    main()
