#!/usr/bin/env python3
"""
02_tge_TE_interesection.py

Find overlaps between TGR loci (file1) and TE hits (file2), and optionally
annotate those TE regions with overlapping genes from a GFF file.

Usage:
    python 02_tge_TE_interesection.py \
        <tgr_file> <te_file> <overlap_output.tsv> \
        [--gff-file <gff_file> \
         [--gene-output <genes_output.tsv>] \
         [--merge-gap N] [--no-extend] \
         [--columns start_col,end_col]]
"""

import argparse
import pandas as pd


def normalize_chrom_names(df: pd.DataFrame, column: str) -> pd.DataFrame:
    df[column] = (
        df[column]
        .astype(str)
        .str.replace(r"v[14]_f1_", "", regex=True)
        .str.replace(r"_p\d+$", "", regex=True)
        .str.replace(
            r"^([A-Za-z0-9]+)_\1_(contig|scaffold|chr)_",
            r"\1_\2_",
            regex=True
        )
    )
    return df


def read_tgr_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        skiprows=1,
        names=["chrom", "start", "end", "value"],
        dtype={"chrom": str, "start": int, "end": int, "value": float},
    )
    if "value" in df.columns and (df["value"] == 1).all():
        df = df.drop(columns=["value"])
    mask = df["start"] > df["end"]
    if mask.any():
        df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values
    return normalize_chrom_names(df, "chrom")


def read_te_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["chrom", "start", "end", "TE_name", "score", "direction"],
        dtype={"chrom": str, "start": int, "end": int, "TE_name": str},
    )
    mask = df["start"] > df["end"]
    if mask.any():
        df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values
    return normalize_chrom_names(df, "chrom")


def find_overlaps(
    tgr_df: pd.DataFrame,
    te_df: pd.DataFrame
) -> pd.DataFrame:
    overlaps: list[list] = []
    for _, row in tgr_df.iterrows():
        chrom = row["chrom"].strip()
        tgr_start, tgr_end = sorted((row["start"], row["end"]))
        subset = te_df[
            (te_df["chrom"].str.strip() == chrom) &
            (te_df["start"] <= tgr_end) &
            (te_df["end"]   >= tgr_start)
        ]
        for _, te in subset.iterrows():
            overlaps.append([
                chrom,
                tgr_start,
                tgr_end,
                te["start"],
                te["end"],
                te["TE_name"],
                te["direction"],
            ])
    return pd.DataFrame(
        overlaps,
        columns=[
            "chrom",
            "tgrb_start",
            "tgrb_end",
            "te_start",
            "te_end",
            "TE_name",
            "direction",
        ],
    )


def find_genes_for_te_regions(
    te_csv: str,
    gff_file: str,
    columns: tuple[str, str]| None,
    merge_gap: int,
    extend: bool
) -> pd.DataFrame:
    gff_df = pd.read_csv(
        gff_file,
        sep="\t",
        header=None,
        names=["chrom", "gene", "start", "end"],
        dtype={"chrom": str},
    )
    gff_df = normalize_chrom_names(gff_df, "chrom")

    te_df = pd.read_csv(te_csv, sep="\t", dtype={"chrom": str})
    if columns is not None and len(columns) != 2:
        raise ValueError("`--columns` must be start_col,end_col")
    start_col, end_col = columns or ("te_start", "te_end")

    genes_list = []
    for _, row in te_df.iterrows():
        chrom = row["chrom"]
        te_start = int(row[start_col])
        te_end   = int(row[end_col])
        name     = row["TE_name"]

        subset = gff_df[
            (gff_df["chrom"] == chrom) &
            (gff_df["start"] <= te_end) &
            (gff_df["end"]   >= te_start)
        ]
        if subset.empty and extend:
            peers = te_df[
                (te_df["chrom"]    == chrom) &
                (te_df["TE_name"]  == name)
            ]
            if len(peers) > 1:
                ext_start = te_start - merge_gap
                ext_end   = te_end   + merge_gap
                subset = gff_df[
                    (gff_df["chrom"] == chrom) &
                    (gff_df["start"] <= ext_end) &
                    (gff_df["end"]   >= ext_start)
                ]
        genes_list.append(subset["gene"].drop_duplicates().tolist())

    te_df["overlapping_genes"] = genes_list
    return te_df


def main():
    parser = argparse.ArgumentParser(
        description="02_tge_TE_interesection.py — find TGR/TE overlaps and optionally annotate with genes."
    )
    parser.add_argument("tgr_file", help="TGR bedGraph file (skip header)")
    parser.add_argument("te_file",  help="TE hits bed file")
    parser.add_argument("overlap_out", help="TSV to write TGR×TE overlaps")
    parser.add_argument(
        "--gff-file",
        help="GFF file to extract overlapping genes; if omitted, gene extraction is skipped"
    )
    parser.add_argument(
        "--gene-output",
        help="Output TSV for TE×genes (default: <overlap_out>_genes.tsv)"
    )
    parser.add_argument(
        "--merge-gap",
        type=int,
        default=1000,
        help="bp to extend TE region if no direct overlap (default: 1000)"
    )
    parser.add_argument(
        "--no-extend",
        action="store_true",
        help="Do not extend TE intervals when no direct overlap"
    )
    parser.add_argument(
        "--columns",
        help="Comma-separated TE interval columns (default: te_start,te_end)"
    )

    args = parser.parse_args()

    # 1) read & normalize
    tgr_df = read_tgr_file(args.tgr_file)
    te_df  = read_te_file(args.te_file)

    # 2) find & save overlaps
    overlap_df = find_overlaps(tgr_df, te_df)
    overlap_df.to_csv(args.overlap_out, sep="\t", index=False)
    print(f"Saved {len(overlap_df)} overlaps → {args.overlap_out}")

    # 3) optionally extract genes
    if args.gff_file:
        cols = tuple(args.columns.split(",", 1)) if args.columns else None
        gene_df = find_genes_for_te_regions(
            args.overlap_out,
            args.gff_file,
            columns=cols,
            merge_gap=args.merge_gap,
            extend=not args.no_extend
        )
        out_genes = args.gene_output or f"{args.overlap_out.rsplit('.',1)[0]}_genes.tsv"
        gene_df.to_csv(out_genes, sep="\t", index=False)
        print(f"Saved {len(gene_df)} TE regions with genes → {out_genes}")


if __name__ == "__main__":
    main()
