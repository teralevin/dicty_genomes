#!/usr/bin/env python3
"""
03_tgr_CSE_intersection.py

Find intervals in file 1 (CRE) that overlap intervals in file 2 (TGR), and report the exact overlap coordinates.

Usage:
    python 03_tgr_CSE_intersection.py <cre_file> <tgr_file> <output.tsv>
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


def read_cre_file(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["contig", "start", "end", "aqi"],
        dtype={"contig": str, "start": int, "end": int, "aqi": float},
    )
    # drop trivial aqi column if it's always 1
    if "aqi" in df.columns and (df["aqi"] == 1).all():
        df = df.drop(columns=["aqi"])
    # fix any reversed coords
    mask = df["start"] > df["end"]
    if mask.any():
        df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values
    return normalize_chrom_names(df, "contig")


def read_tgr_bedgraph(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        skiprows=1,
        names=["contig", "start", "end", "value"],
        dtype={"contig": str, "start": int, "end": int, "value": float},
    )
    # drop trivial value column if it's always 1
    if "value" in df.columns and (df["value"] == 1).all():
        df = df.drop(columns=["value"])
    # fix any reversed coords
    mask = df["start"] > df["end"]
    if mask.any():
        df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values
    return normalize_chrom_names(df, "contig")


def find_overlap(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a DataFrame with columns:
      contig, cre_start, cre_end, aqi,
      tgr_start, tgr_end, overlap_start, overlap_end
    for every interval in df1 that overlaps any interval in df2.
    """
    merged = pd.merge(
        df1, df2,
        on="contig",
        how="inner",
        suffixes=("_cre", "_tgr")
    )
    ov_mask = (
        (merged["end_cre"]   > merged["start_tgr"]) &
        (merged["start_cre"] < merged["end_tgr"])
    )
    ov = merged.loc[ov_mask].copy()
    ov["overlap_start"] = ov[["start_cre", "start_tgr"]].max(axis=1)
    ov["overlap_end"]   = ov[["end_cre",   "end_tgr"  ]].min(axis=1)

    if "aqi" not in ov.columns:
        ov["aqi"] = pd.NA

    return (
        ov[[
            "contig",
            "start_cre", "end_cre", "aqi",
            "start_tgr", "end_tgr",
            "overlap_start", "overlap_end"
        ]]
        .rename(columns={
            "start_cre": "cre_start",
            "end_cre"  : "cre_end",
            "start_tgr": "tgr_start",
            "end_tgr"  : "tgr_end"
        })
        .reset_index(drop=True)
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Find intervals in CRE file that overlap intervals in TGR file, "
            "and output the overlap coordinates."
        )
    )
    parser.add_argument(
        "cre_file",
        help="Path to CRE file (bedGraph-style) with columns: contig, start, end[, aqi]"
    )
    parser.add_argument(
        "tgr_file",
        help="Path to TGR file (bedGraph-style) with columns: contig, start, end[, value]"
    )
    parser.add_argument(
        "output",
        help="Path to write the TSV of overlapping intervals"
    )
    args = parser.parse_args()

    df1 = read_cre_file(args.cre_file)
    df2 = read_tgr_bedgraph(args.tgr_file)
    overlaps = find_overlap(df1, df2)

    overlaps.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(overlaps)} overlapping intervals to {args.output}")


if __name__ == "__main__":
    main()
