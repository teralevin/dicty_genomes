#!/usr/bin/env python3
"""
03_filter_telomeres.py — Filter BED repeats to those near contig ends

Usage:
  03_filter_telomeres.py <BED> <LENGTHS_CSV> [-bp TERMINAL_BP | -f FRACTION] -o OUTPUT

Positional arguments:
  BED             Input BED file (tab-delim columns: chrom, start, end)
  LENGTHS_CSV     CSV of contig lengths (cols “contig” and “size”)

Options (one required):
  -bp, --terminal-bp INT    Distance in bp from contig ends
  -f, --fraction FLOAT      Fraction of contig length (0–1)
  
-o, --output FILE         Output BED file

Examples:
  # Within 5 kb of contig ends:
  python 03_filter_telomeres.py assembly.repeats.bed contig_lengths.csv \
        -bp 5000 \
        -o terminal_5kb.bed

  # Within the last 2 % of each contig’s length:
  python 03_filter_telomeres.py assembly.repeats.bed contig_lengths.csv \
        -f 0.02 \
        -o terminal_2pct.bed
"""
import pandas as pd
import argparse
import os

def normalize_chrom_names(df, column):
    """
    Normalize chromosome names by:
    - Removing 'v1_f1_' or 'v4_f1_'
    - Removing any '_p<digits>' suffix
    - Fixing duplicated prefixes like 'KGL29Av1_KGL29Av1_contig_50' → 'KGL29Av1_contig_50'

    Args:
        df: pandas dataframe with chromosome names
        column: column name with chromosome names

    Returns:
        df: pandas dataframe with normalized chromosome names   
    """
    df[column] = (
        df[column]
        .str.replace(r"v[14]_f1_", "", regex=True)                          # Remove version suffix
        .str.replace(r"_p\d+$", "", regex=True)                            # Remove _p## suffix
        .str.replace(r"^([A-Za-z0-9]+)_\1_(contig|scaffold|chr)_", r"\1_\2_", regex=True)  # Collapse duplicated prefix
    )
    return df

def load_contig_lengths(csv_path):
    """
    Load contig lengths CSV into a DataFrame.
    Expects columns: contig, size, (others ignored)
    """
    df = pd.read_csv(csv_path)
    df = normalize_chrom_names(df, 'contig')
    print(f"Loaded {len(df)} contigs from {csv_path}")
    return df[['contig', 'size']].rename(columns={'contig': 'chrom', 'size': 'contig_length'})


def filter_terminal_repeats(bed_path, lengths_df, output,terminal_bp=None, fraction=None):
    """
    Filter repeats in the BED file to those within terminal regions of contigs.

    Args:
        bed_path (str): path to BED file with at least three columns: chrom, start, end
        lengths_df (pd.DataFrame): DataFrame with 'chrom' and 'contig_length' columns
        terminal_bp (int): absolute bp distance from ends
        fraction (float): fraction of contig length (0-1) defining terminal region

    Returns:
        pd.DataFrame: filtered repeats
    """
    # Load bed
    cols = ['chrom', 'start', 'end']
    bed_df = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0,1,2], names=cols)
    bed_df = normalize_chrom_names(bed_df, 'chrom')
    # Merge lengths
    df = bed_df.merge(lengths_df, on='chrom', how='left')
    missing = df['contig_length'].isna().sum()
    if missing > 0:
        print(f"Warning: {missing} repeats have no matching contig length and will be dropped.")
    df = df.dropna(subset=['contig_length'])
    # Compute windows
    if terminal_bp is not None:
        df['window'] = terminal_bp
    elif fraction is not None:
        df['window'] = (df['contig_length'] * fraction).astype(int)
    else:
        raise ValueError("Either terminal_bp or fraction must be specified")
    # Distance to ends
    df['dist_start'] = df['start']
    df['dist_end'] = df['contig_length'] - df['end']
    # Flag terminal
    df['is_terminal'] = (df['dist_start'] <= df['window']) | (df['dist_end'] <= df['window'])
    # Filter
    filtered = df[df['is_terminal']].copy()
       # Write
    filtered[['chrom', 'start', 'end']].to_csv(
        output,
        sep='\t',
        header=False,
        index=False
    )
    print(f"Filtered {len(filtered)} terminal repeats")
    return filtered


def main():
    parser = argparse.ArgumentParser(
        description="Filter BED repeats to those near contig ends"
    )
    parser.add_argument('bed', help='Input BED file')
    parser.add_argument('lengths_csv', help='CSV of contig lengths')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-bp','--terminal-bp', type=int, help='Distance in bp from end')
    group.add_argument('-f','--fraction', type=float, help='Fraction of contig length')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    args = parser.parse_args()

    lengths_df = load_contig_lengths(args.lengths_csv)
    filtered = filter_terminal_repeats(
        args.bed,
        lengths_df,
        terminal_bp=args.terminal_bp,
        fraction=args.fraction,
        output=args.output
    )
    print(f"Wrote {len(filtered)} terminal repeats to {args.output}")
    
if __name__ == "__main__":
    main()