# Dictyostelium Genome Annotation, Assembly & tgr Locus Analysis

This repository contains a collection of scripts and workflows to:

1. **Annotate** new _Dictyostelium_ genome assemblies (tandem repeats, telomeres, transposable elements, rRNA, and mtDNA)  
2. **Assemble** and **polish** long‐read genome assemblies  
3. **Quality-control** variant calls and structural errors  
4. **Extract** and **analyze** the _tgrBC_ self/non-self recognition locus, including its intersections with transposable elements and structural errors  

Assembled genomes and annotations have been deposited in NCBI under:

- **BioProject**: [PRJNA1300491](http://www.ncbi.nlm.nih.gov/bioproject/1300491) 
- **Locus tag prefix**: `ACTA71` 

---


## Table of Contents
* [Installation](#installation)
    * [Requirements](#requirements)
    * [Steps](#steps)
* [Repository Structure](#repository-structure)
    * [annotation](#annotation)
    * [assembly](#assembly)
    * [quality-control](#quality-control)
    * [tgr\_locus](#tgr_locus)
* [Citation](#citation)
---
## Installation

### Requirements

To use this repository, ensure the following requirements are met:

- **Operating System**: Linux
- **Python Version**: Python 3.10+ (tested on 3.12)
- **Python Packages**: pandas, numpy
- **External Tools**: SAMtools, Bedtools, Minimap2, BWA, BLAST+, TRF, CRAQ, BCFtools

### Steps

1. Clone this repository:

    ```bash
    git clone https://github.com/MayarMAhmed/dicty_genome.git
    cd dicty_genome
    ```

2. Install Python dependencies (recommended using conda):

    ```bash
    conda create -n dicty_env python=3.10 
    conda activate dicty_env
    conda install numpy pandas
    ```

3. Install external tools via conda or your package manager.

---

## Repository Structure

```text
├── annotation
│   ├── 01_run_trf.sh             # Tandem Repeat Finder wrapper
│   ├── 02_extract_telomere_motif.sh
│   ├── 03_filter_telomeres.py
│   ├── 04_find_TE.sh             # RepeatMasker/BLAST TE discovery
│   ├── 05_generate_TE_bedgraph.py
│   ├── 06_run_rrna_blast.sh
│   ├── 07_map_rrna_summary.py
│   ├── 08_search_mtdna_dna.sh
│   ├── 09_search_mtdna_prot.sh
│   ├── 10_parse_mtdna_prot.py
│   └── contig_mapping.txt  #contig mapping for the AX4 reference species and the TNS-C-14 
├── assembly                     # De novo assembly & polishing workflows
├── quality-control
│   ├── 01_variant_calling.sh    # Bcf variant calling
│   └── 02_run_craq.sh           # CRAQ structural error detection
└── tgr_locus
    ├── 01_tgr_extraction.py     # Extract tgrBC syntenic locus from GFF/bedgraph
    ├── 02_tgr_TE_intersection.py# Intersect TEs with tgr loci, optional gene annotation
    ├── 03_tgr_CSE_intersection.py# Intersect tgr loci with conserved sequence elements
    └── tgrBClocus_wacA2chdB.bedgraph.txt #tgrBC syntenic locus across all species
```

Each subdirectory contains a driver script and helper scripts. You can execute steps individually or chain them into a pipeline.

---
## Citation
If you use these scripts or the assemblies/annotations generated with them, please cite:

> Holland M., Ahmed M., Young J.M., McFadyen S., Drurey J.R., Ostrowski E.A., Levin T.C. (2025)
> Hypermutable hotspot enables the rapid evolution of self/non-self recognition genes in *Dictyostelium*.
> *bioRxiv*. [https://doi.org/10.1101/2025.08.01.668227](https://doi.org/10.1101/2025.08.01.668227)

