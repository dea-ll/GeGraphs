# GeGraphs

**Gene-centric visualization and structural variant exploration**

GeGraphs is an R package developed to facilitate the visual inspection and exploration of structural variants (SVs) by combining **GENCODE gene annotations**, **BAM sequencing coverage**, and **SV calls from BED/VCF files** in one framework.

Developed during a Master's project in **Bioinformatics and Data Analysis in Biology at the University of Geneva (UNIGE)**.

> **Version:** 0.2.9  
> **Status:** research/development package. GeGraphs supports manual review and interpretation; it is not a standalone clinically validated diagnostic tool.

## Main functions

| Function | Purpose |
|---|---|
| `plot_gene()` | Gene-centred visualization with coverage and SV tracks |
| `plot_region()` | Visualization of a user-defined genomic interval |
| `get_genes_under_sv()` | Identify genes overlapped by structural variants |
| `get_gene_context()` | Explore genomic context around SVs |
| `get_gene_context_for_gene()` | Retrieve neighbouring genes |
| `panels_overlap_df()` | Intersect a gene panel with SV calls and optionally generate plots |

## Example outputs

### Region / trio comparison

![Region/trio example](docs/images/trio_example.png)

### SV caller comparison

![SV caller comparison](docs/images/tools_example.png)

> Example figures must be anonymised before publication.

## Installation

From unpacked package source:

```r
devtools::install("/path/to/GeGraphs")
library(GeGraphs)
packageVersion("GeGraphs")
```

From the packaged source archive:

```bash
R CMD INSTALL GeGraphs_0.2.9.tar.gz
```

The development environment is provided in:

```text
environment/dev_gegraphs_env.yml
```

and can be recreated with:

```bash
conda env create -f environment/dev_gegraphs_env.yml
conda activate dev_gegraphs_env
```

## Quick start

```r
library(GeGraphs)

sample_bam <- "path/to/sample.bam"

sv_files <- list(
  Delly = "path/to/sample.delly.vcf",
  Dysgu = "path/to/sample.dysgu.vcf"
)

plot_gene(
  gene = "ERGIC1",
  bam_file = sample_bam,
  vcf_list = sv_files,
  genome = "hg19",
  flank = 5000,
  show_gene = TRUE,
  show_cov = TRUE,
  show_sv = TRUE,
  bin_size = 2000,
  out = "ERGIC1_review.pdf"
)
```

See [`examples/quickstart.R`](examples/quickstart.R) for a longer example.

## Documentation

The complete tutorial is available here:

[`docs/tutorial/GeGraphs_v0.2.9_Tutorial_FINAL_UNIGE.pdf`](docs/tutorial/GeGraphs_v0.2.9_Tutorial_FINAL_UNIGE.pdf)

It documents installation, inputs, `plot_gene()`, `plot_region()`, reading the output, coverage bin size, genes under SVs, panel overlap, troubleshooting and session information.

## Repository structure

```text
GeGraphs/
├── README.md
├── DESCRIPTION
├── NAMESPACE
├── NEWS.md
├── CITATION.cff
├── .gitignore
├── R/
├── man/
├── tests/
├── vignettes/
├── examples/
│   └── quickstart.R
├── environment/
│   └── dev_gegraphs_env.yml
├── docs/
│   ├── images/
│   │   ├── trio_example.png
│   │   └── tools_example.png
│   └── tutorial/
│       └── GeGraphs_v0.2.9_Tutorial_FINAL_UNIGE.pdf
├── inst/
│   └── extdata/
└── release/
    └── README.md
```

## Scope

GeGraphs is intended for visualization, exploration and prioritization support. It does not call SVs, determine pathogenicity, replace orthogonal validation, or replace expert biological/clinical interpretation.

## Confidentiality

Do not upload patient-derived BAM/VCF/BED/FASTQ data, sample identifiers, clinical information, institutional usernames/emails, credentials, or internal server paths such as `/scratch/...`, `/data/...` or `/home/...`.

## Citation

Llugiqi, D. (2026). *GeGraphs: Gene-Centric Genomic Visualization Tools*. R package version 0.2.9.

## Author

**Dea Llugiqi**  
Master in Bioinformatics and Data Analysis in Biology (BIADB)  
University of Geneva (UNIGE)  
2026

## Licence

No open-source licence is currently assigned. Ownership and redistribution conditions should be confirmed before public reuse.
