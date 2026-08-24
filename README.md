# GeGraphs

**Gene-centric visualization and exploration of structural variants from short-read sequencing data**

GeGraphs is an R package developed to make structural-variant (SV) review easier by bringing together **gene annotations**, **read coverage**, and **SV calls** in a single, interpretable view.

The package was developed in the context of a Master's project in Bioinformatics and Data Analysis in Biology at the University of Geneva, in collaboration with Geneva University Hospitals (HUG).

> **Development status:** GeGraphs v0.2.9 is a research/development package. It is not a standalone clinically validated diagnostic tool.

## What GeGraphs does

GeGraphs integrates three main data sources:

- **GENCODE gene annotations**, including UTR and coding features;
- **BAM alignment coverage** for a proband and optional controls;
- **structural variants** supplied as BED or VCF files.

The package supports two complementary visualization modes:

- `plot_gene()` — gene-centred visualization;
- `plot_region()` — visualization of a user-defined genomic interval.

It also provides utilities to:

- identify genes overlapped by structural variants with `get_genes_under_sv()`;
- inspect neighbouring genes with `get_gene_context()` and `get_gene_context_for_gene()`;
- intersect structural variants with a gene panel using `panels_overlap_df()`.

## Why GeGraphs?

Short-read SV callers can produce many candidate events and different callers may describe the same region differently. GeGraphs was designed as a **visual review and exploration layer** after SV detection.

A typical use case is:

```text
BAM + BED/VCF + genome build
           │
           ▼
       GeGraphs
           │
     ┌─────┼──────────────┐
     ▼     ▼              ▼
 Gene view Region view   SV/gene tables
     │     │              │
     └─────┴──────┬───────┘
                  ▼
        Candidate interpretation
```

The goal is not to replace SV calling or clinical interpretation, but to make candidate review more structured and reproducible.

## Main functions

| Function | Purpose |
|---|---|
| `plot_gene()` | Multi-track figure centred on a gene |
| `plot_region()` | Coverage/SV visualization over an arbitrary genomic interval |
| `get_genes_under_sv()` | Identify genes overlapped by SVs |
| `get_gene_context()` | Explore genomic context around an SV |
| `get_gene_context_for_gene()` | Retrieve neighbouring genes around a selected gene |
| `panels_overlap_df()` | Intersect a clinical/research gene panel with SV calls and optionally produce plots |

### `plot_gene()`

`plot_gene()` can combine:

1. a genomic coordinate axis;
2. GENCODE gene annotation;
3. raw and normalized BAM coverage;
4. optional control BAMs;
5. one or more SV tracks from BED or VCF files.

Minimal example:

```r
library(GeGraphs)

GeGraphs::plot_gene(
  gene = "GENE1",
  bam_file = "path/to/proband.bam",
  genome = "hg38",
  vcf_list = list(
    CallerA = "path/to/callerA.vcf"
  ),
  flank = 5000,
  show_gene = TRUE,
  show_cov = TRUE,
  show_sv = TRUE,
  bin_size = 2000,
  title = "GENE1 structural-variant review"
)
```

Controls can be added with a named vector or list:

```r
controls <- c(
  Control1 = "path/to/control1.bam",
  Control2 = "path/to/control2.bam"
)

GeGraphs::plot_gene(
  gene = "GENE1",
  bam_file = "path/to/proband.bam",
  control_bams = controls,
  genome = "hg38",
  show_gene = TRUE,
  show_cov = TRUE,
  show_sv = FALSE,
  bin_size = 5000
)
```

### `plot_region()`

For a broader interval:

```r
GeGraphs::plot_region(
  chrom = "5",
  start = 172260000,
  end = 172360000,
  bam_file = "path/to/proband.bam",
  genome = "hg38",
  show_cov = TRUE,
  show_sv = FALSE,
  title = "Regional coverage overview"
)
```

### Genes under SVs

```r
genes_under <- GeGraphs::get_genes_under_sv(
  bed_file = "path/to/structural_variants.bed",
  genome = "hg38"
)

head(genes_under)
```

### Gene-panel overlap

```r
panel <- read.delim("path/to/gene_panel.tsv", stringsAsFactors = FALSE)

overlap <- GeGraphs::panels_overlap_df(
  paneldf = panel,
  genecol = "GeneSymbol",
  bedfile = "path/to/structural_variants.bed",
  vcffile = NULL,
  svtype = "DEL",
  genome = "hg38",
  makeplots = FALSE
)

head(overlap)
```

## Installation

### Option 1 — install from a local source package

```bash
R CMD INSTALL GeGraphs_0.2.9.tar.gz
```

Then in R:

```r
library(GeGraphs)
packageVersion("GeGraphs")
```

### Option 2 — reproduce the development environment

A Conda environment definition used during development is provided in:

```text
environment/dev_gegraphs_env.yml
```

Create it with:

```bash
conda env create -f environment/dev_gegraphs_env.yml
conda activate dev_gegraphs_env
```

Then install GeGraphs from the source repository or source archive.

> The Conda file captures the development environment and is intentionally more comprehensive than the minimal package dependencies. The package `DESCRIPTION` file should remain the authoritative source for R package dependencies.

## Repository structure

The final repository should have the standard R-package structure:

```text
GeGraphs/
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── NEWS.md
├── R/
│   └── *.R
├── man/
│   └── *.Rd
├── inst/
│   └── extdata/
├── examples/
│   └── quickstart.R
├── environment/
│   └── dev_gegraphs_env.yml
├── docs/
│   └── ...
├── tests/
│   └── ...
└── vignettes/
    └── ...
```

`R/`, `man/`, `DESCRIPTION`, and `NAMESPACE` must come from the actual GeGraphs source package. This repository scaffold does **not** invent or replace the package implementation.

## Genome annotations and SV handling

GeGraphs v0.2.9 supports `hg19` and `hg38` workflows and uses the corresponding embedded GENCODE annotations. Where available, MANE transcripts are preferred. If a MANE transcript identifier cannot be matched to the GFF annotation, the package can fall back to the available GENCODE annotations for the gene.

BED/VCF parsing utilities focus analyses on standard chromosomes and remove non-standard contigs such as random or unplaced scaffolds.

## Recommended workflow

```text
1. Detect SVs with one or more callers
2. Select a candidate SV / genomic region
3. Load the corresponding BAM (+ BAI)
4. Add BED/VCF calls to GeGraphs
5. Visualize gene or region
6. Compare coverage with controls when available
7. Explore genes under the SV / panel overlap
8. Use the result as supporting evidence for expert review
```

## Reproducibility

This repository is intended to make the package understandable and reproducible for reviewers:

- the R package source documents the implementation;
- the Conda YAML records the development software environment;
- examples show the main public functions;
- package documentation describes arguments and outputs;
- no patient-derived data are required for the repository structure itself.

## Data and confidentiality

**Do not upload patient-derived BAM, VCF, BED, phenotype files, sample lists, internal server paths, usernames, institutional email addresses, or other confidential material to a public repository.**

Examples in the public repository should use synthetic/example file names such as:

```text
proband.bam
control1.bam
sample_sv.vcf
sample_sv.bed
```

If demonstration data are included under `inst/extdata/`, they should be synthetic, public, or explicitly cleared for redistribution.

## Scope and limitations

GeGraphs is intended for **visualization, exploration, and prioritization support**. It does not:

- call structural variants itself;
- establish pathogenicity;
- replace orthogonal validation;
- replace expert biological or clinical interpretation.

Performance can depend on BAM size, genomic window size, coverage binning, and the size of BED/VCF inputs. For very large files, pre-filtering to a chromosome or target region can improve responsiveness.

## Version

Current development version documented here: **0.2.9**.

## Author

**Dea Llugiqi**  
Master's programme in Bioinformatics and Data Analysis in Biology  
University of Geneva  
2026

## Citation

Citation metadata are provided in `CITATION.cff`. A DOI or formal software publication can be added later if the package is archived or published.

## Licence

No open-source licence is assigned in this scaffold.

Before public reuse or redistribution, ownership and licensing conditions for code developed or adapted within the host laboratory should be confirmed.
