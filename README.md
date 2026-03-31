# Cross-Species Coexpression Comparison Pipeline

Nextflow pipeline for comparing mouse vs human transcription regulator (TR) coexpression networks.

## Overview

This pipeline takes pre-computed per-species coexpression matrices (from `mouse_coexpression-nf` or equivalent) and performs cross-species comparison in a shared gene space of ~16,981 high-confidence one-to-one human–mouse orthologs (DIOPT).

### Workflow

```
┌─────────────────────────┐     ┌─────────────────────────┐
│  Mouse coexpr matrices  │     │  Human coexpr matrices  │
│  (*_corAggSparse.h5ad)  │     │  (*_corAggSparse.h5ad)  │
└───────────┬─────────────┘     └───────────┬─────────────┘
            │                               │
            ▼                               ▼
   ┌────────────────────────────────────────────────┐
   │  1. EXTRACT_TR_PROFILES                        │
   │     Per dataset: extract TR profiles in        │
   │     ortholog gene space (DIOPT 1-to-1)         │
   └────────────────────┬───────────────────────────┘
                        │
                        ▼
   ┌────────────────────────────────────────────────┐
   │  2. AGGREGATE_PROFILES                         │
   │     Per species + cell type + TR:              │
   │     - Impute missing to median                 │
   │     - Rank-standardize per dataset             │
   │     - Average across datasets                  │
   │     => One aggregate profile per TR per species│
   └────────────────────┬───────────────────────────┘
                        │
                        ▼
   ┌────────────────────────────────────────────────┐
   │  3. COMPARE_SPECIES                            │
   │     For each TR with ortholog & ≥5 datasets:   │
   │     - Spearman correlation                     │
   │     - Top-K overlap                            │
   │     - Bottom-K overlap                         │
   │     - Ortholog retrieval score (bidirectional)  │
   └────────────────────────────────────────────────┘
```

## Inputs

| Parameter | Description |
|-----------|-------------|
| `mouse_coexpr_dir` | Directory with mouse coexpression h5ad files (`*_corAggSparse.h5ad`) |
| `human_coexpr_dir` | Directory with human coexpression h5ad files (`*_corAggSparse.h5ad`) |
| `ortholog_file` | DIOPT one-to-one ortholog table (TSV with `Symbol_hg`, `Symbol_mm`, `ID`) |
| `mouse_tf_file` | Mouse TR list from AnimalTFDB |
| `human_tf_file` | Human TR list from AnimalTFDB |

Input h5ad filenames are expected in the format: `{CellType}_{StudyID}_corAggSparse.h5ad`

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `top_k` | 200 | K for Top-K and Bottom-K overlap analysis |
| `min_datasets` | 5 | Minimum datasets in both species for a TR to be compared |

## Usage

```bash
nextflow run main.nf \
    --mouse_coexpr_dir /path/to/mouse/coexpr \
    --human_coexpr_dir /path/to/human/coexpr \
    -profile slurm
```

## Output

```
results_YYYYMMDD_HHMMSS/
├── 01_tr_profiles/          # Per-dataset TR profiles in ortholog space
│   ├── mouse/               # Mouse TR profiles (parquet)
│   └── human/               # Human TR profiles (parquet)
├── 02_aggregated/           # Aggregated profiles per TR per species
│   ├── mouse/               # One file per cell type
│   └── human/
├── 03_comparison/           # Cross-species comparison results
│   ├── {CellType}_cross_species_results.csv   # Per-TR metrics
│   └── {CellType}_summary.csv                 # Summary statistics
└── pipeline_info/           # Nextflow execution reports
```

### Output columns (`*_cross_species_results.csv`)

| Column | Description |
|--------|-------------|
| `human_tr` / `mouse_tr` | Orthologous TR pair |
| `spearman_rho` | Spearman correlation of aggregate profiles |
| `top200_overlap` | Fraction of top-200 coexpressed genes shared |
| `bottom200_overlap` | Fraction of bottom-200 anti-coexpressed genes shared |
| `retrieval_human_in_mouse` | Quantile of mouse ortholog among all mouse TRs ranked by similarity to human TR profile |
| `retrieval_mouse_in_human` | Reciprocal retrieval score |

## Metrics

- **Spearman correlation**: Compares every orthologous gene's rank in human vs mouse aggregate profiles
- **Top-K overlap**: The K most positively coexpressed genes in human vs mouse
- **Bottom-K overlap**: The K most anti-coexpressed genes in each species
- **Ortholog retrieval score**: A score of 1 means the TR's ortholog shared more top coexpressed partners than any other TR in the other species. Computed bidirectionally.
