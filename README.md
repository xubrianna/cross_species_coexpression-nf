# Cross-Species Coexpression Comparison Pipeline

Nextflow (DSL2) pipeline for comparing mouse vs human transcription regulator (TR) coexpression networks.

## Overview

This pipeline takes pre-computed, rank-normalized per-dataset coexpression matrices (`.h5ad`) from both mouse and human, aggregates them into consensus networks per cell type and species using the rank-sum-rerank approach of Crow et al., and evaluates cross-species conservation of TR coexpression profiles in a shared gene space of ~16,981 high-confidence one-to-one human–mouse orthologs (DIOPT v9). TR gene lists are sourced from AnimalTFDB.

### Workflow

```
┌─────────────────────────┐     ┌─────────────────────────┐
│  Mouse coexpr matrices  │     │  Human coexpr matrices  │
│  (*_corAggSparse.h5ad)  │     │  (*_corAggSparse.h5ad)  │
└───────────┬─────────────┘     └───────────┬─────────────┘
            │                               │
            ▼                               ▼
   ┌────────────────────────────────────────────────┐
   │  0. AGGREGATE_NETWORKS                         │
   │     Per species + cell type:                   │
   │     - Align to union gene space (impute 0.5)   │
   │     - Rank each dataset (min-ties)             │
   │     - Sum ranks → re-rank → standardize [0,1]  │
   └──────────┬─────────────────────┬───────────────┘
              │                     │
              ▼                     ▼
   ┌──────────────────────┐  ┌──────────────────────────┐
   │ 1. EXTRACT_AGGREGATE │  │ 2. EXTRACT_TR_PROFILES   │
   │    _PROFILES          │  │    Per dataset: extract  │
   │    From aggregate:    │  │    TR profiles in ortho- │
   │    TR profiles in     │  │    log space (CSV)       │
   │    ortholog space     │  └───────────┬──────────────┘
   │    (JSON)             │              │
   └──────────┬───────────┘              ▼
              │              ┌──────────────────────────┐
              │              │ 3. REPRODUCIBILITY       │
              │              │    Within-species Top-K  │
              │              │    consistency + null    │
              │              └──────────────────────────┘
              ▼
   ┌────────────────────────────────────────────────┐
   │  4. COMPARE_SPECIES                            │
   │     For each TR with ortholog & ≥5 datasets:   │
   │     - Spearman correlation                     │
   │     - Top-K overlap                            │
   │     - Ortholog retrieval score (bidirectional)  │
   │     - Gene-centric conservation                │
   │     - Empirical p-values (1000 permutations)   │
   └────────────────────────────────────────────────┘
```

## Inputs

| Parameter | Description |
|-----------|-------------|
| `mouse_coexpr_dir` | Directory with mouse coexpression h5ad files (`*_corAggSparse.h5ad`) |
| `human_coexpr_dir` | Directory with human coexpression h5ad files (`*_corAggSparse.h5ad`). Optional; omit for mouse-only mode |
| `ortholog_file` | DIOPT one-to-one ortholog table (TSV with `Symbol_hg`, `Symbol_mm`, `ID`) |
| `mouse_tf_file` | Mouse TR list from AnimalTFDB |
| `human_tf_file` | Human TR list from AnimalTFDB |
| `mouse_gene_mapping` | Ensembl ID → gene symbol mapping for mouse (`mouse_gene_meta.csv`) |
| `ribosomal_genes` | List of ribosomal gene symbols (TSV), used for flagging in gene conservation output |

Input h5ad filenames are expected in the format: `{CellType}_{StudyID}_corAggSparse.h5ad`. Human files containing `_AD_` are excluded.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `top_k` | 200 | K for Top-K overlap analysis |
| `min_datasets` | 5 | Minimum datasets in both species for a TR to be compared |
| `n_null` | 1000 | Number of permutation iterations for empirical p-values |

## Usage

```bash
nextflow run main.nf \
    --mouse_coexpr_dir /path/to/mouse/coexpr \
    --human_coexpr_dir /path/to/human/coexpr \
    -profile slurm
```

Mouse-only mode (runs aggregation and reproducibility only):

```bash
nextflow run main.nf \
    --mouse_coexpr_dir /path/to/mouse/coexpr \
    -profile slurm
```

## Output

```
results_YYYYMMDD_HHMMSS/
├── 00_aggregate_networks/       # Consensus coexpression networks (Crow et al.)
│   ├── mouse/                   # {CellType}_mouse_aggregate.h5ad
│   └── human/                   # {CellType}_human_aggregate.h5ad
├── 01_tr_profiles/              # Per-dataset TR profiles in ortholog space
│   ├── mouse/profiles/          # CSV per dataset
│   └── human/profiles/
├── 01b_reproducibility/         # Within-species Top-K consistency
│   ├── mouse/                   # {species}_{CellType}_reproducibility.csv
│   └── human/                   #   + {species}_{CellType}_reprod_null.csv
├── 02_aggregated/               # TR profiles from aggregate networks (JSON)
│   ├── mouse_{CellType}_aggregated.json
│   └── human_{CellType}_aggregated.json
├── 03_comparison/{CellType}/    # Cross-species comparison results
│   ├── {CellType}_cross_species_results.csv
│   ├── {CellType}_summary.csv
│   ├── {CellType}_cross_species_null.csv
│   ├── {CellType}_gene_conservation.csv
│   └── {CellType}_gene_conservation_fullnet.csv
└── pipeline_info/               # Nextflow execution reports
```

### Output columns (`*_cross_species_results.csv`)

| Column | Description |
|--------|-------------|
| `human_tr` / `mouse_tr` | Orthologous TR pair |
| `human_n_datasets` / `mouse_n_datasets` | Number of datasets measuring the TR per species |
| `spearman_rho` | Spearman correlation of aggregate coexpression profiles |
| `spearman_pval` | Spearman p-value |
| `spearman_empirical_pval` | Empirical p-value from permutation null |
| `top200_overlap` | Fraction of top-200 coexpressed genes shared between species |
| `top200_empirical_pval` | Empirical p-value for top-200 overlap |
| `retrieval_human_in_mouse` | Quantile of mouse ortholog among all mouse TRs ranked by Top-K overlap with the human TR |
| `retrieval_mouse_in_human` | Reciprocal retrieval score |
| `n_shared_genes` | Number of shared ortholog genes used for comparison |

## Processes

### 0. AGGREGATE_NETWORKS

Aggregates per-dataset coexpression matrices into a consensus network per species + cell type using the rank-sum-rerank method. Datasets are aligned to a union gene space (unmeasured pairs imputed to 0.5), each matrix is ranked using minimum-ties ranking, ranks are summed across datasets, and the result is re-ranked and standardized to [0, 1]. For mouse, Ensembl IDs are converted to gene symbols via the gene mapping file (many-to-one: average; one-to-many: duplicate).

### 1. EXTRACT_AGGREGATE_PROFILES

Extracts TR coexpression profiles from each aggregate network, restricted to the shared ortholog gene space. Profiles are indexed by ortholog pair ID, ranked in descending order (rank 1 = highest coexpression), and saved as JSON.

### 2. EXTRACT_TR_PROFILES

Extracts TR coexpression profiles from each individual per-dataset matrix in the ortholog gene space. Same gene mapping and ortholog filtering as the aggregate extraction. Output is CSV for use in reproducibility analysis.

### 3. REPRODUCIBILITY

For each TR measured in ≥2 datasets, computes mean pairwise Top-K overlap across all dataset pairs as a within-species consistency score. An empirical null is generated by randomly selecting one TR per dataset across 1,000 iterations.

### 4. COMPARE_SPECIES

For each TR with a one-to-one ortholog and ≥ `min_datasets` datasets in both species, computes:

- **Spearman correlation** of aggregate profiles across shared ortholog genes (requires ≥10 genes, non-zero variance)
- **Top-K overlap**: fraction of the K most positively coexpressed genes shared between species
- **Ortholog retrieval score**: quantile position of the true ortholog's Top-K overlap among all TRs in the opposing species (1.0 = best). Computed bidirectionally
- **Gene-centric conservation**: per-gene Spearman and Top-K overlap across species, computed for both TR-restricted and full-network profiles
- **Empirical p-values**: derived from 1,000 random cross-species TR pairings; p = (count of null ≥ observed + 1) / (n_null + 1)
