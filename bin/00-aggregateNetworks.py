#!/usr/bin/env python3
"""
Aggregate scRNA-seq coexpression networks across datasets.

Follows the rank-sum-rerank approach (Crow et al.), applied across
per-dataset coexpression matrices to produce a single consensus network:

  1. Align all datasets to the union gene space
  2. For each dataset's coexpression matrix:
     - Impute unmeasured gene pairs to 0.5
     - Extract upper triangle (prevent double-ranking symmetric elements)
     - Rank all elements jointly (minimum ties method)
  3. Sum rank matrices across datasets
  4. Re-rank the summed matrix
  5. Standardize to [0, 1] by dividing by the maximum rank

The 0.5 imputation (Step 2) ensures the ranking procedure includes
non-measured genes, placing them in between positive and negative
coexpression — equivalent to Crow et al.'s zero-imputation in correlation
space applied to the [0,1] rank-standardized matrices.

Input:  Directory of per-dataset h5ad coexpression matrices
Output: Single aggregate h5ad coexpression matrix with metadata:
        - var['n_datasets_measured']: per-gene count of contributing datasets
        - varm['dataset_presence']: boolean (n_genes x n_datasets) presence matrix
        - uns['n_datasets']: total number of datasets aggregated
        - uns['dataset_ids']: list of dataset identifiers
"""

import argparse
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path
from scipy.stats import rankdata
from scipy import sparse


# Neutral midpoint for [0,1] rank-standardized coexpression matrices.
NEUTRAL_VALUE = 0.5


def load_and_map_genes(h5ad_path, gene_mapping_df=None):
    """
    Load h5ad coexpression matrix, optionally mapping Ensembl IDs to symbols.

    Args:
        h5ad_path: path to h5ad file
        gene_mapping_df: DataFrame with 'gene.ensembl' and 'gene.symbol' columns

    Returns:
        X: dense numpy array (gene x gene, float64)
        genes: list of gene symbol strings
        symbol_to_ensembl: dict mapping gene symbol -> list of source Ensembl IDs
                           (empty dict when no mapping is applied)
    """
    adata = ad.read_h5ad(h5ad_path)
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)
    genes = adata.var_names.tolist()
    symbol_to_ensembl = {}  # gene symbol -> [ensembl_id, ...]

    if gene_mapping_df is not None:
        ensembl_to_symbols = (
            gene_mapping_df.groupby('gene.ensembl')['gene.symbol']
            .apply(list).to_dict()
        )

        new_names = []
        new_ensembl = []  # parallel list: source Ensembl ID for each new name
        src_indices = []
        for i, name in enumerate(genes):
            symbols = ensembl_to_symbols.get(name)
            if symbols:
                for sym in symbols:
                    new_names.append(sym)
                    new_ensembl.append(name)
                    src_indices.append(i)
            elif not name.startswith('ENSMUSG'):
                # Already a symbol
                new_names.append(name)
                new_ensembl.append(name)  # no Ensembl ID, use symbol itself
                src_indices.append(i)

        # Build symbol -> [ensembl_ids] mapping (before dedup)
        for sym, ens in zip(new_names, new_ensembl):
            symbol_to_ensembl.setdefault(sym, []).append(ens)
        # Deduplicate per symbol
        symbol_to_ensembl = {
            sym: sorted(set(ids)) for sym, ids in symbol_to_ensembl.items()
        }

        n_expanded = len(set(new_names))
        n_one_to_many = sum(
            1 for name in genes
            if name in ensembl_to_symbols and len(ensembl_to_symbols[name]) > 1
        )

        src_indices = np.array(src_indices)
        expanded = X[np.ix_(src_indices, src_indices)]

        # Average rows/cols with same symbol (many-Ensembl-to-1-symbol)
        df_matrix = pd.DataFrame(expanded, index=new_names, columns=new_names)
        if df_matrix.index.duplicated().any():
            df_matrix = df_matrix.groupby(level=0).mean()
            df_matrix = df_matrix.T.groupby(level=0).mean().T

        n_orig = len(adata.var_names)
        genes = df_matrix.index.tolist()
        X = df_matrix.values.astype(np.float64)
        print(f"    Gene mapping: {n_orig} Ensembl -> {len(genes)} symbols "
              f"({n_one_to_many} Ensembl with multiple symbols)")

    return X, genes, symbol_to_ensembl


def aggregate_networks(h5ad_files, gene_mapping_df=None):
    """
    Aggregate coexpression matrices across datasets using rank-sum-rerank.

    Args:
        h5ad_files: list of paths to h5ad coexpression matrices
        gene_mapping_df: optional DataFrame for Ensembl->symbol conversion

    Returns:
        AnnData with aggregate coexpression matrix and metadata
    """
    n_datasets = len(h5ad_files)
    print(f"Aggregating {n_datasets} datasets")

    # Step 1: Load all matrices and determine union gene space
    loaded = []
    dataset_ids = []
    all_symbol_to_ensembl = {}  # accumulated across datasets
    for f in h5ad_files:
        print(f"  Loading {Path(f).name}...")
        X, genes, sym_to_ens = load_and_map_genes(f, gene_mapping_df)
        loaded.append((X, set(genes), genes))
        dataset_id = Path(f).stem.replace('_corAggSparse', '').split('_', 1)[-1]
        dataset_ids.append(dataset_id)
        print(f"    {dataset_id}: {len(genes)} genes, shape {X.shape}")
        # Merge ensembl mappings
        for sym, ens_ids in sym_to_ens.items():
            all_symbol_to_ensembl.setdefault(sym, set()).update(ens_ids)

    # Union gene space
    all_genes = sorted(set(g for _, _, gl in loaded for g in gl))
    n_genes = len(all_genes)
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    print(f"  Union gene space: {n_genes} genes")

    # Track per-gene dataset presence (for downstream N_comeasured)
    gene_presence = np.zeros((n_genes, n_datasets), dtype=bool)
    for d_idx, (_, gene_set, _) in enumerate(loaded):
        for g in gene_set:
            gene_presence[gene_to_idx[g], d_idx] = True
    gene_dataset_counts = gene_presence.sum(axis=1).astype(int)

    print(f"  Per-gene dataset coverage: "
          f"median={np.median(gene_dataset_counts):.0f}, "
          f"min={gene_dataset_counts.min()}, "
          f"max={gene_dataset_counts.max()}")

    # Upper triangle indices (excluding diagonal — no self-coexpression)
    triu_r, triu_c = np.triu_indices(n_genes, k=1)
    n_upper = len(triu_r)
    print(f"  Upper triangle elements: {n_upper:,}")

    # Step 2–3: Rank each dataset and accumulate rank sum
    rank_sum = np.zeros(n_upper, dtype=np.float64)

    for d_idx, (X_orig, _, genes_orig) in enumerate(loaded):
        print(f"  Ranking dataset {d_idx + 1}/{n_datasets}: {dataset_ids[d_idx]}")

        # Align to union gene space.
        # Unmeasured gene pairs → NEUTRAL_VALUE (0.5, midpoint of [0,1]),
        # placing them between positive and negative coexpression when ranked.
        aligned = np.full((n_genes, n_genes), NEUTRAL_VALUE, dtype=np.float64)
        orig_indices = np.array([gene_to_idx[g] for g in genes_orig])
        aligned[np.ix_(orig_indices, orig_indices)] = X_orig

        # Extract upper triangle and rank (minimum ties method)
        upper_vals = aligned[triu_r, triu_c]
        del aligned
        ranks = rankdata(upper_vals, method='min').astype(np.float64)
        del upper_vals

        rank_sum += ranks
        del ranks

        # Free the original matrix
        loaded[d_idx] = None

    # Step 4: Re-rank the summed matrix
    print("  Re-ranking summed matrix...")
    final_ranks = rankdata(rank_sum, method='min')
    del rank_sum

    # Step 5: Standardize to [0, 1]
    max_rank = final_ranks.max()
    standardized = (final_ranks / max_rank).astype(np.float32)
    del final_ranks

    # Reconstruct symmetric matrix
    print("  Reconstructing symmetric matrix...")
    result = np.zeros((n_genes, n_genes), dtype=np.float32)
    result[triu_r, triu_c] = standardized
    result += result.T  # make symmetric
    np.fill_diagonal(result, 1.0)
    del standardized

    # Build AnnData with metadata
    # Map gene symbols back to source Ensembl IDs (semicolon-separated if multiple)
    ensembl_col = [
        ';'.join(sorted(all_symbol_to_ensembl.get(g, set())))
        if g in all_symbol_to_ensembl else ''
        for g in all_genes
    ]
    var_df = pd.DataFrame(
        {
            'n_datasets_measured': gene_dataset_counts,
            'ensembl_id': ensembl_col,
        },
        index=all_genes,
    )

    adata = ad.AnnData(
        X=result,
        obs=pd.DataFrame(index=all_genes),
        var=var_df,
    )
    adata.uns['n_datasets'] = n_datasets
    adata.uns['dataset_ids'] = dataset_ids
    adata.varm['dataset_presence'] = gene_presence

    print(f"  Output: {n_genes} x {n_genes} aggregate matrix from {n_datasets} datasets")
    return adata


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate coexpression networks across datasets '
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory of h5ad coexpression matrices')
    parser.add_argument('--species', required=True, choices=['mouse', 'human'])
    parser.add_argument('--cell_type', required=True,
                        help='Cell type label (for logging)')
    parser.add_argument('--gene_mapping', default=None,
                        help='CSV with gene.ensembl and gene.symbol columns '
                             '(for mouse Ensembl ID conversion)')
    parser.add_argument('--output', required=True, help='Output h5ad path')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    h5ad_files = sorted(input_dir.glob('*.h5ad'))
    print(f"Found {len(h5ad_files)} h5ad files in {input_dir}")
    print(f"Species: {args.species}, Cell type: {args.cell_type}")

    if not h5ad_files:
        print("ERROR: No input files found!")
        return

    gene_mapping_df = None
    if args.species == 'mouse' and args.gene_mapping:
        gene_mapping_df = pd.read_csv(args.gene_mapping)

    adata = aggregate_networks(h5ad_files, gene_mapping_df=gene_mapping_df)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"Saved aggregate network to {output_path}")


if __name__ == '__main__':
    main()
