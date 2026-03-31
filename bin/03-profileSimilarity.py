#!/usr/bin/env python3
"""
Within-species pairwise profile similarity and null distribution.

For each TR measured in >= min_datasets:
  1. From each dataset, extract its Top-K and Bottom-K coexpressed partners.
  2. Compute pairwise overlaps (set intersection size) across all dataset pairs.
  3. Average pairwise overlap = consistency score for that TR.

Null comparison:
  For 1000 iterations, randomly select one TR per dataset to create a set of
  shuffled profiles, compute pairwise Top-K overlaps, and average. This yields
  a null distribution. A TR with avg similarity > all 1000 nulls has empirical
  p-value < 0.001.

Input:  directory of parquet files (per-dataset TR profiles from extractTRProfiles.py)
Output: CSV with per-TR consistency metrics and empirical p-values
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations


def get_topk_genes(profile_df, tr, dataset, k, bottom=False):
    """Get top-K (or bottom-K) genes for a TR in a dataset by coexpr_value."""
    sub = profile_df[(profile_df['tr_symbol'] == tr) &
                     (profile_df['dataset_id'] == dataset)]
    if len(sub) == 0:
        return set()
    if bottom:
        sub = sub.nsmallest(k, 'coexpr_value')
    else:
        sub = sub.nlargest(k, 'coexpr_value')
    return set(sub['ortholog_id'].tolist())


def compute_pairwise_overlap(gene_sets):
    """Compute average pairwise intersection size among a list of gene sets."""
    if len(gene_sets) < 2:
        return np.nan
    overlaps = [len(a & b) for a, b in combinations(gene_sets, 2)]
    return np.mean(overlaps)


def compute_tr_consistency(profiles_df, min_datasets, top_k):
    """
    For each TR, compute average pairwise Top-K and Bottom-K overlaps
    across its datasets.

    Returns DataFrame with columns:
      [tr_symbol, n_datasets, mean_topk_overlap, mean_bottomk_overlap]
    """
    # Count datasets per TR
    tr_datasets = profiles_df.groupby('tr_symbol')['dataset_id'].nunique()
    qualifying_trs = tr_datasets[tr_datasets >= min_datasets].index.tolist()

    results = []
    for tr in qualifying_trs:
        tr_data = profiles_df[profiles_df['tr_symbol'] == tr]
        datasets = tr_data['dataset_id'].unique()

        topk_sets = [get_topk_genes(tr_data, tr, ds, top_k) for ds in datasets]
        bottomk_sets = [get_topk_genes(tr_data, tr, ds, top_k, bottom=True) for ds in datasets]

        # Filter out empty sets
        topk_sets = [s for s in topk_sets if len(s) > 0]
        bottomk_sets = [s for s in bottomk_sets if len(s) > 0]

        results.append({
            'tr_symbol': tr,
            'n_datasets': len(datasets),
            'mean_topk_overlap': compute_pairwise_overlap(topk_sets),
            'mean_bottomk_overlap': compute_pairwise_overlap(bottomk_sets),
        })

    return pd.DataFrame(results)


def compute_null_distribution(profiles_df, n_iterations, top_k):
    """
    Generate null distribution by shuffling TR identities across datasets.

    For each iteration:
      - Select one random TR from each dataset
      - Collect their Top-K gene sets
      - Compute average pairwise overlap
    """
    datasets = profiles_df['dataset_id'].unique()
    # Pre-compute available TRs per dataset
    trs_per_dataset = {
        ds: profiles_df[profiles_df['dataset_id'] == ds]['tr_symbol'].unique()
        for ds in datasets
    }

    rng = np.random.default_rng(42)
    null_values = np.empty(n_iterations)

    for i in range(n_iterations):
        gene_sets = []
        for ds in datasets:
            available_trs = trs_per_dataset[ds]
            if len(available_trs) == 0:
                continue
            random_tr = rng.choice(available_trs)
            gene_set = get_topk_genes(profiles_df, random_tr, ds, top_k)
            if len(gene_set) > 0:
                gene_sets.append(gene_set)

        null_values[i] = compute_pairwise_overlap(gene_sets)

    return null_values


def main():
    parser = argparse.ArgumentParser(
        description='Within-species pairwise profile similarity and null testing'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing per-dataset profile parquets')
    parser.add_argument('--species', required=True, choices=['mouse', 'human'])
    parser.add_argument('--cell_type', required=True, help='Cell type to process')
    parser.add_argument('--top_k', type=int, default=200,
                        help='K for top-K / bottom-K overlap')
    parser.add_argument('--min_datasets', type=int, default=5,
                        help='Minimum datasets for a TR to be included')
    parser.add_argument('--n_null', type=int, default=1000,
                        help='Number of null iterations')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--null_output', required=True,
                        help='Output CSV for null distribution')
    args = parser.parse_args()

    # Load all profile files
    input_dir = Path(args.input_dir)
    parquet_files = sorted(input_dir.glob('*.parquet'))
    print(f"Found {len(parquet_files)} profile files in {input_dir}")

    if not parquet_files:
        print("No input files found!")
        pd.DataFrame(columns=['tr_symbol', 'n_datasets', 'mean_topk_overlap',
                               'mean_bottomk_overlap', 'empirical_pval']
                      ).to_csv(args.output, index=False)
        pd.DataFrame(columns=['iteration', 'null_avg_overlap']
                      ).to_csv(args.null_output, index=False)
        return

    dfs = [pd.read_parquet(f) for f in parquet_files]
    dfs = [df for df in dfs if len(df) > 0]
    if not dfs:
        print("No non-empty profile files!")
        return

    profiles_df = pd.concat(dfs, ignore_index=True)
    profiles_df = profiles_df[profiles_df['cell_type'] == args.cell_type]
    print(f"Loaded {len(profiles_df)} rows, "
          f"{profiles_df['dataset_id'].nunique()} datasets, "
          f"{profiles_df['tr_symbol'].nunique()} TRs for {args.cell_type}")

    # 1. Compute per-TR consistency
    print(f"Computing per-TR pairwise consistency (top_k={args.top_k})...")
    consistency_df = compute_tr_consistency(profiles_df, args.min_datasets, args.top_k)
    print(f"  {len(consistency_df)} TRs with >= {args.min_datasets} datasets")

    # 2. Compute null distribution
    print(f"Computing null distribution ({args.n_null} iterations)...")
    null_values = compute_null_distribution(profiles_df, args.n_null, args.top_k)

    # 3. Empirical p-values: fraction of null >= observed
    consistency_df['empirical_pval'] = consistency_df['mean_topk_overlap'].apply(
        lambda obs: (np.sum(null_values >= obs) + 1) / (len(null_values) + 1)
        if not np.isnan(obs) else np.nan
    )

    # Add metadata
    consistency_df['species'] = args.species
    consistency_df['cell_type'] = args.cell_type
    consistency_df['top_k'] = args.top_k

    # Save
    consistency_df.to_csv(args.output, index=False)
    print(f"Saved TR consistency to {args.output}")

    null_df = pd.DataFrame({
        'iteration': range(1, args.n_null + 1),
        'null_avg_overlap': null_values
    })
    null_df.to_csv(args.null_output, index=False)
    print(f"Saved null distribution to {args.null_output}")

    # Summary
    print(f"\n{'='*60}")
    print(f"Species: {args.species}, Cell type: {args.cell_type}")
    print(f"TRs analysed: {len(consistency_df)}")
    print(f"Median Top-{args.top_k} pairwise overlap: "
          f"{consistency_df['mean_topk_overlap'].median():.2f}")
    print(f"Null mean: {np.nanmean(null_values):.2f}, "
          f"Null max: {np.nanmax(null_values):.2f}")
    n_sig = (consistency_df['empirical_pval'] < 0.001).sum()
    print(f"TRs with empirical p < 0.001: {n_sig}")


if __name__ == '__main__':
    main()
