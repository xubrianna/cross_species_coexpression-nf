#!/usr/bin/env python3
"""
Compute TR Top-K coexpression reproducibility across datasets. For each TR (and ribosomal gene), measure
how consistently its top-K coexpressed partners are replicated across datasets
using mean pairwise overlap of top-K sets.

Also computes a null distribution by randomly sampling different TRs across
datasets and measuring their pairwise top-K overlap.

Input:  Directory of per-dataset TR profile CSVs (from 01-extractTRProfiles.py)
        Each has columns: [tr_symbol, ortholog_id, coexpr_value, dataset_id]
Output:
  - Reproducibility CSV: per-TR mean pairwise top-K overlap across datasets
  - Null CSV: null distribution of random TR overlap (n_null iterations)

The per-dataset TR profile CSVs contain coexpression values in the ortholog
gene space, so the overlap metric is directly comparable across species.
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor


def load_profiles(input_dir):
    """Load all per-dataset TR profile CSVs into a single DataFrame."""
    csv_files = sorted(Path(input_dir).glob('*.csv'))
    if not csv_files:
        return pd.DataFrame()

    usecols = ['tr_symbol', 'ortholog_id', 'coexpr_value', 'dataset_id']

    def _load(f):
        return pd.read_csv(f, usecols=usecols)

    with ThreadPoolExecutor() as pool:
        dfs = [df for df in pool.map(_load, csv_files) if len(df) > 0]

    if not dfs:
        return pd.DataFrame()

    profiles = pd.concat(dfs, ignore_index=True)
    for col in ['tr_symbol', 'ortholog_id', 'dataset_id']:
        profiles[col] = profiles[col].astype('category')
    return profiles


def get_topk_set(df, tr, dataset, k):
    """Get the top-K ortholog_ids for a TR in a specific dataset."""
    subset = df[(df['tr_symbol'] == tr) & (df['dataset_id'] == dataset)]
    if len(subset) < k:
        return set(subset.nlargest(len(subset), 'coexpr_value')['ortholog_id'])
    return set(subset.nlargest(k, 'coexpr_value')['ortholog_id'])


def compute_reproducibility(profiles_df, top_k=200):
    """
    Compute mean pairwise top-K overlap for each TR across datasets.

    For each TR:
      - Extract its top-K coexpressed partners in each dataset
      - Compute all pairwise overlaps (|intersection|) between datasets
      - Report the mean pairwise overlap

    Args:
        profiles_df: DataFrame with [tr_symbol, ortholog_id, coexpr_value, dataset_id]
        top_k: number of top coexpressed partners

    Returns:
        DataFrame with columns [tr_symbol, mean_topk_overlap, n_datasets, n_pairs]
    """
    all_trs = profiles_df['tr_symbol'].cat.categories.tolist()
    datasets = profiles_df['dataset_id'].cat.categories.tolist()

    # Pre-build top-K sets for all (TR, dataset) pairs
    print(f"  Building top-{top_k} sets for {len(all_trs)} TRs x {len(datasets)} datasets...")
    topk_sets = {}
    for tr in all_trs:
        tr_data = profiles_df[profiles_df['tr_symbol'] == tr]
        tr_datasets = tr_data['dataset_id'].cat.remove_unused_categories().categories.tolist()
        for ds in tr_datasets:
            ds_data = tr_data[tr_data['dataset_id'] == ds]
            if len(ds_data) >= top_k:
                topk_ids = set(
                    ds_data.nlargest(top_k, 'coexpr_value')['ortholog_id'].tolist()
                )
            else:
                topk_ids = set(ds_data['ortholog_id'].tolist())
            topk_sets[(tr, ds)] = topk_ids

    # Compute mean pairwise overlap per TR
    print(f"  Computing pairwise overlaps...")
    results = []
    for tr_i, tr in enumerate(all_trs):
        if (tr_i + 1) % 200 == 0:
            print(f"    TR {tr_i + 1}/{len(all_trs)}: {tr}")

        # Datasets that have this TR
        tr_ds_keys = [(tr, ds) for ds in datasets if (tr, ds) in topk_sets]
        n_ds = len(tr_ds_keys)

        if n_ds < 2:
            results.append({
                'tr_symbol': tr,
                'mean_topk_overlap': np.nan,
                'n_datasets': n_ds,
                'n_pairs': 0,
            })
            continue

        # Pairwise overlaps
        overlaps = []
        for (_, ds1), (_, ds2) in combinations(tr_ds_keys, 2):
            s1 = topk_sets[(tr, ds1)]
            s2 = topk_sets[(tr, ds2)]
            overlaps.append(len(s1 & s2))

        results.append({
            'tr_symbol': tr,
            'mean_topk_overlap': float(np.mean(overlaps)),
            'n_datasets': n_ds,
            'n_pairs': len(overlaps),
        })

    return pd.DataFrame(results)


def compute_null_distribution(profiles_df, top_k=200, n_null=1000):
    """
    Null distribution: randomly sample one TR per dataset, compute mean pairwise
    top-K overlap. Repeat n_null times.

    Args:
        profiles_df: DataFrame with [tr_symbol, ortholog_id, coexpr_value, dataset_id]
        top_k: number of top coexpressed partners
        n_null: number of null iterations

    Returns:
        DataFrame with columns [iteration, avg_overlap]
    """
    datasets = profiles_df['dataset_id'].cat.categories.tolist()
    if len(datasets) < 2:
        return pd.DataFrame({'iteration': [], 'avg_overlap': []})

    # Pre-compute per-dataset TR lists
    ds_trs = {}
    for ds in datasets:
        ds_data = profiles_df[profiles_df['dataset_id'] == ds]
        ds_trs[ds] = ds_data['tr_symbol'].cat.remove_unused_categories().categories.tolist()

    # Pre-build all top-K sets (indexed by (tr, dataset))
    topk_cache = {}
    for ds in datasets:
        ds_data = profiles_df[profiles_df['dataset_id'] == ds]
        for tr in ds_trs[ds]:
            tr_data = ds_data[ds_data['tr_symbol'] == tr]
            if len(tr_data) >= top_k:
                topk_cache[(tr, ds)] = set(
                    tr_data.nlargest(top_k, 'coexpr_value')['ortholog_id'].tolist()
                )
            else:
                topk_cache[(tr, ds)] = set(tr_data['ortholog_id'].tolist())

    rng = np.random.default_rng(42)
    null_results = []

    print(f"  Computing null distribution ({n_null} iterations)...")
    for i in range(n_null):
        if (i + 1) % 200 == 0:
            print(f"    Iteration {i + 1}/{n_null}")

        # Sample one random TR per dataset
        sampled_sets = []
        for ds in datasets:
            trs = ds_trs[ds]
            if not trs:
                continue
            tr = trs[rng.integers(len(trs))]
            s = topk_cache.get((tr, ds))
            if s is not None:
                sampled_sets.append(s)

        if len(sampled_sets) < 2:
            null_results.append(np.nan)
            continue

        # Pairwise overlaps
        overlaps = []
        for s1, s2 in combinations(sampled_sets, 2):
            overlaps.append(len(s1 & s2))
        null_results.append(float(np.mean(overlaps)))

    return pd.DataFrame({
        'iteration': range(1, n_null + 1),
        'avg_overlap': null_results,
    })


def main():
    parser = argparse.ArgumentParser(
        description='Compute TR Top-K coexpression reproducibility across datasets'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory of per-dataset TR profile CSVs')
    parser.add_argument('--species', required=True, choices=['mouse', 'human'])
    parser.add_argument('--cell_type', required=True, help='Cell type label')
    parser.add_argument('--top_k', type=int, default=200,
                        help='K for top-K overlap')
    parser.add_argument('--n_null', type=int, default=1000,
                        help='Number of null iterations')
    parser.add_argument('--output', required=True,
                        help='Output CSV: per-TR reproducibility')
    parser.add_argument('--null_output', required=True,
                        help='Output CSV: null distribution')
    args = parser.parse_args()

    profiles_df = load_profiles(args.input_dir)
    if profiles_df.empty:
        print("No profiles found!")
        pd.DataFrame().to_csv(args.output, index=False)
        pd.DataFrame().to_csv(args.null_output, index=False)
        return

    n_trs = profiles_df['tr_symbol'].nunique()
    n_ds = profiles_df['dataset_id'].nunique()
    print(f"Loaded {len(profiles_df)} rows: {n_trs} TRs x {n_ds} datasets")
    print(f"Species: {args.species}, Cell type: {args.cell_type}")

    # Reproducibility
    print(f"\nComputing top-{args.top_k} reproducibility...")
    reprod_df = compute_reproducibility(profiles_df, top_k=args.top_k)
    reprod_df['species'] = args.species
    reprod_df['cell_type'] = args.cell_type
    reprod_df['top_k'] = args.top_k
    reprod_df.to_csv(args.output, index=False)
    print(f"Saved {len(reprod_df)} TR reproducibility scores to {args.output}")

    # Null distribution
    print(f"\nComputing null distribution...")
    null_df = compute_null_distribution(profiles_df, top_k=args.top_k, n_null=args.n_null)
    null_df.to_csv(args.null_output, index=False)
    print(f"Saved {len(null_df)} null iterations to {args.null_output}")

    # Quick summary
    valid = reprod_df.dropna(subset=['mean_topk_overlap'])
    if len(valid) > 0:
        null_mean = null_df['avg_overlap'].mean()
        print(f"\nSummary:")
        print(f"  TRs with >= 2 datasets: {len(valid)}")
        print(f"  Mean top-{args.top_k} overlap: {valid['mean_topk_overlap'].mean():.1f}")
        print(f"  Median top-{args.top_k} overlap: {valid['mean_topk_overlap'].median():.1f}")
        print(f"  Null mean overlap: {null_mean:.1f}")


if __name__ == '__main__':
    main()
