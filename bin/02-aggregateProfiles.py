#!/usr/bin/env python3
"""
Aggregate TR coexpression profiles across datasets for a given species + cell type.

For each TR:
  1. Collect profiles from all datasets
  2. Missing/tied values (from unmeasured genes) are imputed to the median
  3. Rank-standardize each per-dataset profile
  4. Average the rank-standardized profiles across datasets
  => One aggregate coexpression profile per TR

Input:  directory of parquet files (per-dataset TR profiles from extractTRProfiles.py)
Output: pickle dict {tr_symbol: DataFrame} with columns:
        [ortholog_id, Rank_aggr_coexpr, N_comeasured, Avg_aggr_coexpr,
         Rank_single_best, Top200_count]
"""

import argparse
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import rankdata


def rank_standardize(values):
    """Rank-standardize a vector, imputing NaN to the median rank before ranking."""
    arr = np.array(values, dtype=np.float64)
    # Impute NaN/missing to median of non-NaN values
    valid = arr[~np.isnan(arr)]
    if len(valid) == 0:
        return arr
    median_val = np.median(valid)
    arr[np.isnan(arr)] = median_val
    # Rank (average ties)
    ranks = rankdata(arr, method='average')
    # Standardize to [0, 1]
    n = len(ranks)
    if n > 1:
        ranks = (ranks - 1) / (n - 1)
    else:
        ranks = np.array([0.5])
    return ranks


def aggregate_profiles(profiles_df, cell_type):
    """
    Aggregate TR profiles across datasets for one species + cell type.

    Args:
        profiles_df: DataFrame with columns [tr_symbol, ortholog_id, coexpr_value, dataset_id, ...]
        cell_type: cell type being processed

    Returns:
        dict {tr_symbol: DataFrame} where each DataFrame has columns:
        [ortholog_id, Rank_aggr_coexpr, N_comeasured, Avg_aggr_coexpr,
         Rank_single_best, Top200_count]
    """
    # Filter to target cell type
    df = profiles_df[profiles_df['cell_type'] == cell_type].copy()
    if len(df) == 0:
        print(f"  No data for cell type: {cell_type}")
        return {}

    datasets = df['dataset_id'].unique()
    n_datasets = len(datasets)
    print(f"  Cell type: {cell_type}, datasets: {n_datasets}")

    # Get all ortholog IDs across all datasets
    all_ortho_ids = sorted(df['ortholog_id'].unique())
    n_ortho = len(all_ortho_ids)
    ortho_idx = {oid: i for i, oid in enumerate(all_ortho_ids)}

    # Get all TRs
    all_trs = sorted(df['tr_symbol'].unique())
    print(f"  TRs: {len(all_trs)}, Ortholog genes: {n_ortho}")

    tf_dict = {}

    for tr in all_trs:
        tr_data = df[df['tr_symbol'] == tr]
        tr_datasets = tr_data['dataset_id'].unique()
        n_tr_datasets = len(tr_datasets)

        # Build raw matrix: datasets x ortholog_genes
        raw_matrix = np.full((n_tr_datasets, n_ortho), np.nan, dtype=np.float64)

        for d_idx, dataset in enumerate(tr_datasets):
            ds_data = tr_data[tr_data['dataset_id'] == dataset]
            for _, row in ds_data.iterrows():
                o_idx = ortho_idx[row['ortholog_id']]
                raw_matrix[d_idx, o_idx] = row['coexpr_value']

        # Track which values were originally measured (before imputation)
        measured_mask = ~np.isnan(raw_matrix)

        # Rank-standardize each dataset's profile (imputes NaN to median)
        rs_matrix = np.copy(raw_matrix)
        for d_idx in range(n_tr_datasets):
            rs_matrix[d_idx, :] = rank_standardize(rs_matrix[d_idx, :])

        # Average rank-standardized profiles across datasets
        agg_profile = np.nanmean(rs_matrix, axis=0)

        # Per-dataset gene ranks from raw values (measured genes only, rank 1 = highest coexpr)
        per_dataset_ranks = np.full_like(raw_matrix, np.nan)
        for d_idx in range(n_tr_datasets):
            measured = measured_mask[d_idx, :]
            if measured.sum() > 0:
                vals = raw_matrix[d_idx, measured]
                per_dataset_ranks[d_idx, measured] = rankdata(-vals, method='average')

        # Build per-gene summary
        records = []
        for o_idx, ortho_id in enumerate(all_ortho_ids):
            n_comeasured = int(measured_mask[:, o_idx].sum())
            avg_coexpr = float(agg_profile[o_idx])

            gene_ranks = per_dataset_ranks[:, o_idx]
            valid_ranks = gene_ranks[~np.isnan(gene_ranks)]
            if len(valid_ranks) > 0:
                rank_single_best = int(np.min(valid_ranks))
                top200_count = int(np.sum(valid_ranks <= 200))
            else:
                rank_single_best = np.nan
                top200_count = 0

            records.append({
                'ortholog_id': ortho_id,
                'N_comeasured': n_comeasured,
                'Avg_aggr_coexpr': avg_coexpr,
                'Rank_single_best': rank_single_best,
                'Top200_count': top200_count,
            })

        gene_df = pd.DataFrame(records)
        gene_df = gene_df.sort_values('Avg_aggr_coexpr', ascending=False).reset_index(drop=True)
        gene_df['Rank_aggr_coexpr'] = range(1, len(gene_df) + 1)
        gene_df = gene_df[['ortholog_id', 'Rank_aggr_coexpr', 'N_comeasured',
                           'Avg_aggr_coexpr', 'Rank_single_best', 'Top200_count']]
        tf_dict[tr] = gene_df

    print(f"  Output: {len(tf_dict)} TRs as dict of DataFrames")
    return tf_dict


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate TR coexpression profiles across datasets'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing per-dataset profile parquets')
    parser.add_argument('--species', required=True, choices=['mouse', 'human'])
    parser.add_argument('--cell_type', required=True, help='Cell type to aggregate')
    parser.add_argument('--output', required=True, help='Output pickle file')
    args = parser.parse_args()

    # Load all profile files
    input_dir = Path(args.input_dir)
    parquet_files = sorted(input_dir.glob('*.parquet'))
    print(f"Found {len(parquet_files)} profile files in {input_dir}")

    if not parquet_files:
        print("No input files found!")
        with open(args.output, 'wb') as f:
            pickle.dump({}, f)
        return

    # Concatenate all datasets
    dfs = []
    for f in parquet_files:
        df = pd.read_parquet(f)
        if len(df) > 0:
            dfs.append(df)

    if not dfs:
        print("No non-empty profile files found!")
        with open(args.output, 'wb') as f:
            pickle.dump({}, f)
        return

    profiles_df = pd.concat(dfs, ignore_index=True)
    print(f"Combined: {len(profiles_df)} rows from {profiles_df['dataset_id'].nunique()} datasets")

    # Aggregate
    tf_dict = aggregate_profiles(profiles_df, args.cell_type)

    if not tf_dict:
        with open(args.output, 'wb') as f:
            pickle.dump({}, f)
        return

    # Save as pickle: dict of {tr_symbol: DataFrame}
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'wb') as f:
        pickle.dump(tf_dict, f)
    print(f"Saved aggregated profiles ({len(tf_dict)} TRs) to {output_path}")


if __name__ == '__main__':
    main()
