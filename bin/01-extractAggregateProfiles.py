#!/usr/bin/env python3
"""
Extract TR coexpression profiles from an aggregate coexpression network
and format as JSON.

For each TR in the aggregate network:
  - Extract its coexpression profile restricted to ortholog genes
  - Map gene symbols to ortholog pair IDs for cross-species alignment
  - Compute summary metrics: rank, N_comeasured (from aggregate metadata)
  - Format as JSON

Input:  Aggregate h5ad from 00-aggregateNetworks.py
Output: JSON with per-TR profiles
        Schema: {tr_symbol: [{ortholog_id, Rank_aggr_coexpr, N_comeasured,
                              Avg_aggr_coexpr, Rank_single_best, Top200_count}]}
"""

import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path
from scipy import sparse


def load_orthologs(ortholog_file):
    """Load DIOPT one-to-one ortholog table."""
    return pd.read_csv(ortholog_file, sep='\t')


def load_tfs(tf_file):
    """Load TF/TR list from AnimalTFDB."""
    df = pd.read_csv(tf_file, sep='\t')
    df = df[df['Symbol'].notna() & (df['Symbol'] != '')]
    return df[['Symbol', 'Ensembl', 'Family']].drop_duplicates()


def main():
    parser = argparse.ArgumentParser(
        description='Extract TR profiles from aggregate network '
                    'for cross-species comparison'
    )
    parser.add_argument('--input', required=True,
                        help='Aggregate h5ad from 00-aggregateNetworks.py')
    parser.add_argument('--ortholog_file', required=True,
                        help='DIOPT ortholog table')
    parser.add_argument('--mouse_tf_file', required=True,
                        help='Mouse TF/TR list')
    parser.add_argument('--human_tf_file', required=True,
                        help='Human TF/TR list')
    parser.add_argument('--species', required=True, choices=['mouse', 'human'],
                        help='Species of the aggregate network')
    parser.add_argument('--output', required=True, help='Output JSON file')
    args = parser.parse_args()

    # Load aggregate network
    print(f"Loading aggregate network: {args.input}")
    adata = ad.read_h5ad(args.input)
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)

    gene_names = adata.var_names.tolist()
    gene_idx = {g: i for i, g in enumerate(gene_names)}
    n_datasets = adata.uns.get('n_datasets', 1)
    gene_presence = adata.varm.get('dataset_presence', None)

    print(f"  Matrix: {X.shape[0]} genes, aggregated from {n_datasets} datasets")

    # Load reference data
    ortho_df = load_orthologs(args.ortholog_file)
    mouse_tfs = load_tfs(args.mouse_tf_file)
    human_tfs = load_tfs(args.human_tf_file)

    if args.species == 'mouse':
        species_col = 'Symbol_mm'
        tr_symbols = set(mouse_tfs['Symbol'].tolist())
    else:
        species_col = 'Symbol_hg'
        tr_symbols = set(human_tfs['Symbol'].tolist())

    gene_to_orthoid = dict(zip(ortho_df[species_col], ortho_df['ID']))
    ortho_genes = set(ortho_df[species_col].tolist())

    # Find TRs and orthologs present in the aggregate matrix
    trs_present = [g for g in gene_names if g in tr_symbols]
    ortho_present = [g for g in gene_names if g in ortho_genes]
    ortho_indices = np.array([gene_idx[g] for g in ortho_present])
    ortho_ids = [gene_to_orthoid[g] for g in ortho_present]

    print(f"  TRs in reference: {len(tr_symbols)}")
    print(f"  TRs in matrix: {len(trs_present)}")
    print(f"  Ortholog genes in matrix: {len(ortho_present)}")

    output_data = {}

    for tr_i, tr in enumerate(trs_present):
        if (tr_i + 1) % 200 == 0:
            print(f"  Processing TR {tr_i + 1}/{len(trs_present)}: {tr}")

        tr_idx_val = gene_idx[tr]

        # Extract coexpression values with ortholog genes
        values = X[tr_idx_val, ortho_indices].flatten()

        # Exclude self (TR matching its own ortholog ID)
        mask = np.array([g != tr for g in ortho_present])
        filtered_vals = values[mask]
        filtered_ids = [oid for oid, m in zip(ortho_ids, mask) if m]
        filtered_genes = [g for g, m in zip(ortho_present, mask) if m]

        if len(filtered_vals) == 0:
            continue

        n_ortho = len(filtered_vals)

        # Rank: 1 = highest coexpression value
        order = np.argsort(-filtered_vals, kind='mergesort')
        ranks = np.empty(n_ortho, dtype=int)
        ranks[order] = np.arange(1, n_ortho + 1)

        # N_comeasured: exact count of datasets that measured BOTH this TR and
        # each partner gene (from the dataset_presence matrix)
        if gene_presence is not None:
            tr_presence = gene_presence[tr_idx_val]
            n_comeasured = np.array([
                int((tr_presence & gene_presence[gene_idx[g]]).sum())
                for g in filtered_genes
            ])
        else:
            n_comeasured = np.full(n_ortho, n_datasets, dtype=int)

        # Build records in the schema expected by 04-compareSpecies.py
        records = []
        for j in range(n_ortho):
            records.append({
                'ortholog_id': filtered_ids[j],
                'Rank_aggr_coexpr': int(ranks[j]),
                'N_comeasured': int(n_comeasured[j]),
                'Avg_aggr_coexpr': float(filtered_vals[j]),
                'Rank_single_best': int(ranks[j]),
                'Top200_count': int(ranks[j] <= 200),
            })

        output_data[tr] = records

    # Save JSON
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output_data, f)
    print(f"Saved {len(output_data)} TR profiles to {output_path}")


if __name__ == '__main__':
    main()
