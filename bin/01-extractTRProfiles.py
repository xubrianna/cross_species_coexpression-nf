#!/usr/bin/env python3
"""
Extract TR coexpression profiles from aggregated coexpression matrices,
restricted to the shared ortholog gene space (DIOPT one-to-one orthologs).

For each TR found in the matrix:
  - Extract the TR's row from the coexpression matrix
  - Restrict to genes in the ortholog gene space
  - Map gene symbols to ortholog pair IDs for cross-species alignment
  - Output a ranked profile

Input:  aggregated coexpression h5ad (genes x genes, rank-normalized)
Output: parquet with columns [tr_symbol, ortholog_id, coexpr_value, n_genes, dataset_id]
"""

import argparse
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path


def load_orthologs(ortholog_file):
    """Load DIOPT one-to-one ortholog table."""
    df = pd.read_csv(ortholog_file, sep='\t')
    return df


def load_tfs(tf_file):
    """Load TF/TR list from AnimalTFDB."""
    df = pd.read_csv(tf_file, sep='\t')
    df = df[df['Symbol'].notna() & (df['Symbol'] != '')]
    return df[['Symbol', 'Ensembl', 'Family']].drop_duplicates()


def parse_file_metadata(filename):
    """Parse cell type and study ID from filename.
    
    Expected: {CellType}_{StudyID}_corAggSparse.h5ad
    """
    parts = filename.replace('_corAggSparse', '').replace('.h5ad', '').split('_', 1)
    if len(parts) >= 2:
        return parts[0], parts[1]
    return parts[0], 'unknown'


def extract_profiles(adata, tr_symbols, ortho_genes, gene_to_orthoid, dataset_id):
    """
    Extract TR coexpression profiles restricted to ortholog gene space.

    Args:
        adata: AnnData coexpression matrix (genes x genes)
        tr_symbols: set of TR gene symbols to extract
        ortho_genes: set of ortholog gene symbols in this species
        gene_to_orthoid: dict mapping species gene symbol -> ortholog pair ID
        dataset_id: identifier for this dataset

    Returns:
        DataFrame with columns [tr_symbol, ortholog_id, coexpr_value, dataset_id]
    """
    gene_names = adata.var_names.tolist()
    gene_idx = {g: i for i, g in enumerate(gene_names)}

    # TRs present in this matrix
    trs_present = [g for g in gene_names if g in tr_symbols]
    # Ortholog genes present in this matrix
    ortho_present = [g for g in gene_names if g in ortho_genes]
    ortho_indices = np.array([gene_idx[g] for g in ortho_present])

    if not trs_present or len(ortho_present) == 0:
        return pd.DataFrame()

    coexpr_matrix = adata.X
    records = []

    for tr in trs_present:
        tr_idx = gene_idx[tr]
        # Get coexpression values with all ortholog genes
        values = coexpr_matrix[tr_idx, ortho_indices].flatten()

        for gene, val in zip(ortho_present, values):
            if gene == tr:
                continue  # skip self
            ortho_id = gene_to_orthoid.get(gene)
            if ortho_id is not None:
                records.append({
                    'tr_symbol': tr,
                    'ortholog_id': ortho_id,
                    'coexpr_value': float(val),
                    'dataset_id': dataset_id
                })

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(
        description='Extract TR coexpression profiles in ortholog gene space'
    )
    parser.add_argument('--input', required=True, help='Input h5ad coexpression matrix')
    parser.add_argument('--ortholog_file', required=True, help='DIOPT ortholog table')
    parser.add_argument('--mouse_tf_file', required=True, help='Mouse TF/TR list')
    parser.add_argument('--human_tf_file', required=True, help='Human TF/TR list')
    parser.add_argument('--species', required=True, choices=['mouse', 'human'],
                        help='Species of the input data')
    parser.add_argument('--gene_mapping', default=None,
                        help='CSV mapping file with gene.ensembl and gene.symbol columns '
                             '(for converting Ensembl IDs to symbols)')
    parser.add_argument('--output', required=True, help='Output parquet file')
    args = parser.parse_args()

    # Load reference data
    ortho_df = load_orthologs(args.ortholog_file)
    mouse_tfs = load_tfs(args.mouse_tf_file)
    human_tfs = load_tfs(args.human_tf_file)

    # Build ortholog mappings
    if args.species == 'mouse':
        species_col = 'Symbol_mm'
        tr_symbols = set(mouse_tfs['Symbol'].tolist())
    else:
        species_col = 'Symbol_hg'
        tr_symbols = set(human_tfs['Symbol'].tolist())

    # Map species gene symbol -> ortholog pair ID
    gene_to_orthoid = dict(zip(ortho_df[species_col], ortho_df['ID']))
    ortho_genes = set(ortho_df[species_col].tolist())

    # Load coexpression matrix
    input_path = Path(args.input)
    adata = ad.read_h5ad(input_path)
    cell_type, study_id = parse_file_metadata(input_path.name)
    dataset_id = f"{cell_type}_{study_id}"

    # Convert Ensembl IDs to gene symbols if mapping provided
    if args.gene_mapping:
        mapping_df = pd.read_csv(args.gene_mapping)
        ensmug_to_symbol = dict(zip(mapping_df['gene.ensembl'], mapping_df['gene.symbol']))
        original_names = adata.var_names.tolist()
        new_names = [ensmug_to_symbol.get(g, g) for g in original_names]
        # Drop genes that couldn't be mapped (still have ENSMUSG prefix)
        mapped_mask = [not n.startswith('ENSMUSG') for n in new_names]
        n_mapped = sum(mapped_mask)
        n_unmapped = len(mapped_mask) - n_mapped
        adata = adata[:, mapped_mask]
        adata.var_names = [n for n, m in zip(new_names, mapped_mask) if m]
        # Handle any duplicate symbols after mapping (keep first)
        if adata.var_names.duplicated().any():
            adata = adata[:, ~adata.var_names.duplicated()]
        print(f"  Gene mapping: {n_mapped} mapped, {n_unmapped} unmapped Ensembl IDs dropped")

    print(f"Processing {args.species} | {dataset_id} | matrix shape: {adata.shape}")
    print(f"  TRs in reference: {len(tr_symbols)}")
    print(f"  Ortholog genes in reference: {len(ortho_genes)}")

    gene_names = set(adata.var_names.tolist())
    print(f"  TRs in matrix: {len(tr_symbols & gene_names)}")
    print(f"  Ortholog genes in matrix: {len(ortho_genes & gene_names)}")

    # Extract profiles
    profiles_df = extract_profiles(
        adata, tr_symbols, ortho_genes, gene_to_orthoid, dataset_id
    )

    if len(profiles_df) == 0:
        print("  WARNING: No TR profiles extracted!")
        # Write empty parquet with correct schema
        profiles_df = pd.DataFrame(
            columns=['tr_symbol', 'ortholog_id', 'coexpr_value', 'dataset_id']
        )

    # Add metadata
    profiles_df['cell_type'] = cell_type
    profiles_df['species'] = args.species
    profiles_df['n_ortho_genes'] = len(ortho_genes & gene_names)

    # Save
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profiles_df.to_parquet(output_path, index=False)

    n_trs = profiles_df['tr_symbol'].nunique()
    print(f"  Saved {len(profiles_df)} rows for {n_trs} TRs -> {output_path}")


if __name__ == '__main__':
    main()
