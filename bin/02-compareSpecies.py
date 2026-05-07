#!/usr/bin/env python3
"""
Cross-species coexpression comparison: Mouse vs Human.

For each orthologous gene present in both species' aggregate coexpression
networks (and measured in >= min_datasets per species):
  1. Spearman correlation: compare coexpression profiles between human and mouse
  2. Top-K overlap: most positively coexpressed genes shared between species
  3. Ortholog retrieval score (Top-K overlap based):
       - For a human gene, compute Top-K overlap of its profile with
         every mouse gene's profile
       - The quantile at which the true mouse ortholog falls = retrieval score
       - Score of 1 = ortholog shares more top coexpressed partners than any
         other gene in the other species
       - Repeated reciprocally (mouse in human)
  4. Null comparison: similarities between randomly paired genes across species

Input:  aggregate h5ad coexpression networks for human and mouse (one cell type)
Output: CSV with per-gene comparison metrics (is_TR column flags TRs)
"""

import argparse
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr
from scipy import sparse


def load_tfs(tf_file):
    """Load TF/TR list from AnimalTFDB."""
    df = pd.read_csv(tf_file, sep='\t')
    df = df[df['Symbol'].notna() & (df['Symbol'] != '')]
    return set(df['Symbol'].tolist())


def compute_spearman(h_vec, m_vec):
    """Spearman correlation between aligned profiles (NaNs masked)."""
    mask = ~(np.isnan(h_vec) | np.isnan(m_vec))
    if mask.sum() < 10:
        return np.nan, np.nan
    h_valid, m_valid = h_vec[mask], m_vec[mask]
    if np.ptp(h_valid) == 0 or np.ptp(m_valid) == 0:
        return np.nan, np.nan
    rho, pval = spearmanr(h_valid, m_valid)
    return rho, pval


def compute_topk_overlap(h_vec, m_vec, k):
    """Fraction of top-K most positively coexpressed genes shared."""
    mask = ~(np.isnan(h_vec) | np.isnan(m_vec))
    if mask.sum() < k:
        return np.nan
    h_valid = h_vec[mask]
    m_valid = m_vec[mask]
    idx_valid = np.where(mask)[0]
    top_h = set(idx_valid[np.argsort(h_valid)[-k:]])
    top_m = set(idx_valid[np.argsort(m_valid)[-k:]])
    return len(top_h & top_m) / k


def build_topk_sparse(mat, k):
    """
    Build sparse CSR indicator matrix for top-K genes per row.
    Diagonal should already be NaN (set before calling).
    """
    n = mat.shape[0]
    rows, cols = [], []
    for i in range(n):
        profile = mat[i, :]
        valid_mask = ~np.isnan(profile)
        if valid_mask.sum() < k:
            continue
        valid_indices = np.where(valid_mask)[0]
        valid_vals = profile[valid_indices]
        topk_local = np.argsort(valid_vals)[-k:]
        topk_global = valid_indices[topk_local]
        rows.extend([i] * k)
        cols.extend(topk_global)
    data = np.ones(len(rows), dtype=np.float32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n))


def compute_retrieval_scores(h_mat, m_mat, top_k):
    """
    Compute ortholog retrieval scores for all genes using sparse Top-K overlap.

    Since genes are aligned (gene i in human = ortholog of gene i in mouse),
    the retrieval score for gene i is the quantile of overlap[i, i]
    among all overlap[i, :].
    """
    n = h_mat.shape[0]

    print("  Building Top-K sparse indicators...")
    H_topk = build_topk_sparse(h_mat, top_k)
    M_topk = build_topk_sparse(m_mat, top_k)

    h_valid = np.array(H_topk.sum(axis=1)).flatten() > 0
    m_valid = np.array(M_topk.sum(axis=1)).flatten() > 0

    h_in_m_scores = np.full(n, np.nan)
    m_in_h_scores = np.full(n, np.nan)

    print("  Computing human-in-mouse retrieval...")
    M_topk_T = M_topk.T.tocsc()
    for i in range(n):
        if not h_valid[i] or not m_valid[i]:
            continue
        overlap_vec = (H_topk[i] @ M_topk_T).toarray().flatten()
        valid_overlaps = overlap_vec[m_valid]
        ortholog_overlap = overlap_vec[i]
        h_in_m_scores[i] = float(np.mean(valid_overlaps <= ortholog_overlap))

    print("  Computing mouse-in-human retrieval...")
    H_topk_T = H_topk.T.tocsc()
    for i in range(n):
        if not m_valid[i] or not h_valid[i]:
            continue
        overlap_vec = (M_topk[i] @ H_topk_T).toarray().flatten()
        valid_overlaps = overlap_vec[h_valid]
        ortholog_overlap = overlap_vec[i]
        m_in_h_scores[i] = float(np.mean(valid_overlaps <= ortholog_overlap))

    return h_in_m_scores, m_in_h_scores


def compute_null(h_mat, m_mat, top_k, n_iterations):
    """
    Null distribution: random gene pairings across species.
    Diagonal is already NaN, so self-coexpression is automatically excluded.
    """
    n = h_mat.shape[0]
    rng = np.random.default_rng(42)
    null_spearman = np.empty(n_iterations)
    null_topk = np.empty(n_iterations)
    null_bottomk = np.empty(n_iterations)

    for it in range(n_iterations):
        h_idx = rng.integers(n)
        m_idx = rng.integers(n)
        h_vec = h_mat[h_idx, :]
        m_vec = m_mat[m_idx, :]

        mask = ~(np.isnan(h_vec) | np.isnan(m_vec))
        if mask.sum() < max(top_k, 10):
            null_spearman[it] = np.nan
            null_topk[it] = np.nan
            null_bottomk[it] = np.nan
            continue

        h_valid, m_valid = h_vec[mask], m_vec[mask]
        if np.ptp(h_valid) == 0 or np.ptp(m_valid) == 0:
            null_spearman[it] = np.nan
        else:
            null_spearman[it], _ = spearmanr(h_valid, m_valid)

        idx_valid = np.where(mask)[0]
        top_h = set(idx_valid[np.argsort(h_valid)[-top_k:]])
        top_m = set(idx_valid[np.argsort(m_valid)[-top_k:]])
        null_topk[it] = len(top_h & top_m) / top_k

        bot_h = set(idx_valid[np.argsort(h_valid)[:top_k]])
        bot_m = set(idx_valid[np.argsort(m_valid)[:top_k]])
        null_bottomk[it] = len(bot_h & bot_m) / top_k

    return null_spearman, null_topk, null_bottomk


def main():
    parser = argparse.ArgumentParser(
        description='Cross-species coexpression comparison (all orthologous genes)'
    )
    parser.add_argument('--human_aggregate_h5ad', required=True,
                        help='Human aggregate coexpression network h5ad')
    parser.add_argument('--mouse_aggregate_h5ad', required=True,
                        help='Mouse aggregate coexpression network h5ad')
    parser.add_argument('--ortholog_file', required=True,
                        help='DIOPT ortholog table')
    parser.add_argument('--human_tf_file', required=True,
                        help='Human TF/TR list (AnimalTFDB)')
    parser.add_argument('--mouse_tf_file', required=True,
                        help='Mouse TF/TR list (AnimalTFDB)')
    parser.add_argument('--cell_type', required=True, help='Cell type label')
    parser.add_argument('--top_k', type=int, default=200,
                        help='K for top-K overlap')
    parser.add_argument('--min_datasets', type=int, default=5,
                        help='Minimum datasets in both species for a gene')
    parser.add_argument('--n_null', type=int, default=1000,
                        help='Number of null iterations')
    parser.add_argument('--output', required=True, help='Output CSV')
    parser.add_argument('--summary', required=True, help='Summary CSV')
    parser.add_argument('--null_output', default=None,
                        help='Output CSV for null distribution')
    parser.add_argument('--ribosomal_genes', default=None,
                        help='File with ribosomal gene symbols (one per line)')
    args = parser.parse_args()

    # Load ortholog mapping
    ortho_df = pd.read_csv(args.ortholog_file, sep='\t')
    hg2mm = dict(zip(ortho_df['Symbol_hg'], ortho_df['Symbol_mm']))
    hg2id = dict(zip(ortho_df['Symbol_hg'], ortho_df['ID']))

    # Load TF/TR lists
    human_trs = load_tfs(args.human_tf_file)
    mouse_trs = load_tfs(args.mouse_tf_file)

    # Load aggregate networks
    print("Loading human aggregate network...")
    h_adata = ad.read_h5ad(args.human_aggregate_h5ad)
    print(f"  Human: {h_adata.shape[0]} genes")

    print("Loading mouse aggregate network...")
    m_adata = ad.read_h5ad(args.mouse_aggregate_h5ad)
    print(f"  Mouse: {m_adata.shape[0]} genes")

    # Find shared orthologs present in both networks
    human_genes = set(h_adata.var_names)
    mouse_genes = set(m_adata.var_names)
    shared_hg = sorted([g for g in human_genes if g in hg2mm and hg2mm[g] in mouse_genes])
    shared_mm = [hg2mm[g] for g in shared_hg]
    shared_ids = [hg2id[g] for g in shared_hg]

    print(f"Shared ortholog genes: {len(shared_hg)}")

    # Get n_datasets per gene for both species
    h_gene_presence = h_adata.varm.get('dataset_presence', None)
    m_gene_presence = m_adata.varm.get('dataset_presence', None)
    h_gene_idx = {g: i for i, g in enumerate(h_adata.var_names)}
    m_gene_idx = {g: i for i, g in enumerate(m_adata.var_names)}
    h_n_datasets_total = h_adata.uns.get('n_datasets', 1)
    m_n_datasets_total = m_adata.uns.get('n_datasets', 1)

    h_n_ds = {}
    m_n_ds = {}
    for hg, mm in zip(shared_hg, shared_mm):
        if h_gene_presence is not None:
            h_n_ds[hg] = int(h_gene_presence[h_gene_idx[hg]].sum())
        else:
            h_n_ds[hg] = h_n_datasets_total
        if m_gene_presence is not None:
            m_n_ds[mm] = int(m_gene_presence[m_gene_idx[mm]].sum())
        else:
            m_n_ds[mm] = m_n_datasets_total

    # Filter by min_datasets in both species
    qualifying = [
        (hg, mm, sid) for hg, mm, sid in zip(shared_hg, shared_mm, shared_ids)
        if h_n_ds[hg] >= args.min_datasets and m_n_ds[mm] >= args.min_datasets
    ]

    print(f"Qualifying genes (>= {args.min_datasets} datasets both): {len(qualifying)}")

    if not qualifying:
        print("No qualifying ortholog genes found!")
        pd.DataFrame().to_csv(args.output, index=False)
        pd.DataFrame().to_csv(args.summary, index=False)
        if args.null_output:
            pd.DataFrame().to_csv(args.null_output, index=False)
        return

    q_hg = [x[0] for x in qualifying]
    q_mm = [x[1] for x in qualifying]
    q_ids = [x[2] for x in qualifying]

    # Extract aligned submatrices
    print("Extracting aligned submatrices...")
    h_mat = h_adata[q_hg, q_hg].X
    m_mat = m_adata[q_mm, q_mm].X
    if hasattr(h_mat, 'toarray'):
        h_mat = h_mat.toarray()
    if hasattr(m_mat, 'toarray'):
        m_mat = m_mat.toarray()
    h_mat = h_mat.astype(np.float64)
    m_mat = m_mat.astype(np.float64)

    # Set diagonal to NaN (self-coexpression) so all functions exclude it via NaN masking
    np.fill_diagonal(h_mat, np.nan)
    np.fill_diagonal(m_mat, np.nan)

    n_genes = len(qualifying)

    # Compute null distribution
    print(f"Computing null distribution ({args.n_null} iterations)...")
    null_spearman, null_topk, null_bottomk = compute_null(
        h_mat, m_mat, args.top_k, args.n_null
    )
    if args.null_output:
        null_df = pd.DataFrame({
            'iteration': range(1, args.n_null + 1),
            'null_spearman': null_spearman,
            'null_topk_overlap': null_topk,
            'null_bottomk_overlap': null_bottomk
        })
        null_df.to_csv(args.null_output, index=False)
        print(f"Saved null distribution to {args.null_output}")

    # Compute retrieval scores
    print("Computing ortholog retrieval scores...")
    h_in_m_scores, m_in_h_scores = compute_retrieval_scores(
        h_mat, m_mat, args.top_k
    )

    # Compute per-gene metrics
    print("Computing per-gene comparison metrics...")
    results = []
    for i in range(n_genes):
        hg, mm, oid = q_hg[i], q_mm[i], q_ids[i]
        h_profile = h_mat[i, :]
        m_profile = m_mat[i, :]

        rho, pval = compute_spearman(h_profile, m_profile)
        topk = compute_topk_overlap(h_profile, m_profile, args.top_k)

        # Empirical p-values vs null
        spearman_emp_p = (
            (np.sum(null_spearman[~np.isnan(null_spearman)] >= rho) + 1) /
            (np.sum(~np.isnan(null_spearman)) + 1)
            if not np.isnan(rho) else np.nan
        )
        topk_emp_p = (
            (np.sum(null_topk[~np.isnan(null_topk)] >= topk) + 1) /
            (np.sum(~np.isnan(null_topk)) + 1)
            if not np.isnan(topk) else np.nan
        )

        results.append({
            'cell_type': args.cell_type,
            'human_gene': hg,
            'mouse_gene': mm,
            'ortholog_id': oid,
            'is_TR': hg in human_trs or mm in mouse_trs,
            'human_n_datasets': h_n_ds[hg],
            'mouse_n_datasets': m_n_ds[mm],
            'spearman_rho': rho,
            'spearman_pval': pval,
            'spearman_empirical_pval': spearman_emp_p,
            f'top{args.top_k}_overlap': topk,
            f'top{args.top_k}_empirical_pval': topk_emp_p,
            'retrieval_human_in_mouse': h_in_m_scores[i],
            'retrieval_mouse_in_human': m_in_h_scores[i],
            'n_shared_genes': n_genes - 1,
        })

    results_df = pd.DataFrame(results)

    # Flag ribosomal genes
    if args.ribosomal_genes:
        ribo_df = pd.read_csv(args.ribosomal_genes, sep='\t')
        ribo_genes = set(ribo_df.iloc[:, 0].dropna().str.strip())
        results_df['is_ribosomal'] = results_df['human_gene'].isin(ribo_genes)
    else:
        results_df['is_ribosomal'] = False

    results_df.to_csv(args.output, index=False)
    print(f"Saved {len(results_df)} gene comparisons to {args.output}")

    # Summary statistics
    n_tr = int(results_df['is_TR'].sum())
    n_ribo = int(results_df['is_ribosomal'].sum())
    tr_df = results_df[results_df['is_TR']]

    summary = pd.DataFrame([{
        'cell_type': args.cell_type,
        'n_gene_pairs': len(results_df),
        'n_tr_pairs': n_tr,
        'n_ribosomal': n_ribo,
        'n_shared_genes': n_genes,
        'median_spearman': results_df['spearman_rho'].median(),
        'mean_spearman': results_df['spearman_rho'].mean(),
        f'median_top{args.top_k}_overlap': results_df[f'top{args.top_k}_overlap'].median(),
        'median_spearman_tr': tr_df['spearman_rho'].median() if len(tr_df) else np.nan,
        f'median_top{args.top_k}_overlap_tr': tr_df[f'top{args.top_k}_overlap'].median() if len(tr_df) else np.nan,
        'median_retrieval_h_in_m': results_df['retrieval_human_in_mouse'].median(),
        'median_retrieval_m_in_h': results_df['retrieval_mouse_in_human'].median(),
        'min_datasets_threshold': args.min_datasets,
        'top_k': args.top_k,
    }])
    summary.to_csv(args.summary, index=False)
    print(f"Saved summary to {args.summary}")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Cell type: {args.cell_type}")
    print(f"Gene pairs compared: {len(results_df)} ({n_tr} TRs, {n_ribo} ribosomal)")
    print(f"Median Spearman rho: {results_df['spearman_rho'].median():.4f}")
    print(f"Median Top-{args.top_k} overlap: {results_df[f'top{args.top_k}_overlap'].median():.4f}")
    if len(tr_df):
        print(f"TR median Spearman: {tr_df['spearman_rho'].median():.4f}")
        print(f"TR median Top-{args.top_k} overlap: {tr_df[f'top{args.top_k}_overlap'].median():.4f}")
    print(f"Median retrieval (human in mouse): {results_df['retrieval_human_in_mouse'].median():.4f}")
    print(f"Median retrieval (mouse in human): {results_df['retrieval_mouse_in_human'].median():.4f}")


if __name__ == '__main__':
    main()
