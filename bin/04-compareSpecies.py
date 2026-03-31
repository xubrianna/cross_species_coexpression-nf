#!/usr/bin/env python3
"""
Cross-species coexpression comparison: Mouse vs Human.

For each TR with a one-to-one ortholog and >= min_datasets in both species:
  1. Spearman correlation: compare orthologous gene ranks in human vs mouse
  2. Top-K overlap:  most positively coexpressed genes in human vs mouse
  3. Bottom-K overlap: most negatively coexpressed genes (anti-coexpressed)
  4. Rank product consensus: combine human + mouse aggregate profiles via rank product
  5. Ortholog retrieval score (Top-K overlap based):
       - For a human TR, compute Top-K overlap of its aggregate profile with
         every mouse TR's aggregate profile
       - The quantile at which the true mouse ortholog falls = retrieval score
       - Score of 1 = ortholog shares more top coexpressed partners than any other TR
       - Repeated reciprocally (mouse in human)
  6. Null comparison: similarities between randomly shuffled aggregate profiles
     across species over 1000 iterations

Input:  aggregated profile pickles for human and mouse (one cell type)
Output: CSV with per-TR comparison metrics
"""

import argparse
import pickle
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, rankdata
from pathlib import Path


def load_profiles_as_matrix(pkl_path):
    """
    Load aggregated profiles from pickle (dict of TR DataFrames)
    and convert to TR x ortholog_gene matrix.

    Returns:
        profiles: DataFrame with TRs as index, ortholog_ids as columns, Avg_aggr_coexpr as values
        n_datasets: dict mapping tr_symbol -> max N_comeasured (approx datasets)
    """
    with open(pkl_path, 'rb') as f:
        tf_dict = pickle.load(f)

    profiles_data = {}
    n_datasets = {}
    for tr, gene_df in tf_dict.items():
        profiles_data[tr] = gene_df.set_index('ortholog_id')['Avg_aggr_coexpr']
        n_datasets[tr] = int(gene_df['N_comeasured'].max())

    profiles = pd.DataFrame(profiles_data).T
    return profiles, n_datasets


def compute_spearman(human_profile, mouse_profile):
    """Spearman correlation between aligned human and mouse profiles."""
    mask = ~(np.isnan(human_profile) | np.isnan(mouse_profile))
    if mask.sum() < 10:
        return np.nan, np.nan
    rho, pval = spearmanr(human_profile[mask], mouse_profile[mask])
    return rho, pval


def compute_topk_overlap(human_profile, mouse_profile, k):
    """
    Fraction of top-K most positively coexpressed genes shared between species.

    Args:
        human_profile: 1D array of coexpression values (ortholog-aligned)
        mouse_profile: 1D array of coexpression values (ortholog-aligned)
        k: number of top genes

    Returns:
        overlap fraction (intersection / k)
    """
    mask = ~(np.isnan(human_profile) | np.isnan(mouse_profile))
    if mask.sum() < k:
        return np.nan

    h_valid = human_profile[mask]
    m_valid = mouse_profile[mask]
    idx_valid = np.where(mask)[0]

    top_h = set(idx_valid[np.argsort(h_valid)[-k:]])
    top_m = set(idx_valid[np.argsort(m_valid)[-k:]])

    return len(top_h & top_m) / k


def compute_bottomk_overlap(human_profile, mouse_profile, k):
    """
    Fraction of bottom-K most anti-coexpressed genes shared between species.
    """
    mask = ~(np.isnan(human_profile) | np.isnan(mouse_profile))
    if mask.sum() < k:
        return np.nan

    h_valid = human_profile[mask]
    m_valid = mouse_profile[mask]
    idx_valid = np.where(mask)[0]

    bot_h = set(idx_valid[np.argsort(h_valid)[:k]])
    bot_m = set(idx_valid[np.argsort(m_valid)[:k]])

    return len(bot_h & bot_m) / k


def compute_rank_product(human_profiles, mouse_profiles, shared_genes):
    """
    Compute rank product consensus profiles between human and mouse.

    For each TR present in both species, rank each species' aggregate profile,
    then take the geometric mean of the ranks (rank product).

    Returns:
        DataFrame with columns [tr_symbol, ortholog_id, rank_product]
    """
    common_trs = set(human_profiles.index) & set(mouse_profiles.index)
    results = []

    for tr in common_trs:
        h_vals = human_profiles.loc[tr, shared_genes].values.astype(np.float64)
        m_vals = mouse_profiles.loc[tr, shared_genes].values.astype(np.float64)

        mask = ~(np.isnan(h_vals) | np.isnan(m_vals))
        if mask.sum() < 10:
            continue

        h_ranks = rankdata(h_vals[mask], method='average')
        m_ranks = rankdata(m_vals[mask], method='average')
        rp = np.sqrt(h_ranks * m_ranks)

        valid_genes = np.array(shared_genes)[mask]
        for gene, val in zip(valid_genes, rp):
            results.append({
                'tr_symbol': tr,
                'ortholog_id': gene,
                'rank_product': float(val)
            })

    return pd.DataFrame(results)


def _get_topk_set(profile_vec, gene_ids, k):
    """Return set of top-K gene IDs by value from a profile vector."""
    mask = ~np.isnan(profile_vec)
    valid_vals = profile_vec[mask]
    valid_genes = np.array(gene_ids)[mask]
    if len(valid_vals) < k:
        return set()
    topk_idx = np.argsort(valid_vals)[-k:]
    return set(valid_genes[topk_idx])


def compute_ortholog_retrieval_scores(human_profiles, mouse_profiles,
                                       ortholog_tr_map, shared_genes, top_k=200):
    """
    Compute ortholog retrieval scores using Top-K overlap.

    For each human TR:
      - Get its Top-K coexpressed partners (gene set)
      - Compute Top-K overlap with EVERY mouse TR's Top-K gene set
      - The quantile at which the true mouse ortholog falls among all
        mouse TRs = "human in mouse" retrieval score
      - Score of 1 = ortholog shares more top coexpressed partners than
        any other TR in the other species

    Reciprocally for each mouse TR (mouse in human).

    Args:
        human_profiles: DataFrame (human TRs x ortholog genes)
        mouse_profiles: DataFrame (mouse TRs x ortholog genes)
        ortholog_tr_map: dict mapping human_tr -> mouse_tr
        shared_genes: list of ortholog gene columns shared by both
        top_k: number of top genes for overlap

    Returns:
        dict: {(human_tr, mouse_tr): (human_in_mouse_score, mouse_in_human_score)}
    """
    results = {}

    human_trs = human_profiles.index.tolist()
    mouse_trs = mouse_profiles.index.tolist()
    human_mat = human_profiles[shared_genes].values
    mouse_mat = mouse_profiles[shared_genes].values
    human_tr_idx = {tr: i for i, tr in enumerate(human_trs)}
    mouse_tr_idx = {tr: i for i, tr in enumerate(mouse_trs)}
    gene_ids = list(shared_genes)

    inverse_map = {v: k for k, v in ortholog_tr_map.items()}

    # Pre-compute Top-K gene sets for all TRs in both species
    human_topk_sets = {
        tr: _get_topk_set(human_mat[i], gene_ids, top_k)
        for tr, i in human_tr_idx.items()
    }
    mouse_topk_sets = {
        tr: _get_topk_set(mouse_mat[i], gene_ids, top_k)
        for tr, i in mouse_tr_idx.items()
    }

    # Human-in-mouse: for each human TR, how well does its mouse ortholog rank?
    for h_tr, m_tr in ortholog_tr_map.items():
        if h_tr not in human_tr_idx or m_tr not in mouse_tr_idx:
            continue

        h_set = human_topk_sets[h_tr]
        if len(h_set) == 0:
            results[(h_tr, m_tr)] = [np.nan, np.nan]
            continue

        # Compute Top-K overlap with every mouse TR
        overlaps = []
        for m_candidate in mouse_trs:
            m_set = mouse_topk_sets[m_candidate]
            if len(m_set) == 0:
                overlaps.append(np.nan)
            else:
                overlaps.append(len(h_set & m_set))

        overlaps = np.array(overlaps, dtype=np.float64)
        valid = overlaps[~np.isnan(overlaps)]
        ortholog_overlap = overlaps[mouse_tr_idx[m_tr]]

        if len(valid) == 0 or np.isnan(ortholog_overlap):
            h_in_m_score = np.nan
        else:
            h_in_m_score = float(np.mean(valid <= ortholog_overlap))

        results[(h_tr, m_tr)] = [h_in_m_score, np.nan]

    # Mouse-in-human: for each mouse TR, how well does its human ortholog rank?
    for m_tr, h_tr in inverse_map.items():
        if m_tr not in mouse_tr_idx or h_tr not in human_tr_idx:
            continue

        m_set = mouse_topk_sets[m_tr]
        if len(m_set) == 0:
            key = (h_tr, m_tr)
            if key in results:
                results[key][1] = np.nan
            else:
                results[key] = [np.nan, np.nan]
            continue

        overlaps = []
        for h_candidate in human_trs:
            h_set = human_topk_sets[h_candidate]
            if len(h_set) == 0:
                overlaps.append(np.nan)
            else:
                overlaps.append(len(m_set & h_set))

        overlaps = np.array(overlaps, dtype=np.float64)
        valid = overlaps[~np.isnan(overlaps)]
        ortholog_overlap = overlaps[human_tr_idx[h_tr]]

        if len(valid) == 0 or np.isnan(ortholog_overlap):
            m_in_h_score = np.nan
        else:
            m_in_h_score = float(np.mean(valid <= ortholog_overlap))

        key = (h_tr, m_tr)
        if key in results:
            results[key][1] = m_in_h_score
        else:
            results[key] = [np.nan, m_in_h_score]

    return results


def compute_cross_species_null(human_profiles, mouse_profiles, shared_genes,
                                top_k=200, n_iterations=1000):
    """
    Generate null distribution for cross-species comparison.

    For each iteration, randomly shuffle aggregate profiles between species
    (i.e., randomly pair a human TR with a mouse TR) and compute Spearman
    correlation and Top-K overlap. Returns distributions of null similarities.
    """
    human_trs = human_profiles.index.tolist()
    mouse_trs = mouse_profiles.index.tolist()
    human_mat = human_profiles[shared_genes].values
    mouse_mat = mouse_profiles[shared_genes].values
    gene_ids = list(shared_genes)

    rng = np.random.default_rng(42)
    null_spearman = np.empty(n_iterations)
    null_topk = np.empty(n_iterations)
    null_bottomk = np.empty(n_iterations)

    for i in range(n_iterations):
        h_idx = rng.integers(len(human_trs))
        m_idx = rng.integers(len(mouse_trs))
        h_vec = human_mat[h_idx]
        m_vec = mouse_mat[m_idx]

        mask = ~(np.isnan(h_vec) | np.isnan(m_vec))
        if mask.sum() < max(top_k, 10):
            null_spearman[i] = np.nan
            null_topk[i] = np.nan
            null_bottomk[i] = np.nan
            continue

        rho, _ = spearmanr(h_vec[mask], m_vec[mask])
        null_spearman[i] = rho

        h_valid = h_vec[mask]
        m_valid = m_vec[mask]
        idx_valid = np.where(mask)[0]

        top_h = set(idx_valid[np.argsort(h_valid)[-top_k:]])
        top_m = set(idx_valid[np.argsort(m_valid)[-top_k:]])
        null_topk[i] = len(top_h & top_m) / top_k

        bot_h = set(idx_valid[np.argsort(h_valid)[:top_k]])
        bot_m = set(idx_valid[np.argsort(m_valid)[:top_k]])
        null_bottomk[i] = len(bot_h & bot_m) / top_k

    return null_spearman, null_topk, null_bottomk


def main():
    parser = argparse.ArgumentParser(
        description='Cross-species coexpression comparison'
    )
    parser.add_argument('--human_profiles', required=True,
                        help='Human aggregated profiles pickle')
    parser.add_argument('--mouse_profiles', required=True,
                        help='Mouse aggregated profiles pickle')
    parser.add_argument('--ortholog_file', required=True,
                        help='DIOPT ortholog table')
    parser.add_argument('--cell_type', required=True, help='Cell type label')
    parser.add_argument('--top_k', type=int, default=200,
                        help='K for top-K / bottom-K overlap')
    parser.add_argument('--min_datasets', type=int, default=5,
                        help='Minimum datasets in both species for a TR')
    parser.add_argument('--n_null', type=int, default=1000,
                        help='Number of null iterations for cross-species comparison')
    parser.add_argument('--output', required=True, help='Output CSV')
    parser.add_argument('--summary', required=True, help='Summary CSV')
    parser.add_argument('--rank_product_output', default=None,
                        help='Output parquet for rank product consensus profiles')
    parser.add_argument('--null_output', default=None,
                        help='Output CSV for cross-species null distribution')
    args = parser.parse_args()

    # Load ortholog mapping
    ortho_df = pd.read_csv(args.ortholog_file, sep='\t')
    human_to_mouse = dict(zip(ortho_df['Symbol_hg'], ortho_df['Symbol_mm']))
    mouse_to_human = dict(zip(ortho_df['Symbol_mm'], ortho_df['Symbol_hg']))
    ortho_pair_ids = dict(zip(ortho_df['Symbol_hg'], ortho_df['ID']))

    # Load aggregated profiles
    print("Loading human profiles...")
    human_profiles, human_n_datasets = load_profiles_as_matrix(args.human_profiles)
    print(f"  Human: {human_profiles.shape[0]} TRs x {human_profiles.shape[1]} genes")

    print("Loading mouse profiles...")
    mouse_profiles, mouse_n_datasets = load_profiles_as_matrix(args.mouse_profiles)
    print(f"  Mouse: {mouse_profiles.shape[0]} TRs x {mouse_profiles.shape[1]} genes")

    # Find shared ortholog gene columns
    shared_genes = sorted(set(human_profiles.columns) & set(mouse_profiles.columns))
    print(f"Shared ortholog genes: {len(shared_genes)}")

    # Find TR pairs: human TR with mouse ortholog, both present, >= min_datasets each
    tr_pairs = []
    for h_tr in human_profiles.index:
        m_tr = human_to_mouse.get(h_tr)
        if m_tr is None or m_tr not in mouse_profiles.index:
            continue
        h_n = human_n_datasets.get(h_tr, 0)
        m_n = mouse_n_datasets.get(m_tr, 0)
        if h_n >= args.min_datasets and m_n >= args.min_datasets:
            tr_pairs.append((h_tr, m_tr, h_n, m_n))

    print(f"TR pairs (>= {args.min_datasets} datasets both): {len(tr_pairs)}")

    if not tr_pairs:
        print("No qualifying TR pairs found!")
        pd.DataFrame().to_csv(args.output, index=False)
        pd.DataFrame().to_csv(args.summary, index=False)
        return

    # Build ortholog TR map for retrieval score (only qualifying TRs)
    ortholog_tr_map = {h: m for h, m, _, _ in tr_pairs}

    # Compute rank product consensus profiles
    print("Computing rank product consensus profiles...")
    rp_df = compute_rank_product(human_profiles, mouse_profiles, shared_genes)
    if args.rank_product_output and len(rp_df) > 0:
        rp_df['cell_type'] = args.cell_type
        rp_out = Path(args.rank_product_output)
        rp_out.parent.mkdir(parents=True, exist_ok=True)
        rp_df.to_parquet(rp_out, index=False)
        print(f"Saved rank product consensus profiles ({rp_df['tr_symbol'].nunique()} TRs) to {rp_out}")

    # Compute ortholog retrieval scores (using Top-K overlap)
    print("Computing ortholog retrieval scores (Top-K overlap)...")
    retrieval_scores = compute_ortholog_retrieval_scores(
        human_profiles, mouse_profiles, ortholog_tr_map, shared_genes,
        top_k=args.top_k
    )

    # Compute cross-species null distribution
    print(f"Computing cross-species null distribution ({args.n_null} iterations)...")
    null_spearman, null_topk, null_bottomk = compute_cross_species_null(
        human_profiles, mouse_profiles, shared_genes,
        top_k=args.top_k, n_iterations=args.n_null
    )
    if args.null_output:
        null_df = pd.DataFrame({
            'iteration': range(1, args.n_null + 1),
            'null_spearman': null_spearman,
            'null_topk_overlap': null_topk,
            'null_bottomk_overlap': null_bottomk
        })
        null_df.to_csv(args.null_output, index=False)
        print(f"Saved cross-species null distribution to {args.null_output}")

    # Compute per-TR metrics
    print("Computing per-TR comparison metrics...")
    results = []
    for h_tr, m_tr, h_n, m_n in tr_pairs:
        h_profile = human_profiles.loc[h_tr, shared_genes].values.astype(np.float64)
        m_profile = mouse_profiles.loc[m_tr, shared_genes].values.astype(np.float64)

        # Spearman
        rho, pval = compute_spearman(h_profile, m_profile)

        # Top-K overlap
        topk = compute_topk_overlap(h_profile, m_profile, args.top_k)

        # Bottom-K overlap
        bottomk = compute_bottomk_overlap(h_profile, m_profile, args.top_k)

        # Retrieval scores
        key = (h_tr, m_tr)
        h_in_m_score, m_in_h_score = retrieval_scores.get(key, [np.nan, np.nan])

        ortho_id = ortho_pair_ids.get(h_tr, f"{h_tr}_{m_tr}")

        # Empirical p-values vs null
        spearman_empirical_p = (
            (np.sum(null_spearman[~np.isnan(null_spearman)] >= rho) + 1) /
            (np.sum(~np.isnan(null_spearman)) + 1)
            if not np.isnan(rho) else np.nan
        )
        topk_empirical_p = (
            (np.sum(null_topk[~np.isnan(null_topk)] >= topk) + 1) /
            (np.sum(~np.isnan(null_topk)) + 1)
            if not np.isnan(topk) else np.nan
        )

        results.append({
            'cell_type': args.cell_type,
            'human_tr': h_tr,
            'mouse_tr': m_tr,
            'ortholog_id': ortho_id,
            'human_n_datasets': h_n,
            'mouse_n_datasets': m_n,
            'spearman_rho': rho,
            'spearman_pval': pval,
            'spearman_empirical_pval': spearman_empirical_p,
            f'top{args.top_k}_overlap': topk,
            f'top{args.top_k}_empirical_pval': topk_empirical_p,
            f'bottom{args.top_k}_overlap': bottomk,
            'retrieval_human_in_mouse': h_in_m_score,
            'retrieval_mouse_in_human': m_in_h_score,
            'n_shared_genes': len(shared_genes)
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, index=False)
    print(f"Saved {len(results_df)} TR comparisons to {args.output}")

    # Summary statistics
    summary = pd.DataFrame([{
        'cell_type': args.cell_type,
        'n_tr_pairs': len(tr_pairs),
        'n_shared_genes': len(shared_genes),
        'median_spearman': results_df['spearman_rho'].median(),
        'mean_spearman': results_df['spearman_rho'].mean(),
        f'median_top{args.top_k}_overlap': results_df[f'top{args.top_k}_overlap'].median(),
        f'median_bottom{args.top_k}_overlap': results_df[f'bottom{args.top_k}_overlap'].median(),
        'median_retrieval_h_in_m': results_df['retrieval_human_in_mouse'].median(),
        'median_retrieval_m_in_h': results_df['retrieval_mouse_in_human'].median(),
        'n_human_trs_total': human_profiles.shape[0],
        'n_mouse_trs_total': mouse_profiles.shape[0],
        'min_datasets_threshold': args.min_datasets,
        'top_k': args.top_k
    }])
    summary.to_csv(args.summary, index=False)
    print(f"Saved summary to {args.summary}")

    # Print summary
    print(f"\n{'='*60}")
    print(f"Cell type: {args.cell_type}")
    print(f"TR pairs compared: {len(tr_pairs)}")
    print(f"Shared ortholog genes: {len(shared_genes)}")
    print(f"Median Spearman rho: {results_df['spearman_rho'].median():.4f}")
    print(f"Median Top-{args.top_k} overlap: {results_df[f'top{args.top_k}_overlap'].median():.4f}")
    print(f"Median Bottom-{args.top_k} overlap: {results_df[f'bottom{args.top_k}_overlap'].median():.4f}")
    print(f"Median retrieval (human in mouse): {results_df['retrieval_human_in_mouse'].median():.4f}")
    print(f"Median retrieval (mouse in human): {results_df['retrieval_mouse_in_human'].median():.4f}")


if __name__ == '__main__':
    main()
