#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Cross-Species Coexpression Comparison Pipeline
 * Mouse vs Human TR coexpression network analysis
 *
 * Workflow:
 *   1. Extract TR coexpression profiles in shared ortholog gene space
 *   2. Aggregate profiles across datasets per TR per species
 *   3. Compare mouse vs human (Spearman, Top-K, Bottom-K, Ortholog retrieval)
 */

log.info """\
    CROSS-SPECIES COEXPRESSION PIPELINE
    ====================================
    mouse_coexpr_dir   : ${params.mouse_coexpr_dir}
    human_coexpr_dir   : ${params.human_coexpr_dir ?: 'not provided (mouse-only mode)'}
    ortholog_file      : ${params.ortholog_file}
    mouse_tf_file      : ${params.mouse_tf_file}
    human_tf_file      : ${params.human_tf_file}
    mouse_gene_mapping : ${params.mouse_gene_mapping}
    top_k              : ${params.top_k}
    min_datasets       : ${params.min_datasets}
    outdir             : ${params.outdir}
    """
    .stripIndent()

/*
 * Process 1: Extract TR coexpression profiles in the ortholog gene space
 *
 * For each coexpression matrix (h5ad), extract one profile per TR:
 *   a ranked list of coexpression values restricted to ortholog genes.
 */
process EXTRACT_TR_PROFILES {
    tag "${species}/${h5ad.baseName}"
    publishDir "${params.outdir}/01_tr_profiles/${species}", mode: 'copy'

    conda "${params.conda_env}"

    input:
    tuple val(species), path(h5ad)
    path ortholog_file
    path mouse_tf_file
    path human_tf_file
    path mouse_gene_mapping

    output:
    tuple val(species), path("profiles/${h5ad.baseName}_tr_profiles.pkl"), emit: profiles

    script:
    def mapping_arg = species == 'mouse' ? "--gene_mapping ${mouse_gene_mapping}" : ''
    """
    python ${projectDir}/bin/01-extractTRProfiles.py \\
        --input ${h5ad} \\
        --ortholog_file ${ortholog_file} \\
        --mouse_tf_file ${mouse_tf_file} \\
        --human_tf_file ${human_tf_file} \\
        --species ${species} \\
        ${mapping_arg} \\
        --output profiles/${h5ad.baseName}_tr_profiles.pkl
    """
}

/*
 * Process 2: Aggregate TR profiles across datasets
 *
 * Per species + cell type + TR:
 *   - Missing/tied values imputed to median
 *   - Rank-standardized profiles averaged
 *   => One aggregate coexpression profile per TR per species
 */
process AGGREGATE_PROFILES {
    tag "${species}/${cell_type}"
    publishDir "${params.outdir}/02_aggregated/${species}", mode: 'copy'
    memory '16 GB'

    conda "${params.conda_env}"

    input:
    tuple val(species), val(cell_type), path("profiles/*")

    output:
    tuple val(species), val(cell_type), path("${cell_type}_aggregated.pkl"), emit: aggregated

    script:
    """
    python ${projectDir}/bin/02-aggregateProfiles.py \\
        --input_dir profiles \\
        --species ${species} \\
        --cell_type ${cell_type} \\
        --output ${cell_type}_aggregated.pkl
    """
}

/*
 * Process 2b: Within-species pairwise profile similarity + null distribution
 *
 * For each TR measured in >= min_datasets:
 *   - Compute pairwise Top-K and Bottom-K overlaps across datasets
 *   - Average = consistency score
 *   - Generate null distribution (1000 iterations) for empirical p-values
 */
process PROFILE_SIMILARITY {
    tag "${species}/${cell_type}"
    publishDir "${params.outdir}/02b_profile_similarity/${species}", mode: 'copy'
    memory '16 GB'

    conda "${params.conda_env}"

    input:
    tuple val(species), val(cell_type), path("profiles/*")

    output:
    tuple val(species), val(cell_type), path("${cell_type}_consistency.csv"), emit: consistency
    tuple val(species), val(cell_type), path("${cell_type}_null_distribution.csv"), emit: null_dist

    script:
    """
    python ${projectDir}/bin/03-profileSimilarity.py \\
        --input_dir profiles \\
        --species ${species} \\
        --cell_type ${cell_type} \\
        --top_k ${params.top_k} \\
        --min_datasets ${params.min_datasets} \\
        --n_null ${params.n_null} \\
        --output ${cell_type}_consistency.csv \\
        --null_output ${cell_type}_null_distribution.csv
    """
}

/*
 * Process 3: Cross-species comparison
 *
 * For each TR with a one-to-one ortholog and >= min_datasets in both species:
 *   - Spearman correlation of aggregate profiles
 *   - Top-K overlap (most positively coexpressed)
 *   - Bottom-K overlap (most negatively coexpressed)
 *   - Rank product consensus profile
 *   - Ortholog retrieval score (Top-K overlap based)
 *   - Null distribution (shuffled aggregate profiles across species)
 */
process COMPARE_SPECIES {
    tag "${cell_type}"
    publishDir "${params.outdir}/03_comparison", mode: 'copy'
    memory '16 GB'

    conda "${params.conda_env}"

    input:
    tuple val(cell_type), path(human_agg), path(mouse_agg)
    path ortholog_file

    output:
    path "${cell_type}_cross_species_results.csv", emit: results
    path "${cell_type}_summary.csv", emit: summary
    path "${cell_type}_rank_product.pkl", emit: rank_product
    path "${cell_type}_cross_species_null.csv", emit: null_dist

    script:
    """
    python ${projectDir}/bin/04-compareSpecies.py \\
        --human_profiles ${human_agg} \\
        --mouse_profiles ${mouse_agg} \\
        --ortholog_file ${ortholog_file} \\
        --cell_type ${cell_type} \\
        --top_k ${params.top_k} \\
        --min_datasets ${params.min_datasets} \\
        --n_null ${params.n_null} \\
        --output ${cell_type}_cross_species_results.csv \\
        --summary ${cell_type}_summary.csv \\
        --rank_product_output ${cell_type}_rank_product.pkl \\
        --null_output ${cell_type}_cross_species_null.csv
    """
}

/*
 * Main workflow
 */
workflow {
    // Reference files
    ortholog_file = file(params.ortholog_file)
    mouse_tf_file = file(params.mouse_tf_file)
    human_tf_file = file(params.human_tf_file)
    mouse_gene_mapping = file(params.mouse_gene_mapping)

    // Create channel for mouse coexpression matrices
    mouse_h5ad_ch = Channel
        .fromPath("${params.mouse_coexpr_dir}/*_corAggSparse.h5ad")
        .map { file -> tuple('mouse', file) }

    // Optionally add human matrices if provided
    if (params.human_coexpr_dir) {
        human_h5ad_ch = Channel
            .fromPath("${params.human_coexpr_dir}/*_corAggSparse.h5ad")
            .map { file -> tuple('human', file) }
        all_h5ad_ch = mouse_h5ad_ch.mix(human_h5ad_ch)
    } else {
        all_h5ad_ch = mouse_h5ad_ch
    }

    // Step 1: Extract TR profiles in ortholog gene space
    tr_profiles = EXTRACT_TR_PROFILES(
        all_h5ad_ch,
        ortholog_file,
        mouse_tf_file,
        human_tf_file,
        mouse_gene_mapping
    )

    // Step 2: Group profiles by species + cell_type, then aggregate
    // Parse cell_type from filename: e.g. "Ast_GSE124952_tr_profiles.pkl" -> "Ast"
    grouped_profiles = tr_profiles.profiles
        .map { species, profile ->
            def cell_type = profile.baseName.toString().split('_')[0]
            tuple(species, cell_type, profile)
        }
        .groupTuple(by: [0, 1])  // group by (species, cell_type)
        .map { species, cell_type, profiles ->
            tuple(species, cell_type, profiles)
        }

    aggregated = AGGREGATE_PROFILES(grouped_profiles)

    // Step 2b: Within-species pairwise similarity and null testing
    PROFILE_SIMILARITY(grouped_profiles)

    // Step 3: Cross-species comparison (only if both species provided)
    if (params.human_coexpr_dir) {
        human_agg = aggregated.aggregated
            .filter { species, cell_type, file -> species == 'human' }
            .map { species, cell_type, file -> tuple(cell_type, file) }

        mouse_agg = aggregated.aggregated
            .filter { species, cell_type, file -> species == 'mouse' }
            .map { species, cell_type, file -> tuple(cell_type, file) }

        // Join on cell_type to pair human + mouse for comparison
        paired_ch = human_agg.join(mouse_agg)
            .map { cell_type, human_file, mouse_file ->
                tuple(cell_type, human_file, mouse_file)
            }

        COMPARE_SPECIES(paired_ch, ortholog_file)
    }
}

workflow.onComplete {
    log.info """\
        Pipeline completed!
        ==================
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        Results   : ${params.outdir}

        Output Structure:
        - 01_tr_profiles/    : Per-dataset TR profiles in ortholog space
        - 02_aggregated/     : Aggregated profiles per TR per species per cell type
        - 03_comparison/     : Cross-species comparison results
        """
        .stripIndent()
}
