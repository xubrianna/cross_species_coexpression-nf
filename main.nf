#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Cross-Species Coexpression Comparison Pipeline
 * Mouse vs Human TR coexpression network analysis
 *
 * Workflow:
 *   0. Aggregate coexpression networks across datasets (Crow et al. rank-sum-rerank)
 *   1. Extract TR profiles from aggregate networks in ortholog gene space
 *   2. Cross-species comparison (Spearman, Top-K, Ortholog retrieval)
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
 * Process 0: Aggregate coexpression networks across datasets
 *
 * For each species + cell type, aggregate per-dataset coexpression matrices
 * into a single consensus network using the Crow et al. rank-sum-rerank approach:
 *   - Align datasets to union gene space (unmeasured pairs → neutral midpoint)
 *   - Rank each dataset's matrix jointly (minimum ties)
 *   - Sum ranks across datasets → re-rank → standardize to [0, 1]
 */
process AGGREGATE_NETWORKS {
    tag "${species}/${cell_type}"
    publishDir "${params.outdir}/00_aggregate_networks/${species}", mode: 'copy'

    conda "${params.conda_env}"

    input:
    tuple val(species), val(cell_type), path("matrices/*")
    path mouse_gene_mapping

    output:
    tuple val(species), val(cell_type), path("${cell_type}_${species}_aggregate.h5ad"), emit: aggregate

    script:
    def mapping_arg = species == 'mouse' ? "--gene_mapping ${mouse_gene_mapping}" : ''
    """
    python ${projectDir}/bin/00-aggregateNetworks.py \\
        --input_dir matrices \\
        --species ${species} \\
        --cell_type ${cell_type} \\
        ${mapping_arg} \\
        --output ${cell_type}_${species}_aggregate.h5ad
    """
}

/*
 * Process 1: Extract TR profiles from aggregate networks
 *
 * Extract TR coexpression profiles from the aggregate network in the shared
 * ortholog gene space, formatted for cross-species comparison.
 */
process EXTRACT_AGGREGATE_PROFILES {
    tag "${species}/${cell_type}"
    publishDir "${params.outdir}/01_aggregated", mode: 'copy'

    conda "${params.conda_env}"

    input:
    tuple val(species), val(cell_type), path(aggregate_h5ad)
    path ortholog_file
    path mouse_tf_file
    path human_tf_file

    output:
    tuple val(species), val(cell_type), path("${species}_${cell_type}_aggregated.json"), emit: aggregated

    script:
    """
    python ${projectDir}/bin/01-extractAggregateProfiles.py \\
        --input ${aggregate_h5ad} \\
        --ortholog_file ${ortholog_file} \\
        --mouse_tf_file ${mouse_tf_file} \\
        --human_tf_file ${human_tf_file} \\
        --species ${species} \\
        --output ${species}_${cell_type}_aggregated.json
    """
}

/*
 * Process 2: Cross-species comparison
 *
 * For each orthologous gene present in both species' aggregate networks
 * and measured in >= min_datasets per species:
 *   - Spearman correlation of coexpression profiles
 *   - Top-K overlap (most positively coexpressed)
 *   - Ortholog retrieval score (Top-K overlap based)
 *   - Null distribution (randomly paired genes across species)
 *   - is_TR flag for transcription regulators
 */
process COMPARE_SPECIES {
    tag "${cell_type}"
    publishDir "${params.outdir}/02_comparison/${cell_type}", mode: 'copy'
    memory '16 GB'

    conda "${params.conda_env}"

    input:
    tuple val(cell_type), path(human_agg_h5ad, stageAs: 'human_aggregate.h5ad'), path(mouse_agg_h5ad, stageAs: 'mouse_aggregate.h5ad')
    path ortholog_file
    path ribosomal_genes
    path human_tf_file
    path mouse_tf_file

    output:
    path "${cell_type}_cross_species_results.csv", emit: results
    path "${cell_type}_summary.csv", emit: summary
    path "${cell_type}_cross_species_null.csv", emit: null_dist

    script:
    def ribo_arg = ribosomal_genes.name != 'NO_FILE' ? "--ribosomal_genes ${ribosomal_genes}" : ''
    """
    python ${projectDir}/bin/02-compareSpecies.py \\
        --human_aggregate_h5ad human_aggregate.h5ad \\
        --mouse_aggregate_h5ad mouse_aggregate.h5ad \\
        --ortholog_file ${ortholog_file} \\
        --human_tf_file ${human_tf_file} \\
        --mouse_tf_file ${mouse_tf_file} \\
        --cell_type ${cell_type} \\
        --top_k ${params.top_k} \\
        --min_datasets ${params.min_datasets} \\
        --n_null ${params.n_null} \\
        --output ${cell_type}_cross_species_results.csv \\
        --summary ${cell_type}_summary.csv \\
        --null_output ${cell_type}_cross_species_null.csv \\
        ${ribo_arg}
    """
}

/*
 * Main workflow
 */
workflow {
    // Reference files
    ortholog_file = file(params.ortholog_file)
    ribosomal_genes = params.ribosomal_genes ? file(params.ribosomal_genes) : file('NO_FILE')
    mouse_tf_file = file(params.mouse_tf_file)
    human_tf_file = file(params.human_tf_file)
    mouse_gene_mapping = file(params.mouse_gene_mapping)

    // Create channel for mouse coexpression matrices
    mouse_h5ad_ch = Channel
        .fromPath("${params.mouse_coexpr_dir}/*_corAggSparse.h5ad")
        .map { file -> tuple('mouse', file) }

    // Optionally add human matrices if provided
    // Priority: use ALL if available, otherwise CTL; never AD
    if (params.human_coexpr_dir) {
        human_h5ad_ch = Channel
            .fromPath("${params.human_coexpr_dir}/*_corAggSparse.h5ad")
            .filter { file -> !file.name.contains('_AD_') }
            .map { file ->
                // Filename: {CellType}_{StudyID}_{Condition}_corAggSparse.h5ad
                def parts = file.baseName.toString().replace('_corAggSparse', '').split('_')
                def condition = parts[-1]   // ALL or CTL
                def study_id = parts[1]     // e.g. Lau-2020
                def cell_type = parts[0]
                tuple("${cell_type}_${study_id}", condition, file)
            }
            .groupTuple(by: 0)  // group by cell_type + study_id
            .map { key, conditions, files ->
                // Pick ALL if available, otherwise CTL
                def idx = conditions.indexOf('ALL')
                if (idx == -1) idx = conditions.indexOf('CTL')
                if (idx == -1) return null
                tuple('human', files[idx])
            }
            .filter { it != null }
        all_h5ad_ch = mouse_h5ad_ch.mix(human_h5ad_ch)
    } else {
        all_h5ad_ch = mouse_h5ad_ch
    }

    // Process 0: Group by species + cell_type, aggregate networks across datasets
    // Parse cell_type from filename: e.g. "Ast_GSE124952_corAggSparse.h5ad" -> "Ast"
    grouped_h5ad = all_h5ad_ch
        .map { species, h5ad ->
            def cell_type = h5ad.baseName.toString().split('_')[0]
            tuple(species, cell_type, h5ad)
        }
        .groupTuple(by: [0, 1])  // group by (species, cell_type)

    aggregate_networks = AGGREGATE_NETWORKS(
        grouped_h5ad,
        mouse_gene_mapping
    )

    // Process 1: Extract TR profiles from aggregate networks in ortholog space
    aggregated = EXTRACT_AGGREGATE_PROFILES(
        aggregate_networks.aggregate,
        ortholog_file,
        mouse_tf_file,
        human_tf_file
    )

    // Process 2: Cross-species comparison (only if both species provided)
    if (params.human_coexpr_dir) {
        // Pair aggregate h5ad networks by cell_type
        human_agg_h5ad = aggregate_networks.aggregate
            .filter { species, cell_type, file -> species == 'human' }
            .map { species, cell_type, file -> tuple(cell_type, file) }

        mouse_agg_h5ad = aggregate_networks.aggregate
            .filter { species, cell_type, file -> species == 'mouse' }
            .map { species, cell_type, file -> tuple(cell_type, file) }

        paired_h5ad_ch = human_agg_h5ad.join(mouse_agg_h5ad)
            .map { cell_type, human_file, mouse_file ->
                tuple(cell_type, human_file, mouse_file)
            }

        COMPARE_SPECIES(paired_h5ad_ch, ortholog_file, ribosomal_genes, human_tf_file, mouse_tf_file)
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
        - 00_aggregate_networks/{species}/    : Aggregate coexpression matrices (Crow et al.)
        - 01_aggregated/                      : TR profiles from aggregate ({species}_{cell_type}_aggregated.json)
        - 02_comparison/{cell_type}/          : Cross-species comparison results
        """
        .stripIndent()
}
