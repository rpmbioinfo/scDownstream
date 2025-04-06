process CLUSTERING {
    tag "Performing clustering..."
    stageInMode 'copy'

    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val integrate_datasets
    val sketch_cells
    val rna_normalization_method


    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_clustered.RDS", emit:seurat
    path "_freeze", emit: freeze

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P integrate_datasets:${integrate_datasets} \
                    -P sketch_cells:"${sketch_cells}" \
                    -P rna_normalization_method:"${rna_normalization_method}"

    bash chapter_package.sh "${rmd.baseName}"
    """

}


process SET_CLUSTERS {
    tag "Comparing clusters..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_clustered.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val clustering2_res
    val integrate_datasets
    val outcomes
    val sketch_cells
    val de_method
    val de_latent_vars
    val de_min_pct
    val de_logfc
    val markers_rna
    val markers_adt
    val rna_normalization_method

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_clustered.RDS", emit:seurat
    path "_freeze", emit: freeze

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P clustering2_res:${clustering2_res} \
                    -P integrate_datasets:${integrate_datasets} \
                    -P outcomes:"${outcomes}" \
                    -P sketch_cells:"${sketch_cells}" \
                    -P de_method:"${de_method}" \
                    -P de_latent_vars:"${de_latent_vars}" \
                    -P de_min_pct:"${de_min_pct}" \
                    -P de_logfc:"${de_logfc}" \
                    -P markers_rna:"${markers_rna}" \
                    -P markers_adt:"${markers_adt}" \
                    -P rna_normalization_method:"${rna_normalization_method}"

    bash chapter_package.sh "${rmd.baseName}"
    """

}


