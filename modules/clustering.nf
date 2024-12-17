process CLUSTERING {
    tag "Performing clustering..."
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

    bash chapter_package.sh "${rmd.baseName}"
    """

}

