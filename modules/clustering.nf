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
                    -P outcomes:"${outcomes}"


    bash chapter_package.sh "${rmd.baseName}"
    """

}

