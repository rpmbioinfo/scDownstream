process CELL_INTERACTOME {
    tag "Computing cell-cell interactome..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_cellchat.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val integrate_datasets
    val outcomes
    val species
    val cc_avg_method
    val cc_pathways
    

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_chromvar.RDS", emit:seurat
    path "_freeze", emit: freeze

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P integrate_datasets:${integrate_datasets} \
                    -P cpus:"${task.cpus}" \
                    -P species:"${species}" 

    bash chapter_package.sh "${rmd.baseName}"
    """

}
