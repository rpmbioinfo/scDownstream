process MULTIMODAL_INTEGRATION {
    tag "Performing multimodal integration..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_integrated_wnn.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"
    
    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val umap2_ndims
    val first_lsi_pc
    val integrate_datasets
    val sketch_cells
 

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_integrated_wnn.RDS", emit:seurat


    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P umap2_ndims:${umap2_ndims} \
                    -P first_lsi_pc:${first_lsi_pc} \
                    -P integrate_datasets:${integrate_datasets} \
                    -P sketch_cells:"${sketch_cells}" 

    bash chapter_package.sh "${rmd.baseName}"
    """

}