process INTEGRATE_DATASETS {
    tag "Performing dataset integration..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_integrated.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val integrate_datasets
    val integration_method
    val integrate_by
    val view_batch
    val umap2_ndims
    val first_lsi_pc
    val sketch_cells
    val rna_normalization_method

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_integrated.RDS", emit:seurat
    path rmd, emit:script

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P integrate_datasets:"${integrate_datasets}" \
                    -P integration_method:"${integration_method}" \
                    -P integrate_by:"${integrate_by}" \
                    -P view_batch:"${view_batch}" \
                    -P umap2_ndims:${umap2_ndims} \
                    -P first_lsi_pc:${first_lsi_pc} \
                    -P rna_normalization_method:"${rna_normalization_method}" \
                    -P sketch_cells:"${sketch_cells}" 

    bash chapter_package.sh "${rmd.baseName}"
    """

}

