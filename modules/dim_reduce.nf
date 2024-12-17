process DIMENSION_REDUCTION {
    tag "Performing dimension reduction..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_dim_reduced.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"

    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val umap2_ndims
    val first_lsi_pc
    val rna_normalization_method
    val sketch_cells
    val sketch_n
 

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_dim_reduced.RDS", emit:seurat
    path rmd, emit:script

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P umap2_ndims:${umap2_ndims} \
                    -P first_lsi_pc:${first_lsi_pc} \
                    -P rna_normalization_method:"${rna_normalization_method}" \
                    -P sketch_cells:"${sketch_cells}" \
                    -P sketch_n:"${sketch_n}"

    bash chapter_package.sh "${rmd.baseName}"
    """

}