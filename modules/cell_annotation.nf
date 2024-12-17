process CELL_ANNOTATION{
    tag "Performing cell annotation ..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_annotated.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"

    memory { run_scimilarity == true ? 96.GB : 24.GB }
    cpus { run_scimilarity == true ? 12 : 3 }


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val integrate_datasets
    val annotation_resolution
    val run_azimuth
    val azimuth_reference
    val run_singler
    val singler_reference
    val run_scimilarity
    path scimilarity_reference
    val selected_method
    val markers_rna
    val markers_adt
    val outcomes
    val sketch_cells
    path celltype_dictionary

 

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_annotated.RDS", emit:seurat
    path rmd, emit:script
    path "_freeze", emit: freeze

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P annotation_resolution:"${annotation_resolution}" \
                    -P integrate_datasets:"${integrate_datasets}" \
                    -P run_azimuth:"${run_azimuth}" \
                    -P azimuth_reference:"${azimuth_reference}" \
                    -P run_singler:"${run_singler}" \
                    -P singler_reference:"${singler_reference}" \
                    -P run_scimilarity:"${run_scimilarity}" \
                    -P scimilarity_reference:"${scimilarity_reference}" \
                    -P selected_method:"${selected_method}" \
                    -P markers_rna:"${markers_rna}" \
                    -P markers_adt:"${markers_adt}" \
                    -P selected_method:"${selected_method}" \
                    -P outcomes:"${outcomes}" \
                    -P sketch_cells:"${sketch_cells}" \
                    -P celltype_dictionary:"${celltype_dictionary}"


    bash chapter_package.sh "${rmd.baseName}"
    """

}



