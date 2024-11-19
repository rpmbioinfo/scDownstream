
process CONVERT_MATRIX {
    tag "Bitpacking counts for (${sample}"

    stageInMode 'copy'
    publishDir "$params.outdir/seurat", mode:'copy', pattern: "*.RDS"
    publishDir "$params.outdir/fragments", mode:'copy', pattern: "*_sample_fragments.tsv*"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.Rmd"



    memory '16 GB'
    cpus 2

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix), path(fragments_index), path(fragments), path(barcode_metrics)
    path annotation
    val pipeline
    path rmd


    output:
    path "*_raw_counts_seurat.RDS", emit: raw
    path "*_filtered_counts_seurat.RDS", emit: filt
    path "*_sample_fragments.tsv.gz", emit: fragments
    path "*_sample_fragments.tsv.gz.tbi", emit: fragment_index
    path rmd, emit:script


    script:
    """
    mv atac_fragments.tsv.gz ${sample}_sample_fragments.tsv.gz
    mv atac_fragments.tsv.gz.tbi ${sample}_sample_fragments.tsv.gz.tbi
    Rscript -e 'rmarkdown::render("${rmd}", params = list(annotation = "${annotation}", pipeline = "${pipeline}", sample = "${sample}"))'
    """

}


process CREATE_SEURAT_MULTIMODAL {
    tag "Making Seurat object"

    stageInMode 'copy'
    publishDir "$params.outdir/seurat", mode:'copy', pattern: "*.RDS"
    publishDir "$params.outdir/fragments", mode:'copy', pattern: "*_sample_fragments.tsv*"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.Rmd"



    memory '16 GB'
    cpus 2

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix), path(fragments_index), path(fragments), path(barcode_metrics)
    path annotation
    val pipeline
    path rmd


    output:
    path "*_raw_counts_seurat.RDS", emit: raw
    path "*_filtered_counts_seurat.RDS", emit: filt
    path "*_sample_fragments.tsv.gz", emit: fragments
    path "*_sample_fragments.tsv.gz.tbi", emit: fragment_index
    path rmd, emit:script


    script:
    """
    mv atac_fragments.tsv.gz ${sample}_sample_fragments.tsv.gz
    mv atac_fragments.tsv.gz.tbi ${sample}_sample_fragments.tsv.gz.tbi
    Rscript -e 'rmarkdown::render("${rmd}", params = list(annotation = "${annotation}", pipeline = "${pipeline}", sample = "${sample}"))'
    """

}


process CREATE_SEURAT_GEX {
    tag "Making Seurat object"
    stageInMode 'copy'

    // publishDir "$params.outdir", mode:'copy', pattern: "*.html"


    memory '8 GB'
    cpus 1

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix)
    val pipeline
    path rmd


    output:
    path "*_raw_counts_seurat.RDS", emit: raw
    path "*_filtered_counts_seurat.RDS", emit: filt

    script:
    """
    Rscript -e 'rmarkdown::render("${rmd}", params = list(pipeline = "${pipeline}", sample = "${sample}"))'
    """

}

process CREATE_SEURAT_ATAC {
    tag "Making Seurat object"
    stageInMode 'copy'

    // publishDir "$params.outdir", mode:'copy', pattern: "*.html"


    memory '16 GB'
    cpus 2

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix), path(fragments_index), path(fragments)
    path(annotation)
    path rmd


    output:
    path "*_raw_counts_seurat.RDS", emit: raw
    path "*_filtered_counts_seurat.RDS", emit: filt
    path "*_sample_fragments.tsv.gz", emit: fragments
    path "*_sample_fragments.tsv.gz.tbi", emit: fragment_index

    script:
    """
    mv atac_fragments.tsv.gz ${sample}_sample_fragments.tsv.gz
    mv atac_fragments.tsv.gz.tbi ${sample}_sample_fragments.tsv.gz.tbi
    Rscript -e 'rmarkdown::render("${rmd}", params = list(annotation = "${annotation}", pipeline = "${pipeline}", sample = "${sample}"))'
    """

}
