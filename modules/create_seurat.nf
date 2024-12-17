
process CONVERT_MATRIX {
    tag "Bitpacking counts for ${sample}"

    stageInMode 'copy'
    publishDir "$params.outdir/BPcells", mode:'copy', pattern: "*matrix_BP.tar.gz"

    memory '16 GB'
    cpus 2

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix), path(fragments_index), path(fragments), path(barcode_metrics)
    val pipeline
    val matrix_sel
    path rmd
    path annotation


    output:
    path "*_matrix_BP.tar.gz", emit: bp_mat
    path rmd, emit:script


    script:
    """
    mkdir ${sample}
    mv atac_fragments.tsv.gz ${sample}/${sample}_sample_fragments.tsv.gz
    mv atac_fragments.tsv.gz.tbi ${sample}/${sample}_sample_fragments.tsv.gz.tbi
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
    Rscript -e 'rmarkdown::render("${rmd}", params = list(pipeline = "${pipeline}", sample = "${sample}", matrix_sel = "${matrix_sel}", annotation = "${annotation}"))'
    """

}



process CREATE_SEURAT_ATAC {
    tag "Making Seurat object"

    stageInMode 'copy'
    publishDir "$params.outdir/seurat", mode:'copy', pattern: "*.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.Rmd"



    memory '16 GB'
    cpus 2

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix), path(fragments_index), path(fragments), path(barcode_metrics)
    val pipeline
    val matrix_sel
    path rmd
    path annotation
    path book_assets

    output:
    path "*_seurat.RDS", emit: seurat
    path "*_matrix_BP.tar.gz", emit: count_pkg
    path rmd, emit:script
    path "*_freeze.zip", emit: quarto



    script:
    """
    mkdir ${sample}
    mv atac_fragments.tsv.gz ${sample}/${sample}_sample_fragments.tsv.gz
    mv atac_fragments.tsv.gz.tbi ${sample}/${sample}_sample_fragments.tsv.gz.tbi
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P pipeline:"${pipeline}" \
                    -P sample:"${sample}" \
                    -P matrix_sel:"${matrix_sel}" \
                    -P annotation:"${annotation}" 

    bash chapter_package.sh "${rmd.baseName}"
    """

}


process CREATE_SEURAT {
    tag "Making Seurat object"
    stageInMode 'copy'

    // publishDir "$params.outdir", mode:'copy', pattern: "*.html"


    memory '8 GB'
    cpus 1

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix)
    val pipeline
    val matrix_sel
    path rmd
    path book_assets


    output:
    path "*_seurat.RDS", emit: seurat
    path "*_matrix_BP.tar.gz", emit: count_pkg
    path rmd, emit:script
    path "*_freeze.zip", emit: quarto

    script:
    """
    mkdir ${sample}
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P pipeline:"${pipeline}" \
                    -P sample:"${sample}" \
                    -P matrix_sel:"${matrix_sel}" 

    bash chapter_package.sh "${rmd.baseName}"

    """

}
