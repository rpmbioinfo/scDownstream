process QUALITY_CONTROL {
    tag "Performing Quality Control..."

    stageInMode 'copy'
    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_after_qc.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "_quarto_full.yml"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "index.qmd"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "run_funcs.R"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "favicon-32x32.png"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "rpm_en_logo_lowres.jpg"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.basename}/*"


    input:
    path seurat_solo
    path fragment_files
    path rmd
    path book_assets
    val pipeline
    path metadata
    val genome
    val mito_regex
    val ribo_regex
    val percent_mt
    val percent_ribo
    val nCount_RNA_min
    val nCount_RNA_max
    val nFeature_RNA_min
    val nFeature_RNA_max
    val nCount_ATAC_min
    val nCount_ATAC_max
    val nucleosome_signal_min
    val nucleosome_signal_max
    val TSS_enrichment
    val atac_peak_region_fragments_min
    val atac_peak_region_fragments_max
    val pct_reads_in_peaks_min
    val blacklist_fraction
    val sample_exclusion
    val umap1_ndims
    val clustering1_res
    val doublet_removal
    val doublet_confidence

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_after_qc.RDS", emit:seurat
    path rmd, emit:script
    path book_assets, emit: assets


    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P pipeline:"${pipeline}" \
                    -P metadata:"${metadata}" \
                    -P genome:"${genome}" \
                    -P mito_regex:"${mito_regex}" \
                    -P ribo_regex:"${ribo_regex}" \
                    -P percent_mt:"${percent_mt}" \
                    -P percent_ribo:"${percent_ribo}" \
                    -P nCount_RNA_min:"${nCount_RNA_min}" \
                    -P nCount_RNA_max:"${nCount_RNA_max}" \
                    -P nFeature_RNA_min:"${nFeature_RNA_min}" \
                    -P nFeature_RNA_max:"${nFeature_RNA_max}" \
                    -P nCount_ATAC_min:"${nCount_ATAC_min}" \
                    -P nCount_ATAC_max:"${nCount_ATAC_max}" \
                    -P nucleosome_signal_min:"${nucleosome_signal_min}" \
                    -P nucleosome_signal_max:"${nucleosome_signal_max}" \
                    -P TSS_enrichment:"${TSS_enrichment}" \
                    -P atac_peak_region_fragments_min:"${atac_peak_region_fragments_min}" \
                    -P atac_peak_region_fragments_max:"${atac_peak_region_fragments_max}" \
                    -P pct_reads_in_peaks_min:"${pct_reads_in_peaks_min}" \
                    -P blacklist_fraction:"${blacklist_fraction}" \
                    -P sample_exclusion:"${sample_exclusion}" \
                    -P umap1_ndims:"${umap1_ndims}" \
                    -P clustering1_res:"${clustering1_res}" \
                    -P doublet_removal:"${doublet_removal}" \
                    -P doublet_confidence:"${doublet_confidence}"

    bash chapter_package.sh "${rmd.baseName}"
    """

}
