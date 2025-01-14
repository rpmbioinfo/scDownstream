#!/usr/bin/env nextflow
import groovy.json.JsonOutput

//include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

nextflow.enable.dsl = 2

params.help      = false





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INITIALISATION } from './subworkflows/local/pipeline_init'

include { QUALITY_CONTROL } from './modules/qc.nf'
include { CALL_PEAKS } from './modules/qc.nf'
include { CREATE_SEURAT_ATAC } from './modules/create_seurat.nf'
include { CREATE_SEURAT } from './modules/create_seurat.nf'
include { CONVERT_MATRIX } from './modules/create_seurat.nf'
include { GET_ENSDB_ANNOTATION } from './modules/get_ensdb_annot.nf'
include { DIMENSION_REDUCTION } from './modules/dim_reduce.nf'
include { INTEGRATE_DATASETS } from './modules/integrate_datasets.nf'
include { MULTIMODAL_INTEGRATION } from './modules/multimodal_integration.nf'
include { CLUSTERING } from './modules/clustering.nf'
include { CELL_ANNOTATION } from './modules/cell_annotation.nf'
include { BOOK_RENDER } from './modules/render.nf'




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/




workflow {



    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.pipeline
    )


    if (params.rna_normalization_method == 'sct' && params.run_scimilarity) {
      error "Error: SCimilarity cell annotation can only be performed on log transformed data."
    }

    if (params.run_scimilarity) {
      params.scimilarity_reference = params.annotations["scimilarity"]["annotation_path"]
    } else {
        params.scimilarity_reference = "."

    }


    book_assets = channel.fromPath( 'assets/book_assets/*')
                .mix(channel.fromPath(params.logo))
                .mix(channel.fromPath(params.quarto_index))
                .mix(channel.fromPath(params.quarto_yml))
                .collect()


    scripts_ch = channel.fromPath("bin/*")
                 .collect()


    

    params.annotation = getGenomeAttribute('annotation')
    params.mito_regex = getGenomeAttribute('mito_regex')
    params.ribo_regex = getGenomeAttribute('ribo_regex')
    params.bsgenome = getGenomeAttribute('bs_genome')
    params.scimilarity_dictionary = params.annotations["scimilarity"]["celltype_dictionary"]

    params.name = getGenomeAttribute('name')

    if(params.input_type == "counts") {
        Channel.fromPath(  "${params.input}/**/outs/filtered_*bc_matrix.h5", checkIfExists : true )
        .map { it -> tuple( it.parent.parent.name, it)}
        .set { filt_counts_ch }

        Channel.fromPath(  "${params.input}/**/outs/atac_fragments.tsv*", checkIfExists : true )
        | map { it -> tuple( it.parent.parent.name, it)}
        | groupTuple()
        | map { sample, files -> [sample , files[0], files[1]]}
        | set { fragments_ch }

        Channel.fromPath(  "${params.input}/**/outs/per_barcode_metrics.csv", checkIfExists : true )
        | map { it -> tuple( it.parent.parent.name, it)}
        | set { barcode_ch }

        Channel.fromPath(  "${params.input}/**/outs/raw_*bc_matrix.h5", checkIfExists : true )
        | map { it -> tuple( it.parent.parent.name, it)}
        | join(filt_counts_ch)
        | join(fragments_ch)
        | join(barcode_ch)
        | set { counts_ch }


        

        GET_ENSDB_ANNOTATION(params.annotation, params.get_annot, params.genome)
        | set { ens_annot_ch }
        
        if (['multiome', 'teaseq', 'atac'].contains(params.pipeline)) {
            CREATE_SEURAT_ATAC(counts_ch, params.pipeline, params.matrix, params.create_seurat,ens_annot_ch.annotation, book_assets)
            | set { seurat_ch }
            //CREATE_SEURAT_ATAC(conv_mat_ch.bp_mat.collect(), ens_annot_ch.annotation, params.pipeline, params.create_seurat)
            //| set { seurat_ch }
        } else {
            CREATE_SEURAT(counts_ch, params.pipeline, params.matrix, params.create_seurat, book_assets)
            | set { seurat_ch }
            //CREATE_SEURAT(conv_mat_ch.bp_mat.collect(), params.pipeline, params.create_seurat)
        }
        
    }

    count_pkg_ch = seurat_ch.count_pkg.collect()

    quarto_ch = seurat_ch.quarto.first()
   

    QUALITY_CONTROL(seurat_ch.seurat.collect(),
                    count_pkg_ch,
                    params.qc_script,
                    book_assets,
                    params.pipeline,
                    params.metadata,
                    params.genome,
                    params.mito_regex,
                    params.ribo_regex,
                    params.percent_mt,
                    params.percent_ribo,
                    params.nCount_RNA_min,
                    params.nCount_RNA_max,
                    params.nFeature_RNA_min,
                    params.nFeature_RNA_max,
                    params.nCount_ATAC_min,
                    params.nCount_ATAC_max,
                    params.nucleosome_signal_min,
                    params.nucleosome_signal_max,
                    params.TSS_enrichment,
                    params.atac_peak_region_fragments_min,
                    params.atac_peak_region_fragments_max,
                    params.pct_reads_in_peaks_min,
                    params.blacklist_fraction,
                    params.sample_exclusion,
                    params.umap1_ndims,
                    params.clustering1_res,
                    params.doublet_removal,
                    params.doublet_confidence,
                    params.sketch_cells,
                    params.sketch_n).set {qc_ch}
    quarto_ch = quarto_ch.concat(qc_ch.quarto)
 

    if (params.pipeline in ["multiome", 'teaseq', 'atac'] & params.recall_peaks) {
        CALL_PEAKS(qc_ch.seurat,
                    params.pipeline,
                    params.peak_call,
                    book_assets)
        | set {post_peaks_ch }
        post_qc_seurat = post_peaks_ch.seurat
        quarto_ch = quarto_ch.concat(post_peaks_ch.quarto)
    } else {
        post_qc_seurat = qc_ch.seurat
    }
    

    if (!['QC'].contains(params.stop_after)) {
        DIMENSION_REDUCTION(post_qc_seurat,
                           params.dimreduc_script,
                           book_assets,
                           params.pipeline,
                           params.umap2_ndims,
                           params.first_lsi_pc,
                           params.rna_normalization_method,
                           params.sketch_cells,
                           params.sketch_n)
        | set {umap_ch }
        quarto_ch = quarto_ch.concat(umap_ch.quarto)
    }
   
    if (!['QC', "DR"].contains(params.stop_after)) {
        if (params.integrate_datasets ) {
            INTEGRATE_DATASETS(umap_ch.seurat,
                            params.int_script,
                            book_assets,
                            params.pipeline,
                            params.integrate_datasets,
                            params.integration_method,
                            params.integrate_by,
                            params.view_batch,
                            params.umap2_ndims,
                            params.first_lsi_pc,
                            params.sketch_cells,
                            params.rna_normalization_method)
            | set { int_ch }
            quarto_ch = quarto_ch.concat(int_ch.quarto)
            
        } else {
            int_ch = umap_ch
        }
    }
    
    if (!['QC', "DR", 'INT'].contains(params.stop_after)) {
        if (params.pipeline in ["multiome", 'teaseq', 'cite'] & params.integrate_modalities) {
            MULTIMODAL_INTEGRATION(int_ch.seurat,
                                params.int_multimod_script,
                                book_assets,
                                params.pipeline,
                                params.umap2_ndims,
                                params.first_lsi_pc,
                                params.integrate_datasets,
                                params.sketch_cells)
            | set { intwnn_ch }
            quarto_ch = quarto_ch.concat(intwnn_ch.quarto)
            
        } else if (params.integrate_datasets){
            intwnn_ch = int_ch
        } else {
            intwnn_ch = umap_ch
        }
        
        CLUSTERING(intwnn_ch.seurat,
                params.clustering_script,
                book_assets,
                params.pipeline,
                params.clustering2_res,
                params.integrate_datasets,
                params.outcomes,
                params.sketch_cells,
                params.de_method,
                params.de_latent_vars,
                params.de_min_pct,
                params.de_logfc)
        | set { cluster_ch }

        quarto_ch = quarto_ch.concat(cluster_ch.quarto)
        
    }




    if (!['QC', "DR", 'INT', "CLUST"].contains(params.stop_after) ) {
        CELL_ANNOTATION(cluster_ch.seurat,
                        params.annotation_script,
                        book_assets,
                        params.pipeline,
                        params.integrate_datasets,
                        params.annotation_resolution,
                        params.run_azimuth,
                        params.azimuth_reference,
                        params.run_singler,
                        params.singler_reference,
                        params.run_scimilarity,
                        params.scimilarity_reference,
                        params.selected_method,
                        params.markers_rna,
                        params.markers_adt,
                        params.outcomes,
                        params.sketch_cells,
                        params.scimilarity_dictionary)
        | set { cell_ann_ch }

        quarto_ch = quarto_ch.concat(cell_ann_ch.quarto)
    }

    BOOK_RENDER(scripts_ch, 
                     book_assets,
                     quarto_ch.collect())
    
    
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}



def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

