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
include { CREATE_SEURAT_MULTIMODAL } from './modules/create_seurat.nf'
include { CREATE_SEURAT_GEX } from './modules/create_seurat.nf'
include { CREATE_SEURAT_ATAC } from './modules/create_seurat.nf'
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



    book_assets = channel.fromPath( 'assets/book_assets/*/*')
                .collect()
                .mix(channel.fromPath(params.logo))

    book_assets.view()

    //book_assets.view()
    scripts_ch = channel.fromPath("bin/*")
                 .collect()


    

    params.annotation = getGenomeAttribute('annotation')
    params.mito_regex = getGenomeAttribute('mito_regex')
    params.ribo_regex = getGenomeAttribute('ribo_regex')
    params.bsgenome = getGenomeAttribute('bs_genome')
    //params.genbank = getGenomeAttribute('genbank')
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


        CREATE_SEURAT_MULTIMODAL(counts_ch, ens_annot_ch.annotation, params.pipeline, params.create_seurat)
        | set { seurat_ch }

        if (params.matrix == "filtered") {
        seurat_ch.filt
        | collect 
        | set { seurat_combined_ch}
        
        } else {
            seurat_ch.raw
            | collect 
            | set { seurat_combined_ch}
        }

        seurat_ch.fragments
        | mix(seurat_ch.fragment_index)
        | collect
        | set { combined_fragments_ch }

    } else {
        if(params.matrix == "filtered") {
            Channel.fromPath(  "${params.input}/**/outs/filtered_matrix_seurat.rds", checkIfExists : true )
            .set { seurat_combined_ch }
        } else {
            Channel.fromPath(  "${params.input}/**/outs/raw_matrix_seurat.rds", checkIfExists : true )
            .set { seurat_combined_ch }
        }
    }
    

    /*
    QUALITY_CONTROL(seurat_combined_ch,
                    combined_fragments_ch,
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
                    params.doublet_confidence).set {qc_ch}
 
    quarto_ch = qc_ch.quarto


    if (!['QC'].contains(params.stop_after)) {
        DIMENSION_REDUCTION(qc_ch.seurat,
                           params.dimreduc_script,
                           book_assets,
                           params.pipeline,
                           params.umap2_ndims,
                           params.first_lsi_pc,
                           params.rna_normalization_method)
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
                            params.first_lsi_pc)
            | set { int_ch }
            quarto_ch = quarto_ch.concat(int_ch.quarto)
            
        } else {
            int_ch = umap_ch
        }
    }

    if (!['QC', "DR", 'INT'].contains(params.stop_after)) {
        if (params.pipeline in ["multiome", 'teaseq', 'cite', 'cite_crispr'] & params.integrate_modalities) {
            MULTIMODAL_INTEGRATION(int_ch.seurat,
                                params.int_multimod_script,
                                book_assets,
                                params.pipeline,
                                params.umap2_ndims,
                                params.first_lsi_pc,
                                params.integrate_datasets)
            | set { intwnn_ch }
            quarto_ch = quarto_ch.concat(intwnn_ch.quarto)
            
        } else {
            intwnn_ch = int_ch
        }

        CLUSTERING(intwnn_ch.seurat,
                params.clustering_script,
                book_assets,
                params.pipeline,
                params.clustering2_res,
                params.integrate_datasets,
                params.outcomes)
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
                        params.outcomes)
        | set { cell_ann_ch }

        quarto_ch = quarto_ch.concat(cell_ann_ch.quarto)
    }

    BOOK_RENDER(scripts_ch, 
                     book_assets,
                     quarto_ch.collect())
    */
    
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