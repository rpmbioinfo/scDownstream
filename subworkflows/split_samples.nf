


def create_sample_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.assay        = row.Assay

    def fastq_meta = []


    fastq_meta = [meta.id]
    
    return fastq_meta
}

def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.assay        = row.Assay

    def fastq_meta = []

    
    fastq_meta = [meta.id, meta.assay, file(row.file) ]
    
    return fastq_meta
}

def sample_file(LinkedHashMap i) {
    // create meta map
    def meta = [:]
    meta.sample        = tuple[0]
    meta.assay        = tuple[1]
    meta.files          = tuple[2]
    
    
    
    return meta
}

def criteria = multiMapCriteria {
    gex: [it[1], it[0]['Assay'] == "Gene Expression"]
    atac: [it[1], it[0]['Assay'] == "Chromatin Accessbility"]
}

workflow SPLIT_SAMPLES {
    take: 
    samplesheet_ch


    main:
        samplesheet_ch
        | splitCsv(header:true, sep:',')
        | map { row -> 
            meta = row.subMap('sample', 'Assay')
            [meta , 
                file(row.file, checkIfExists: true)]

        }
        | groupTuple()
        | branch { meta, files -> 
            atac: meta.Assay == "Chromatin Accessibility"
            rna: meta.Assay == "Gene Expression"
            adt: meta.Assay == "Antibody Capture"
            image: meta.Assay == "Image"
        }
        | set { samples }


        //.map { meta, files -> [sample : meta.sample, Assay: meta.Assay, files: files]}
        

        /*

        .map{  create_fastq_channel ( it) }

        . map { row -> 
            meta = row.subMap('sample', 'Assay')
            [meta , 
                file(row.file, checkIfExists: true)]

        }
        .groupTuple(by:  0)
        .map { meta, files - > ['sample' : meta.sample, ["Assay" : meta.Assay, [files]]]}
        .view( it -> it[0])
      
        .set { ch1 }
        ch1.gex.view()

        
        //.map { create_sample_channel( it ) }
        //.unique()
        //.set { samplels_ch }


        samplesheet_ch
        .splitCsv(header:true, sep:',')
        .map { create_fastq_channel( it ) }
        .groupTuple()
        .multiMap(criteria).set { ch1 }
        
        ch1.gex.view()
        */

        //.groupTuple(by: [0,1])
        //.map { sample, assay, files -> ['sample': sample, 'assay': assay, "files":files] }
        //.groupTuple(by:  [0,1], sort: true)
        //.groupTuple(by:0 )
        //.map { sample, assay, files -> [ sample, assay.flatten(), files.flatten()]}
        //.map { sample, assay, files -> [ sample , ['assay' : assay, 'files' : files]]}
        //.groupTuple(by:0)
        .set { reads}

    emit:
    samples

}




