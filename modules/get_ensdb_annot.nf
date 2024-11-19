process GET_ENSDB_ANNOTATION {
    tag "Fetching EnsDB annotation for ID... (${annotation})"


    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.Rmd"



    memory '8 GB'
    cpus 1

    input:
    val annotation
    path rmd
    val genome


    output:
    path "signac_annotation.RDS", emit:annotation
    path rmd, emit:script



    script:
    """
    Rscript -e 'rmarkdown::render("${rmd}", params = list(AH_ID = "${annotation}", genome = "${genome}"))'
    """

}