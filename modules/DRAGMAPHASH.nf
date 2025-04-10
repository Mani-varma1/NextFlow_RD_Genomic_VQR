process DRAGMAP_HASHTABLE {
    label 'process_high'

    container 'quay.io/biocontainers/dragmap:1.3.0--h5ca1c30_5'


    publishDir("$params.outdir/GENOME_IDX/", mode: "copy")
    
    input:
    path genomeFasta

    output:
    path "DRAG"


    script:
    """
    mkdir DRAG
    # Split the table into chunks (uses less memory) and limit memory usage
    dragen-os \\
        --build-hash-table true \\
        --ht-reference ${genomeFasta} \\
        --output-directory DRAG \\
        --num-threads ${task.cpus} \\
        --ht-max-table-chunks 8

    ## Copy the reference FASTA into the hash directory if it's not already there
    if [ ! -f "DRAG/${genomeFasta}" ]; then
        cp ${genomeFasta} DRAG/${genomeFasta}
    fi
    """


}
