process DRAGMAP_HASHTABLE {
    label 'process_high'

    container 'gambalab/dragmap'

    publishDir("$params.outdir/dragmap", mode: "copy")
    input:
    path genomeFasta

    output:
    path "dragmap/*"


    script:

    """
    mkdir dragmap
    dragen-os \\
        --build-hash-table true \\
        --ht-reference $genomeFasta \\
        --output-directory dragmap \\
        --num-threads 8
    """


}
