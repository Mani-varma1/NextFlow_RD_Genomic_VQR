process DRAGMAP_ALIGN {
    label 'process_single'

    container 'quay.io/biocontainers/dragmap:1.3.0--h91baf5a_3'

    input:
    path(reads)
    path(hashmap)

    output:
    path("*.sam") 


    script:

    """
    dragen-os \\
        -r $hashmap \\
        "$reads[1]" "$reads[2]" > sample.sam\\
    """

}