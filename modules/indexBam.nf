/*
 * Index the BAM files
 */
process indexBam {

    if (params.platform == 'local') {
        label 'process_medium'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/indexgenome:1.1.0'
 
    tag "$bamFile"

    // Publish indexed BAM files to the specified directory
    publishDir("$params.outdir/BAM", mode: "copy")

    input:
    tuple val(sample_id), file(bamFile)

    output:
    tuple val(sample_id), file("${bamFile}"), file("${bamFile}.bai")

    script:
    """
    echo "Running Index Bam"
 
    # Use samtools to create the index for the input BAM file
    samtools index -@ $task.cpus ${bamFile} ${bamFile}.bai

    echo "Bam Indexing Complete"
    """
}

/*
 * Index the BAM files
 */
process indexMapDamageBam {

    label 'process_low'
    container 'variantvalidator/indexgenome:1.1.0'

    tag "$bamFile"

    // Publish indexed BAM files to the specified directory
    publishDir("$params.outdir/BAM", mode: "copy")

    input:
    tuple val(sample_id), file(bamFile)

    output:
    tuple val(sample_id), file("${bamFile}"), file("${bamFile}.bai")

    script:
    """
    echo "Running Index Bam"

    # Use samtools to create the index for the input BAM file
    samtools index ${bamFile} ${bamFile}.bai

    echo "Bam Indexing Complete"
    """
}
