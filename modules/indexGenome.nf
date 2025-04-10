/*
 * Define the indexGenome process that creates a BWA index
 * given the genome fasta file
 */
process indexGenome {

    if (params.platform == 'local') {
        label 'process_high'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/indexgenome:1.1.0'


    // Publish indexed files to the specified directory
    publishDir("$params.outdir/GENOME_IDX/BWA/", mode: "copy")

    input:
    path genomeFasta

    output:
    tuple path(genomeFasta), path("${genomeFasta}.*")

    script:
    """
    echo "Running Index Genome"

    # Generate BWA index
    bwa index "${genomeFasta}"

    # Generate samtools faidx
    samtools faidx "${genomeFasta}"

    # Generate Fasta dict
    picard CreateSequenceDictionary R="${genomeFasta}" O="${genomeFasta}.dict"

    echo "Genome Indexing complete."
    """
}
