/*
 * Align reads using DRAGMAP and process to BAM with read groups
 */

// First process: DRAGMAP alignment to produce SAM
process DRAGMAP_ALIGN_SAM {
    label 'process_high'
    container 'quay.io/biocontainers/dragmap:1.3.0--h5ca1c30_5'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)
    path dragmap_index
    
    output:
    tuple val(sample_id), path("${sample_id}.sam")
    
    script:
    """
    # Check if the input FASTQ files exist
    if [ -f "${reads[0]}" ]; then
        if [ \$(ls ${reads} | wc -l) -eq 2 ] && [ -f "${reads[1]}" ]; then
            # Paired-end mode
            dragen-os \\
                -r $dragmap_index \\
                -1 ${reads[0]} -2 ${reads[1]} > ${sample_id}.sam
        else
            # Single FASTQ mode
            dragen-os \\
                -r $dragmap_index \\
                -1 ${reads[0]} > ${sample_id}.sam
        fi
    else
        echo "Error: Read file ${reads[0]} does not exist for sample ${sample_id}."
        exit 1
    fi
    """
}

// Second process: Convert SAM to BAM and add read groups
process SAM_TO_BAM {
    label 'process_medium'
    container 'staphb/samtools:latest'  // Container with samtools
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(sam_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    samtools view -b $sam_file | \\
    samtools addreplacerg -r "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:illumina" - > ${sample_id}.bam
    """
}

// Combined workflow for DRAGMAP alignment
workflow DRAGMAP_ALIGN {
    take:
    read_pairs_ch
    dragmap_index_ch
    
    main:
    sam_ch = DRAGMAP_ALIGN_SAM(read_pairs_ch, dragmap_index_ch)
    bam_ch = SAM_TO_BAM(sam_ch)
    
    emit:
    bam_ch
}