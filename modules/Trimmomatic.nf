/*
 * ReadTrimming
 */

 /*
 * Run fastq on the read fastq files
 */
process TRIMMOMATIC{

    label 'process_high'

    container 'staphb/trimmomatic'

    // Add a tag to identify the process
    tag "$sample_id"

    // Specify the output directory for the FASTQC results
    publishDir("$params.outdir/Trimmomatic", pattern: "${sample_id}_trimmed_*.fq.gz", mode: "copy")

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_*.fq.gz")

    script:
    """
    echo "Running Trimmomatic"

    # Check the number of files in reads and run fastqc accordingly
    # Paired ends
    if [ -f "${reads[0]}" ] && [ -f "${reads[1]}" ]; then
        trimmomatic PE ${reads[0]} ${reads[1]} \\
        ${sample_id}_trimmed_P_R1.fq.gz ${sample_id}_trimmed_U_R1.fq.gz \\
        ${sample_id}_trimmed_P_R2.fq.gz ${sample_id}_trimmed_U_R2.fq.gz \\
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \\
        2> trimming.log
    # Single end
    elif [ -f "${reads[0]}" ]; then
        trimmomatic PE ${reads[0]} trimmed_${sample_id}.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 2> trimming.log
    else
        echo "No valid read files found for sample ${sample_id}"
        exit 1
    fi

    echo "Trimming complete Complete"
    """
}