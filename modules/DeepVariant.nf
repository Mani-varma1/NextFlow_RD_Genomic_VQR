nextflow.enable.dsl = 2

process deepvariant {
    if (params.platform == 'local') {
        label 'process_high'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    }
    container 'google/deepvariant'

    input:
    tuple val(sample_id), path(bam_files)
    path reference

    output:
    path("*")

    script:
    def bam_file = bam_files[0]
    def bai_file = bam_files[1]
    """
    outputVcf="${sample_id}.vcf"
    outputgVcf="${sample_id}.gvcf"

    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --ref=${reference} \
        --reads=${bam_file} \
        --output_vcf=\${outputVcf} \
        --output_gvcf=\${outputgVcf} \
        --vcf_stats_report=true

    echo "Sample: ${sample_id} VCF: \${outputVcf}"
    echo "Variant Calling Complete"
    """
}

