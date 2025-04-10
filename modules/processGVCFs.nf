process combineGVCFs {
    if (params.platform == 'local') {
        label 'process_high'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/gatk4:4.3.0.0'
    tag "${sample_ids.join('_')}" // Add a tag based on the sample IDs

    input:
    tuple val(sample_ids), path(gvcf_files), path(gvcf_index_files)
    // path indexFiles
    path genomeFasta
    path genomeFai
    path genomeDict

    output:
    tuple val("${sample_ids.join('_')}"), file("*_combined.vcf"), file("*_combined.vcf.idx")

    script:
    def merged_sample_id = "${sample_ids.join('_')}"
    def gvcf_files_args = gvcf_files.collect { file -> "-V ${file}" }.join(' ')

    """
    echo "Combining GVCFs for samples: ${gvcf_files.collect { it.baseName }.join(', ')}"

    if [[ -n ${params.genome_file} ]]; then
        genomeFasta=${genomeFasta}
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"
    echo "Genom Fai : ${genomeFai}"
    echo "Geome dict : ${genomeDict}"

    gatk CombineGVCFs -R "\${genomeFasta}"\
        ${gvcf_files_args} \
        -O ${merged_sample_id}_combined.vcf
    """
}

process genotypeGVCFs {
    if (params.platform == 'local') {
        label 'process_high'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'variantvalidator/gatk4:4.3.0.0'
    tag "$combined_sample_id"

    input:
    tuple val(combined_sample_id), file(combined_gvcf), file(combined_gvcf_idx)
    // path indexFiles
    path genomeFasta
    path genomeFai
    path genomeDict

    output:
    tuple val(combined_sample_id), file("*_genotyped.vcf"), file("*_genotyped.vcf.idx")

    script:
    def merged_sample_id = combined_gvcf.baseName

    """
    echo "Genotyping combined GVCF: ${combined_gvcf.baseName}"

    if [[ -n ${params.genome_file} ]]; then
        genomeFasta=${genomeFasta}
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"
    echo "Genom Fai : ${genomeFai}"
    echo "Geome dict : ${genomeDict}"

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    gatk GenotypeGVCFs -R "\${genomeFasta}" \
        -V ${combined_gvcf} \
        -O ${merged_sample_id}_genotyped.vcf

    """
}