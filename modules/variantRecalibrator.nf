process variantRecalibrator {

    if (params.platform == 'local') {
        label 'process_high'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    }
    container 'broadinstitute/gatk:4.6.1.0'

    tag "$vcf"

    publishDir("$params.outdir/VCF", mode: "copy")

    input:
    tuple val(sample_id), file(vcf), file(vcfIndex)
    val knownSitesArgs
    // path indexFiles
    path genomeFasta
    path genomeFai
    path genomeDict
    path qsrc_vcf

    output:
    tuple val(sample_id), file("${vcf.baseName}.recalibrated.vcf")

    script:
    def knownSitesArgsStr = knownSitesArgs.join(' ')
    def degradedDna = params.degraded_dna == "true"

    """
    echo "Running VQSR"

    if [[ -n ${params.genome_file} ]]; then
        genomeFasta=${genomeFasta}
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"
    echo "Genom Fai : ${genomeFai}"
    echo "Geome dict : ${genomeDict}"

    if ${degradedDna}; then
        echo "Running VQSR for degraded DNA (1x coverage)"
        # relaxed parameters for SNPs and INDELs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an FS -an SOR \
            -mode SNP \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -O ${vcf.baseName}.recalibrated_SNP.recal
        # Apply VQSR for SNPs
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -recal-file ${vcf.baseName}.recalibrated_SNP.recal \
            -mode SNP \
            -O ${vcf.baseName}.output_SNP.vcf
        # relaxed parameters for INDELs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an FS -an SOR \
            -mode INDEL \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -O ${vcf.baseName}.recalibrated_INDEL.recal
        # Apply VQSR for INDELs
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -recal-file ${vcf.baseName}.recalibrated_INDEL.recal \
            -mode INDEL \
            -O ${vcf.baseName}.output_INDEL.vcf
    else
        echo "Running VQSR for standard DNA (10x+ coverage)"
        # stricter parameters for SNPs and INDELs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            -mode SNP \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -O ${vcf.baseName}.recalibrated_SNP.recal
        # Apply VQSR for SNPs
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_SNP.tranches \
            -recal-file ${vcf.baseName}.recalibrated_SNP.recal \
            -mode SNP \
            -O ${vcf.baseName}.output_SNP.vcf
        # stricter parameters for INDELs
        gatk VariantRecalibrator \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            ${knownSitesArgsStr} \
            -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
            -mode INDEL \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -O ${vcf.baseName}.recalibrated_INDEL.recal
        # Apply VQSR for INDELs
        gatk ApplyVQSR \
            -R "\${genomeFasta}" \
            -V ${vcf} \
            --truth-sensitivity-filter-level 99.0 \
            -tranches-file ${vcf.baseName}.recalibrated_INDEL.tranches \
            -recal-file ${vcf.baseName}.recalibrated_INDEL.recal \
            -mode INDEL \
            -O ${vcf.baseName}.output_INDEL.vcf
    fi

    gatk MergeVcfs \
        -I ${vcf.baseName}.output_SNP.vcf \
        -I ${vcf.baseName}.output_INDEL.vcf \
        -O ${vcf.baseName}.recalibrated.vcf

    echo "VQSR Complete"
    """
}
