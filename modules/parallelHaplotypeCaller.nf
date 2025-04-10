nextflow.enable.dsl = 2

/**
 * Create intervals for each chromosome and everything else.
 * - One .interval_list per standard chromosome (1..22, X, Y, M)
 * - One .interval_list named "extra.interval_list" for all other contigs
 */
process createChromosomeIntervals {
    if (params.platform == 'local') {
        label 'process_low'
    } else if (params.platform == 'cloud') {
        label 'process_low'
    } else {
        // For slurm
        label 'process_low'
    }
    container 'variantvalidator/gatk4:4.3.0.0'

    input:
    path reference_files // e.g. your .fasta (plus .fai / .dict if they exist)

    output:
    path "splitted_by_chrom/*.interval_list"

    script:
    """
    # 1) Identify reference FASTA
    genomeFasta=\$(find -L . -name '*.fasta')

    # 2) Ensure .fai index exists
    if [ ! -f "\${genomeFasta}.fai" ]; then
        echo "Creating FASTA index"
        samtools faidx \${genomeFasta}
    fi

    # 3) Ensure .dict exists
    dictFile=\${genomeFasta%.fasta}.dict
    if [ ! -f "\$dictFile" ]; then
        echo "Creating dictionary"
        gatk CreateSequenceDictionary -R \${genomeFasta} -O \$dictFile
    fi

    # 4) Make an output dir for chromosome-based interval lists
    mkdir splitted_by_chrom

    # We'll treat these as 'main' chromosomes (adjust as you like).
    # E.g. you can remove 'chrM' or add 'chrUn' etc. if desired.
    declare -a MAIN_CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \\
                            chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \\
                            chr20 chr21 chr22 chrX chrY chrM)

    # 5) Parse the .dict to make one interval_list per main chromosome,
    #    and an "extra.interval_list" for all other contigs.
    #    Each line @SQ SN:chrN LN:length ...
    #    We'll create a single line:  "chrN 1 length + chrN"
    #    (Picard-style interval_list).
    while read -r line; do
        # Only process lines starting with '@SQ'
        [[ ! "\$line" =~ ^@SQ ]] && continue

        # Example line: "@SQ  SN:chr1  LN:248956422 ..."
        # We'll extract SN=chrN and LN=length
        contig=\$(echo "\$line" | awk '{for(i=1;i<=NF;i++){ if(\$i ~ /^SN:/){sub("SN:","",\$i);print \$i} } }')
        length=\$(echo "\$line" | awk '{for(i=1;i<=NF;i++){ if(\$i ~ /^LN:/){sub("LN:","",\$i);print \$i} } }')

        # Output line in Picard interval_list format: contig <TAB> 1 <TAB> length <TAB> + <TAB> contig
        # Decide if contig is a "main" chromosome or "extra"
        if [[ " \${MAIN_CHROMS[@]} " =~ " \$contig " ]]; then
            echo -e "\$contig\\t1\\t\$length\\t+\\t\$contig" >> splitted_by_chrom/\${contig}.interval_list
        else
            # everything else lumps into one file
            echo -e "\$contig\\t1\\t\$length\\t+\\t\$contig" >> splitted_by_chrom/extra.interval_list
        fi
    done < "\$dictFile"
    """
}


/**
 * Call variants on each chromosome-based interval list
 */
process callVariantsOnInterval {
    if (params.platform == 'local') {
        label 'process_medium'
    } else if (params.platform == 'cloud') {
        label 'process_medium'
    } else {
        // For slurm
        label 'process_medium'
    }
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$sample_id - $interval.baseName"

    input:
    tuple val(sample_id), file(bamFile), file(bamIndex), path(interval)
    path reference_files

    output:
    tuple val(sample_id),
          file("${interval.baseName}.vcf"),
          file("${interval.baseName}.vcf.idx")

    script:
    """
    # Find the reference fasta
    genomeFasta=\$(find -L . -name '*.fasta')

    echo "Running HaplotypeCaller for Sample: ${sample_id} on interval: ${interval.baseName}"
    gatk --java-options "-Xmx${task.memory.toGiga()-2}g" HaplotypeCaller \\
        -R \${genomeFasta} \\
        -I ${bamFile} \\
        -O ${interval.baseName}.vcf \\
        -ERC GVCF \\
        -L ${interval} \\
        -A BaseQuality -A DepthPerSampleHC -A MappingQuality -A QualByDepth \\
        -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A StrandOddsRatio \\
        -A MappingQualityZero -A InbreedingCoeff -A BaseQualityRankSumTest -A HaplotypeFilteringAnnotation
    """
}


/**
 * Merge the output VCFs for each sample
 */
process mergeVcfs {
    if (params.platform == 'local') {
        label 'process_high'
    } else if (params.platform == 'cloud') {
        label 'process_high'
    } else {
        // For slurm
        label 'process_high'
    }
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcfs), path(indices)

    output:
    tuple val(sample_id), file("${sample_id}.vcf"), file("${sample_id}.vcf.idx")

    script:
    def vcfInputs = vcfs.collect { "-I $it" }.join(' ')
    """
    echo "Merging per-chromosome VCFs for Sample: ${sample_id}"
    gatk MergeVcfs \\
        ${vcfInputs} \\
        -O ${sample_id}.vcf
    """
}


/**
 * Top-level workflow
 */
workflow parallelHaplotypeCaller {
    take:
    bam_channel
    reference_files

    main:
    // 1) Create intervals for each chromosome + one 'extra'
    chr_intervals = createChromosomeIntervals(reference_files)

    // 2) Combine each BAM with each interval
    bam_intervals = bam_channel.combine(chr_intervals.flatten())

    // 3) HaplotypeCaller per chromosome/extra
    interval_vcfs = callVariantsOnInterval(bam_intervals, reference_files)

    // 4) Group all VCFs by sample
    grouped_vcfs = interval_vcfs.groupTuple()

    // 5) Merge per-sample VCF
    merged_vcfs = mergeVcfs(grouped_vcfs)

    emit:
    vcfs = merged_vcfs
}
