#!/usr/bin/env python3
"""
Simple VCF Statistics Calculator
-------------------------------
Extracts Ti/Tv ratios and variant counts from VCF files.

Usage:
    Edit the file paths below and run: python simple_vcf_stats.py

Requirements:
    - cyvcf2 (or PyVCF)
"""

###############################
# EDIT THESE FILE PATHS
###############################
VCF_FILES = [
    "/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/BWA_VCF_NOQSR/NA12878_CCP19_combined_genotyped_filtered.vcf",      # BWA without VQSR/BQSR
    "/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/BWA_VCF_QSR/NA12878_CCP19_combined_genotyped_recalibrated_dedup.vcf",    # BWA with VQSR/BQSR
    "/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/DRAG_VCF_NOQSR/NA12878_CCP19_combined_genotyped_filtered.vcf",  # DRAGMAP without VQSR/BQSR
    "/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/DRAG_VCF_QSR/NA12878_CCP19_combined_genotyped_recalibrated_dedup.vcf" # DRAGMAP with VQSR/BQSR
]

# Labels for each VCF file
LABELS = [
    "BWA with VQSR/BQSR OFF",
    "BWA with VQSR/BQSR ON",
    "DRAGMAP with VQSR/BQSR OFF",
    "DRAGMAP with VQSR/BQSR ON"
]
###############################

import sys
import cyvcf2  # Much faster than PyVCF for large files

# Define nucleotide changes
TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

def is_transition(ref, alt):
    """Check if a variant is a transition (vs transversion)"""
    if len(ref) != 1 or len(alt) != 1:
        return None  # Not a simple SNP
    return (ref, alt) in TRANSITIONS

def analyze_vcf(vcf_file):
    """Calculate basic statistics for a VCF file"""
    # Initialize counters
    total_variants = 0
    snps = 0
    indels = 0
    transitions = 0
    transversions = 0
    
    # Process variants
    print(f"Analyzing {vcf_file}...")
    vcf_reader = cyvcf2.VCF(vcf_file)
    
    for variant in vcf_reader:
        total_variants += 1
        
        # Get type information
        if variant.is_snp:
            snps += 1
            
            # Only consider first alternative allele for simplicity
            ref = variant.REF
            alt = variant.ALT[0]
            
            if is_transition(ref, alt):
                transitions += 1
            else:
                transversions += 1
                
        elif variant.is_indel:
            indels += 1
    
    # Calculate Ti/Tv ratio
    titv_ratio = transitions / transversions if transversions > 0 else float('nan')
    
    return {
        'total_variants': total_variants,
        'snps': snps,
        'indels': indels,
        'transitions': transitions,
        'transversions': transversions,
        'titv_ratio': titv_ratio
    }

def main():
    print("\n{:<25} {:<15} {:<15} {:<15} {:<15}".format(
        "Pipeline Configuration", "Total Variants", "SNPs", "Indels", "Ti/Tv Ratio"))
    print("-" * 85)
    
    for i, vcf_file in enumerate(VCF_FILES):
        label = LABELS[i] if i < len(LABELS) else f"VCF {i+1}"
        
        try:
            stats = analyze_vcf(vcf_file)
            
            print("{:<25} {:<15,} {:<15,} {:<15,} {:<15.2f}".format(
                label,
                stats['total_variants'],
                stats['snps'],
                stats['indels'],
                stats['titv_ratio']
            ))
        except Exception as e:
            print(f"{label}: Error analyzing file - {str(e)}")
    
    print("\n")
    print("Note: For whole exome data, a Ti/Tv ratio of ~2.0-2.1 is expected.")
    print("Significant deviations may indicate issues with variant calling.")

if __name__ == "__main__":
    main()