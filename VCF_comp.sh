#!/bin/bash
#SBATCH --job-name=NF_GATK
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --output=haplog.txt
#SBATCH --error=haperror.txt




REF="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/data/genome/hg38.fa"
FALSE_POSTIVES="/mnt/scratch1/projects/rudram_home/home_space/nf_compare/data/GIAB/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
TRUTH_VCF="/mnt/scratch1/projects/rudram_home/home_space/nf_compare/data/GIAB/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"


BEDFILE="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/data/bed/CCP19_sortmerged.bed"

#BWA VQSR
QUERY_VCF="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/BWA_VCF_QSR/NA12878_CCP19_combined_genotyped_recalibrated_dedup.vcf"

#DRAGMAP VQSR
# QUERY_VCF="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/DRAG_VCF_QSR/NA12878_CCP19_combined_genotyped_recalibrated_dedup.vcf"

#BWA
# BEDFILE="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/analysis/bwa_noqsr_hap/lowcov_bwa_noqsr/merged_filtered_exome.bed"

# DRAGMAP
#BEDFILE="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/analysis/dragmap_noqsr_hap/low_dragmap_qsr/merged_filtered_exome.bed"
apptainer exec \
  --bind /mnt/scratch1/projects/ \
  hap.py_latest.sif /opt/hap.py/bin/hap.py $TRUTH_VCF $QUERY_VCF \
  --false-positives $FALSE_POSTIVES \
  --restrict-regions $BEDFILE \
  --reference $REF \
  --report-prefix happy.output \
  --engine=vcfeval



# apptainer exec \
#   --bind /mnt/scratch1/projects/ \
#   hap.py_latest.sif /opt/hap.py/bin/hap.py $TRUTH_VCF $QUERY_VCF \
#   --restrict-regions $BEDFILE \
#   --reference $REF \
#   --report-prefix happy.output \
#   --engine=vcfeval





  # samtools depth -a -b /home/rudram/nf_panel/WRGLpipeline-files/WRGL5_hg19_v02.bed_hg38lift.bed NA12878C.recal.bam > coverage.txt
  # awk '$3 < 20' coverage.txt > low_coverage_regions.txt
  # ls
  # awk '{print $1"\t"$2-1"\t"$2}' low_coverage_regions.txt > low_coverage.bed
  # bedtools subtract -a /home/rudram/nf_panel/WRGLpipeline-files/WRGL5_hg19_v02.bed_hg38lift.bed -b low_coverage.bed > filtered_exome.bed
  # bedtools merge -i filtered_exome.bed > merged_filtered_exome.bed