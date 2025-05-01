
#!/bin/bash

# Define your BAM files from two different aligners
ALIGNER1_BAM="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/BWA_BAM_NOQSR/NA12878_CCP19_filtered_sorted_dedup.bam"
ALIGNER2_BAM="/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/DRAG_BAM_NOQSR/NA12878_CCP19_filtered_sorted_dedup.bam"

echo "===== Simple Alignment Statistics Comparison ====="
echo "Metric                 | ALIGNER 1      | ALIGNER 2"
echo "-----------------------|----------------|----------------"

# Total reads
TOTAL1=$(samtools flagstat $ALIGNER1_BAM | grep 'in total' | cut -d' ' -f1)
TOTAL2=$(samtools flagstat $ALIGNER2_BAM | grep 'in total' | cut -d' ' -f1)
echo "Total reads            | $TOTAL1 | $TOTAL2"

# Mapped reads
MAPPED1=$(samtools flagstat $ALIGNER1_BAM | grep 'mapped (' | head -n1 | cut -d' ' -f1)
MAPPED2=$(samtools flagstat $ALIGNER2_BAM | grep 'mapped (' | head -n1 | cut -d' ' -f1)
MAPPED_PCT1=$(samtools flagstat $ALIGNER1_BAM | grep 'mapped (' | head -n1 | grep -o '[0-9]*\.[0-9]*%')
MAPPED_PCT2=$(samtools flagstat $ALIGNER2_BAM | grep 'mapped (' | head -n1 | grep -o '[0-9]*\.[0-9]*%')
echo "Mapped reads           | $MAPPED1 ($MAPPED_PCT1) | $MAPPED2 ($MAPPED_PCT2)"

# Primary mapped reads (mapped - secondary - supplementary)
SECONDARY1=$(samtools flagstat $ALIGNER1_BAM | grep 'secondary' | cut -d' ' -f1)
SECONDARY2=$(samtools flagstat $ALIGNER2_BAM | grep 'secondary' | cut -d' ' -f1)
SUPPLEMENTARY1=$(samtools flagstat $ALIGNER1_BAM | grep 'supplementary' | cut -d' ' -f1)
SUPPLEMENTARY2=$(samtools flagstat $ALIGNER2_BAM | grep 'supplementary' | cut -d' ' -f1)
PRIMARY1=$((MAPPED1 - SECONDARY1 - SUPPLEMENTARY1))
PRIMARY2=$((MAPPED2 - SECONDARY2 - SUPPLEMENTARY2))
PRIMARY_PCT1=$(awk "BEGIN {printf \"%.2f%%\", 100*$PRIMARY1/$TOTAL1}")
PRIMARY_PCT2=$(awk "BEGIN {printf \"%.2f%%\", 100*$PRIMARY2/$TOTAL2}")
echo "Primary alignments     | $PRIMARY1 ($PRIMARY_PCT1) | $PRIMARY2 ($PRIMARY_PCT2)"

# Secondary alignments
SECONDARY_PCT1=$(awk "BEGIN {printf \"%.2f%%\", 100*$SECONDARY1/$TOTAL1}")
SECONDARY_PCT2=$(awk "BEGIN {printf \"%.2f%%\", 100*$SECONDARY2/$TOTAL2}")
echo "Secondary alignments   | $SECONDARY1 ($SECONDARY_PCT1) | $SECONDARY2 ($SECONDARY_PCT2)"

# Properly paired
PROPER1=$(samtools flagstat $ALIGNER1_BAM | grep 'properly paired' | cut -d' ' -f1)
PROPER2=$(samtools flagstat $ALIGNER2_BAM | grep 'properly paired' | cut -d' ' -f1)
PROPER_PCT1=$(samtools flagstat $ALIGNER1_BAM | grep 'properly paired' | grep -o '[0-9]*\.[0-9]*%')
PROPER_PCT2=$(samtools flagstat $ALIGNER2_BAM | grep 'properly paired' | grep -o '[0-9]*\.[0-9]*%')
echo "Properly paired        | $PROPER1 ($PROPER_PCT1) | $PROPER2 ($PROPER_PCT2)"