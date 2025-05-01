#!/usr/bin/env python3
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import random
from tqdm import tqdm
from collections import defaultdict
import matplotlib.ticker as ticker

# Set the aesthetic style for plots
plt.style.use('ggplot')
sns.set_theme(style="whitegrid", font_scale=1.2)
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'DejaVu Sans'

def extract_data_from_bam(bam_file, num_reads=100000):
    """
    Extract quality scores, mapping quality, and edit distances from BAM file
    """
    qualities = []
    mapping_quals = []
    edit_distances = []
    
    print(f"Processing {bam_file}...")
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Sample random reads from the file
        try:
            all_reads = list(bam.fetch(until_eof=True))
        except ValueError:
            print(f"BAM file {bam_file} does not have an index. Creating temporary index...")
            pysam.index(bam_file)
            with pysam.AlignmentFile(bam_file, "rb") as indexed_bam:
                all_reads = list(indexed_bam.fetch(until_eof=True))
        
        if len(all_reads) > num_reads:
            sampled_reads = random.sample(all_reads, num_reads)
        else:
            sampled_reads = all_reads
            print(f"Warning: Only {len(all_reads)} reads found in {bam_file}")
        
        for read in tqdm(sampled_reads, desc=f"Extracting from {os.path.basename(bam_file)}"):
            # Skip unmapped reads
            if read.is_unmapped:
                continue
                
            # Extract quality scores for each base in the read
            read_quals = read.query_qualities
            if read_quals:
                qualities.extend(read_quals)
            
            # Extract mapping quality
            mapping_quals.append(read.mapping_quality)
            
            # Extract edit distance (NM tag)
            if read.has_tag('NM'):
                edit_distances.append(read.get_tag('NM'))
    
    return qualities, mapping_quals, edit_distances

def create_dataframes(data_dict):
    """
    Create DataFrames for plotting from the extracted data
    """
    quality_data = []
    mapq_data = []
    edit_data = []
    
    for source, data in data_dict.items():
        qualities, mapping_quals, edit_distances = data
        
        # Add quality scores to DataFrame
        for qual in qualities:
            quality_data.append({
                'Base Quality': qual,
                'Source': source
            })
        
        # Add mapping qualities to DataFrame
        for mapq in mapping_quals:
            mapq_data.append({
                'Mapping Quality': mapq,
                'Source': source
            })
        
        # Add edit distances to DataFrame
        for edit_dist in edit_distances:
            edit_data.append({
                'Edit Distance': edit_dist,
                'Source': source
            })
    
    # Create DataFrames
    quality_df = pd.DataFrame(quality_data)
    mapq_df = pd.DataFrame(mapq_data)
    edit_df = pd.DataFrame(edit_data)
    
    # Add Aligner and BQSR columns
    for df in [quality_df, mapq_df, edit_df]:
        df['Aligner'] = df['Source'].apply(lambda x: x.split('_')[0])
        df['BQSR'] = df['Source'].apply(lambda x: 'After BQSR' if 'after' in x.lower() else 'Before BQSR')
    
    # Create quality range column
    quality_df['Quality Range'] = pd.cut(
        quality_df['Base Quality'], 
        bins=[0, 10, 20, 30, 40, 100], 
        labels=['0-10', '11-20', '21-30', '31-40', '41+']
    )
    
    return quality_df, mapq_df, edit_df

def plot_quality_violin(quality_df, output_dir):
    """
    Create split violin plot of base quality scores
    """
    plt.figure(figsize=(14, 9))
    ax = sns.violinplot(
        data=quality_df, 
        x='Aligner', 
        y='Base Quality', 
        hue='BQSR',
        split=True, 
        palette='viridis',
        inner='quartile'
    )
    plt.title('Distribution of Base Quality Scores Before and After BQSR', fontsize=18)
    plt.xlabel('Alignment Method', fontsize=16)
    plt.ylabel('Base Quality Score (Phred)', fontsize=16)
    plt.legend(title='BQSR Status', fontsize=14, title_fontsize=14)
    plt.grid(True, alpha=0.3)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=20))
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'quality_violin.png'))
    plt.close()

def plot_quality_ranges(quality_df, output_dir):
    """
    Create bar plot of quality score ranges
    """
    # Calculate percentage in each range
    range_stats = quality_df.groupby(['Aligner', 'BQSR', 'Quality Range']).size().reset_index(name='Count')
    total_counts = range_stats.groupby(['Aligner', 'BQSR'])['Count'].sum().reset_index(name='Total')
    range_stats = range_stats.merge(total_counts, on=['Aligner', 'BQSR'])
    range_stats['Percentage'] = (range_stats['Count'] / range_stats['Total']) * 100
    
    # Plot
    plt.figure(figsize=(16, 10))
    g = sns.catplot(
        data=range_stats,
        x='Quality Range',
        y='Percentage',
        hue='BQSR',
        col='Aligner',
        kind='bar',
        height=8,
        aspect=1.2,
        palette=['#2a9d8f', '#e76f51'],
        alpha=0.8
    )
    g.set_axis_labels('Quality Score Range', 'Percentage of Bases (%)')
    g.set_titles('{col_name}')
    g.figure.suptitle('Distribution of Quality Scores Before and After BQSR', fontsize=18, y=0.98)
    
    # Add percentage labels on bars
    for ax in g.axes.flat:
        for container in ax.containers:
            ax.bar_label(container, fmt='%.1f%%', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'quality_ranges.png'))
    plt.close()

def plot_mapping_quality_violin(mapq_df, output_dir):
    """
    Create split violin plot of mapping qualities
    """
    plt.figure(figsize=(14, 9))
    ax = sns.violinplot(
        data=mapq_df, 
        x='Aligner', 
        y='Mapping Quality', 
        hue='BQSR',
        split=True, 
        palette='magma',
        inner='quartile'
    )
    plt.title('Distribution of Mapping Quality Before and After BQSR', fontsize=18)
    plt.xlabel('Alignment Method', fontsize=16)
    plt.ylabel('Mapping Quality', fontsize=16)
    plt.legend(title='BQSR Status', fontsize=14, title_fontsize=14)
    plt.grid(True, alpha=0.3)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mapping_quality_violin.png'))
    plt.close()

def plot_edit_distance_histogram(edit_df, output_dir):
    """
    Create histogram plot of edit distances
    """
    plt.figure(figsize=(16, 10))
    g = sns.FacetGrid(edit_df, col='Aligner', row='BQSR', height=4.5, aspect=1.5)
    g.map_dataframe(sns.histplot, x='Edit Distance', binwidth=1, kde=False, color='#5a189a')
    
    # Customize the plot
    g.set_axis_labels('Edit Distance (NM)', 'Count')
    g.set_titles('{row_name} - {col_name}')
    g.figure.suptitle('Distribution of Edit Distances Before and After BQSR', fontsize=18, y=0.98)
    
    # Ensure integer x-axis ticks
    for ax in g.axes.flat:
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'edit_distance_histogram.png'))
    plt.close()

def create_quality_stats_table(quality_df, output_dir):
    """
    Create a table with mean Phred scores and 25-75 quantiles
    """
    # Calculate statistics
    stats = quality_df.groupby(['Aligner', 'BQSR'])['Base Quality'].agg([
        ('Mean', 'mean'),
        ('25th Percentile', lambda x: np.percentile(x, 25)),
        ('75th Percentile', lambda x: np.percentile(x, 75))
    ]).reset_index()
    
    # Round values to 2 decimal places
    for col in ['Mean', '25th Percentile', '75th Percentile']:
        stats[col] = stats[col].round(2)
    
    # Save as CSV
    stats.to_csv(os.path.join(output_dir, 'quality_stats_table.csv'), index=False)
    
    # Print table for immediate view
    print("\nBase Quality Statistics:")
    print(stats)
    
    return stats

def main():
    # Define BAM files
    bam_files = {
        'BWA_before': '/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/BWA_BAM_NOQSR/NA12878_CCP19_filtered_sorted_dedup.bam',
        'BWA_after': '/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/BWA_BAM_QSR/NA12878_CCP19_filtered_sorted_dedup_recalibrated.bam',
        'Dragmap_before': '/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/DRAG_BAM_QSR/NA12878_CCP19_filtered_sorted_dedup.bam',
        'Dragmap_after': '/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/DRAG_BAM_QSR/NA12878_CCP19_filtered_sorted_dedup_recalibrated.bam'
    }
    
    # Number of reads to sample
    num_reads = 100000
    
    # Output directory for plots
    output_dir = '/mnt/scratch1/projects/rudram_home/home_space/NextFlow_RD_Genomic_VQR/results/analysis/bqsr_claude'
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract data from each BAM file
    extracted_data = {}
    for label, file_path in bam_files.items():
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} not found. Please provide the correct path.")
            continue
            
        qualities, mapping_quals, edit_distances = extract_data_from_bam(file_path, num_reads)
        extracted_data[label] = (qualities, mapping_quals, edit_distances)
    
    # Create DataFrames
    quality_df, mapq_df, edit_df = create_dataframes(extracted_data)
    
    # Generate plots
    print("Creating plots...")
    plot_quality_violin(quality_df, output_dir)
    plot_quality_ranges(quality_df, output_dir)
    plot_mapping_quality_violin(mapq_df, output_dir)
    plot_edit_distance_histogram(edit_df, output_dir)
    
    # Generate quality statistics table
    create_quality_stats_table(quality_df, output_dir)
    
    print(f"Analysis complete! All visualizations saved to '{output_dir}' directory.")

if __name__ == "__main__":
    main()