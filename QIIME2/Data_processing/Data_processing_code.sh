#!/bin/bash

# Create a new directory for team 6
(qiime2-amplicon-2025.4) root@stu-1130:/data# mkdir team5_parkinsons_data

# Go to /data/team5_parkinsons_data

# Importing and demultiplex data using manifest (detached screen named parkinsons-import)
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/parkinsons/parkinsons_manifest.tsv \
  --output-path /data/team5_parkinsons_data/team5_parkinsons_demux_seqs.qza

# Creating a visualization of the demultiplexed samples
qiime demux summarize \
  --i-data team5_parkinsons_demux_seqs.qza \
  --o-visualization team5_parkinsons_demux_seqs.qzv 

# Transferred team5_parkinsons_demux_seqs.qzv onto local computer and visualized on view.QIIME2.org

# Determine ASVs with DADA2 (detached screen named denoising)
# Creating representative sequences and table
qiime dada2 denoise-single \
  --i-demultiplexed-seqs team5_parkinsons_demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 251 \
  --o-representative-sequences parkinsons_rep-seqs.qza \
  --o-table parkinsons_table.qza \
  --o-denoising-stats denoising-stats.qza

# Visualizing ASVs stats
qiime feature-table summarize \
  --i-table parkinsons_table.qza \
  --o-visualization parkinsons_table.qzv \
  --m-sample-metadata-file /datasets/project_2/parkinsons/parkinsons_metadata.txt

qiime feature-table tabulate-seqs \
  --i-data parkinsons_rep-seqs.qza \
  --o-visualization parkinsons_rep-seqs.qzv

# Data processing for GC dataset

# Create a new directory for the gastric cancer dataset within data folder and navigate to it
(qiime2-amplicon-2025.4) root@stu-1130:/data# mkdir team5_gc_data
(qiime2-amplicon-2025.4) root@stu-1130:/data# cd team5_gc_data

# Import and demultiplex gastric cancer sequencing data using a manifest file
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/gastric_cancer/gastric_cancer_manifest.tsv \
  --output-path team5_gc_demux_seqs.qza

# Create a visualization of the demultiplexed gastric cancer sequences
qiime demux summarize \
  --i-data team5_gc_demux_seqs.qza \
  --o-visualization team5_gc_demux_seqs.qzv

# Transferred team5_gc_demux_seqs.qzv to local computer and visualized quality metrics using view.qiime2.org

# Based on high quality scores, determined that trimming is not required, confirmed with Ritu, and ran denoising:
qiime dada2 denoise-single \
  --i-demultiplexed-seqs team5_gc_demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --o-representative-sequences gc_rep-seqs.qza \
  --o-table gc_table.qza \
  --o-denoising-stats gc-denoising-stats.qza

# Visualizing ASVs stats
qiime feature-table summarize \
  --i-table gc_table.qza \
  --o-visualization gc_table.qzv \
  --m-sample-metadata-file /datasets/project_2/gastric_cancer/gastric_cancer_metadata.tsv

  qiime feature-table tabulate-seqs \
  --i-data gc_rep-seqs.qza \
  --o-visualization gc_rep-seqs.qzv


