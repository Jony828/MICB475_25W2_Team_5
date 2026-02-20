#!/bin/bash

#### Data processing for Parkinson dataset ####

# Create a new directory for the gastric cancer dataset within data folder and navigate to it
(qiime2-amplicon-2025.4) root@stu-1130:/data# mkdir team5_parkinsons_data
(qiime2-amplicon-2025.4) root@stu-1130:/data# cd team5_parkinsons_data

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

# Taxonomic analysis (detached screen named taxonomy)
qiime feature-classifier classify-sklearn \
  --i-classifier /datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads parkinsons_rep-seqs.qza \
  --o-classification parkinsons_taxonomy.qza

qiime metadata tabulate \
  --m-input-file parkinsons_taxonomy.qza \
  --o-visualization parkinsons_taxonomy.qzv

# Removing mitochondria and chloroplast sequences
qiime taxa filter-table \
  --i-table parkinsons_table.qza \
  --i-taxonomy parkinsons_taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table parkinsons_table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
  --i-table parkinsons_table-no-mitochondria-no-chloroplast.qza \
  --o-visualization parkinsons_table-no-mitochondria-no-chloroplast.qzv \
  --m-sample-metadata-file /datasets/project_2/parkinsons/parkinsons_metadata.txt

# Generating a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences parkinsons_rep-seqs.qza \
  --o-alignment parkinsons_aligned-rep-seqs.qza \
  --o-masked-alignment parkinsons_masked-aligned-rep-seqs.qza \
  --o-tree parkinsons_unrooted-tree.qza \
  --o-rooted-tree parkinsons_rooted-tree.qza

# Alpha-rarefaction (detached screen named alpha_rarefaction)
qiime diversity alpha-rarefaction \
  --i-table  parkinsons_table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny parkinsons_rooted-tree.qza \
  --p-max-depth 20000 \
  --m-metadata-file /datasets/project_2/parkinsons/parkinsons_metadata.txt \
  --o-visualization parkinsons_alpha-rarefaction-no-mitochondria-no-chloroplast.qzv

# Calculate alpha- and beta-diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny parkinsons_rooted-tree.qza \
  --i-table parkinsons_table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 10232 \
  --m-metadata-file /datasets/project_2/parkinsons/parkinsons_metadata.txt \
  --output-dir core-metrics-results

# Calculate alpha-group-significance
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file /datasets/project_2/parkinsons/parkinsons_metadata.txt \
  --o-visualization core-metrics-results/parkinsons_shannon_group_significance.qzv

# Transferred parkinsons_shannon_group_significance.qzv to local computer

#### Data processing for GC dataset ####

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


