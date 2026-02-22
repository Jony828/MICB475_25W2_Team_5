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

# Transferred arkinsons_table-no-mitochondria-no-chloroplast.qzv onto local computer and visualized on view.QIIME2.org

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

# Transferred parkinsons_alpha-rarefaction-no-mitochondria-no-chloroplast.qzv onto local computer and visualized on view.QIIME2.org

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

# Create a new directory for parkinsons export and navigate to it

(qiime2-amplicon-2025.4) root@stu-1130:/data/team5_parkinsons_data# mkdir parkinsons_export
(qiime2-amplicon-2025.4) root@stu-1130:/data/team5_parkinsons_data# cd parkinsons_export

# Export parkinsons_table-no-mitochondria-no-chloroplast.qza
qiime tools export \
  --input-path ../parkinsons_table-no-mitochondria-no-chloroplast.qza \
  --output-path parkinsons_table-no-mitochondria-no-chloroplast_export

# Navigate to parkinsons_table-no-mitochondria-no-chloroplast_export and convert .biom file to .txt file
(qiime2-amplicon-2025.4) root@stu-1130:/data/team5_parkinsons_data/parkinsons_export# cd parkinsons_table-no-mitochondria-no-chloroplast_export/

biom convert -i feature-table.biom --to-tsv -o feature-table.txt

# Export parkinsons_rooted-tree.qza
(qiime2-amplicon-2025.4) root@stu-1130:/data/team5_parkinsons_data# cd parkinsons_export

qiime tools export \
  --input-path ../parkinsons_rooted-tree.qza \
  --output-path parkinsons_rooted-tree_export

# Export parkinsons_taxonomy.qza
(qiime2-amplicon-2025.4) root@stu-1130:/data/team5_parkinsons_data# cd parkinsons_export

qiime tools export \
  --input-path ../parkinsons_taxonomy.qza \
  --output-path parkinsons_taxonomy_export

# Transferred parkinsons_export to local computer
# Transferred parkinsons_metadata.txt to local computer

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

    # Taxonomic analysis 
qiime feature-classifier classify-sklearn \
  --i-classifier /datasets/classifiers/silva-138-99-515-806-nb-classifier.qza \
  --i-reads gc_rep-seqs.qza \
  --o-classification gc_taxonomy.qza

qiime metadata tabulate \
  --m-input-file gc_taxonomy.qza \
  --o-visualization gc_taxonomy.qzv

# Removing mitochondria and chloroplast sequences
qiime taxa filter-table \
  --i-table gc_table.qza \
  --i-taxonomy gc_taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table gc_table-no-mitochondria-no-chloroplast.qza

  # Removing samples below age 40 and over age 80
  qiime feature-table filter-samples \
  --i-table gc_table-no-mitochondria-no-chloroplast.qza \
  --m-metadata-file /datasets/project_2/gastric_cancer/gastric_cancer_metadata.tsv \
  --p-where "Age >= 40 AND Age <= 80" \
  --o-filtered-table gc_table_age40-80.qza
  
  qiime feature-table summarize \
  --i-table gc_table_age40-80.qza \
  --o-visualization gc_table_age40-80.qzv \
  --m-sample-metadata-file /datasets/project_2/gastric_cancer/gastric_cancer_metadata.tsv

  # Transferred gc_table_age40-80.qzv onto local computer and visualized on view.QIIME2.org

  # Generating a tree for phylogenetic diversity analyses
  qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences gc_rep-seqs.qza \
  --o-alignment gc_aligned-rep-seqs.qza \
  --o-masked-alignment gc_masked-aligned-rep-seqs.qza \
  --o-tree gc_unrooted-tree.qza \
  --o-rooted-tree gc_rooted-tree.qza

# Alpha-rarefaction
qiime diversity alpha-rarefaction \
  --i-table  parkinsons_table-no-mitochondria-no-chloroplast.qza \
  --i-phylogeny parkinsons_rooted-tree.qza \
  --p-max-depth 20000 \
  --m-metadata-file /datasets/project_2/parkinsons/parkinsons_metadata.txt \
  --o-visualization parkinsons_alpha-rarefaction-no-mitochondria-no-chloroplast.qzv

qiime diversity alpha-rarefaction \
  --i-table gc_table_age40-80.qza \
  --i-phylogeny gc_rooted-tree.qza \
  --p-max-depth 180000 \
  --m-metadata-file /datasets/project_2/gastric_cancer/gastric_cancer_metadata.tsv \
  --o-visualization gc_alpha-rarefaction_age40-80.qzv

  ## Confirm sampling depth then run below codes

  qiime diversity core-metrics-phylogenetic \
  --i-phylogeny gc_rooted-tree.qza \
  --i-table gc_table_age40-80.qza \
  --p-sampling-depth 26244 \
  --m-metadata-file /datasets/project_2/gastric_cancer/gastric_cancer_metadata.tsv \
  --output-dir gc_core-metrics-results

  qiime diversity alpha-group-significance \
  --i-alpha-diversity gc_core-metrics-results/shannon_vector.qza \
  --m-metadata-file /datasets/project_2/gastric_cancer/gastric_cancer_metadata.tsv \
  --o-visualization gc_shannon_group_significance.qzv

mkdir gc_export
cd gc_export

qiime tools export \
  --input-path ../gc_table_age40-80.qza \
  --output-path gc_table_age40-80_export

cd gc_table_age40-80_export
biom convert -i feature-table.biom --to-tsv -o feature-table.txt

cd ../
qiime tools export \
  --input-path ../gc_rooted-tree.qza \
  --output-path gc_rooted-tree_export

  qiime tools export \
  --input-path ../gc_taxonomy.qza \
  --output-path gc_taxonomy_export


