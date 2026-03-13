#!/usr/bin/env Rscript
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

#### Load data ####
PD_metadata <- read_delim(file="parkinsons_metadata.txt", delim = "\t")
PD_otu <- read_delim(file="PD_feature-table.txt", delim = "\t", skip=1)
PD_tax <- read_delim(file = "PD_taxonomy.tsv", delim="\t")
PD_tree <- read.tree("PD_tree.nwk")

GC_metadata <- read_delim(file = "gastric_cancer_metadata.txt", delim = "\t")
GC_otu <- read_delim(file = "gc_feature-table.txt", delim = "\t", skip = 1)
GC_tax <- read_delim(file = "gc_taxonomy.tsv", delim = "\t")
GC_tree <- read.tree("gc_tree.nwk")

#### Format OTU table ####
# Save everything but first column into a matrix
PD_otu_mat <- as.matrix(PD_otu[,-1])

GC_otu_mat <- as.matrix(GC_otu[,-1])

# Make #OTU ID the rownames of the new matrix
rownames(PD_otu_mat) <- PD_otu$`#OTU ID`

rownames(GC_otu_mat) <- GC_otu$`#OTU ID`

# Make an OTU table
PD_OTU <- otu_table(PD_otu_mat, taxa_are_rows = TRUE)

GC_OTU <- otu_table (GC_otu_mat, taxa_are_rows = TRUE)

#### Format metadata ####
# Save everything but sampleid as new data frame
PD_df <- as.data.frame(PD_metadata[,-1])

GC_df <- as.data.frame(GC_metadata[,-1])
# Make SampleID the rownames
rownames(PD_df) <- PD_metadata$`#SampleID`

rownames(GC_df) <- GC_metadata$`sample-id`
#Make phyloseq sample data
PD_SAMP <- sample_data(PD_df)

GC_SAMP <- sample_data(GC_df)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
PD_tax_mat <- PD_tax %>% select(-Confidence) %>%
  separate(col=Taxon, sep="; "
           , into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()

GC_tax_mat <- GC_tax %>% select(-Confidence) %>%
  separate(col=Taxon, sep=";"
           , into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()

# Save everything except feature IDs
PD_tax_mat <- PD_tax_mat[,-1]

GC_tax_mat <- GC_tax_mat[, -1]

# Make Feature ID the rownames
rownames(PD_tax_mat) <- PD_tax$`Feature ID`

rownames(GC_tax_mat) <- GC_tax$`Feature ID`

# Make taxa table
PD_TAX <- tax_table(PD_tax_mat)

GC_TAX <- tax_table(GC_tax_mat)

#### Create phyloseq object ####
PD <- phyloseq(PD_OTU, PD_SAMP, PD_TAX, PD_tree)

GC <- phyloseq(GC_OTU, GC_SAMP, GC_TAX, GC_tree)
