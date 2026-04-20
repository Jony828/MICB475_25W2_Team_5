#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### load data ####
load("GC_PD_taxonomy.RData")
load("GC_PD_phyloseq(no_tree).RData")
load("GC_PD_Otu.RData")

#Attain relative abundance
GCPD_RA <- transform_sample_counts(GC_PD_phy, fun=function(x) x/sum(x))

#Subset data into treatment and control groups
GC_stat <- subset_samples(GCPD_RA, `Disease`=="GC")
PD_stat <- subset_samples(GCPD_RA, `Disease`=="PD")
control_statGC <- subset_samples(GCPD_RA, `Disease`=="HC")
control_statPD <- subset_samples(GCPD_RA, `Disease`=="Control")


#Set prevalence and detection threshold. Comparing ASVs across GC, PD, and Healthy Control 
GC_ASVs <- core_members(GC_stat, detection=0.001, prevalence = 0.1)
PD_ASVs <- core_members(PD_stat, detection=0.001, prevalence = 0.1)

control_ASVsGC <- core_members(control_statGC, detection=0.001, prevalence = 0.1)
control_ASVsPD <- core_members(control_statPD, detection=0.001, prevalence = 0.1)

ASV_listGC <- list(GC = GC_ASVs, control = control_ASVsGC)
ASV_listPD <- list(PD = PD_ASVs, control = control_ASVsPD)
ASV_listxyz <- list(PD = PD_ASVs, GC = GC_ASVs)


#### Save the Venn Diagrem as a png ####
max_count <- 220

#GC vs Healthy Control
GCPD_venndiagramGC <- ggVennDiagram(ASV_listGC) +
  scale_fill_gradient(limits = c(0, 220))
ggsave("GCPDVennDiagramGC.png", GCPD_venndiagramGC)

#PD vs Healthy Control
GCPD_venndiagramPD <- ggVennDiagram(ASV_listPD) +
  scale_fill_gradient(limits = c(0, 220))
ggsave("GCPDVennDiagramPD.png", GCPD_venndiagramPD)

#GC vs PD
GCPD_venndiagramxyz <- ggVennDiagram(ASV_listxyz) +
  scale_fill_gradient(limits = c(0, 220))
ggsave("GCPDVennDiagramxyz.png", GCPD_venndiagramxyz)

