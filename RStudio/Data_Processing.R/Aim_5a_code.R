#!/usr/bin/env Rscript

# MICB 475 Team 5 - Aim 5A: DESeq2 Differential Abundance Analysis
# Script: Aim_5a_code.R
# GitHub: MICB475_25W2_Team_5

# SECTION 1: Load Required Libraries

library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(patchwork)

# Install and load "ggrepel"
install.packages("ggrepel") 
library(ggrepel)


# SECTION 2: Load Phyloseq Objects

load("RStudio/Phyloseq objects/PD_final.RData")
load("RStudio/Phyloseq objects/GC_rare.RData")


# SECTION 3: Subset, Pseudocount, Convert to DESeq2

GC_subset <- subset_samples(GC_rare, Group %in% c("Gastric cancer (GC)", "Healthy control (HC)"))

PD_plus1 <- transform_sample_counts(PD_final, function(x) x + 1)
GC_plus1 <- transform_sample_counts(GC_subset, function(x) x + 1)

PD_deseq_obj <- phyloseq_to_deseq2(PD_plus1, ~ Disease)
GC_deseq_obj <- phyloseq_to_deseq2(GC_plus1, ~ Group)


# SECTION 4: Run DESeq2

PD_DESeq <- DESeq(PD_deseq_obj)
GC_DESeq <- DESeq(GC_deseq_obj)


# SECTION 5: CHECK GC CONTRAST DIRECTION

# IMPORTANT: DESeq2 uses alphabetical order to set the reference level by
# default if no contrast is specified. We explicitly set:
#   numerator   = "Gastric cancer (GC)"   -> positive LFC = enriched in GC
#   denominator = "Healthy control (HC)"  -> the baseline/reference
# This means:
#   POSITIVE log2FoldChange = more abundant in GC patients (e.g. Helicobacter)
#   NEGATIVE log2FoldChange = less abundant in GC patients (depleted)
# We can verify this is correct by checking where Helicobacter lands.

GC_res <- results(GC_DESeq,
                  tidy = TRUE,
                  contrast = c("Group",
                               "Gastric cancer (GC)",    # numerator (disease)
                               "Healthy control (HC)"))  # denominator (reference)

PD_res <- results(PD_DESeq,
                  tidy = TRUE,
                  contrast = c("Disease",
                               "PD",       # numerator (disease)
                               "Control")) # denominator (reference)

# GC DIRECTION SANITY CHECK
# Find any ASVs that are classified as Helicobacter in the full GC results
# (before filtering) and check their log2FoldChange sign.
# If the contrast is correct, Helicobacter should have POSITIVE LFC (enriched in GC).

# Step 1: Get full taxonomy from GC phyloseq
GC_full_tax <- tax_table(GC_subset) %>%
  as.data.frame() %>%
  rownames_to_column(var = "row")

# Step 2: Join with the full (unfiltered) DESeq2 results
GC_res_with_tax <- left_join(GC_res, GC_full_tax, by = "row")

# Step 3: Find Helicobacter rows
helicobacter_check <- GC_res_with_tax %>%
  filter(grepl("Helicobacter", Genus, ignore.case = TRUE)) %>%
  select(row, Genus, log2FoldChange, padj)

cat("\n=== HELICOBACTER DIRECTION CHECK ===\n")
print(helicobacter_check)
cat("\nIf log2FoldChange is POSITIVE: contrast is correct (Helicobacter enriched in GC)\n")
cat("If log2FoldChange is NEGATIVE: contrast is FLIPPED - swap numerator/denominator\n\n")

# Helicobacter shows NEGATIVE LFC, so replace the GC_res line above with:

GC_res <- results(GC_DESeq,
                  tidy = TRUE,
                  contrast = c("Group",
                               "Healthy control (HC)",   # swapped
                               "Gastric cancer (GC)"))   # swapped


# SECTION 6: HIGH THRESHOLD BAR PLOTS

# High threshold filtering
PD_sigASVs_high <- PD_res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

GC_sigASVs_high <- GC_res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 6) %>%
  dplyr::rename(ASV = row)

# Attach taxonomy - PD high threshold
PD_pruned_high <- prune_taxa(PD_sigASVs_high %>% pull(ASV), PD_final)
PD_merged_high <- tax_table(PD_pruned_high) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(PD_sigASVs_high, by = "ASV") %>%
  filter(!is.na(Genus)) %>%
  filter(!grepl("Incertae_Sedis", Genus)) %>%
  filter(!grepl("^NA", Genus)) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# Attach taxonomy - GC high threshold
GC_pruned_high <- prune_taxa(GC_sigASVs_high %>% pull(ASV), GC_subset)
GC_merged_high <- tax_table(GC_pruned_high) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(GC_sigASVs_high, by = "ASV") %>%
  filter(!is.na(Genus)) %>%
  filter(!grepl("Incertae_Sedis", Genus)) %>%
  filter(!grepl("^NA", Genus)) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# PD High Threshold Bar Plot
PD_deseq_plot_high <- ggplot(PD_merged_high) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = log2FoldChange > 0),
           stat = "identity") +
  geom_errorbar(aes(x = Genus,
                    ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#E07B7B", "FALSE" = "#7BB2E0"),
                    labels = c("TRUE" = "Enriched in PD", "FALSE" = "Depleted in PD"),
                    name = "") +
  labs(title = "Differentially Abundant Taxa: PD vs. Control (High Threshold)",
       x = "Genus", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title  = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

PD_deseq_plot_high
ggsave("Figures/PD_DESeq2_barplot_high_threshold.png",
       plot = PD_deseq_plot_high, width = 10, height = 6, dpi = 300)

# GC High Threshold Bar Plot
GC_deseq_plot_high <- ggplot(GC_merged_high) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = log2FoldChange > 0),
           stat = "identity") +
  geom_errorbar(aes(x = Genus,
                    ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#E07B7B", "FALSE" = "#7BB2E0"),
                    labels = c("TRUE" = "Enriched in GC", "FALSE" = "Depleted in GC"),
                    name = "") +
  labs(title = "Differentially Abundant Taxa: GC vs. Healthy Control (High Threshold)",
       x = "Genus", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title   = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

GC_deseq_plot_high
ggsave("Figures/GC_DESeq2_barplot_high_threshold.png",
       plot = GC_deseq_plot_high, width = 10, height = 6, dpi = 300)


# SECTION 7: VOLCANO PLOTS

# A volcano plot shows ALL ASVs (not just significant ones), plotting:
#   X axis = log2FoldChange (effect size - how different between groups)
#   Y axis = -log10(padj)  (significance - higher = more significant)
# Points in the upper-right  = enriched in disease AND significant
# Points in the upper-left   = depleted in disease AND significant
# Points near the bottom     = not significant (regardless of fold change)
# We colour points by significance category and label the top hits by genus.
# This gives a much better overview of the full results than a bar plot alone.

# Helper: attach taxonomy to full (unfiltered) results for labelling

# PD: join full results with full taxonomy
PD_full_tax <- tax_table(PD_final) %>%
  as.data.frame() %>%
  rownames_to_column(var = "row")

PD_volcano_data <- left_join(PD_res, PD_full_tax, by = "row") %>%
  filter(!is.na(padj)) %>%                        # remove ASVs with no p-value
  mutate(
    neg_log10_padj = -log10(padj),                # transform p-value for y axis
    significance = case_when(
      padj < 0.01 & log2FoldChange > 2  ~ "Enriched in PD",   # sig & enriched
      padj < 0.01 & log2FoldChange < -2 ~ "Depleted in PD",   # sig & depleted
      TRUE                               ~ "Not significant"   # everything else
    ),
    # Genus only - no family fallback (consistent with bar chart)
    label = case_when(
      padj < 0.01 & abs(log2FoldChange) > 2 &
        !is.na(Genus) &
        !grepl("Incertae_Sedis", Genus) &
        !grepl("^NA", Genus) ~ Genus,
      TRUE ~ NA_character_
    )
  )

# GC: join full results with full taxonomy
GC_volcano_data <- left_join(GC_res, GC_full_tax, by = "row") %>%
  filter(!is.na(padj)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = case_when(
      padj < 0.01 & log2FoldChange > 2  ~ "Enriched in GC",
      padj < 0.01 & log2FoldChange < -2 ~ "Depleted in GC",
      TRUE                               ~ "Not significant"
    ),
    # Genus only - no family fallback (consistent with bar chart)
    label = case_when(
      padj < 0.01 & abs(log2FoldChange) > 6 &
        !is.na(Genus) &
        !grepl("Incertae_Sedis", Genus) &
        !grepl("^NA", Genus) ~ Genus,
      TRUE ~ NA_character_
    )
  )

# PD Volcano Plot
PD_volcano_plot <- ggplot(PD_volcano_data,
                          aes(x = log2FoldChange,
                              y = neg_log10_padj,
                              colour = significance)) +
  # plot all non-significant points first (grey, small)
  geom_point(data = filter(PD_volcano_data, significance == "Not significant"),
             alpha = 0.4, size = 1.2) +
  # plot significant points on top (coloured, larger)
  geom_point(data = filter(PD_volcano_data, significance != "Not significant"),
             size = 2, alpha = 0.8) +
  # add genus labels for top hits only (ggrepel avoids overlap)
  geom_label_repel(aes(label = label),
                   na.rm = TRUE,
                   size = 2.8,
                   max.overlaps = 20,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   show.legend = FALSE) +
  # threshold lines
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", colour = "grey40") +
  # colours: enriched = red, depleted = blue, ns = grey
  scale_colour_manual(values = c("Enriched in PD"  = "#E07B7B",
                                 "Depleted in PD"  = "#7BB2E0",
                                 "Not significant" = "grey70"),
                      name = "") +
  labs(title = "Volcano Plot: PD vs. Control",
       subtitle = "Dashed lines: padj = 0.01, |log2FC| = 2",
       x = "Log2 Fold Change (PD / Control)",
       y = "-log10(adjusted p-value)") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title   = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

PD_volcano_plot
ggsave("Figures/PD_DESeq2_volcano.png",
       plot = PD_volcano_plot, width = 9, height = 7, dpi = 300)

# GC Volcano Plot
GC_volcano_plot <- ggplot(GC_volcano_data,
                          aes(x = log2FoldChange,
                              y = neg_log10_padj,
                              colour = significance)) +
  geom_point(data = filter(GC_volcano_data, significance == "Not significant"),
             alpha = 0.4, size = 1.2) +
  geom_point(data = filter(GC_volcano_data, significance != "Not significant"),
             size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = label),
                   na.rm = TRUE,
                   size = 2.8,
                   max.overlaps = 20,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Enriched in GC"  = "#E07B7B",
                                 "Depleted in GC"  = "#7BB2E0",
                                 "Not significant" = "grey70"),
                      name = "") +
  labs(title = "Volcano Plot: GC vs. Healthy Control",
       subtitle = "Dashed lines: padj = 0.01, |log2FC| = 2  |  Labels shown for |log2FC| > 6",
       x = "Log2 Fold Change (GC / Healthy Control)",
       y = "-log10(adjusted p-value)") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title    = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title.x  = element_text(size = 13),
        axis.title.y  = element_text(size = 13))

GC_volcano_plot
ggsave("Figures/GC_DESeq2_volcano.png",
       plot = GC_volcano_plot, width = 9, height = 7, dpi = 300)


# SECTION 8: LOW THRESHOLD RESULTS + BAR PLOTS (for overlap analysis)

# Re-filter with relaxed thresholds to find cross-disease overlap

PD_sigASVs <- PD_res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
  dplyr::rename(ASV = row)

GC_sigASVs <- GC_res %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)

nrow(PD_sigASVs)
nrow(GC_sigASVs)

# Attach taxonomy - PD
PD_pruned <- prune_taxa(PD_sigASVs %>% pull(ASV), PD_final)
PD_merged_results <- tax_table(PD_pruned) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(PD_sigASVs, by = "ASV") %>%
  filter(!is.na(Genus)) %>%
  filter(!grepl("Incertae_Sedis", Genus)) %>%
  filter(!grepl("^NA", Genus)) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# Attach taxonomy - GC
GC_pruned <- prune_taxa(GC_sigASVs %>% pull(ASV), GC_subset)
GC_merged_results <- tax_table(GC_pruned) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(GC_sigASVs, by = "ASV") %>%
  filter(!is.na(Genus)) %>%
  filter(!grepl("Incertae_Sedis", Genus)) %>%
  filter(!grepl("^NA", Genus)) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# Save CSVs
write.csv(PD_merged_results, "RStudio/DESeq2_objects/PD_DESeq2_results.csv", row.names = FALSE)
write.csv(GC_merged_results, "RStudio/DESeq2_objects/GC_DESeq2_results.csv", row.names = FALSE)

# PD Bar Plot (low threshold)
PD_deseq_plot <- ggplot(PD_merged_results) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = log2FoldChange > 0),
           stat = "identity") +
  geom_errorbar(aes(x = Genus,
                    ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#E07B7B", "FALSE" = "#7BB2E0"),
                    labels = c("TRUE" = "Enriched in PD", "FALSE" = "Depleted in PD"),
                    name = "") +
  labs(title = "Differentially Abundant Taxa: PD vs. Control",
       x = "Genus", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title   = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

PD_deseq_plot
ggsave("Figures/PD_DESeq2_barplot.png",
       plot = PD_deseq_plot, width = 10, height = 6, dpi = 300)

# GC Bar Plot (low threshold)
GC_deseq_plot <- ggplot(GC_merged_results) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = log2FoldChange > 0),
           stat = "identity") +
  geom_errorbar(aes(x = Genus,
                    ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#E07B7B", "FALSE" = "#7BB2E0"),
                    labels = c("TRUE" = "Enriched in GC", "FALSE" = "Depleted in GC"),
                    name = "") +
  labs(title = "Differentially Abundant Taxa: GC vs. Healthy Control",
       x = "Genus", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title   = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

GC_deseq_plot
ggsave("Figures/GC_DESeq2_barplot.png",
       plot = GC_deseq_plot, width = 10, height = 6, dpi = 300)


# SECTION 9: OVERLAP ANALYSIS

PD_sig_genera <- PD_merged_results %>% pull(Genus) %>% as.character()
GC_sig_genera <- GC_merged_results %>% pull(Genus) %>% as.character()

shared_genera <- intersect(PD_sig_genera, GC_sig_genera)
print(shared_genera)

PD_unique_genera <- setdiff(PD_sig_genera, GC_sig_genera)
GC_unique_genera <- setdiff(GC_sig_genera, PD_sig_genera)

cat("Genera significant in BOTH PD and GC:", length(shared_genera), "\n")
cat("Genera significant in PD only:", length(PD_unique_genera), "\n")
cat("Genera significant in GC only:", length(GC_unique_genera), "\n")

if (length(shared_genera) > 0) {
  PD_shared_direction <- PD_merged_results %>%
    filter(as.character(Genus) %in% shared_genera) %>%
    select(Genus, log2FoldChange) %>%
    dplyr::rename(PD_log2FC = log2FoldChange)
  
  GC_shared_direction <- GC_merged_results %>%
    filter(as.character(Genus) %in% shared_genera) %>%
    select(Genus, log2FoldChange) %>%
    dplyr::rename(GC_log2FC = log2FoldChange)
  
  shared_comparison <- left_join(PD_shared_direction, GC_shared_direction, by = "Genus") %>%
    mutate(
      PD_direction = ifelse(PD_log2FC > 0, "Enriched in PD", "Depleted in PD"),
      GC_direction = ifelse(GC_log2FC > 0, "Enriched in GC", "Depleted in GC"),
      same_direction = (PD_log2FC > 0) == (GC_log2FC > 0)
    )
  
  View(shared_comparison)
  write.csv(shared_comparison, "RStudio/DESeq2_objects/Shared_taxa_PD_GC.csv", row.names = FALSE)
  
# Shared overlap bar plot
  shared_plot_data <- bind_rows(
    PD_merged_results %>%
      filter(as.character(Genus) %in% shared_genera) %>%
      mutate(Disease = "PD"),
    GC_merged_results %>%
      filter(as.character(Genus) %in% shared_genera) %>%
      mutate(Disease = "GC")
  )
  
  shared_overlap_plot <- ggplot(shared_plot_data,
                                aes(x = as.character(Genus),
                                    y = log2FoldChange,
                                    fill = Disease)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = log2FoldChange - lfcSE,
                      ymax = log2FoldChange + lfcSE),
                  position = position_dodge(0.9), width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("PD" = "#9B59B6", "GC" = "#E67E22")) +
    labs(title = "Shared Differentially Abundant Taxa: PD and GC vs. Controls",
         x = "Genus", y = "Log2 Fold Change") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title   = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13))
  
  shared_overlap_plot
  ggsave("Figures/Shared_taxa_PD_GC_plot.png",
         plot = shared_overlap_plot, width = 10, height = 6, dpi = 300)
}

#SECTION 10: COMBINED 4-PANEL FIGURE (Figure 3)
#Panels A & B = Volcano plots (PD and GC)
#Panels C & D = High threshold bar plots (PD and GC)

# Rebuild all 4 plots without titles (the panel labels will serve as identifiers)
# and with consistent theming for the combined figure

#Panel A: PD Volcano
panel_A <- ggplot(PD_volcano_data,
                  aes(x = log2FoldChange,
                      y = neg_log10_padj,
                      colour = significance)) +
  geom_point(data = filter(PD_volcano_data, significance == "Not significant"),
             alpha = 0.4, size = 1.2) +
  geom_point(data = filter(PD_volcano_data, significance != "Not significant"),
             size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = label),
                   na.rm = TRUE,
                   size = 2.8,
                   max.overlaps = 20,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Enriched in PD"  = "#E07B7B",
                                 "Depleted in PD"  = "#7BB2E0",
                                 "Not significant" = "grey70"),
                      name = "") +
  labs(title = "PD vs. Control",
       x = "Log2 Fold Change (PD / Control)",
       y = "-log10(adjusted p-value)") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title   = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

#Panel B: GC Volcano
panel_B <- ggplot(GC_volcano_data,
                  aes(x = log2FoldChange,
                      y = neg_log10_padj,
                      colour = significance)) +
  geom_point(data = filter(GC_volcano_data, significance == "Not significant"),
             alpha = 0.4, size = 1.2) +
  geom_point(data = filter(GC_volcano_data, significance != "Not significant"),
             size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = label),
                   na.rm = TRUE,
                   size = 2.8,
                   max.overlaps = 20,
                   box.padding = 0.4,
                   label.padding = 0.2,
                   show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Enriched in GC"  = "#E07B7B",
                                 "Depleted in GC"  = "#7BB2E0",
                                 "Not significant" = "grey70"),
                      name = "") +
  labs(title = "GC vs. Healthy Control",
       x = "Log2 Fold Change (GC / Healthy Control)",
       y = "-log10(adjusted p-value)") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title   = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

#Panel C: PD High Threshold Bar Plot
panel_C <- ggplot(PD_merged_high) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = log2FoldChange > 0),
           stat = "identity") +
  geom_errorbar(aes(x = Genus,
                    ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#E07B7B", "FALSE" = "#7BB2E0"),
                    labels = c("TRUE" = "Enriched in PD", "FALSE" = "Depleted in PD"),
                    name = "") +
  labs(title = "PD vs. Control",
       x = "Genus", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title   = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

#Panel D: GC High Threshold Bar Plot
panel_D <- ggplot(GC_merged_high) +
  geom_bar(aes(x = Genus, y = log2FoldChange, fill = log2FoldChange > 0),
           stat = "identity") +
  geom_errorbar(aes(x = Genus,
                    ymin = log2FoldChange - lfcSE,
                    ymax = log2FoldChange + lfcSE),
                width = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#E07B7B", "FALSE" = "#7BB2E0"),
                    labels = c("TRUE" = "Enriched in GC", "FALSE" = "Depleted in GC"),
                    name = "") +
  labs(title = "GC vs. Healthy Control",
       x = "Genus", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title   = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13))

#Combine into 2x2 grid with panel labels A, B, C, D
#Layout: A (PD volcano) | B (GC volcano)
#        C (PD bar)     | D (GC bar)
combined_figure <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_annotation(
    title = "",
    tag_levels = "A",          # adds A, B, C, D labels automatically
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

combined_figure

ggsave("Figures/Aim5a_Figure_combined_DESeq2.png",
       plot = combined_figure,
       width = 16, height = 14, dpi = 300)  

