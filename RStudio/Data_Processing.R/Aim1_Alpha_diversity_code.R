#!/usr/bin/env Rscript
library(phyloseq)
library(ggplot2)
library(ggsignif)

#### Load in RData ####
load("GC_rare.RData")
load("PD_final.RData")

#### Shannon diversity ######
# Keep only GC and healthy control
GC_two <- subset_samples(
  GC_rare,
  Group %in% c("Healthy control (HC)", "Gastric cancer (GC)")
)

sample_data(GC_two)$Group <- factor(
  sample_data(GC_two)$Group,
  levels = c("Healthy control (HC)", "Gastric cancer (GC)"),
  labels = c("Control", "GC")
)

GC_gg_richness_two <- plot_richness(GC_two, x = "Group", measures = "Shannon") +
  xlab("Group") +
  geom_boxplot()

ggsave(filename = "GC_plot_richness_two.png"
       , GC_gg_richness_two
       , height=4, width=6)

PD_gg_richness <- plot_richness(PD_final, x = "Disease", measures = c("Shannon")) +
  xlab("Group") +
  geom_boxplot()

ggsave(filename = "PD_plot_richness.png"
       , PD_gg_richness
       , height=4, width=6)

GC_alphadiv <- estimate_richness(GC_two)
GC_samp_dat <- sample_data(GC_two)
GC_samp_dat_wdiv <- data.frame(GC_samp_dat, GC_alphadiv)

PD_alphadiv <- estimate_richness(PD_final)
PD_samp_dat <- sample_data(PD_final)
PD_samp_dat_wdiv <- data.frame(PD_samp_dat, PD_alphadiv)

#### Significance for GC ####
GC_wilcox_shannon <- wilcox.test(Shannon ~ Group, data = GC_samp_dat_wdiv)
GC_wilcox_shannon

GC_plot_richness_sig <- plot_richness(GC_two, x = "Group", measures = c("Shannon")) +
  xlab("Group") +
  geom_boxplot() + 
  geom_signif(
    comparisons = list(
      c("Control", "GC")
    ),
    y_position = c(6),
    annotations = c("****")
  )+
  theme_bw()

ggsave(filename = "GC_plot_richness_two_sig.png"
       , GC_plot_richness_sig
       , height=5, width=6)

#### Significance for PD ####
wilcox.test(Shannon ~ Disease, data = PD_samp_dat_wdiv)

PD_plot_richness_sig <- plot_richness(PD_final, x = "Disease", measures = "Shannon") +
  xlab("Group") +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("Control", "PD")),
    y_position = 5,
    annotations = "n.s."
    ) +
  theme_bw()

ggsave(filename = "PD_plot_richness_sig.png"
       , PD_plot_richness_sig
       , height=5, width=6)
