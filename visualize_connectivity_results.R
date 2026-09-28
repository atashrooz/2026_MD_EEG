# ============================================================
# Connectivity Results Visualization
#
# Two figures, matching the visual style already used for the
# spectral-power figures (Okabe-Ito colors, theme_classic base):
#
#   Figure 1: global_wpli by group (raw data, box + jitter)
#   Figure 2: 6x6 ROI heatmap of MD-nonMD wPLI differences,
#             with FDR-significant regions marked
#
# Inputs (all already produced by earlier scripts):
#   - eeg_connectivity_features_by_band_v3_6ROI.csv
#   - connectivity_group_emmeans_v3_secondary.csv
#
# Requires: tidyverse, stringr
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# File paths - edit these to point at your local copies
# ------------------------------------------------------------

connectivity_file       <- "eeg_connectivity_features_by_band_v3_6ROI.csv"
group_emmeans_secondary <- "connectivity_group_emmeans_v3_secondary.csv"

plot_directory <- "plots"
dir.create(plot_directory, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Shared style (matches the spectral-power figures)
# ------------------------------------------------------------

group_colors <- c(nonMD = "#0072B2", MD = "#D55E00")
group_labels <- c(nonMD = "non MDer", MD = "Probable MDer")

fig_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none"
  )


# ============================================================
# FIGURE 1: global_wpli by group
# ============================================================

conn <- read_csv(connectivity_file, show_col_types = FALSE) %>%
  mutate(group = factor(md_status, levels = c("nonMD", "MD"))) %>%
  filter(band %in% c("theta", "alpha", "beta"))

# One value per subject: global_wpli averaged across the three
# a priori bands and both conditions, matching what the marginal
# group contrast (Section 5 of the analysis script) tested.
global_by_subject <- conn %>%
  group_by(id, group) %>%
  summarise(global_wpli = mean(global_wpli, na.rm = TRUE), .groups = "drop")

fig1_global_wpli <- ggplot(
  global_by_subject,
  aes(x = group, y = global_wpli, fill = group)
) +
  geom_boxplot(width = 0.5, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.7, aes(color = group)) +
  scale_x_discrete(labels = group_labels) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Global wPLI Connectivity by Group",
    x = NULL,
    y = "Global wPLI (subject mean, theta/alpha/beta)"
  ) +
  fig_theme

print(fig1_global_wpli)

ggsave(
  filename = file.path(plot_directory, "fig_global_wpli_by_group.tiff"),
  plot = fig1_global_wpli,
  width = 5, height = 5, units = "in", dpi = 300, bg = "white"
)


# ============================================================
# FIGURE 2: 6x6 ROI heatmap of MD - nonMD wPLI differences
# ============================================================

roi_order <- c("frontal", "frontocentral", "central",
               "centroparietal", "parietal", "occipital")

roi_labels <- c(
  frontal        = "Frontal",
  frontocentral  = "Fronto-\ncentral",
  central        = "Central",
  centroparietal = "Centro-\nparietal",
  parietal       = "Parietal",
  occipital      = "Occipital"
)

emmeans_secondary <- read_csv(group_emmeans_secondary, show_col_types = FALSE)

# Parse outcome names ("within_frontal_wpli", "frontal_central_wpli", ...)
# into roi1 / roi2 pairs. Within-ROI outcomes get roi1 == roi2.
roi_pairs <- emmeans_secondary %>%
  mutate(outcome_stub = str_remove(outcome, "_wpli$")) %>%
  mutate(
    is_within = str_starts(outcome_stub, "within_"),
    roi1 = if_else(
      is_within,
      str_remove(outcome_stub, "^within_"),
      str_extract(outcome_stub, paste0("^(", paste(roi_order, collapse = "|"), ")"))
    )
  ) %>%
  mutate(
    roi2 = if_else(
      is_within,
      roi1,
      str_remove(outcome_stub, paste0("^", roi1, "_"))
    )
  )

# Sanity check: every roi1/roi2 should be one of the six known ROI names
stopifnot(all(roi_pairs$roi1 %in% roi_order))
stopifnot(all(roi_pairs$roi2 %in% roi_order))

# Build a full symmetric matrix (mirror off-diagonal entries) so the
# heatmap reads the same regardless of which ROI is listed first.
roi_matrix_data <- bind_rows(
  roi_pairs %>% select(roi1, roi2, estimate, p_fdr_across_regions),
  roi_pairs %>% filter(roi1 != roi2) %>%
    transmute(roi1 = roi2, roi2 = roi1, estimate, p_fdr_across_regions)
) %>%
  mutate(
    roi1 = factor(roi1, levels = roi_order),
    roi2 = factor(roi2, levels = rev(roi_order)),
    significant = p_fdr_across_regions < 0.05,
    sig_label = if_else(significant, "*", "")
  )

fig2_roi_heatmap <- ggplot(
  roi_matrix_data,
  aes(x = roi1, y = roi2, fill = estimate)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_label), size = 6, vjust = 0.75) +
  scale_x_discrete(labels = roi_labels, position = "top") +
  scale_y_discrete(labels = roi_labels) +
  scale_fill_gradient2(
    name = "MD \u2212 nonMD\n(wPLI)",
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0
  ) +
  coord_fixed() +
  labs(
    title = "Regional wPLI Connectivity: MD \u2212 nonMD",
    subtitle = "* FDR-corrected across all 21 regions, p < .05",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(color = "black", size = 9),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9)
  )

print(fig2_roi_heatmap)

ggsave(
  filename = file.path(plot_directory, "fig_roi_connectivity_heatmap.tiff"),
  plot = fig2_roi_heatmap,
  width = 6.5, height = 5.5, units = "in", dpi = 300, bg = "white"
)
