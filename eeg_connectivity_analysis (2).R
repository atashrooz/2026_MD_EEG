# ============================================================
# Resting-State EEG Connectivity Analysis (v3, 6-ROI scheme)
#
# Tests whether debiased weighted phase-lag index (wPLI)
# connectivity differs between groups (MD vs. nonMD), eye
# conditions (EO vs. EC), and frequency bands - and whether
# connectivity scales with MD symptom severity.
#
# Connectivity outcomes are split into two families, mirroring
# the primary/secondary split already used for the alpha-
# asymmetry analysis:
#
#   PRIMARY   (1 outcome):  global_wpli
#             whole-scalp connectivity - the confirmatory test
#
#   SECONDARY (21 outcomes): six within-ROI + fifteen between-ROI
#             wPLI measures, from the six-ROI scheme (frontal,
#             frontocentral, central, centroparietal, parietal,
#             occipital) described in the manuscript Methods.
#             Treated as exploratory: in addition to the
#             within-model correction (across terms, same as
#             the primary model), p-values for each fixed-effect
#             term are ALSO FDR-corrected across all 21 regional
#             outcomes, since this family involves testing the
#             same hypothesis at many spatial locations.
#
# Each outcome is modeled with a linear mixed model
# (group/severity x condition x band, random intercept per
# subject).
#
# Input:
#   - eeg_connectivity_features_by_band_v3_6ROI.csv
#       one row per subject x condition x band, with wPLI/
#       coherence connectivity measures (6-ROI scheme) and
#       demographic/severity columns
#
# Requires: tidyverse, lme4, lmerTest, broom.mixed
# ============================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)


# Sum-to-zero contrasts, so fixed-effect main effects are marginal
# (averaged across other factors) rather than simple effects at
# reference levels — matches v1/v2 and is required for the results
# to be comparable across versions.
options(contrasts = c("contr.sum", "contr.poly"))
# ------------------------------------------------------------
# File paths - edit these to point at your local copies
# ------------------------------------------------------------

connectivity_file <- "eeg_connectivity_features_by_band_v3_6ROI.csv"

group_primary_file    <- "connectivity_group_models_v3_primary.csv"
group_secondary_file  <- "connectivity_group_models_v3_secondary.csv"
mds_primary_file      <- "connectivity_mds_models_v3_primary.csv"
mds_secondary_file    <- "connectivity_mds_models_v3_secondary.csv"


# ============================================================
# 1. LOAD AND PREPARE DATA
# ============================================================

conn <- read_csv(connectivity_file, show_col_types = FALSE)

glimpse(conn)

# Cell sizes per group x condition x band
conn %>%
  count(md_status, condition, band)

# All wPLI outcome columns present in the 6-ROI file
wpli_columns <- names(conn)[str_detect(names(conn), "_wpli$")]
wpli_columns

# Missingness on each wPLI outcome
conn %>%
  summarise(across(all_of(wpli_columns), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "outcome", values_to = "n_missing")

conn <- conn %>%
  mutate(
    id        = factor(id),
    group     = factor(md_status, levels = c("nonMD", "MD")),
    condition = factor(condition, levels = c("EO", "EC")),
    band      = factor(band, levels = c("delta", "theta", "alpha", "beta", "gamma")),
    gender    = factor(gender),
    age_z     = as.numeric(scale(age)),
    mds_value = as.numeric(str_replace(mds_total, "/", ".")),
    mds_z     = as.numeric(scale(mds_value))
  )

# Restrict to the three a priori bands of interest
conn_main <- conn %>%
  filter(band %in% c("theta", "alpha", "beta"))

# ------------------------------------------------------------
# Outcome families
# ------------------------------------------------------------

primary_outcome <- "global_wpli"

secondary_outcomes <- setdiff(wpli_columns, primary_outcome)

length(secondary_outcomes)  # should be 21: 6 within-ROI + 15 between-ROI
secondary_outcomes


# ============================================================
# 2. GROUP MODELS (MD vs. nonMD)
#
# For each outcome, fits:
#   outcome ~ group * condition * band + age_z + gender + (1 | id)
# ============================================================

fit_group_model <- function(outcome_name, data) {
  model_formula <- as.formula(
    paste0(outcome_name, " ~ group * condition * band + age_z + gender + (1 | id)")
  )

  model <- lmer(model_formula, data = data, REML = FALSE)

  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(outcome = outcome_name)
}

# --- Primary: global connectivity, one confirmatory test -----

group_primary_results <- fit_group_model(primary_outcome, conn_main) %>%
  group_by(outcome) %>%
  mutate(p_fdr_within_model = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

write_csv(group_primary_results, group_primary_file)
group_primary_results

# --- Secondary: 21 regional wPLI outcomes ---------------------

group_secondary_results <- map_dfr(secondary_outcomes, fit_group_model, data = conn_main)

group_secondary_results_fdr <- group_secondary_results %>%
  # (a) within-model correction: across fixed-effect terms,
  #     same logic used for the primary model and for v1/v2
  group_by(outcome) %>%
  mutate(p_fdr_within_model = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  # (b) cross-outcome correction: for each fixed-effect term,
  #     across all 21 regional outcomes - this is the stricter,
  #     spatially-aware correction to lead with when reporting
  #     the secondary/exploratory family
  group_by(term) %>%
  mutate(p_fdr_across_regions = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

write_csv(group_secondary_results_fdr, group_secondary_file)
group_secondary_results_fdr


# ============================================================
# 4. QUICK SUMMARY: which secondary regions survive correction?
#
# Focuses on the group1 / mds_z terms specifically, since those
# are the terms of substantive interest.
# ============================================================

group_secondary_results_fdr %>%
  filter(term %in% c("group1")) %>%
  arrange(p_fdr_across_regions) %>%
  select(outcome, term, estimate, p.value, p_fdr_within_model, p_fdr_across_regions)

