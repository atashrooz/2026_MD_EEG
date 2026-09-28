# ============================================================
# RESTING-STATE ALPHA ASYMMETRY ANALYSIS
# ============================================================

rm(list = ls())

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

# Needed for interpretable Type III tests
options(contrasts = c("contr.sum", "contr.poly"))

# ------------------------------------------------------------
# 1. Read data
# ------------------------------------------------------------

eeg <- read_csv(
  "eeg_band_power_features.csv",
  show_col_types = FALSE
) %>%
  mutate(
    subject_id = as.integer(subject_id),
    group = factor(group, levels = c("nonMD", "MD")),
    condition = factor(condition, levels = c("EO", "EC"))
  )

phase2 <- read_excel("phase2_data_mod2 - Comp.xlsx")

# Alpha power must be positive before log transformation
if (any(eeg$alpha <= 0, na.rm = TRUE)) {
  stop("Some alpha-power values are zero or negative.")
}

# ------------------------------------------------------------
# 2. Prepare participant metadata
# ------------------------------------------------------------

metadata <- phase2 %>%
  transmute(
    id = as.character(id),
    age = readr::parse_number(as.character(age)),
    gender = factor(gender),
    mds_total = readr::parse_number(as.character(mds_total))
  )

eeg_subjects <- eeg %>%
  distinct(subject_id, group) %>%
  mutate(id = paste0("Sub", subject_id))

metadata <- eeg_subjects %>%
  left_join(metadata, by = "id") %>%
  mutate(
    group = factor(group, levels = c("nonMD", "MD"))
  ) %>%
  select(subject_id, group, age, gender, mds_total)

# Check missing metadata
metadata %>%
  summarise(
    missing_age = sum(is.na(age)),
    missing_gender = sum(is.na(gender)),
    missing_mds = sum(is.na(mds_total))
  )

# ------------------------------------------------------------
# 3. Define homologous electrode pairs
# ------------------------------------------------------------

alpha_pairs <- tribble(
  ~pair,       ~left,  ~right, ~analysis_family,
  "F4-F3",     "F3",   "F4",   "primary",
  "F8-F7",     "F7",   "F8",   "primary",
  "FC2-FC1",   "FC1",  "FC2",  "secondary",
  "FC4-FC3",   "FC3",  "FC4",  "secondary",
  "P4-P3",     "P3",   "P4",   "secondary",
  "PO4-PO3",   "PO3",  "PO4",  "secondary",
  "O2-O1",     "O1",   "O2",   "secondary"
)

# Check whether required channels exist
required_channels <- unique(c(alpha_pairs$left, alpha_pairs$right))

missing_channels <- setdiff(
  required_channels,
  unique(eeg$channel)
)

missing_channels

if (length(missing_channels) > 0) {
  stop(
    paste(
      "Missing channels:",
      paste(missing_channels, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 4. Calculate alpha asymmetry
#
# Asymmetry = ln(right alpha) - ln(left alpha)
#
# Positive value:
# greater right alpha power
# conventionally interpreted as relatively greater left activity
# ------------------------------------------------------------

alpha_wide <- eeg %>%
  filter(channel %in% required_channels) %>%
  mutate(log_alpha = log(alpha)) %>%
  select(
    subject_id,
    group,
    condition,
    channel,
    log_alpha
  ) %>%
  pivot_wider(
    names_from = channel,
    values_from = log_alpha
  )

asymmetry_data <- alpha_pairs %>%
  pmap_dfr(
    function(pair, left, right, analysis_family) {
      
      alpha_wide %>%
        transmute(
          subject_id,
          group,
          condition,
          pair = pair,
          analysis_family = analysis_family,
          left_log_alpha = .data[[left]],
          right_log_alpha = .data[[right]],
          asymmetry = right_log_alpha - left_log_alpha
        )
    }
  ) %>%
  left_join(
    metadata %>%
      select(subject_id, age, gender, mds_total),
    by = "subject_id"
  ) %>%
  mutate(
    pair = factor(pair, levels = alpha_pairs$pair),
    group = factor(group, levels = c("nonMD", "MD")),
    condition = factor(condition, levels = c("EO", "EC")),
    gender = factor(gender),
    age_z = as.numeric(scale(age)),
    mds_z = as.numeric(scale(mds_total))
  )

# Check sample sizes and missingness
asymmetry_data %>%
  group_by(pair, group, condition) %>%
  summarise(
    n = sum(!is.na(asymmetry)),
    missing = sum(is.na(asymmetry)),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5. Descriptive statistics
# ------------------------------------------------------------

asymmetry_descriptives <- asymmetry_data %>%
  group_by(
    analysis_family,
    pair,
    group,
    condition
  ) %>%
  summarise(
    n = sum(!is.na(asymmetry)),
    mean = mean(asymmetry, na.rm = TRUE),
    sd = sd(asymmetry, na.rm = TRUE),
    se = sd / sqrt(n),
    median = median(asymmetry, na.rm = TRUE),
    .groups = "drop"
  )

asymmetry_descriptives

# ============================================================
# PRIMARY FRONTAL ANALYSIS: F4-F3 AND F8-F7
# ============================================================

faa_primary <- asymmetry_data %>%
  filter(analysis_family == "primary") %>%
  droplevels() %>%
  drop_na(
    asymmetry,
    group,
    condition,
    pair,
    age_z,
    gender
  )

# ------------------------------------------------------------
# 6. Unadjusted group model
# ------------------------------------------------------------

model_group_unadjusted <- lmer(
  asymmetry ~ group * condition * pair +
    (1 | subject_id),
  data = faa_primary,
  REML = FALSE
)

anova(
  model_group_unadjusted,
  type = 3,
  ddf = "Satterthwaite"
)

summary(model_group_unadjusted)

# ------------------------------------------------------------
# 7. Group model adjusted for age and gender
# ------------------------------------------------------------

model_group_adjusted <- lmer(
  asymmetry ~ group * condition * pair +
    age_z + gender +
    (1 | subject_id),
  data = faa_primary,
  REML = FALSE
)

anova(
  model_group_adjusted,
  type = 3,
  ddf = "Satterthwaite"
)

summary(model_group_adjusted)
confint(model_group_adjusted, method = "Wald")

# ------------------------------------------------------------
# 8. Group comparisons within each condition and pair
#
# Positive contrast means MD > nonMD
# ------------------------------------------------------------

emm_group <- emmeans(
  model_group_adjusted,
  ~ group | condition * pair
)

group_contrasts <- contrast(
  emm_group,
  method = "revpairwise",
  adjust = "BH"
)

group_contrasts

# ------------------------------------------------------------
# 9. EO versus EC comparisons within each group and pair
#
# With EO as the first factor level,
# revpairwise gives EC - EO
# ------------------------------------------------------------

emm_condition <- emmeans(
  model_group_adjusted,
  ~ condition | group * pair
)

condition_contrasts <- contrast(
  emm_condition,
  method = "revpairwise",
  adjust = "BH"
)

condition_contrasts

# ------------------------------------------------------------
# 10. Estimated asymmetry values
# ------------------------------------------------------------

estimated_asymmetry <- emmeans(
  model_group_adjusted,
  ~ group * condition * pair
)

estimated_asymmetry

# Test whether asymmetry differs from zero
asymmetry_vs_zero <- test(
  estimated_asymmetry,
  null = 0,
  adjust = "BH"
)

asymmetry_vs_zero

# ============================================================
# DIMENSIONAL MDS-SEVERITY ANALYSIS
# ============================================================

faa_primary_mds <- asymmetry_data %>%
  filter(analysis_family == "primary") %>%
  droplevels() %>%
  drop_na(
    asymmetry,
    mds_z,
    condition,
    pair,
    age_z,
    gender
  )

# ------------------------------------------------------------
# 11. MDS severity × condition × pair model
# ------------------------------------------------------------

model_mds <- lmer(
  asymmetry ~ mds_z * condition * pair +
    age_z + gender +
    (1 | subject_id),
  data = faa_primary_mds,
  REML = FALSE
)

anova(
  model_mds,
  type = 3,
  ddf = "Satterthwaite"
)

summary(model_mds)
confint(model_mds, method = "Wald")

# ------------------------------------------------------------
# 12. MDS slopes within each condition and electrode pair
# ------------------------------------------------------------

mds_slopes <- emtrends(
  model_mds,
  ~ condition * pair,
  var = "mds_z"
)

mds_slopes

# Test each MDS slope against zero
test(
  mds_slopes,
  null = 0,
  adjust = "BH"
)

# Compare EO and EC MDS slopes within each pair
contrast(
  mds_slopes,
  method = "revpairwise",
  by = "pair",
  adjust = "BH"
)

# ============================================================
# SECONDARY ELECTRODE-PAIR ANALYSIS
# ============================================================

faa_secondary <- asymmetry_data %>%
  filter(analysis_family == "secondary") %>%
  droplevels() %>%
  drop_na(
    asymmetry,
    group,
    condition,
    pair,
    age_z,
    gender
  )

# ------------------------------------------------------------
# 13. Secondary group model
# ------------------------------------------------------------

model_secondary <- lmer(
  asymmetry ~ group * condition * pair +
    age_z + gender +
    (1 | subject_id),
  data = faa_secondary,
  REML = FALSE
)

anova(
  model_secondary,
  type = 3,
  ddf = "Satterthwaite"
)

summary(model_secondary)

# Group differences for each secondary pair and condition
secondary_group_emm <- emmeans(
  model_secondary,
  ~ group | condition * pair
)

contrast(
  secondary_group_emm,
  method = "revpairwise",
  adjust = "BH"
)

# Condition differences for each secondary pair and group
secondary_condition_emm <- emmeans(
  model_secondary,
  ~ condition | group * pair
)

contrast(
  secondary_condition_emm,
  method = "revpairwise",
  adjust = "BH"
)

# ============================================================
# SIDE-SPECIFIC ALPHA-POWER MODELS
#
# These determine whether an asymmetry effect is caused by
# changes in the left electrode, right electrode, or both.
# ============================================================

# ------------------------------------------------------------
# 14. F3/F4 model
# ------------------------------------------------------------

f43_side_data <- eeg %>%
  filter(channel %in% c("F3", "F4")) %>%
  mutate(
    hemisphere = case_when(
      channel == "F3" ~ "Left",
      channel == "F4" ~ "Right"
    ),
    hemisphere = factor(
      hemisphere,
      levels = c("Left", "Right")
    ),
    log_alpha = log(alpha)
  ) %>%
  left_join(
    metadata %>%
      select(subject_id, age, gender, mds_total),
    by = "subject_id"
  ) %>%
  mutate(
    age_z = as.numeric(scale(age)),
    mds_z = as.numeric(scale(mds_total))
  ) %>%
  drop_na(log_alpha, group, condition, hemisphere, age_z, gender)

model_f43_sides <- lmer(
  log_alpha ~ group * condition * hemisphere +
    age_z + gender +
    (1 | subject_id),
  data = f43_side_data,
  REML = FALSE
)

anova(
  model_f43_sides,
  type = 3,
  ddf = "Satterthwaite"
)

summary(model_f43_sides)

emmeans(
  model_f43_sides,
  ~ hemisphere | group * condition
)

# ------------------------------------------------------------
# 15. F7/F8 model
# ------------------------------------------------------------

f87_side_data <- eeg %>%
  filter(channel %in% c("F7", "F8")) %>%
  mutate(
    hemisphere = case_when(
      channel == "F7" ~ "Left",
      channel == "F8" ~ "Right"
    ),
    hemisphere = factor(
      hemisphere,
      levels = c("Left", "Right")
    ),
    log_alpha = log(alpha)
  ) %>%
  left_join(
    metadata %>%
      select(subject_id, age, gender, mds_total),
    by = "subject_id"
  ) %>%
  mutate(
    age_z = as.numeric(scale(age)),
    mds_z = as.numeric(scale(mds_total))
  ) %>%
  drop_na(log_alpha, group, condition, hemisphere, age_z, gender)

model_f87_sides <- lmer(
  log_alpha ~ group * condition * hemisphere +
    age_z + gender +
    (1 | subject_id),
  data = f87_side_data,
  REML = FALSE
)

anova(
  model_f87_sides,
  type = 3,
  ddf = "Satterthwaite"
)

summary(model_f87_sides)

emmeans(
  model_f87_sides,
  ~ hemisphere | group * condition
)

# ============================================================
# PLOTS
# ============================================================

# ------------------------------------------------------------
# 16. Primary frontal asymmetry plot
# ------------------------------------------------------------

ggplot(
  faa_primary,
  aes(
    x = condition,
    y = asymmetry,
    group = subject_id
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(alpha = 0.18) +
  geom_point(alpha = 0.45) +
  stat_summary(
    aes(group = group),
    fun = mean,
    geom = "line",
    linewidth = 1.1
  ) +
  stat_summary(
    aes(group = group),
    fun = mean,
    geom = "point",
    size = 3
  ) +
  facet_grid(group ~ pair) +
  theme_classic() +
  labs(
    x = "Resting-state condition",
    y = "Alpha asymmetry: ln(right) - ln(left)",
    title = "Frontal Alpha Asymmetry"
  )

# ------------------------------------------------------------
# 17. Group-level boxplot
# ------------------------------------------------------------

ggplot(
  faa_primary,
  aes(
    x = condition,
    y = asymmetry,
    fill = group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_boxplot(
    alpha = 0.60,
    outlier.shape = NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.10,
    alpha = 0.65
  ) +
  facet_wrap(~ pair) +
  theme_classic() +
  labs(
    x = "Resting-state condition",
    y = "Alpha asymmetry: ln(right) - ln(left)",
    title = "Alpha Asymmetry by Group and Condition"
  )

# ------------------------------------------------------------
# 18. MDS severity scatterplots
# ------------------------------------------------------------

ggplot(
  faa_primary_mds,
  aes(
    x = mds_total,
    y = asymmetry
  )
) +
  geom_point(alpha = 0.70) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  facet_grid(condition ~ pair) +
  theme_classic() +
  labs(
    x = "MDS-16 total score",
    y = "Alpha asymmetry: ln(right) - ln(left)",
    title = "MD Severity and Frontal Alpha Asymmetry"
  )

# ============================================================
# BASIC MODEL DIAGNOSTICS
# ============================================================

# Group model residuals
par(mfrow = c(1, 2))

plot(
  fitted(model_group_adjusted),
  resid(model_group_adjusted),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals versus fitted"
)

abline(h = 0, lty = 2)

qqnorm(
  resid(model_group_adjusted),
  main = "Residual Q-Q plot"
)

qqline(resid(model_group_adjusted))

par(mfrow = c(1, 1))

# Check singularity
isSingular(model_group_adjusted)
isSingular(model_mds)
isSingular(model_secondary)

# ============================================================
# VISUAL INSPECTION OF ALPHA-ASYMMETRY RESULTS
# ============================================================

library(tidyverse)
library(emmeans)

# ------------------------------------------------------------
# 1. Raw primary asymmetry distributions
#    F4-F3 and F8-F7
# ------------------------------------------------------------

p_primary_distribution <- ggplot(
  faa_primary,
  aes(
    x = condition,
    y = asymmetry,
    fill = group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.60,
    alpha = 0.45,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = group),
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.75
    ),
    alpha = 0.65,
    size = 2
  ) +
  stat_summary(
    aes(group = group),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.75),
    shape = 23,
    size = 3.5,
    fill = "white"
  ) +
  facet_wrap(~ pair) +
  theme_classic(base_size = 13) +
  labs(
    title = "Primary Frontal Alpha-Asymmetry Distributions",
    subtitle = "Positive values indicate greater right than left alpha power",
    x = "Condition",
    y = "Alpha asymmetry: ln(right) - ln(left)",
    fill = "Group",
    color = "Group"
  )

print(p_primary_distribution)
