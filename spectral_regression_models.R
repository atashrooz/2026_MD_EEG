# ============================================================
# Spectral Regression Models
# EEG MD Project
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(readxl)
library(effectsize)
library(broom)
library(emmeans)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggridges)
library(grid)

options(contrasts = c("contr.sum", "contr.poly"),
        scipen = 999)
# ------------------------------------------------------------
# 1. Read files
# ------------------------------------------------------------

eeg <- read_csv("eeg_band_power_features.csv")

phase2 <- read_excel("md_vviq_suis_sorted_2.xlsx")

# ------------------------------------------------------------
# 2. Create participant-level metadata
# ------------------------------------------------------------

eeg_subjects <- eeg %>%
  distinct(subject_id, group) %>%
  mutate(
    subject_id = as.integer(subject_id),
    id = paste0("Sub", subject_id),
    group = factor(group, levels = c("nonMD", "MD"))
  ) %>%
  arrange(subject_id)

# Define questionnaire items
mds_items  <- paste0("mds", 1:16)
vviq_items <- paste0("vviq", 1:32)
suis_items <- paste0("suis", 1:12)
dass_items <- paste0("dass", 1:21)

# DASS-21 standard subscale item mapping
dass_depression_items <- paste0("dass", c(3, 5, 10, 13, 16, 17, 21))
dass_anxiety_items    <- paste0("dass", c(2, 4, 7, 9, 15, 19, 20))
dass_stress_items     <- paste0("dass", c(1, 6, 8, 11, 12, 14, 18))

# Clean Phase 2 data
phase2_clean <- phase2 %>%
  mutate(
    id = as.character(id),
    md_status = as.character(md_status),
    age = readr::parse_number(as.character(age)),
    gender = as.character(gender),
    marital_status = as.character(marital_status),
    education = as.character(education),
    across(
      all_of(c("mds_total", mds_items, vviq_items, suis_items, dass_items)),
      ~ readr::parse_number(as.character(.x))
    )
  )

# Helper function: complete-case sum score
score_sum_complete <- function(data, items, multiplier = 1) {
  item_data <- data %>%
    select(all_of(items))
  
  ifelse(
    rowSums(is.na(item_data)) == 0,
    rowSums(item_data, na.rm = TRUE) * multiplier,
    NA_real_
  )
}

# Merge metadata
metadata_raw <- eeg_subjects %>%
  left_join(phase2_clean, by = "id") %>%
  mutate(
    md_status_clean = recode(
      md_status,
      "non_MD" = "nonMD",
      "nonMD" = "nonMD",
      "MD" = "MD"
    ),
    group = factor(group, levels = c("nonMD", "MD")),
    md_status_clean = factor(md_status_clean, levels = c("nonMD", "MD")),
    gender = factor(gender),
    marital_status = factor(marital_status),
    education = factor(education)
  )

# Check group matching
table(metadata_raw$group, metadata_raw$md_status_clean)

metadata_raw %>%
  filter(as.character(group) != as.character(md_status_clean)) %>%
  select(subject_id, id, group, md_status, md_status_clean)

# Create final metadata
metadata <- metadata_raw

metadata$mds_raw_sum <- score_sum_complete(metadata, mds_items)

metadata$mds_total_from_items <- score_sum_complete(
  metadata,
  mds_items,
  multiplier = 10 / 16
)

metadata$vviq_total <- score_sum_complete(metadata, vviq_items)
metadata$suis_total <- score_sum_complete(metadata, suis_items)


metadata$dass_depression <- score_sum_complete(
  metadata,
  dass_depression_items,
  multiplier = 2
)

metadata$dass_anxiety <- score_sum_complete(
  metadata,
  dass_anxiety_items,
  multiplier = 2
)

metadata$dass_stress <- score_sum_complete(
  metadata,
  dass_stress_items,
  multiplier = 2
)


metadata <- metadata %>%
  mutate(
    has_age = !is.na(age),
    has_mds = !is.na(mds_total),
    has_vviq = !is.na(vviq_total),
    has_suis = !is.na(suis_total),
    has_dass_depression = !is.na(dass_depression),
    has_dass_anxiety = !is.na(dass_anxiety),
    has_dass_stress = !is.na(dass_stress),
    complete_core = complete.cases(age, gender, mds_total)
  )


eeg_metadata <- metadata %>%
  select(
    subject_id,
    id,
    group,
    md_status_clean,
    age,
    gender,
    marital_status,
    education,
    mds_total,
    mds_raw_sum,
    mds_total_from_items,
    vviq_total,
    suis_total,
    dass_depression,
    dass_anxiety,
    dass_stress,
    has_age,
    has_mds,
    has_vviq,
    has_suis,
    has_dass_depression,
    has_dass_anxiety,
    has_dass_stress,
    complete_core
  ) %>%
  arrange(subject_id)


# Save reusable metadata
# write_csv(eeg_metadata, "eeg_participant_metadata.csv")
# writexl::write_xlsx(eeg_metadata, "eeg_participant_metadata.xlsx")
# writexl::write_xlsx(metadata, "eeg_metadata.xlsx")


# =========================================
#           Reliability Assessment
# ===========================================


# Use the EEG participant sample
reliability_data <- metadata_raw

# Ensure all questionnaire items are numeric
all_scale_items <- c(
  mds_items,
  vviq_items,
  suis_items,
  dass_items
)

reliability_data <- reliability_data %>%
  mutate(
    across(
      all_of(all_scale_items),
      ~ readr::parse_number(as.character(.x))
    )
  )

# ------------------------------------------------------------
# 1. MDS-16
# ------------------------------------------------------------

mds_reliability_data <- reliability_data %>%
  select(all_of(mds_items)) %>%
  drop_na()

mds_alpha <- psych::alpha(
  mds_reliability_data,
  check.keys = FALSE,
  warnings = FALSE
)

mds_omega <- psych::omega(
  mds_reliability_data,
  nfactors = 1,
  fm = "minres",
  plot = FALSE,
  flip = FALSE
)

summary(mds_alpha)
summary(mds_omega)

# ------------------------------------------------------------
# 2. VVIQ-32
# ------------------------------------------------------------

vviq_reliability_data <- reliability_data %>%
  select(all_of(vviq_items)) %>%
  drop_na()

vviq_alpha <- psych::alpha(
  vviq_reliability_data,
  check.keys = FALSE,
  warnings = FALSE
)

vviq_omega <- psych::omega(
  vviq_reliability_data,
  nfactors = 1,
  fm = "minres",
  plot = FALSE,
  flip = FALSE
)

summary(vviq_alpha)
summary(vviq_omega)

# ------------------------------------------------------------
# 3. SUIS-12
# ------------------------------------------------------------

suis_reliability_data <- reliability_data %>%
  select(all_of(suis_items)) %>%
  drop_na()

suis_alpha <- psych::alpha(
  suis_reliability_data,
  check.keys = FALSE,
  warnings = FALSE
)

suis_omega <- psych::omega(
  suis_reliability_data,
  nfactors = 1,
  fm = "minres",
  plot = FALSE,
  flip = FALSE
)

summary(suis_alpha)
summary(suis_omega)

# ------------------------------------------------------------
# 4. DASS-21 Depression
# ------------------------------------------------------------

dass_depression_reliability_data <- reliability_data %>%
  select(all_of(dass_depression_items)) %>%
  drop_na()

dass_depression_alpha <- psych::alpha(
  dass_depression_reliability_data,
  check.keys = FALSE,
  warnings = FALSE
)

dass_depression_omega <- psych::omega(
  dass_depression_reliability_data,
  nfactors = 1,
  fm = "minres",
  plot = FALSE,
  flip = FALSE
)

summary(dass_depression_alpha)
summary(dass_depression_omega)

# ------------------------------------------------------------
# 5. DASS-21 Anxiety
# ------------------------------------------------------------

dass_anxiety_reliability_data <- reliability_data %>%
  select(all_of(dass_anxiety_items)) %>%
  drop_na()

dass_anxiety_alpha <- psych::alpha(
  dass_anxiety_reliability_data,
  check.keys = FALSE,
  warnings = FALSE
)

dass_anxiety_omega <- psych::omega(
  dass_anxiety_reliability_data,
  nfactors = 1,
  fm = "minres",
  plot = FALSE,
  flip = FALSE
)

summary(dass_anxiety_alpha)
summary(dass_anxiety_omega)

# ------------------------------------------------------------
# 6. DASS-21 Stress
# ------------------------------------------------------------

dass_stress_reliability_data <- reliability_data %>%
  select(all_of(dass_stress_items)) %>%
  drop_na()

dass_stress_alpha <- psych::alpha(
  dass_stress_reliability_data,
  check.keys = FALSE,
  warnings = FALSE
)

dass_stress_omega <- psych::omega(
  dass_stress_reliability_data,
  nfactors = 1,
  fm = "minres",
  plot = FALSE,
  flip = FALSE
)

summary(dass_stress_alpha)
summary(dass_stress_omega)


# ------------------------------------------------------------
# 3. Define ROIs
# ------------------------------------------------------------

posterior_roi <- c(
  "O1", "Oz", "O2",
  "PO3", "POz", "PO4",
  "P3", "Pz", "P4"
)

frontal_roi <- c(
  "F7", "F3", "Fz", "F4", "F8",
  "FC1", "FCz", "FC2",
  "FC3", "FC4"
)

# Check missing channels
setdiff(posterior_roi, unique(eeg$channel))
setdiff(frontal_roi, unique(eeg$channel))

# ------------------------------------------------------------
# 4. Create ROI-level spectral variables
# ------------------------------------------------------------

eeg_roi <- eeg %>%
  mutate(
    subject_id = as.integer(subject_id),
    condition = factor(condition, levels = c("EO", "EC")),
    group = factor(group, levels = c("nonMD", "MD")),
    
    delta_db = 10 * log10(delta),
    theta_db = 10 * log10(theta),
    alpha_db = 10 * log10(alpha),
    beta_db  = 10 * log10(beta),
    gamma_db = 10 * log10(gamma)
  )

# Posterior alpha
posterior_alpha <- eeg_roi %>%
  filter(channel %in% posterior_roi) %>%
  group_by(subject_id, group, condition) %>%
  summarise(
    posterior_alpha_db = mean(alpha_db, na.rm = TRUE),
    .groups = "drop"
  )

# Frontal theta, beta, and theta/beta ratio
frontal_spectral <- eeg_roi %>%
  filter(channel %in% frontal_roi) %>%
  group_by(subject_id, group, condition) %>%
  summarise(
    frontal_theta_db = mean(theta_db, na.rm = TRUE),
    frontal_beta_db  = mean(beta_db, na.rm = TRUE),
    
    # dB theta/beta ratio:
    # 10log10(theta/beta) = theta_db - beta_db
    frontal_tbr_db = frontal_theta_db - frontal_beta_db,
    
    .groups = "drop"
  )

# Combine ROI outcomes
spectral_roi <- posterior_alpha %>%
  left_join(
    frontal_spectral,
    by = c("subject_id", "group", "condition")
  ) %>%
  left_join(
    eeg_metadata,
    by = c("subject_id", "group")
  ) %>%
  mutate(
    condition = factor(condition, levels = c("EO", "EC")),
    group = factor(group, levels = c("nonMD", "MD")),
    gender = factor(gender)
  )

spectral_roi <- spectral_roi %>%
  mutate(
    subject_id = factor(subject_id),
    condition = factor(condition, levels = c("EO", "EC")),
    group = factor(group, levels = c("nonMD", "MD")),
    gender = factor(gender),
    
    age_z = as.numeric(scale(age)),
    mds_z = as.numeric(scale(mds_total)),
    vviq_z = as.numeric(scale(vviq_total)),
    suis_z = as.numeric(scale(suis_total)),
    dass_depression_z = as.numeric(scale(dass_depression)),
    dass_anxiety_z = as.numeric(scale(dass_anxiety)),
    dass_stress_z = as.numeric(scale(dass_stress))
  )

# Save ROI-level spectral dataset
# write_csv(spectral_roi, "spectral_roi_with_metadata.csv")

# ------------------------------------------------------------
# 5. Check final spectral dataset
# ------------------------------------------------------------

spectral_roi %>%
  count(group, condition)

spectral_roi %>%
  summarise(
    n_subjects = n_distinct(subject_id),
    n_rows = n(),
    n_missing_age = sum(is.na(age)),
    n_missing_mds = sum(is.na(mds_total)),
    n_missing_vviq = sum(is.na(vviq_total)),
    n_missing_suis = sum(is.na(suis_total))
  )

# ------------------------------------------------------------
# 6. Spectral outcomes
# ------------------------------------------------------------

spectral_outcomes <- c(
  "posterior_alpha_db",
  "frontal_theta_db",
  "frontal_beta_db",
  "frontal_tbr_db"
)

# ------------------------------------------------------------
# 7. Descriptive statistics for spectral ROI outcomes
# ------------------------------------------------------------

spectral_descriptives <- spectral_roi %>%
  select(subject_id, group, condition, all_of(spectral_outcomes)) %>%
  pivot_longer(
    cols = all_of(spectral_outcomes),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  group_by(outcome, group, condition) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

spectral_descriptives

# write_csv(spectral_descriptives, "spectral_roi_descriptives.csv")

# ------------------------------------------------------------
# 8. Helper function for mixed models
# ------------------------------------------------------------

run_lmer_model <- function(data, outcome, predictor_formula, model_label) {
  
  formula_text <- paste0(
    outcome,
    " ~ ",
    predictor_formula,
    " + (1 | subject_id)"
  )
  
  model <- lmer(
    as.formula(formula_text),
    data = data,
    REML = TRUE
  )
  
  coefficient_table <- broom.mixed::tidy(
    model,
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(
      outcome = .env$outcome,
      model = .env$model_label,
      n_subjects = n_distinct(data$subject_id),
      n_observations = nrow(data)
    ) %>%
    select(
      model,
      outcome,
      n_subjects,
      n_observations,
      term,
      estimate,
      std.error,
      statistic,
      df,
      p.value,
      conf.low,
      conf.high
    )
  
  omnibus_table <- anova(
    model,
    type = 3,
    ddf = "Satterthwaite"
  ) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    rename(
      df_num = NumDF,
      df_den = DenDF,
      statistic = `F value`,
      p.value = `Pr(>F)`
    ) %>%
    mutate(
      partial_eta2 =
        (statistic * df_num) /
        (statistic * df_num + df_den),
      
      outcome = .env$outcome,
      model = .env$model_label,
      n_subjects = n_distinct(data$subject_id),
      n_observations = nrow(data)
    ) %>%
    select(
      model,
      outcome,
      n_subjects,
      n_observations,
      term,
      df_num,
      df_den,
      statistic,
      partial_eta2,
      p.value
    )
  
  list(
    model = model,
    coefficients = coefficient_table,
    omnibus = omnibus_table
  )
}

# ------------------------------------------------------------
# 9. Primary analysis:
#  group × condition models
# ------------------------------------------------------------

primary_group_results <- map(
  spectral_outcomes,
  ~ run_lmer_model(
    data = spectral_roi %>%
      drop_na(.data[[.x]], group, condition),
    outcome = .x,
    predictor_formula = "group * condition",
    model_label = "primary_group_condition"
  )
)

primary_group_coefficients <- map_dfr(
  primary_group_results,
  ~ .x$coefficients
)

primary_group_omnibus <- map_dfr(
  primary_group_results,
  ~ .x$omnibus
)

# write_csv(
#   primary_group_coefficients,
#   "primary_group_condition_coefficients.csv"
# )
# 
# write_csv(
#   primary_group_omnibus,
#   "primary_group_condition_omnibus.csv"
# )

# ------------------------------------------------------------
# 10. Sensitivity analysis:
# Group × condition adjusted for age and gender
# ------------------------------------------------------------

adjusted_group_results <- map(
  spectral_outcomes,
  ~ run_lmer_model(
    data = spectral_roi %>%
      drop_na(
        .data[[.x]],
        group,
        condition,
        age_z,
        gender
      ),
    outcome = .x,
    predictor_formula =
      "group * condition + age_z + gender",
    model_label = "sensitivity_group_age_gender"
  )
)

adjusted_group_coefficients <- map_dfr(
  adjusted_group_results,
  ~ .x$coefficients
)

adjusted_group_omnibus <- map_dfr(
  adjusted_group_results,
  ~ .x$omnibus
)

# 
# write_csv(
#   adjusted_group_coefficients,
#   "sensitivity_group_age_gender_coefficients.csv"
# )
# 
# write_csv(
#   adjusted_group_omnibus,
#   "sensitivity_group_age_gender_omnibus.csv"
# )




spectral_extra <- spectral_roi



spectral_extra %>%
  count(group, condition)

# Helper for exploratory predictor models
run_exploratory_predictor <- function(data, outcome, predictor, model_label) {
  
  predictor_formula <- paste0(
    predictor,
    " * condition + age_z + gender"
  )
  
  run_lmer_model(
    data = data %>%
      drop_na(.data[[outcome]], .data[[predictor]], condition, age_z, gender),
    outcome = outcome,
    predictor_formula = predictor_formula,
    model_label = model_label
  )
}

exploratory_predictors <- c(
  "vviq_z",
  "suis_z",
  "dass_depression_z",
  "dass_anxiety_z",
  "dass_stress_z"
)

exploratory_results <- list()

for (pred in exploratory_predictors) {
  for (out in spectral_outcomes) {
    
    model_label <- paste0(pred, "_condition_age_gender")
    
    exploratory_results[[paste(out, pred, sep = "_")]] <-
      run_exploratory_predictor(
        data = spectral_extra,
        outcome = out,
        predictor = pred,
        model_label = model_label
      )
  }
}

exploratory_tables <- map_dfr(
  exploratory_results,
  ~ .x$coefficients
)

exploratory_tables_fdr <- exploratory_tables %>%
  group_by(model) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "fdr")
  ) %>%
  ungroup()

exploratory_tables_fdr

# write_csv(
#   exploratory_tables_fdr,
#   "spectral_regression_exploratory_vviq_suis_dass_models.csv"
# )

run_covariate_sensitivity <- function(data, outcomes, covariates) {
  
  results <- list()
  
  for (covariate in covariates) {
    for (outcome in outcomes) {
      
      results[[paste(outcome, covariate, sep = "_")]] <-
        run_lmer_model(
          data = data %>%
            drop_na(all_of(c(
              outcome, "group", "condition",
              "age_z", "gender", covariate
            ))),
          outcome = outcome,
          predictor_formula = paste0(
            "group * condition + age_z + gender + ",
            covariate
          ),
          model_label = paste0("sensitivity_", covariate)
        )
    }
  }
  
  list(
    coefficients = map_dfr(results, ~ .x$coefficients),
    omnibus = map_dfr(results, ~ .x$omnibus)
  )
}

covariate_sensitivity <- run_covariate_sensitivity(
  data = spectral_roi,
  outcomes = spectral_outcomes,
  covariates = c(
    "vviq_z",
    "suis_z",
    "dass_depression_z",
    "dass_anxiety_z",
    "dass_stress_z"
  )
)

covariate_sensitivity$omnibus
covariate_sensitivity$coefficients

# write_csv(
#   covariate_sensitivity$coefficients,
#   "covariate_sensitivity_coefficients.csv"
# )
# 
# write_csv(
#   covariate_sensitivity$omnibus,
#   "covariate_sensitivity_omnibus.csv"
# )



adjust_fdr <- function(table, effects) {
  table %>%
    filter(term %in% effects) %>%
    group_by(term) %>%
    mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

primary_group_omnibus_fdr <- adjust_fdr(
  primary_group_omnibus,
  c("group", "condition", "group:condition")
)

adjusted_group_omnibus_fdr <- adjust_fdr(
  adjusted_group_omnibus,
  c("group", "condition", "group:condition")
)


# write_csv(primary_group_omnibus_fdr,
#           "primary_group_condition_omnibus_fdr.csv")
# 
# write_csv(adjusted_group_omnibus_fdr,
#           "sensitivity_group_age_gender_omnibus_fdr.csv")


primary_models <- set_names(
  primary_group_results,
  spectral_outcomes
)

interaction_outcomes <- c(
  "frontal_theta_db",
  "frontal_tbr_db"
)

group_simple_effects <- map_dfr(
  interaction_outcomes,
  function(outcome_name) {
    
    emmeans(
      primary_models[[outcome_name]]$model,
      ~ group | condition
    ) %>%
      contrast(
        method = list("MD - nonMD" = c(-1, 1))
      ) %>%
      summary(infer = c(TRUE, TRUE), adjust = "none") %>%
      as_tibble() %>%
      mutate(outcome = outcome_name, .before = 1)
  }
) %>%
  group_by(outcome) %>%
  mutate(p_holm = p.adjust(p.value, method = "holm")) %>%
  ungroup()

condition_simple_effects <- map_dfr(
  interaction_outcomes,
  function(outcome_name) {
    
    emmeans(
      primary_models[[outcome_name]]$model,
      ~ condition | group
    ) %>%
      contrast(
        method = list("EC - EO" = c(-1, 1))
      ) %>%
      summary(infer = c(TRUE, TRUE), adjust = "none") %>%
      as_tibble() %>%
      mutate(outcome = outcome_name, .before = 1)
  }
) %>%
  group_by(outcome) %>%
  mutate(p_holm = p.adjust(p.value, method = "holm")) %>%
  ungroup()

# write_csv(
#   group_simple_effects,
#   "primary_group_simple_effects.csv"
# )
# 
# write_csv(
#   condition_simple_effects,
#   "primary_condition_simple_effects.csv"
# )


# ------------------------------------------------------------
# 17B. Relative spectral power:
# calculation, mixed-model analysis, and visualization
# ------------------------------------------------------------

relative_bands <- c("delta", "theta", "alpha", "beta", "gamma")

missing_relative_bands <- setdiff(relative_bands, names(eeg))

if (length(missing_relative_bands) > 0) {
  stop(
    "These band-power columns are missing from eeg: ",
    paste(missing_relative_bands, collapse = ", ")
  )
}

if (any(as.matrix(eeg[relative_bands]) < 0, na.rm = TRUE)) {
  stop("Negative power values were detected. Relative power requires non-negative linear power values.")
}

# Relative power = PSD_band / PSD_total, where PSD_total is the sum
# of the five bands' mean PSD values 

relative_power_channel <- eeg %>%
  mutate(
    subject_id = as.integer(subject_id),
    group = factor(group, levels = c("nonMD", "MD")),
    condition = factor(condition, levels = c("EO", "EC")),
    across(all_of(relative_bands), as.numeric)
  ) %>%
  mutate(
    total_analyzed_power =
      delta + theta + alpha + beta + gamma,

    valid_total_power =
      is.finite(total_analyzed_power) & total_analyzed_power > 0,

    delta_rel = if_else(
      valid_total_power,
      delta / total_analyzed_power,
      NA_real_
    ),
    theta_rel = if_else(
      valid_total_power,
      theta / total_analyzed_power,
      NA_real_
    ),
    alpha_rel = if_else(
      valid_total_power,
      alpha / total_analyzed_power,
      NA_real_
    ),
    beta_rel = if_else(
      valid_total_power,
      beta / total_analyzed_power,
      NA_real_
    ),
    gamma_rel = if_else(
      valid_total_power,
      gamma / total_analyzed_power,
      NA_real_
    )
  )

# Convert channel-level relative power to long format and average within ROI.
relative_power_roi <- relative_power_channel %>%
  select(
    subject_id,
    group,
    condition,
    channel,
    all_of(paste0(relative_bands, "_rel"))
  ) %>%
  pivot_longer(
    cols = all_of(paste0(relative_bands, "_rel")),
    names_to = "band",
    names_pattern = "^(.*)_rel$",
    values_to = "relative_power"
  ) %>%
  mutate(
    region = case_when(
      channel %in% posterior_roi ~ "Posterior",
      channel %in% frontal_roi ~ "Frontal",
      TRUE ~ NA_character_
    ),
    band = factor(band, levels = relative_bands)
  ) %>%
  filter(!is.na(region)) %>%
  group_by(subject_id, group, condition, region, band) %>%
  summarise(
    n_channels_used = sum(!is.na(relative_power)),
    relative_power = mean(relative_power, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    relative_power = if_else(
      is.nan(relative_power),
      NA_real_,
      relative_power
    )
  ) %>%
  left_join(
    eeg_metadata,
    by = c("subject_id", "group")
  )

# Center and scale age using participant-level metadata.
relative_age_mean <- mean(eeg_metadata$age, na.rm = TRUE)

relative_age_sd <- sd(eeg_metadata$age, na.rm = TRUE)

if (!is.finite(relative_age_sd) || relative_age_sd == 0) {
  stop("Age cannot be standardized because its SD is zero or undefined.")
}

# Percentages are easier to interpret and plot.
# Logit-transformed proportions are used as model outcomes because
# untransformed relative power is bounded between 0 and 1.
relative_power_roi <- relative_power_roi %>%
  mutate(
    subject_id = factor(subject_id),
    group = factor(group, levels = c("nonMD", "MD")),
    condition = factor(condition, levels = c("EO", "EC")),
    region = factor(region, levels = c("Posterior", "Frontal")),
    band = factor(band, levels = relative_bands),
    gender = factor(gender),
    age_z = (age - relative_age_mean) / relative_age_sd,
    relative_power_pct = 100 * relative_power,

    # Clamp only for the mathematical logit transformation.
    relative_power_clamped = pmin(
      pmax(relative_power, 0.000001),
      0.999999
    ),
    relative_power_logit = qlogis(relative_power_clamped)
  )

# write_csv(
#   relative_power_roi,
#   "spectral_relative_power_roi_with_metadata.csv"
# )

# Check that the five relative bands sum to approximately 1 within each ROI.
relative_power_sum_check <- relative_power_roi %>%
  group_by(subject_id, group, condition, region) %>%
  summarise(
    total_relative_power = sum(relative_power, na.rm = FALSE),
    .groups = "drop"
  ) %>%
  mutate(
    absolute_deviation_from_one = abs(total_relative_power - 1)
  )

relative_power_sum_check %>%
  summarise(
    maximum_absolute_deviation = max(
      absolute_deviation_from_one,
      na.rm = TRUE
    )
  ) %>%
  print()

# write_csv(
#   relative_power_sum_check,
#   "spectral_relative_power_sum_check.csv"
# )

# ------------------------------------------------------------
# Relative-power descriptive statistics
# ------------------------------------------------------------


relative_power_descriptives <- relative_power_roi %>%
  group_by(region, band, group, condition) %>%
  summarise(
    n = sum(!is.na(relative_power_pct)),
    mean_pct = mean(relative_power_pct, na.rm = TRUE),
    sd_pct = sd(relative_power_pct, na.rm = TRUE),
    median_pct = median(relative_power_pct, na.rm = TRUE),
    min_pct = min(relative_power_pct, na.rm = TRUE),
    max_pct = max(relative_power_pct, na.rm = TRUE),
    .groups = "drop"
  )

relative_power_descriptives

# write_csv(
#   relative_power_descriptives,
#   "spectral_relative_power_descriptives.csv"
# )

# ------------------------------------------------------------
# Relative-power mixed-effects models
# Primary: group x condition
# Sensitivity: group x condition + age + gender
# ------------------------------------------------------------

run_relative_power_model <- function(
    data,
    region_name,
    band_name,
    adjusted = FALSE
) {
  predictor_formula <- if (adjusted) {
    "group * condition + age_z + gender"
  } else {
    "group * condition"
  }

  required_variables <- if (adjusted) {
    c(
      "relative_power_logit",
      "subject_id",
      "group",
      "condition",
      "age_z",
      "gender"
    )
  } else {
    c(
      "relative_power_logit",
      "subject_id",
      "group",
      "condition"
    )
  }

  model_data <- data %>%
    filter(
      region == region_name,
      band == band_name
    ) %>%
    drop_na(all_of(required_variables))

  model_label <- if (adjusted) {
    "relative_power_group_condition_age_gender"
  } else {
    "relative_power_group_condition"
  }

  fitted <- run_lmer_model(
    data = model_data,
    outcome = "relative_power_logit",
    predictor_formula = predictor_formula,
    model_label = model_label
  )

  list(
    model = fitted$model,
    coefficients = fitted$coefficients %>%
      mutate(
        region = region_name,
        band = band_name,
        .before = model
      ),
    omnibus = fitted$omnibus %>%
      mutate(
        region = region_name,
        band = band_name,
        .before = model
      )
  )
}

relative_model_grid <- expand_grid(
  region = c("Posterior", "Frontal"),
  band = relative_bands
)

# Primary models
relative_primary_results <- pmap(
  relative_model_grid,
  function(region, band) {
    run_relative_power_model(
      data = relative_power_roi,
      region_name = region,
      band_name = band,
      adjusted = FALSE
    )
  }
)

relative_primary_coefficients <- map_dfr(
  relative_primary_results,
  ~ .x$coefficients
)

relative_primary_omnibus <- map_dfr(
  relative_primary_results,
  ~ .x$omnibus
)

# FDR is applied separately to each omnibus effect across the 10 outcomes:
# 2 ROIs x 5 bands.
relative_primary_omnibus_fdr <- relative_primary_omnibus %>%
  filter(term %in% c("group", "condition", "group:condition")) %>%
  group_by(term) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup()
# 
# write_csv(
#   relative_primary_coefficients,
#   "relative_power_primary_coefficients.csv"
# )
# 
# write_csv(
#   relative_primary_omnibus,
#   "relative_power_primary_omnibus.csv"
# )
# 
# write_csv(
#   relative_primary_omnibus_fdr,
#   "relative_power_primary_omnibus_fdr.csv"
# )

# Age- and gender-adjusted sensitivity models
relative_adjusted_results <- pmap(
  relative_model_grid,
  function(region, band) {
    run_relative_power_model(
      data = relative_power_roi,
      region_name = region,
      band_name = band,
      adjusted = TRUE
    )
  }
)

relative_adjusted_coefficients <- map_dfr(
  relative_adjusted_results,
  ~ .x$coefficients
)

relative_adjusted_omnibus <- map_dfr(
  relative_adjusted_results,
  ~ .x$omnibus
)

relative_adjusted_omnibus_fdr <- relative_adjusted_omnibus %>%
  filter(term %in% c("group", "condition", "group:condition")) %>%
  group_by(term) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup()

# write_csv(
#   relative_adjusted_coefficients,
#   "relative_power_adjusted_coefficients.csv"
# )
# 
# write_csv(
#   relative_adjusted_omnibus,
#   "relative_power_adjusted_omnibus.csv"
# )
# 
# write_csv(
#   relative_adjusted_omnibus_fdr,
#   "relative_power_adjusted_omnibus_fdr.csv"
# )

# Print the FDR-corrected inferential results.
relative_primary_omnibus_fdr %>%
  arrange(term, p_fdr) %>%
  print(n = Inf)

relative_adjusted_omnibus_fdr %>%
  arrange(term, p_fdr) %>%
  print(n = Inf)

# ------------------------------------------------------------
# Relative-power visualization
# ------------------------------------------------------------

relative_group_colors <- c(
  nonMD = "#0072B2",
  MD = "#D55E00"
)

relative_group_shapes <- c(
  nonMD = 16,
  MD = 17
)

relative_group_linetypes <- c(
  nonMD = "solid",
  MD = "solid"
)

relative_group_labels <- c(
  nonMD = "non MDer",
  MD = "Probable MDer"
)

relative_condition_labels <- c(
  EO = "EO",
  EC = "EC"
)

relative_band_labels <- c(
  delta = "Delta",
  theta = "Theta",
  alpha = "Alpha",
  beta = "Beta",
  gamma = "Gamma"
)

relative_plot_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  )

relative_plot_directory <- "plots"

dir.create(
  relative_plot_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

# Plot 1: mean relative power and 95% CI for every ROI and band.
relative_power_group_condition_plot <- ggplot(
  relative_power_roi,
  aes(
    x = condition,
    y = relative_power_pct,
    group = group,
    color = group,
    shape = group,
    linetype = group
  )
) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.9,
    position = position_dodge(width = 0.08),
    na.rm = TRUE
  ) +
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "pointrange",
    linewidth = 0.65,
    size = 0.8,
    position = position_dodge(width = 0.08),
    na.rm = TRUE
  ) +
  facet_wrap(
    vars(region, band),
    ncol = 5,
    scales = "free_y",
    labeller = labeller(
      band = relative_band_labels
    )
  ) +
  scale_x_discrete(
    limits = c("EO", "EC"),
    labels = relative_condition_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    name = "Group",
    values = relative_group_colors,
    breaks = c("nonMD", "MD"),
    labels = relative_group_labels
  ) +
  scale_shape_manual(
    name = "Group",
    values = relative_group_shapes,
    breaks = c("nonMD", "MD"),
    labels = relative_group_labels
  ) +
  scale_linetype_manual(
    name = "Group",
    values = relative_group_linetypes,
    breaks = c("nonMD", "MD"),
    labels = relative_group_labels
  ) +
  labs(
  # title = "Relative Spectral Power by Group and Condition",
  #  subtitle = "Points represent means; error bars represent 95% confidence intervals",
    x = NULL,
    y = "Relative power (%)"
  ) +
  relative_plot_theme

print(relative_power_group_condition_plot)

# ------------------------------------------------------------
# PSD Visualization
# ------------------------------------------------------------

# Okabe-Ito colorblind-safe colors
group_colors <- c(
  nonMD = "#0072B2",
  MD    = "#D55E00"
)

group_shapes <- c(
  nonMD = 16,
  MD    = 17
)

group_linetypes <- c(
  nonMD = "solid",
  MD    = "solid"
)

group_labels <- c(
  nonMD = "Control",
  MD    = "Probable MDer"
)

condition_labels <- c(
  EO = "EO",
  EC = "EC"
)

spectral_outcome_labels <- c(
  posterior_alpha_db = "Posterior alpha",
  frontal_theta_db   = "Frontal theta",
  frontal_beta_db    = "Frontal beta",
  frontal_tbr_db     = "Frontal theta/beta ratio"
)


predictor_labels <- c(
  vviq_z            = "VVIQ",
  suis_z             = "SUIS",
  dass_depression_z  = "DASS depression",
  dass_anxiety_z     = "DASS anxiety",
  dass_stress_z      = "DASS stress"
)

band_labels <- c(
  delta = "Delta",
  theta = "Theta",
  alpha = "Alpha",
  beta  = "Beta",
  gamma = "Gamma"
)

eeg_plot_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.4, "lines"),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(margin = margin(r = 8)),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 10)
  )

plot_directory <- "plots"

dir.create(
  plot_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
#  Prepare long-format spectral data
# ------------------------------------------------------------

spectral_plot_long <- spectral_roi %>%
  select(
    subject_id,
    group,
    condition,
    all_of(spectral_outcomes)
  ) %>%
  pivot_longer(
    cols = all_of(spectral_outcomes),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  mutate(
    outcome_label = factor(
      spectral_outcome_labels[outcome],
      levels = unname(spectral_outcome_labels)
    ),
    group = factor(group, levels = c("nonMD", "MD")),
    condition = factor(condition, levels = c("EO", "EC"))
  )

spectral_plot_summary <- spectral_plot_long %>%
  group_by(outcome, outcome_label, group, condition) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_multiplier = if_else(n > 1, qt(0.975, df = n - 1), NA_real_),
    ci_low = mean - ci_multiplier * se,
    ci_high = mean + ci_multiplier * se,
    .groups = "drop"
  )

# ------------------------------------------------------------
#  Combined spectral outcomes: group × condition
# ------------------------------------------------------------

plot_all_spectral_group_condition <- ggplot(
  spectral_plot_summary,
  aes(
    x = condition,
    y = mean,
    group = group,
    color = group,
    shape = group,
    linetype = group
  )
) +
  geom_line(
    linewidth = 0.9,
    position = position_dodge(width = 0.08),
    na.rm = TRUE
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.06,
    linewidth = 0.65,
    position = position_dodge(width = 0.08),
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.8,
    position = position_dodge(width = 0.08),
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ outcome_label,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_x_discrete(
    limits = c("EO", "EC"),
    labels = condition_labels,
    drop = FALSE
  ) +
  scale_color_manual(
    name = "Group",
    values = group_colors,
    breaks = c("nonMD", "MD"),
    labels = group_labels
  ) +
  scale_shape_manual(
    name = "Group",
    values = group_shapes,
    breaks = c("nonMD", "MD"),
    labels = group_labels
  ) +
  scale_linetype_manual(
    name = "Group",
    values = group_linetypes,
    breaks = c("nonMD", "MD"),
    labels = group_labels
  ) +
  labs(
    x = NULL,
    y = "Power (dB)"
  ) +
  eeg_plot_theme

print(plot_all_spectral_group_condition)

