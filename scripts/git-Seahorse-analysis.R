
# seahorse data analysis and parameters

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(broom)

file_path <- "P:/your_path/to seahorse data in supplemental table 6"

# Group names
ctrl1_label <- "CTRL1"
control_labels <- c("CTRL1", "CTRL2", "CTRL3")   # pooled controls



# Evaluate timepoints deviating more than 20% from the mean of the 3 timepoints within each well/phase/measure
deviation_threshold <- 0.20

# Measurement numbers for each phase
# Edit if needed
phase_map <- list(
  baseline       = 1:3,
  post_oligo     = 4:6,
  post_fccp      = 7:9,
  post_antimycin = 10:12
)

raw <- read_excel(file_path)

print(names(raw))
print(head(raw))

df <- raw %>%
  rename(
    Well = Well,
    Group = Group,
    Measurement = Measurement,
    OCR = OCR,
    ECAR = ECAR
  ) %>%
  mutate(
    Well = as.character(Well),
    Group = as.character(Group),
    Measurement = as.integer(Measurement),
    OCR = as.numeric(OCR),
    ECAR = as.numeric(ECAR)
  )


# assign assay phases

assign_phase <- function(x) {
  case_when(
    x %in% phase_map$baseline       ~ "baseline",
    x %in% phase_map$post_oligo     ~ "post_oligo",
    x %in% phase_map$post_fccp      ~ "post_fccp",
    x %in% phase_map$post_antimycin ~ "post_antimycin",
    TRUE ~ NA_character_
  )
}

df <- df %>%
  mutate(
    Phase = assign_phase(Measurement),
    Phase = factor(
      Phase,
      levels = c("baseline", "post_oligo", "post_fccp", "post_antimycin")
    )
  ) %>%
  filter(!is.na(Phase))

print(head(df))

# long format
df_long <- df %>%
  pivot_longer(
    cols = c(OCR, ECAR),
    names_to = "MeasureType",
    values_to = "Value"
  )

print(head(df_long))


# flag devaiting time points

df_groom_flagged <- df_long %>%
  group_by(Group, Well, Phase, MeasureType) %>%
  mutate(
    PhaseMean_raw = mean(Value, na.rm = TRUE),
    RelDev = abs(Value - PhaseMean_raw) / ifelse(PhaseMean_raw == 0, NA, abs(PhaseMean_raw)),
    Keep = ifelse(is.na(RelDev), TRUE, RelDev <= deviation_threshold)
  ) %>%
  ungroup()

flagged_points <- df_groom_flagged %>%
  filter(!Keep)

print(head(flagged_points))

df_groomed <- df_groom_flagged %>%
  filter(Keep)

print(head(df_groomed))

# compute means after filtering

phase_means <- df_groomed %>%
  group_by(Group, Well, Phase, MeasureType) %>%
  summarise(
    n_timepoints_kept = sum(!is.na(Value)),
    PhaseMean = mean(Value, na.rm = TRUE),
    PhaseSD = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

print(head(phase_means))

phase_wide <- phase_means %>%
  select(Group, Well, Phase, MeasureType, PhaseMean) %>%
  pivot_wider(
    names_from = c(MeasureType, Phase),
    values_from = PhaseMean
  )

print(head(phase_wide))


# calulate parameters

# Antimycin A used as the non-mitochondrial respiration phase

params <- phase_wide %>%
  mutate(
    OCR_basal = OCR_baseline - OCR_post_antimycin,
    OCR_ATP_linked = OCR_baseline - OCR_post_oligo,
    OCR_proton_leak = OCR_post_oligo - OCR_post_antimycin,
    OCR_maximal = OCR_post_fccp - OCR_post_antimycin,
    OCR_spare_capacity = OCR_maximal - OCR_basal,
    OCR_non_mito = OCR_post_antimycin,
    OCR_coupling_eff = ifelse(OCR_basal == 0, NA, 100 * OCR_ATP_linked / OCR_basal),
    
    ECAR_basal = ECAR_baseline
  )

print(head(params))


# long format
params_long <- params %>%
  select(
    Group, Well,
    OCR_basal,
    OCR_ATP_linked,
    OCR_proton_leak,
    OCR_maximal,
    OCR_spare_capacity,
    OCR_non_mito,
    OCR_coupling_eff,
    ECAR_basal
  ) %>%
  pivot_longer(
    cols = -c(Group, Well),
    names_to = "Parameter",
    values_to = "Value"
  )

print(head(params_long))


#  CTRL1 mean and pooled control mean

reference_means <- params_long %>%
  group_by(Parameter) %>%
  summarise(
    CTRL1_mean = mean(Value[Group == ctrl1_label], na.rm = TRUE),
    pooled_control_mean = mean(Value[Group %in% control_labels], na.rm = TRUE),
    .groups = "drop"
  )

print(head(reference_means))


# values relative to CTRL1 = 100%
# and pooled controls  = 100%

params_norm <- params_long %>%
  left_join(reference_means, by = "Parameter") %>%
  mutate(
    Percent_of_CTRL1 = 100 * Value / CTRL1_mean,
    Percent_of_pooled_controls = 100 * Value / pooled_control_mean
  )

print(head(params_norm))

#summarize 
summary_abs <- params_long %>%
  group_by(Group, Parameter) %>%
  summarise(
    n = sum(!is.na(Value)),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    SEM = SD / sqrt(n),
    .groups = "drop"
  )

summary_ctrl1 <- params_norm %>%
  group_by(Group, Parameter) %>%
  summarise(
    n = sum(!is.na(Percent_of_CTRL1)),
    Mean_percent_CTRL1 = mean(Percent_of_CTRL1, na.rm = TRUE),
    SD_percent_CTRL1 = sd(Percent_of_CTRL1, na.rm = TRUE),
    SEM_percent_CTRL1 = SD_percent_CTRL1 / sqrt(n),
    .groups = "drop"
  )

summary_pooled <- params_norm %>%
  group_by(Group, Parameter) %>%
  summarise(
    n = sum(!is.na(Percent_of_pooled_controls)),
    Mean_percent_pooled = mean(Percent_of_pooled_controls, na.rm = TRUE),
    SD_percent_pooled = sd(Percent_of_pooled_controls, na.rm = TRUE),
    SEM_percent_pooled = SD_percent_pooled / sqrt(n),
    .groups = "drop"
  )

print(head(summary_abs))
print(head(summary_ctrl1))
print(head(summary_pooled))

# test for normality 

normality_by_group <- params_long %>%
  group_by(Parameter, Group) %>%
  summarise(
    n = sum(!is.na(Value)),
    shapiro_p = if (n >= 3) shapiro.test(Value)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    Normality = case_when(
      is.na(shapiro_p) ~ "not tested",
      shapiro_p > 0.05 ~ "approx normal",
      TRUE ~ "non-normal"
    )
  )

print(head(normality_by_group, 20))

test_residual_normality <- function(dat) {
  dat <- dat %>% filter(!is.na(Value), !is.na(Group))
  if (nrow(dat) < 4 || length(unique(dat$Group)) < 2) {
    return(tibble(n = nrow(dat), shapiro_p_resid = NA_real_))
  }
  
  fit <- lm(Value ~ Group, data = dat)
  res <- residuals(fit)
  
  if (length(res) < 3) {
    return(tibble(n = nrow(dat), shapiro_p_resid = NA_real_))
  }
  
  tibble(
    n = nrow(dat),
    shapiro_p_resid = shapiro.test(res)$p.value
  )
}

normality_residuals <- params_long %>%
  group_by(Parameter) %>%
  group_modify(~ test_residual_normality(.x)) %>%
  ungroup() %>%
  mutate(
    Residual_normality = case_when(
      is.na(shapiro_p_resid) ~ "not tested",
      shapiro_p_resid > 0.05 ~ "approx normal",
      TRUE ~ "non-normal"
    )
  )

print(head(normality_residuals))


#suggest test

test_suggestions <- normality_residuals %>%
  mutate(
    Suggested_test = case_when(
      is.na(shapiro_p_resid) ~ "inspect manually / likely nonparametric",
      shapiro_p_resid > 0.05 ~ "ANOVA or t-test possible",
      TRUE ~ "Kruskal-Wallis recommended"
    )
  )

print(head(test_suggestions))


# QQ plots 
# ggplot(params_long %>% filter(Parameter == "OCR_basal"),
#        aes(sample = Value)) +
#   stat_qq() +
#   stat_qq_line() +
#   facet_wrap(~ Group, scales = "free") +
#   theme_bw()


# outputs

cat("\n--- Flagged timepoints ---\n")
print(head(flagged_points))

cat("\n--- Groomed phase means ---\n")
print(head(phase_means))

cat("\n--- Calculated parameters ---\n")
print(head(params))

cat("\n--- Long parameter table ---\n")
print(head(params_long))

cat("\n--- Values normalized to CTRL1 and pooled controls ---\n")
print(head(params_norm))

cat("\n--- Normality by group ---\n")
print(head(normality_by_group, 20))

cat("\n--- Residual normality ---\n")
print(head(normality_residuals))

cat("\n--- Suggested tests ---\n")
print(head(test_suggestions))