# Flux analysis + supplemental figure plots (supplemental table 6 data)


# GLUCOSE TRACE

#  all isotopologues
# normalized to Glucose M+6


library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)


file_path <- "your file path to glucose data in supplemental table 6"

raw <- read_excel(file_path, col_names = FALSE)

header1 <- as.character(unlist(raw[1, ]))
header2 <- as.character(unlist(raw[2, ]))

header1[1] <- "Metabolite"
header2[1] <- "Metabolite"

for (i in seq_along(header1)) {
  if (i > 1 && (is.na(header1[i]) || header1[i] == "")) {
    header1[i] <- header1[i - 1]
  }
}

new_names <- c("Metabolite", paste(header1[-1], header2[-1], sep = "__"))

data_only <- raw[-c(1, 2), ]
col_is_empty <- sapply(data_only, function(x) all(is.na(x) | x == ""))

keep_cols <- c(TRUE, !col_is_empty[-1])

raw2 <- raw[, keep_cols]
new_names2 <- new_names[keep_cols]
colnames(raw2) <- new_names2

df <- raw2[-c(1, 2), ]


# long format

df_long <- df %>%
  pivot_longer(
    cols = -Metabolite,
    names_to = "ConditionSample",
    values_to = "Abundance"
  ) %>%
  separate(
    ConditionSample,
    into = c("ConditionTime", "SampleID"),
    sep = "__",
    remove = TRUE
  ) %>%
  mutate(
    Metabolite = as.character(Metabolite),
    Metabolite = str_replace_all(Metabolite, "\u00A0", " "),
    Metabolite = str_squish(Metabolite),
    
    ConditionTime = str_squish(as.character(ConditionTime)),
    SampleID = str_squish(as.character(SampleID)),
    
    Abundance = as.character(Abundance),
    Abundance = str_replace_all(Abundance, ",", "."),
    Abundance = as.numeric(Abundance)
  ) %>%
  filter(!is.na(Abundance), !is.na(Metabolite), Metabolite != "")


# 4. parse metabolite + isotopologue 

df_long <- df_long %>%
  mutate(
    Metabolite_clean = str_squish(Metabolite),
    
    BaseMetabolite = case_when(
      str_detect(Metabolite_clean, "\\s*\\+\\d+$") ~ str_remove(Metabolite_clean, "\\s*\\+\\d+$"),
      TRUE ~ Metabolite_clean
    ),
    
    iso_num = case_when(
      str_detect(Metabolite_clean, "\\s*\\+\\d+$") ~ as.numeric(str_extract(Metabolite_clean, "\\d+$")),
      TRUE ~ 0
    ),
    
    Isotopologue = paste0("M+", iso_num),
    BaseMetabolite = str_squish(BaseMetabolite)
  )


df_long <- df_long %>%
  mutate(
    BaseMetabolite = recode(
      BaseMetabolite,
      "aKG" = "a-KG",
      "alpha-KG" = "a-KG"
    ),
    BaseMetabolite = str_squish(BaseMetabolite)
  )

# Optional check:
# print(sort(unique(df_long$BaseMetabolite)))


# order groups

condition_order <- c(
  "Control 30min",
  "ERCC6L2 30min",
  "SDS 30min",
  "Control 120min",
  "ERCC6L2 120min",
  "SDS 120min"
)

condition_labels <- c(
  "Control 30min" = "Control 30 min",
  "ERCC6L2 30min" = "ED 30 min",
  "SDS 30min" = "SDS 30 min",
  "Control 120min" = "Control 120 min",
  "ERCC6L2 120min" = "ED 120 min",
  "SDS 120min" = "SDS 120 min"
)

group_order <- c(
  "Control 30 min",
  "ED 30 min",
  "SDS 30 min",
  "Control 120 min",
  "ED 120 min",
  "SDS 120 min"
)

df_long <- df_long %>%
  mutate(
    ConditionTime = factor(ConditionTime, levels = condition_order),
    GroupLabel = recode(as.character(ConditionTime), !!!condition_labels),
    GroupLabel = factor(GroupLabel, levels = group_order)
  )


# isotopologue order

iso_levels_all <- df_long %>%
  distinct(Isotopologue, iso_num) %>%
  arrange(iso_num) %>%
  pull(Isotopologue)

iso_levels_noM0 <- df_long %>%
  filter(Isotopologue != "M+0") %>%
  distinct(Isotopologue, iso_num) %>%
  arrange(iso_num) %>%
  pull(Isotopologue)

df_long$Isotopologue <- factor(df_long$Isotopologue, levels = iso_levels_all)


# colors

palette_all <- c(
  "M+0" = "#4E79A7",
  "M+1" = "#A0CBE8",
  "M+2" = "#59A14F",
  "M+3" = "#8CD17D",
  "M+4" = "#F28E2B",
  "M+5" = "#FFBE7D",
  "M+6" = "#B07AA1"
)

palette_noM0 <- c(
  "M+1" = "#A0CBE8",
  "M+2" = "#59A14F",
  "M+3" = "#8CD17D",
  "M+4" = "#F28E2B",
  "M+5" = "#FFBE7D",
  "M+6" = "#B07AA1"
)

palette_all <- palette_all[names(palette_all) %in% iso_levels_all]
palette_noM0 <- palette_noM0[names(palette_noM0) %in% iso_levels_noM0]


# normalize to glucose M+6

glucose_m6 <- df_long %>%
  filter(BaseMetabolite == "Glucose", Isotopologue == "M+6") %>%
  select(ConditionTime, GroupLabel, SampleID, Glucose_M6 = Abundance)

df_glc_norm <- df_long %>%
  filter(Isotopologue != "M+0") %>%
  left_join(glucose_m6, by = c("ConditionTime", "GroupLabel", "SampleID")) %>%
  mutate(
    NormToGlucoseM6 = Abundance / Glucose_M6
  ) %>%
  filter(!is.na(NormToGlucoseM6), is.finite(NormToGlucoseM6))

df_glc_norm$Isotopologue <- factor(df_glc_norm$Isotopologue, levels = iso_levels_noM0)

# summary

plot_df_abundance <- df_long %>%
  group_by(BaseMetabolite, GroupLabel, Isotopologue) %>%
  summarise(
    MeanAbundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

plot_df_glc_norm <- df_glc_norm %>%
  group_by(BaseMetabolite, GroupLabel, Isotopologue) %>%
  summarise(
    MeanNorm = mean(NormToGlucoseM6, na.rm = TRUE),
    .groups = "drop"
  )


# metabolite order

metabolite_order <- c(
  "Glucose", 
  "Pyruvate", 
  "Lactate", 
  "Serine", 
  "Glycine", 
  "Malate", 
  "Aspartate", 
  "Glutamate", 
  "Glutathione"
)

if (length(metabolite_order) == 1) {
  found_mets <- unique(as.character(plot_df_abundance$BaseMetabolite))
  found_mets <- found_mets[!is.na(found_mets)]
  found_mets <- setdiff(found_mets, "Glucose")
  metabolite_order <- c("Glucose", found_mets)
}

plot_df_abundance <- plot_df_abundance %>%
  filter(BaseMetabolite %in% metabolite_order) %>%
  mutate(
    BaseMetabolite = factor(BaseMetabolite, levels = metabolite_order),
    Isotopologue = factor(Isotopologue, levels = iso_levels_all)
  )

plot_df_glc_norm <- plot_df_glc_norm %>%
  filter(BaseMetabolite %in% metabolite_order) %>%
  mutate(
    BaseMetabolite = factor(BaseMetabolite, levels = metabolite_order),
    Isotopologue = factor(Isotopologue, levels = iso_levels_noM0)
  )


# all isotopologues

p_all_isotopes <- ggplot(
  plot_df_abundance,
  aes(x = GroupLabel, y = MeanAbundance, fill = Isotopologue)
) +
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE),
    width = 0.8,
    color = "black",
    linewidth = 0.2
  ) +
  facet_wrap(~ BaseMetabolite, scales = "free_y", ncol = 4, drop = FALSE) +
  scale_fill_manual(values = palette_all, drop = TRUE) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    title = "Glucose tracer: all isotopologues",
    x = NULL,
    y = "Mean abundance",
    fill = "Isotopologue"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_all_isotopes)

# normalized to glucose M+6

p_norm_glc_m6 <- ggplot(
  plot_df_glc_norm,
  aes(x = GroupLabel, y = MeanNorm, fill = Isotopologue)
) +
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE),
    width = 0.8,
    color = "black",
    linewidth = 0.2
  ) +
  facet_wrap(~ BaseMetabolite, scales = "free_y", ncol = 4, drop = FALSE) +
  scale_fill_manual(values = palette_noM0, drop = TRUE) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    title = "Normalized to Glucose M+6",
    x = NULL,
    y = "Abundance / Glucose M+6",
    fill = "Isotopologue"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_norm_glc_m6)