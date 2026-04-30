# Seahorse data analysis

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(tibble)

infile_inhib <- "your data location for data in supplemental table 7"

# ----------------------------
# OPTIONS
# ----------------------------
inhibitors_keep <- c("UK5099", "BPTES")
samplegroups_keep <- c("CTRL1","CTRL2","CTRL3","ED1","ED2","SDS1","SDS2")

drop_negative_raw     <- FALSE
drop_negative_derived <- TRUE
metrics_to_plot <- c("OCR")
show_points <- TRUE
p_adjust_method <- "BH"
p_label_size <- 5.4

# ----------------------------
# WINDOWS FOR INHIBITOR RUNS
# ----------------------------
windows_inhib <- tibble::tribble(
  ~Window,             ~m_min, ~m_max,
  "Basal",                 1,      3,
  "After_inhibitor",       4,      8,
  "After_OM",              9,     11,
  "Maximal",              12,     14,
  "Non_mito",             15,     17
)

# ----------------------------
# PANEL-SPECIFIC X-AXIS ORDER
# CTRL1 can appear in both panels because
# panel membership is assigned per run
# ----------------------------
group_treat_levels_ED <- c(
  "CTRL1", "CTRL1 + inh",
  "ED1",   "ED1 + inh",
  "ED2",   "ED2 + inh"
)

group_treat_levels_SDS <- c(
  "CTRL1", "CTRL1 + inh",
  "SDS1",  "SDS1 + inh",
  "SDS2",  "SDS2 + inh"
)

# ----------------------------
# COLORS
# ----------------------------
cols_group_treat <- c(
  "CTRL1"      = "forestgreen",
  "CTRL1 + inh"= "darkgreen",
  "ED1"        = "orange2",
  "ED1 + inh"  = "darkorange3",
  "ED2"        = "orange2",
  "ED2 + inh"  = "darkorange3",
  "SDS1"       = "plum3",
  "SDS1 + inh" = "plum4",
  "SDS2"       = "plum3",
  "SDS2 + inh" = "plum4"
)

# ----------------------------
# HELPERS
# ----------------------------
as_num_comma <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s+", "")
  x <- str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

safe_mean <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.nan(m)) NA_real_ else m
}

safe_median <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_two_group_p <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(tibble(p=NA_real_, test_used="too_few"))
  
  out <- tryCatch({
    tt <- t.test(x, y)
    tibble(p = as.numeric(tt$p.value), test_used="welch_t")
  }, error=function(e) NULL)
  
  if (!is.null(out) && is.finite(out$p)) return(out)
  
  out2 <- tryCatch({
    wt <- wilcox.test(x, y, exact = FALSE)
    tibble(p = as.numeric(wt$p.value), test_used="wilcox")
  }, error=function(e) NULL)
  
  if (!is.null(out2) && is.finite(out2$p)) return(out2)
  
  tibble(p = NA_real_, test_used = "failed")
}

# ----------------------------
# READ + PARSE
# ----------------------------
dfi <- readxl::read_excel(infile_inhib) %>%
  mutate(
    Measurement = as.integer(Measurement),
    Well = as.character(Well),
    Group = as.character(Group),
    Replicate = as.character(Replicate),
    File = as.character(File),
    Inhibitor = as.character(Inhibitor),
    Group_raw = as.character(Group_raw),
    sample_id = as.character(sample_id),
    Treatment = as.character(Treatment),
    SampleGroup = as.character(SampleGroup),
    Time = as_num_comma(Time),
    OCR  = as_num_comma(OCR),
    RunID = paste(File, Replicate, sep = "__")
  ) %>%
  filter(SampleGroup %in% samplegroups_keep) %>%
  filter(Inhibitor %in% inhibitors_keep) %>%
  filter(Treatment %in% c("Medium", "Inhibitor"))

if (drop_negative_raw) {
  dfi <- dfi %>%
    mutate(OCR = ifelse(!is.na(OCR) & OCR < 0, NA_real_, OCR))
}

# ----------------------------
# LONG OCR + PER-WELL WINDOW MEANS
# ----------------------------
long_metric_inhib <- dfi %>%
  transmute(
    Inhibitor, File, Replicate, RunID,
    SampleGroup, Treatment, Well, Measurement,
    Metric = "OCR", Value = OCR
  )

per_well_window_inhib <- long_metric_inhib %>%
  crossing(windows_inhib) %>%
  filter(Measurement >= m_min, Measurement <= m_max) %>%
  group_by(Inhibitor, File, Replicate, RunID, SampleGroup, Treatment, Well, Metric, Window) %>%
  summarise(
    WindowMean = safe_mean(Value),
    n_points   = sum(!is.na(Value)),
    .groups = "drop"
  )

wide_inhib <- per_well_window_inhib %>%
  select(Inhibitor, File, Replicate, RunID, SampleGroup, Treatment, Well, Metric, Window, WindowMean) %>%
  pivot_wider(names_from = Window, values_from = WindowMean)

# ----------------------------
# DERIVED OCR PARAMETERS
# ----------------------------
derived_inhib <- wide_inhib %>%
  filter(Metric == "OCR") %>%
  mutate(
    ATP_linked          = Basal - After_OM,
    Proton_leak         = After_OM - Non_mito,
    Reserve_capacity    = Maximal - Basal,
    Max_resp            = Maximal - Non_mito,
    Coupling_efficiency = ifelse(
      !is.na(Basal) & Basal != 0,
      (Basal - Proton_leak) / Basal,
      NA_real_
    )
  ) %>%
  pivot_longer(
    cols = c(
      "Basal","After_inhibitor","After_OM","Maximal","Non_mito",
      "ATP_linked","Proton_leak","Reserve_capacity","Max_resp",
      "Coupling_efficiency"
    ),
    names_to = "Parameter",
    values_to = "Value"
  )

if (drop_negative_derived) {
  derived_inhib <- derived_inhib %>%
    mutate(
      Value = ifelse(
        Parameter != "Coupling_efficiency" & !is.na(Value) & Value < 0,
        NA_real_,
        Value
      )
    )
}

# ----------------------------
# ASSIGN PANEL PER RUN
# If a run contains ED samples -> CTRL_ED
# If a run contains SDS samples -> CTRL_SDS
# ----------------------------
run_panel_map <- derived_inhib %>%
  distinct(Inhibitor, File, Replicate, RunID, SampleGroup) %>%
  group_by(Inhibitor, File, Replicate, RunID) %>%
  summarise(
    has_ED  = any(SampleGroup %in% c("ED1", "ED2")),
    has_SDS = any(SampleGroup %in% c("SDS1", "SDS2")),
    .groups = "drop"
  ) %>%
  mutate(
    Panel = case_when(
      has_ED  & !has_SDS ~ "CTRL_ED",
      has_SDS & !has_ED  ~ "CTRL_SDS",
      has_ED  & has_SDS  ~ "MIXED",
      TRUE ~ NA_character_
    )
  )

# ----------------------------
# NORMALIZE: Medium within same SampleGroup/run = 100
# ----------------------------
medium_ref_run <- derived_inhib %>%
  filter(Treatment == "Medium") %>%
  group_by(Inhibitor, File, Replicate, RunID, SampleGroup, Parameter) %>%
  summarise(
    medium_median_run = safe_median(Value),
    medium_n_run      = sum(!is.na(Value)),
    .groups = "drop"
  )

di <- derived_inhib %>%
  left_join(
    medium_ref_run,
    by = c("Inhibitor","File","Replicate","RunID","SampleGroup","Parameter")
  ) %>%
  left_join(
    run_panel_map %>% select(Inhibitor, File, Replicate, RunID, Panel),
    by = c("Inhibitor","File","Replicate","RunID")
  ) %>%
  mutate(
    Value_pct = ifelse(
      !is.na(Value) & !is.na(medium_median_run) & medium_median_run != 0,
      100 * (Value / medium_median_run),
      NA_real_
    ),
    Value_pct_use = ifelse(!is.na(Value_pct) & Value_pct > 500, NA_real_, Value_pct),
    Treatment_plot = factor(Treatment, levels = c("Medium","Inhibitor")),
    Group_treat = case_when(
      Treatment == "Medium"    ~ SampleGroup,
      Treatment == "Inhibitor" ~ paste0(SampleGroup, " + inh"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Panel %in% c("CTRL_ED", "CTRL_SDS"))

# ----------------------------
# PANEL-SPECIFIC FACTOR ORDER
# ----------------------------
di <- di %>%
  mutate(
    Group_treat = case_when(
      Panel == "CTRL_ED"  ~ as.character(factor(Group_treat, levels = group_treat_levels_ED)),
      Panel == "CTRL_SDS" ~ as.character(factor(Group_treat, levels = group_treat_levels_SDS)),
      TRUE ~ Group_treat
    ),
    Group_treat = factor(Group_treat)
  )

# ----------------------------
# SUMMARY STATS
# ----------------------------
summary_bar_inhib <- di %>%
  group_by(Inhibitor, Panel, Metric, Parameter, SampleGroup, Group_treat) %>%
  summarise(
    mean = mean(Value_pct_use, na.rm = TRUE),
    sd   = sd(Value_pct_use, na.rm = TRUE),
    n_wells = sum(!is.na(Value_pct_use)),
    .groups = "drop"
  ) %>%
  mutate(mean = ifelse(is.nan(mean), NA_real_, mean))

# ----------------------------
# STATS: Medium vs Inhibitor within each SampleGroup
# ----------------------------
pvals_inhib <- di %>%
  filter(Metric %in% metrics_to_plot) %>%
  filter(!is.na(Value_pct_use)) %>%
  group_by(Inhibitor, Panel, Metric, Parameter, SampleGroup) %>%
  group_modify(function(dat, keys) {
    
    x_medium <- dat %>% filter(Treatment_plot == "Medium") %>% pull(Value_pct_use)
    x_inhib  <- dat %>% filter(Treatment_plot == "Inhibitor") %>% pull(Value_pct_use)
    
    if (length(x_medium) < 2 || length(x_inhib) < 2) {
      return(tibble(
        n_Medium = length(x_medium),
        n_Inhibitor = length(x_inhib),
        p = NA_real_,
        test_used = "too_few"
      ))
    }
    
    res <- safe_two_group_p(x_medium, x_inhib)
    
    tibble(
      n_Medium = length(x_medium),
      n_Inhibitor = length(x_inhib),
      p = res$p,
      test_used = res$test_used
    )
  }) %>%
  ungroup() %>%
  group_by(Inhibitor, Panel) %>%
  mutate(
    p_adj = p.adjust(p, method = p_adjust_method),
    label = ifelse(is.na(p_adj), "p=NA", paste0("p=", signif(p_adj, 2)))
  ) %>%
  ungroup()

# ----------------------------
# PLOT FUNCTION
# ----------------------------
plot_one_inhib_panel <- function(inhib, panel_name, met, param) {
  
  df_bar <- summary_bar_inhib %>%
    filter(Inhibitor == inhib, Panel == panel_name, Metric == met, Parameter == param) %>%
    filter(n_wells > 0)
  
  if (nrow(df_bar) == 0) return(invisible(NULL))
  
  if (panel_name == "CTRL_ED") {
    level_order <- group_treat_levels_ED
  } else if (panel_name == "CTRL_SDS") {
    level_order <- group_treat_levels_SDS
  } else {
    level_order <- unique(as.character(df_bar$Group_treat))
  }
  
  df_bar <- df_bar %>%
    mutate(Group_treat = factor(as.character(Group_treat), levels = level_order))
  
  df_pts <- di %>%
    filter(Inhibitor == inhib, Panel == panel_name, Metric == met, Parameter == param) %>%
    filter(!is.na(Value_pct_use)) %>%
    mutate(Group_treat = factor(as.character(Group_treat), levels = level_order))
  
  x_labels <- df_bar %>%
    mutate(label = paste0(as.character(Group_treat), "\n(n=", n_wells, ")")) %>%
    select(Group_treat, label) %>%
    distinct()
  
  df_lab <- pvals_inhib %>%
    filter(Inhibitor == inhib, Panel == panel_name, Metric == met, Parameter == param) %>%
    mutate(
      Group_treat = factor(paste0(as.character(SampleGroup), " + inh"), levels = level_order)
    ) %>%
    left_join(
      df_bar %>% select(Group_treat, mean, sd),
      by = "Group_treat"
    ) %>%
    mutate(
      y = (mean + ifelse(is.na(sd), 0, sd)) * 1.10
    )
  
  ymax_bar <- suppressWarnings(max(df_bar$mean + ifelse(is.na(df_bar$sd), 0, df_bar$sd), na.rm = TRUE))
  ymax_lab <- if (nrow(df_lab) > 0) suppressWarnings(max(df_lab$y, na.rm = TRUE)) else NA_real_
  ymax_all <- max(c(ymax_bar, ymax_lab, 120), na.rm = TRUE)
  
  p <- ggplot(df_bar, aes(x = Group_treat, y = mean, fill = Group_treat)) +
    geom_col(width = 0.75) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    scale_fill_manual(values = cols_group_treat, drop = FALSE) +
    scale_x_discrete(labels = setNames(x_labels$label, x_labels$Group_treat)) +
    geom_hline(yintercept = 100, linetype = "dashed") +
    coord_cartesian(ylim = c(0, ymax_all * 1.10)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(size = 12)
    ) +
    labs(
      title = paste(inhib, "|", panel_name, "|", met, "|", param, "| Medium within run = 100"),
      x = NULL,
      y = "% of Medium (within run)"
    )
  
  if (show_points) {
    p <- p +
      geom_point(
        data = df_pts,
        aes(x = Group_treat, y = Value_pct_use, fill = Group_treat),
        shape = 21,
        color = "black",
        stroke = 0.4,
        size = 2.2,
        alpha = 0.85,
        position = position_jitter(width = 0.12, height = 0),
        inherit.aes = FALSE
      )
  }
  
  if (nrow(df_lab) > 0) {
    p <- p +
      geom_text(
        data = df_lab,
        aes(x = Group_treat, y = y, label = label),
        inherit.aes = FALSE,
        vjust = 0,
        size = p_label_size
      )
  }
  
  print(p)
}

# ----------------------------
# RUN PLOTS
# ----------------------------
inhibs_here <- sort(unique(di$Inhibitor))
panels_here <- c("CTRL_ED", "CTRL_SDS")
params_here <- sort(unique(di$Parameter))
mets_here   <- intersect(metrics_to_plot, unique(di$Metric))

for (inhib in inhibs_here) {
  for (panel_name in panels_here) {
    for (met in mets_here) {
      for (param in params_here) {
        plot_one_inhib_panel(inhib, panel_name, met, param)
      }
    }
  }
}

# ----------------------------
# QC
# ----------------------------
cat("\n--- QC: run panel map ---\n")
print(run_panel_map %>% arrange(Inhibitor, File, Replicate), n = 200)

cat("\n--- QC: non-NA counts after normalization + >500% cutoff ---\n")
print(
  di %>%
    group_by(Inhibitor, Panel, Parameter, Group_treat) %>%
    summarise(n = sum(!is.na(Value_pct_use)), .groups = "drop") %>%
    arrange(Inhibitor, Panel, Parameter, Group_treat),
  n = 400
)

cat("\n--- QC: stats table (Medium vs Inhibitor within SampleGroup) ---\n")
print(
  pvals_inhib %>%
    arrange(Inhibitor, Panel, Metric, Parameter, SampleGroup),
  n = 400
)

cat("\n--- QC: Medium references per run (median) ---\n")
print(
  medium_ref_run %>%
    arrange(Inhibitor, File, Replicate, SampleGroup, Parameter),
  n = 400
)
