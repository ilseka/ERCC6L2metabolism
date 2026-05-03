

#### standard vs low glucose plots and t-test

library(tidyverse)
library(stringr)

base_dir  <- "your base directory"
conf_file <- "file.path(base_dir, confluence data from suppl. table 3)"
cnt_file  <- "file.path(base_dir, confluence data from suppl. table 3)"

conf <- read.delim(conf_file, check.names = FALSE, stringsAsFactors = FALSE)
cnt  <- read.delim(cnt_file,  check.names = FALSE, stringsAsFactors = FALSE)


conditions_keep <- c("Control", "ERCC6L2", "SDS")

media_A <- "Low_glucose"   # solid line
media_B <- "Low_gluc_gln"  # dotted line (+ points + CI)

end_tp <- 48
ttest_tp <- "48"
allowed_tps <- as.character(seq(0, end_tp, by = 6))

pal <- c(
  "Control" = "darkolivegreen",
  "ERCC6L2" = "orange2",
  "SDS"     = "plum3"
)

# helper: numeric time ordering
order_time <- function(tp_chr) {
  tp_num <- suppressWarnings(as.numeric(tp_chr))
  if (all(is.na(tp_num))) sort(unique(tp_chr)) else as.character(sort(unique(tp_num)))
}

# label every 12h
x_breaks_12h <- function(tp_levels) {
  tp_num <- suppressWarnings(as.numeric(tp_levels))
  if (all(is.na(tp_num))) tp_levels else tp_levels[tp_num %% 12 == 0]
}

mean_ci <- function(x) {
  m  <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  ci <- 1.96 * se
  tibble(mean = m, ymin = m - ci, ymax = m + ci)
}


conf <- conf %>%
  mutate(
    Timepoint  = str_squish(as.character(Timepoint)),
    Set        = as.character(Set),
    Condition  = as.character(Condition),
    Media      = str_squish(as.character(Media)),
    Timepoint_num = suppressWarnings(as.numeric(Timepoint))
  )

cnt <- cnt %>%
  mutate(
    Timepoint  = str_squish(as.character(Timepoint)),
    Set        = as.character(Set),
    Condition  = as.character(Condition),
    Media      = str_squish(as.character(Media)),
    Timepoint_num = suppressWarnings(as.numeric(Timepoint))
  )

conf_use <- conf
cnt_use  <- cnt

diag_cnt <- cnt_use %>%
  filter(Condition %in% conditions_keep, Media %in% c(media_A, media_B)) %>%
  distinct(Condition, Media, Set) %>%
  count(Condition, Media, name = "n_sets")

diag_conf <- conf_use %>%
  filter(Condition %in% conditions_keep, Media %in% c(media_A, media_B)) %>%
  distinct(Condition, Media, Set) %>%
  count(Condition, Media, name = "n_sets")

message("Counts: number of unique sets per Condition x Media")
print(diag_cnt)
message("Confluence: number of unique sets per Condition x Media")
print(diag_conf)

summarize_metric <- function(df, media_keep, conditions_keep, value_col, metric_label) {
  
  df %>%
    filter(
      Condition %in% conditions_keep,
      Media %in% media_keep,
      !is.na(Timepoint_num),
      Timepoint %in% allowed_tps   # <-- hard guard (fixes 'after 48h' artifacts)
    ) %>%
    group_by(Condition, Media, Timepoint) %>%
    summarise(mean_ci(.data[[value_col]]), .groups = "drop") %>%
    mutate(Metric = metric_label)
}

# Plot 
plot_low_glucose_vs_gln_all_conditions_faceted <- function(save = TRUE) {
  
  media_keep <- c(media_A, media_B)
  
  conf_sum <- summarize_metric(conf_use, media_keep, conditions_keep, "Norm.Confluence", "Confluence")
  cnt_sum  <- summarize_metric(cnt_use,  media_keep, conditions_keep, "Norm.Counts",     "Nuclei count")
  
  df_plot <- bind_rows(conf_sum, cnt_sum) %>%
    mutate(
      Condition = factor(Condition, levels = conditions_keep),
      Media     = factor(Media, levels = c(media_A, media_B)),
      MediaType = if_else(Media == media_A, "Low glucose", "Low glucose + glutamine")
    )
  
  # enforce x-axis order ONLY over allowed timepoints
  tp_levels <- order_time(allowed_tps)
  df_plot <- df_plot %>% mutate(Timepoint = factor(Timepoint, levels = tp_levels))
  breaks12 <- x_breaks_12h(tp_levels)
  
  df_A <- df_plot %>% filter(Media == media_A)
  df_B <- df_plot %>% filter(Media == media_B)
  
  p <- ggplot() +
    geom_ribbon(
      data = df_B,
      aes(x = Timepoint, ymin = ymin, ymax = ymax, fill = Condition,
          group = interaction(Condition, Metric)),
      alpha = 0.12, color = NA, show.legend = FALSE
    ) +
    geom_line(
      data = df_A,
      aes(x = Timepoint, y = mean, color = Condition,
          linetype = MediaType,
          group = interaction(Condition, Metric)),
      linewidth = 1.4
    ) +
    geom_line(
      data = df_B,
      aes(x = Timepoint, y = mean, color = Condition,
          linetype = MediaType,
          group = interaction(Condition, Metric)),
      linewidth = 1.4
    ) +
    geom_point(
      data = df_B,
      aes(x = Timepoint, y = mean, color = Condition,
          group = interaction(Condition, Metric)),
      size = 3.2,
      show.legend = FALSE
    ) +
    facet_wrap(~ Metric, ncol = 2, scales = "free_y") +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    scale_linetype_manual(
      values = c("Low glucose" = "solid", "Low glucose + glutamine" = "dotted"),
      name = "Media"
    ) +
    scale_x_discrete(breaks = breaks12, drop = FALSE) +
    labs(
      title = paste0(media_A, " vs ", media_B, " (0–", end_tp, "h)"),
      x = "Time (h)",
      y = "Normalized to baseline (0 h = 1)",
      color = "Condition"
    ) +
    theme_bw(base_size = 22) +
    theme(
      plot.title   = element_text(face = "bold", size = 24),
      strip.text   = element_text(face = "bold", size = 22),
      legend.title = element_text(size = 20),
      legend.text  = element_text(size = 20),
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 22),
      axis.text.x  = element_text(size = 20, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 20)
    )
  
  if (save) {
    out_name <- paste0("Incucyte_AllConditions_", media_A, "_vs_", media_B, "_0to", end_tp, "h_FACET_Metric.png")
    ggsave(file.path(base_dir, out_name), p, width = 14, height = 7, dpi = 300)
    message("Saved: ", file.path(base_dir, out_name))
  }
  
  p
}

#Welch t-tests at 48h 
run_welch_test <- function(x, y, label = NULL) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(
      comparison = label, n_x = length(x), n_y = length(y),
      t_statistic = NA_real_, df = NA_real_, p_value = NA_real_,
      conf_low = NA_real_, conf_high = NA_real_,
      mean_x = ifelse(length(x) == 0, NA_real_, mean(x)),
      mean_y = ifelse(length(y) == 0, NA_real_, mean(y))
    ))
  }
  t_res <- t.test(x, y)
  tibble(
    comparison  = label,
    n_x         = length(x),
    n_y         = length(y),
    t_statistic = unname(t_res$statistic),
    df          = unname(t_res$parameter),
    p_value     = t_res$p.value,
    conf_low    = t_res$conf.int[1],
    conf_high   = t_res$conf.int[2],
    mean_x      = unname(t_res$estimate[1]),
    mean_y      = unname(t_res$estimate[2])
  )
}

make_low_gln_ttest_table <- function(timepoint = "48") {
  
  out <- list()
  
  for (cond in conditions_keep) {
    
    gA_c <- conf_use %>% filter(Condition == cond, Media == media_A, Timepoint == timepoint)
    gB_c <- conf_use %>% filter(Condition == cond, Media == media_B, Timepoint == timepoint)
    
    out[[length(out)+1]] <- run_welch_test(
      gA_c$Norm.Confluence, gB_c$Norm.Confluence,
      label = paste(cond, media_A, "vs", media_B, paste0("t=", timepoint, "h"), "- Confluence")
    ) %>% mutate(Condition = cond, Media_A = media_A, Media_B = media_B, Metric = "Confluence", Timepoint = timepoint)
    
    gA_n <- cnt_use %>% filter(Condition == cond, Media == media_A, Timepoint == timepoint)
    gB_n <- cnt_use %>% filter(Condition == cond, Media == media_B, Timepoint == timepoint)
    
    out[[length(out)+1]] <- run_welch_test(
      gA_n$Norm.Counts, gB_n$Norm.Counts,
      label = paste(cond, media_A, "vs", media_B, paste0("t=", timepoint, "h"), "- Counts")
    ) %>% mutate(Condition = cond, Media_A = media_A, Media_B = media_B, Metric = "Counts", Timepoint = timepoint)
  }
  
  bind_rows(out)
}

p <- plot_low_glucose_vs_gln_all_conditions_faceted(save = TRUE)
print(p)

ttest_low_gln <- make_low_gln_ttest_table(timepoint = ttest_tp)

ttest_low_gln
