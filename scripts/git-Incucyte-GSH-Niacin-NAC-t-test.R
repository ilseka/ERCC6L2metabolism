### incucyte supplement plots and t-tests

library(tidyverse)


base_dir  <- "your base directory"
conf_file <- "file.path(base_dir, confluence data from suppl. table 3)"
cnt_file  <- "file.path(base_dir, confluence data from suppl. table 3)"

conf <- read.delim(conf_file, check.names = FALSE, stringsAsFactors = FALSE)
cnt  <- read.delim(cnt_file,  check.names = FALSE, stringsAsFactors = FALSE)


conf <- conf %>% mutate(Timepoint = as.character(Timepoint), Set = as.character(Set),
                        Condition = as.character(Condition), Media = as.character(Media))
cnt  <- cnt  %>% mutate(Timepoint = as.character(Timepoint), Set = as.character(Set),
                        Condition = as.character(Condition), Media = as.character(Media))

conditions_keep <- c("Control", "ERCC6L2", "SDS")
supps <- c("Glutathione", "NAC", "Niacin")

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

# label every 12h, but keep all timepoints present
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


#Summarize 
summarize_metric <- function(df, media_keep, conditions_keep, value_col, metric_label) {
  df %>%
    filter(Condition %in% conditions_keep, Media %in% media_keep) %>%
    group_by(Condition, Media, Timepoint) %>%
    summarise(mean_ci(.data[[value_col]]), .groups = "drop") %>%
    mutate(
      Metric = metric_label,
      mean = 100 * mean,
      ymin = 100 * ymin,
      ymax = 100 * ymax
    )
}

# Plot: standard glucose = solid LINE only; supplement = dotted LINE + POINTS + CI 
plot_supplement_all_conditions_faceted <- function(supplement, save = TRUE) {
  
  media_keep <- c("Standard_glucose", supplement)
  
  conf_sum <- summarize_metric(conf_for_supp, media_keep, conditions_keep, "Norm.Confluence", "Confluence")
  cnt_sum  <- summarize_metric(cnt_for_supp,  media_keep, conditions_keep, "Norm.Counts",     "Nuclei count")
  
  df_plot <- bind_rows(conf_sum, cnt_sum) %>%
    mutate(
      Condition = factor(Condition, levels = conditions_keep),
      Media     = factor(Media, levels = c("Standard_glucose", supplement))
    )
  
  tp_levels <- order_time(df_plot$Timepoint)
  df_plot <- df_plot %>% mutate(Timepoint = factor(Timepoint, levels = tp_levels))
  breaks12 <- x_breaks_12h(tp_levels)
  
  df_std  <- df_plot %>% filter(Media == "Standard_glucose")
  df_supp <- df_plot %>% filter(Media == supplement)
  
  p <- ggplot() +
    # CI ONLY for supplement
    geom_ribbon(
      data = df_supp,
      aes(x = Timepoint, ymin = ymin, ymax = ymax, fill = Condition,
          group = interaction(Condition, Metric)),
      alpha = 0.12, color = NA, show.legend = FALSE
    ) +
    # Standard glucose = line only
    geom_line(
      data = df_std,
      aes(x = Timepoint, y = mean, color = Condition,
          group = interaction(Condition, Metric)),
      linewidth = 1.4, linetype = "solid"
    ) +
    # Supplement = dotted line + points
    geom_line(
      data = df_supp,
      aes(x = Timepoint, y = mean, color = Condition,
          group = interaction(Condition, Metric)),
      linewidth = 1.4, linetype = "dotted"
    ) +
    geom_point(
      data = df_supp,
      aes(x = Timepoint, y = mean, color = Condition,
          group = interaction(Condition, Metric)),
      size = 3.2
    ) +
    facet_wrap(~ Metric, ncol = 1, scales = "free_y") +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    scale_x_discrete(breaks = breaks12) +
    labs(
      title = paste0("Standard glucose vs ", supplement),
      x = "Time (h)",
      y = "Percent of baseline (0 h = 100%)",
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
  
  
  p
}

#T-tests at 72h: supplement vs standard, per condition, for confluence and counts 
run_welch_test <- function(x, y, label = NULL) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(
      comparison = label, n_x = length(x), n_y = length(y),
      t_statistic = NA_real_, df = NA_real_, p_value = NA_real_,
      conf_low = NA_real_, conf_high = NA_real_,
      mean_x = mean(x), mean_y = mean(y)
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

make_supp_ttest_table <- function(timepoint = "72") {
  
  out <- list()
  
  for (supp in supps) {
    for (cond in conditions_keep) {
      
      # Confluence uses conf_for_supp (Set 13 removed for supplement)
      g1c <- conf_for_supp %>% filter(Condition == cond, Media == supp,            Timepoint == timepoint)
      g2c <- conf_for_supp %>% filter(Condition == cond, Media == "Standard_glucose", Timepoint == timepoint)
      
      out[[length(out)+1]] <- run_welch_test(
        g1c$Norm.Confluence, g2c$Norm.Confluence,
        label = paste(cond, supp, "vs Standard_glucose", paste0("t=", timepoint, "h"), "- Confluence")
      ) %>% mutate(Condition = cond, Supplement = supp, Metric = "Confluence", Timepoint = timepoint)
      
      # Counts uses cnt_for_supp (Set 13 INCLUDED)
      g1n <- cnt_for_supp %>% filter(Condition == cond, Media == supp,            Timepoint == timepoint)
      g2n <- cnt_for_supp %>% filter(Condition == cond, Media == "Standard_glucose", Timepoint == timepoint)
      
      out[[length(out)+1]] <- run_welch_test(
        g1n$Norm.Counts, g2n$Norm.Counts,
        label = paste(cond, supp, "vs Standard_glucose", paste0("t=", timepoint, "h"), "- Counts")
      ) %>% mutate(Condition = cond, Supplement = supp, Metric = "Counts", Timepoint = timepoint)
    }
  }
  
  bind_rows(out)
}

# plots + table
plots <- setNames(lapply(supps, function(s) {
  p <- plot_supplement_all_conditions_faceted(s, save = TRUE)
  print(p)
  p
}), supps)

ttest_table <- make_supp_ttest_table(timepoint = "72")


ttest_table
