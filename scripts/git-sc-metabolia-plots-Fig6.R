### scRNAseq Metabolism pathways figure 6

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(forcats)


input_dir <- "your directory to scRNAseq/patwhay_lists"

files <- list.files(
  path = input_dir,
  pattern = "enrichR_Reactome_Pathways_2024_.*significant-pathways_updown.*",
  full.names = TRUE
)


clean_celltype <- function(x) {
  case_when(
    str_detect(x, "HSCs_and_MPPs") ~ "HSC & MPP",
    str_detect(x, "HSC") ~ "HSC & MPP",
    str_detect(x, "EMP") ~ "EMP",
    str_detect(x, "EEP") ~ "EEP",
    str_detect(x, "LEP") ~ "LEP",
    str_detect(x, "Late_eryth_prog") ~ "LEP",
    str_detect(x, "Early_eryth_prog") ~ "EEP",
    TRUE ~ x
  )
}

get_comparison_group <- function(x) {
  case_when(
    str_detect(x, "ED_BMF_vs_HD|ED_BM_normal_vs_HD") ~ "ED vs HC",
    str_detect(x, "SDS_BMF_vs_HD|SDS_BM_normal_vs_HD") ~ "SDS vs HC",
    TRUE ~ NA_character_
  )
}

get_celltype_from_filename <- function(x) {
  out <- str_match(basename(x), "vs_HD_(.*?)_CDR")[, 2]
  clean_celltype(out)
}

rename_terms <- function(df) {
  df %>%
    mutate(
      Term = recode(Term,
                    "Citric Acid Cycle (TCA Cycle)" = "TCA Cycle",
                    "Oxidative Stress Induced Senescence" = "Senescence By Oxidative Stress",
                    "Detoxification of Reactive Oxygen Species" = "ROS Detoxification",
                    "Metabolism of Nucleotides" = "Nucleotide Metabolism",
                    "Maturation of TCA Enzymes and Regulation of TCA Cycle" = "TCA Cycle Regulation",
                    "Integration of Energy Metabolism" = "Energy Metabolism",
                    "Cellular Response to Mitochondrial Stress" = "Mitochondrial Stress Response",
                    "Cellular Responses to Mitochondrial Stress" = "Mitochondrial Stress Response",
                    "NFE2L2 Regulating Anti-Oxidant Detoxification Enzymes" = "NFE2L2 anti-oxidant regulation",
                    "Defects in Vitamin and Cofactor Metabolism" = "Vitamin & Cofactor Metabolism"
      )
    )
}

count_genes <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(0)
  vals <- unlist(str_split(x, ";"))
  vals <- str_trim(vals)
  vals <- vals[vals != ""]
  length(unique(vals))
}

infer_direction_majority <- function(genes_up, genes_down, threshold = 0.6) {
  n_up <- count_genes(genes_up)
  n_down <- count_genes(genes_down)
  total <- n_up + n_down
  
  if (total == 0) return(NA_character_)
  
  prop_up <- n_up / total
  prop_down <- n_down / total
  
  case_when(
    prop_up >= threshold ~ "up",
    prop_down >= threshold ~ "down",
    TRUE ~ "mixed"
  )
}

read_pathway_file <- function(f) {
  raw <- read.delim(
    f,
    header = TRUE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(raw) < 8) {
    stop(paste("Could not parse file:", f))
  }
  
  names(raw)[1:8] <- c(
    "Term", "Overlap", "P-value", "Adjusted P-value",
    "Old P-value", "Old Adjusted P-value",
    "Odds Ratio", "Combined Score"
  )
  
  extra <- raw[, -(1:8), drop = FALSE]
  
  collapse_nonempty <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- str_trim(x)
    x <- x[x != ""]
    paste(x, collapse = ";")
  }
  
  if (ncol(extra) == 0) {
    raw$Genes <- NA_character_
    raw$`Genes Up` <- NA_character_
    raw$`Genes Down` <- NA_character_
  } else if (ncol(extra) == 1) {
    raw$Genes <- apply(extra[, 1, drop = FALSE], 1, collapse_nonempty)
    raw$`Genes Up` <- NA_character_
    raw$`Genes Down` <- NA_character_
  } else {
    genes_col <- extra[, 1, drop = FALSE]
    rest <- extra[, -1, drop = FALSE]
    
    half_point <- ceiling(ncol(rest) / 2)
    
    up_part <- rest[, seq_len(half_point), drop = FALSE]
    down_part <- if (ncol(rest) > half_point) {
      rest[, (half_point + 1):ncol(rest), drop = FALSE]
    } else {
      NULL
    }
    
    raw$Genes <- apply(genes_col, 1, collapse_nonempty)
    raw$`Genes Up` <- apply(up_part, 1, collapse_nonempty)
    raw$`Genes Down` <- if (!is.null(down_part)) {
      apply(down_part, 1, collapse_nonempty)
    } else {
      NA_character_
    }
  }
  
  raw %>%
    transmute(
      Term = Term,
      p.adj = as.numeric(`Adjusted P-value`),
      file = basename(f),
      comparison_group = get_comparison_group(basename(f)),
      celltype = get_celltype_from_filename(basename(f)),
      Genes = Genes,
      `Genes Up` = `Genes Up`,
      `Genes Down` = `Genes Down`,
      direction = purrr::map2_chr(
        `Genes Up`, `Genes Down`,
        ~ infer_direction_majority(.x, .y, threshold = 0.6)
      )
    )
}


# read files

all_pathways <- map_dfr(files, read_pathway_file)

print(unique(all_pathways$comparison_group))
print(unique(all_pathways$celltype))
print(table(all_pathways$direction, useNA = "ifany"))

#patwhays 
ed_terms <- c(
  "Glycogen Metabolism",
  "Integration of Energy Metabolism",
  "Glycolysis",
  "Detoxification of Reactive Oxygen Species",
  "Citric Acid Cycle (TCA Cycle)",
  "Cellular Response to Mitochondrial Stress",
  "Oxidative Stress Induced Senescence",
  "Pyruvate Metabolism",
  "Metabolism of Nucleotides",
  "NFE2L2 Regulating Anti-Oxidant Detoxification Enzymes",
  "Maturation of TCA Enzymes and Regulation of TCA Cycle",
  "Glucose Metabolism"
)

sds_terms <- c(
  "Glycogen Metabolism",
  "Detoxification of Reactive Oxygen Species",
  "Citric Acid Cycle (TCA Cycle)",
  "Cellular Response to Mitochondrial Stress",
  "Oxidative Stress Induced Senescence",
  "Pyruvate Metabolism",
  "Metabolism of Nucleotides",
  "NFE2L2 Regulating Anti-Oxidant Detoxification Enzymes",
  "Maturation of TCA Enzymes and Regulation of TCA Cycle",
  "Inositol Phosphate Metabolism",
  "Defects in Vitamin and Cofactor Metabolism",
  "Glutathione Conjugation"
)


ed_order <- c(
  "Glucose Metabolism",
  "TCA Cycle",
  "TCA Cycle Regulation",
  "Pyruvate Metabolism",
  "Glycolysis",
  "Glycogen Metabolism",
  "Energy Metabolism",
  "Nucleotide Metabolism",
  "Mitochondrial Stress Response",
  "Senescence By Oxidative Stress",
  "ROS Detoxification",
  "NFE2L2 anti-oxidant regulation"
)

sds_order <- c(
  "TCA Cycle",
  "TCA Cycle Regulation",
  "Pyruvate Metabolism",
  "Inositol Phosphate Metabolism",
  "Glycogen Metabolism",
  "Nucleotide Metabolism",
  "Mitochondrial Stress Response",
  "Senescence By Oxidative Stress",
  "ROS Detoxification",
  "NFE2L2 anti-oxidant regulation",
  "Vitamin & Cofactor Metabolism",
  "Glutathione Conjugation"
)

# filter / rename
plot_ED <- all_pathways %>%
  filter(comparison_group == "ED vs HC", Term %in% ed_terms) %>%
  rename_terms() %>%
  filter(Term %in% ed_order) %>%
  group_by(Term, celltype) %>%
  slice_min(order_by = p.adj, n = 1, with_ties = FALSE) %>%
  ungroup()

plot_SDS <- all_pathways %>%
  filter(comparison_group == "SDS vs HC", Term %in% sds_terms) %>%
  rename_terms() %>%
  filter(Term %in% sds_order) %>%
  group_by(Term, celltype) %>%
  slice_min(order_by = p.adj, n = 1, with_ties = FALSE) %>%
  ungroup()


plot_ED <- plot_ED %>%
  mutate(
    p.adj = pmax(p.adj, 1e-300),
    logp = -log10(p.adj),
    Term = factor(Term, levels = rev(ed_order)),
    celltype = factor(celltype, levels = c("HSC & MPP", "EMP", "EEP", "LEP")),
    direction = factor(direction, levels = c("up", "down", "mixed"))
  )

plot_SDS <- plot_SDS %>%
  mutate(
    p.adj = pmax(p.adj, 1e-300),
    logp = -log10(p.adj),
    Term = factor(Term, levels = rev(sds_order)),
    celltype = factor(celltype, levels = c("HSC & MPP", "EMP", "EEP", "LEP")),
    direction = factor(direction, levels = c("up", "down", "mixed"))
  )

#show overlapping celltypes
celltype_offset <- c(
  "HSC & MPP" = -0.18,
  "EMP" = -0.06,
  "EEP" = 0.06,
  "LEP" = 0.18
)

plot_ED <- plot_ED %>%
  mutate(logp_offset = logp + unname(celltype_offset[as.character(celltype)]))

plot_SDS <- plot_SDS %>%
  mutate(logp_offset = logp + unname(celltype_offset[as.character(celltype)]))


# keep log values on axis, use fewer clean ticks

min_logp <- floor(min(c(plot_ED$logp, plot_SDS$logp), na.rm = TRUE))
max_logp <- ceiling(max(c(plot_ED$logp_offset, plot_SDS$logp_offset), na.rm = TRUE))

x_limits <- c(min_logp, max_logp + 0.5)
x_breaks <- seq(ceiling(min_logp), floor(max_logp), by = 1)


# colors and shapes for celltypes

celltype_colors <- c(
  "HSC & MPP" = "#8fb0ff",
  "EMP" = "#004d43",
  "EEP" = "#0000a6",
  "LEP" = "#5a0007"
)

direction_shapes <- c(
  "up" = 24,
  "down" = 25,
  "mixed" = 23
)

#plots
bm_plot_ED <- ggplot(
  plot_ED,
  aes(x = logp_offset, y = Term, fill = celltype, shape = direction)
) +
  geom_segment(
    aes(x = logp, xend = logp_offset, y = Term, yend = Term, color = celltype),
    linewidth = 0.35,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  geom_point(
    size = 5.5,
    alpha = 0.95,
    stroke = 0
  ) +
  scale_x_continuous(
    limits = x_limits,
    breaks = x_breaks,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    name = "Cell type",
    values = celltype_colors,
    breaks = c("HSC & MPP", "EMP", "EEP", "LEP")
  ) +
  scale_color_manual(
    values = celltype_colors,
    guide = "none"
  ) +
  scale_shape_manual(
    name = "Direction",
    values = direction_shapes,
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 5, stroke = 0)
    ),
    shape = guide_legend(
      override.aes = list(fill = "grey60", size = 5, stroke = 0)
    )
  ) +
  labs(
    title = "ED vs HC",
    x = expression(-log[10]("adjusted p-value")),
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold")
  )

bm_plot_SDS <- ggplot(
  plot_SDS,
  aes(x = logp_offset, y = Term, fill = celltype, shape = direction)
) +
  geom_segment(
    aes(x = logp, xend = logp_offset, y = Term, yend = Term, color = celltype),
    linewidth = 0.35,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  geom_point(
    size = 5.5,
    alpha = 0.95,
    stroke = 0
  ) +
  scale_x_continuous(
    limits = x_limits,
    breaks = x_breaks,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    name = "Cell type",
    values = celltype_colors,
    breaks = c("HSC & MPP", "EMP", "EEP", "LEP")
  ) +
  scale_color_manual(
    values = celltype_colors,
    guide = "none"
  ) +
  scale_shape_manual(
    name = "Direction",
    values = direction_shapes,
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 5, stroke = 0)
    ),
    shape = guide_legend(
      override.aes = list(fill = "grey60", size = 5, stroke = 0)
    )
  ) +
  labs(
    title = "SDS vs HC",
    x = expression(-log[10]("adjusted p-value")),
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold")
  )


print(bm_plot_ED)
print(bm_plot_SDS)

